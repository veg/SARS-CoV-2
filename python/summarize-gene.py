"""
Summarize HyPhy selection analyses for a gene from JSON files

    SLAC (required)
    FEL  (required)
    MEME (required)
    PRIME (optional)

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2020-03-23)


"""

import argparse
import sys
import json
import re
import datetime
import os
import math, csv
from   os import  path
from   Bio import SeqIO
import operator
import collections 


arguments = argparse.ArgumentParser(description='Summarize selection analysis results.')

arguments.add_argument('-o', '--output', help = 'Write results here', type = argparse.FileType('w'), default = sys.stdout)
arguments.add_argument('-s', '--slac',   help = 'SLAC results file', required = True, type = argparse.FileType('r'))
arguments.add_argument('-f', '--fel',   help = 'FEL results file', required = True, type = argparse.FileType('r'))
arguments.add_argument('-m', '--meme',   help = 'MEME results file', required = True, type = argparse.FileType('r'))
arguments.add_argument('-p', '--prime',  help = 'PRIME results file', required = False, type = argparse.FileType('r'))
arguments.add_argument('-u', '--fubar',  help = 'FUBAR results file', required = False, type = argparse.FileType('r'))
arguments.add_argument('-P', '--pvalue',  help = 'p-value', required = False, type = float, default = 0.1)
arguments.add_argument('-c', '--coordinates',  help = 'An alignment with reference sequence (assumed to start with NC)', required = True, type = argparse.FileType('r'))

arguments.add_argument('-T', '--epitopes',  help = 'If provided, an epitope map in a JSON format', required = False, type = argparse.FileType('r'))
arguments.add_argument('-D', '--database', help ='Primary database record to extract sequence information from', required = True, type = argparse.FileType('r'))
arguments.add_argument('-d', '--duplicates', help ='The JSON file recording compressed sequence duplicates', required = True, type = argparse.FileType('r'))
arguments.add_argument('-M', '--MAF', help ='Also include sites with hapoltype MAF >= this frequency', required = False, type = float, default = 0.2)
arguments.add_argument('-E', '--evolutionary_annotation', help ='If provided use evolutionary likelihood annotation', required = False, type = argparse.FileType('r'))
arguments.add_argument('-F', '--evolutionary_fragment', help ='Used in conjunction with evolutionary annotation to designate the fragment to look up', required = False, type = str)
arguments.add_argument('-A', '--mafs', help ='If provided, write a CSV file with MAF/p-value tables', required = False, type = str)
arguments.add_argument('-V', '--evolutionary_csv', help ='If provided, write a CSV file with observed/predicted frequncies', required = False, type = str)
arguments.add_argument('-O', '--overall', help ='If provided, write site annotations to this JSON', required = False, type = str)
arguments.add_argument('-S', '--offset',  help ='Nucleotide position of the start of this gene in terms of reference genome', required = False, type = int, default = 0)
#arguments.add_argument('-V', '--variants', help ='If provided, write JSON-like VCF to this file', required = False, type = str)


import_settings = arguments.parse_args()

db = json.load (import_settings.database)
dups = json.load (import_settings.duplicates)

sequences_with_dates = {}
sequences_with_locations = {}
country_to_sub = {}

def get_location (v):
    if 'country' in v['location']:
        country_to_sub[v['location']['country']] = v['location']['subregion']
        return v['location']['country']
    if 'subregion' in v['location']:
        country_to_sub[v['location']['subregion']] = v['location']['subregion']
        return v['location']['subregion']
    return None


now      = datetime.datetime.now()
min_date = now
max_date = datetime.datetime (1900,1,1)

for id, record in db.items():
    try:
        date_check = datetime.datetime.strptime (record['collected'], "%Y%m%d")
        if date_check.year < 2019 or date_check.year == 2019 and date_check.month < 10 or date_check >= now: 
            continue
        if date_check < min_date:
            min_date = date_check
        if date_check > max_date:
            max_date = date_check
        sequences_with_dates[id] = record['collected']
        sequences_with_locations[id] = get_location (record)
    except Exception as e:
        pass
        
date_dups     = {}

maf_writer = None

if import_settings.mafs:
    try:
        maf_file = open (import_settings.mafs, "r+")
        maf_writer = csv.writer (maf_file)
        maf_file.seek (0,2)
    except FileNotFoundError as fnf:
        maf_file = open (import_settings.mafs, "w")
        maf_writer = csv.writer (maf_file)
        maf_writer.writerow (["Gene","Site","aa","count","freq","total"])
    
evo_writer = None  

if import_settings.evolutionary_csv:
    try:
        evo_file = open (import_settings.evolutionary_csv, "r+")
        evo_writer = csv.writer (evo_file)
        evo_file.seek (0,2)
    except FileNotFoundError as fnf:
        evo_file = open (import_settings.evolutionary_csv, "w")
        evo_writer = csv.writer (evo_file)
        evo_writer.writerow (["Gene","Site","Codon","Count","Observed","Predicted","Mostlikely"])
        

evo_annotation = None

if import_settings.evolutionary_annotation:
    evo_annotation = json.load (import_settings.evolutionary_annotation)
    if import_settings.evolutionary_fragment not in evo_annotation:
        evo_annotation = None
    else:
        evo_annotation = evo_annotation[import_settings.evolutionary_fragment]

annotation_json = None

if (import_settings.overall):
    try:
        with open (import_settings.overall,"r") as ann:
            try: 
                annotation_json = json.load (ann)
            except:
                annotation_json = {}
    except FileNotFoundError as fnf:
        annotation_json = {}
        


for seq, copies in dups.items():
    date_collection = {}
    location_collection = []
    for cp in copies.values():
        cpv = "_".join (cp.split ('_')[:3])
        if cpv in sequences_with_dates:
            cdate = sequences_with_dates[cpv]
            location = sequences_with_locations[cpv]
            
            tag = (cdate, location,db[cpv]['age'],db[cpv]['gender'])
 
            if not tag in date_collection:
                date_collection[tag] = 1
            else:
                date_collection[tag] += 1
            
                
    date_dups[seq] = date_collection   
 
slac = json.load (import_settings.slac)
fel  = json.load (import_settings.fel)
meme = json.load (import_settings.meme)


if import_settings.prime:
    prime  = json.load (import_settings.prime)
else:
    prime = None

if import_settings.fubar:
    fubar  = json.load (import_settings.fubar)
else:
    fubar = None


if import_settings.epitopes:
    epitopes  = json.load (import_settings.epitopes)
else:
    epitopes = None

    
ref_seq_map = None
ref_seq_re = re.compile ("^NC")

for seq_record in SeqIO.parse(import_settings.coordinates, "fasta"):
    seq_id   = seq_record.description
    if ref_seq_re.search (seq_id):
        ref_seq = str(seq_record.seq).upper()
        i = 0
        c = 0
        ref_seq_map = []
        while i < len (ref_seq):
            if ref_seq[i:i+3] != '---':
                ref_seq_map.append (c)
                c += 1
            else:
                ref_seq_map.append (-1)

            i+=3
        break
        
if ref_seq_map is None:
    raise Exception ("Misssing reference sequence for coordinate mapping")

    
# compile the list of sites that are under selection by either MEME or FEL

def map_site_to_genome_nuc (site): #0-based nucleotide coordinate for ref genome
    rsc = ref_seq_map[site]
    if rsc >= 0:
        return rsc*3 + import_settings.offset
    return -1

site_list = {}
sequences = slac["input"]["number of sequences"]
sites = slac["input"]["number of sites"]
tree = slac["input"]["trees"]["0"]
branch_lengths = {}

L = 0

variants_by_site       = [{} for k in range (sites)]
aa_variants_by_site    = [{} for k in range (sites)]
counts_by_site         = [{} for k in range (sites)]
aa_counts_by_site      = [{} for k in range (sites)]

def compute_site_MAF (site, source = None):
    variants = source[site] if source else counts_by_site [site]
    if len (variants):
        total = sum (variants.values())
        majority = max (variants.values()) / total
        return 1-majority
    return 0

def compute_site_entropy (site, source = None):
    variants = source[site] if source else counts_by_site [site]
    total = sum (variants.values())
    return -sum ([k/total * math.log (k/total,2) for k in variants.values()])
        
    
aa_letters = set ("ACDEFGHIKLMNPQRSTVWY")

def add_annotation_to_site (site, annotation):
    gs = map_site_to_genome_nuc (site)
    if gs >= 0:
        gs = "%d" % gs
        if not gs in annotation_json:
            annotation_json [gs] = {}
        for k, v in annotation.items():
            annotation_json [gs][k] = v
    else:
        print ("Site %d is unmapped" % site, file = sys.stderr)
            

for b,v in slac["tested"]["0"].items():
    branch_lengths[b] = slac["branch attributes"]["0"][b]["Global MG94xREV"]
    if v == "test":
        L += slac["branch attributes"]["0"][b]["Global MG94xREV"]
    else:
        for k in range (sites):
            codon = slac["branch attributes"]["0"][b]["codon"][0][k]
            if codon != '---':
                if codon not in variants_by_site[k]:
                    variants_by_site[k][codon] = 1
                    counts_by_site[k][codon] = len (dups[b])
                else:
                    variants_by_site[k][codon] += 1
                    counts_by_site[k][codon] += len (dups[b])
                aa = slac["branch attributes"]["0"][b]["amino-acid"][0][k]
                if aa in aa_letters:
                    if aa not in aa_variants_by_site[k]:
                        aa_variants_by_site[k][aa] = 1
                        aa_counts_by_site[k][aa] = len (dups[b])
                    else:
                        aa_variants_by_site[k][aa] += 1
                        aa_counts_by_site[k][aa] += len (dups[b])
                
if annotation_json is not None: # if this is specified, write everything out
    for k,v in enumerate (aa_counts_by_site):
        add_annotation_to_site (k, {'cdn' : counts_by_site[k], 'aa' : v})
            
        
variant_count_total = 0
variant_count_NS    = 0
        
valid_nucs = set (["A","C","G","T"])       
include_in_annotation = {}
    
for i, row in enumerate (fel["MLE"]["content"]["0"]):

    if annotation_json is not None: # if this is specified, write everything out
        add_annotation_to_site (i, {'FEL' : {
                    'a' : row[0],
                    'b' : row[1],
                    'p' : row[4]
            }})
                        
    if row[0] + row[1] > 0:
        maf = compute_site_MAF (i)
        if maf_writer:
            ##print (variants_by_site[i], file = sys.stderr)
            total = sum (aa_counts_by_site[i].values())
            for aa, count in aa_counts_by_site[i].items():
                maf_writer.writerow ([import_settings.evolutionary_fragment, "%d" % (ref_seq_map[i] + 1), aa, "%d" % count, "%g" % (count/total), "%d" % total])
            #maf_writer.writerow ([import_settings.evolutionary_fragment, "%d" % (ref_seq_map[i] + 1), "%g" % maf, "%g" % compute_site_MAF (i, aa_counts_by_site), "%g" % compute_site_entropy (i), "%g" % compute_site_entropy (i,aa_counts_by_site), "%g" % meme["MLE"]["content"]["0"][i][6]])
        
        if evo_writer and evo_annotation:
            check_key = "%d" % ref_seq_map[i]
            if evo_annotation and check_key in evo_annotation:
                total = sum (counts_by_site[i].values())
                for codon, freq in counts_by_site[i].items():
                    if codon[0] in valid_nucs and codon[1] in valid_nucs and codon[2] in valid_nucs:
                        evo_writer.writerow ([import_settings.evolutionary_fragment, "%d" % (ref_seq_map[i] + 1), codon, "%g" % freq, "%g" % (freq/total), 
                                              "%g" % (evo_annotation[check_key][codon] if codon in evo_annotation[check_key] else 1e-8), max(evo_annotation[check_key].items(), key=operator.itemgetter(1))[0] ])
            
            
            
        if row[4] < import_settings.pvalue :
            site_list[i] = {'fel' : row[4], 'kind' : 'positive' if row[1] > row[0] else 'negative', 'MAF' : maf}
        else:
            if maf >= import_settings.MAF:
                site_list[i] = {'fel' : row[4],  'MAF' : maf}
        
    
for i, row in enumerate (meme["MLE"]["content"]["0"]):
    if annotation_json is not None: # if this is specified, write everything out
        add_annotation_to_site (i, {'MEME' : {
                            'p'  : row[6],
                            'a'  : row[0],
                            'b+' : row[3],
                            'w+' : row[4],
                            'b-' : row[1],
                            'w-' : row[2],
                            'br' : row [7]
                        }})

    if row[6] < import_settings.pvalue or i in site_list:
        if i in site_list:
            site_list[i]['meme'] = row[6]
            site_list[i]['meme-fraction'] = row[4]
        else:
            site_list[i] = {'meme' : row[6], 'fel' : fel["MLE"]["content"]["0"][i][4], 'meme-fraction' : row[4], 'MAF' : compute_site_MAF (i)}
            
def newick_parser(nwk_str, bootstrap_values):
    clade_stack = []
    automaton_state = 0
    current_node_name = ""
    current_node_attribute = ""
    current_node_annotation = ""
    quote_delimiter = None
    name_quotes = {
      "'": 1,
      '"': 1
    }
    
    def add_new_tree_level():
      new_level = {
        "name": None
      };
      the_parent = clade_stack[len(clade_stack) - 1]
      if (not "children" in the_parent):
        the_parent["children"] = [];
      
      clade_stack.append (new_level);
      the_parent["children"].append(clade_stack[len(clade_stack) - 1]);
      clade_stack[len(clade_stack)-1]["original_child_order"] = len(the_parent["children"])
    

    def finish_node_definition():
      nonlocal current_node_name
      nonlocal current_node_annotation
      nonlocal current_node_attribute
      
      this_node = clade_stack.pop()
      if (bootstrap_values and "children" in this_node):
        this_node["bootstrap_values"] = current_node_name
      else:
        this_node["name"] = current_node_name
      
      this_node["attribute"] = current_node_attribute
      this_node["annotation"] = current_node_annotation
      current_node_name = ""
      current_node_attribute = ""
      current_node_annotation = ""
    

    def generate_error(location):
      return {
        'json': None,
        'error':
          "Unexpected '" +
          nwk_str[location] +
          "' in '" +
          nwk_str[location - 20 : location + 1] +
          "[ERROR HERE]" +
          nwk_str[location + 1 : location + 20] +
          "'"
      }


    tree_json = {
      "name" : "root"
    }
    
    clade_stack.append(tree_json);

    space = re.compile("\s")

    for char_index in range (len(nwk_str)):
      try:
        current_char = nwk_str[char_index]
        if automaton_state == 0:
           #look for the first opening parenthesis
           if (current_char == "("):
              add_new_tree_level()
              automaton_state = 1
        elif automaton_state == 1 or automaton_state == 3:
            #case 1: // name
            #case 3: { // branch length
            #reading name
            if (current_char == ":"):
              automaton_state = 3;
            elif current_char == "," or current_char == ")":
              try:
                finish_node_definition()
                automaton_state = 1
                if (current_char == ","):
                  add_new_tree_level()
              except Exception as e:
                return generate_error(char_index)
              
            elif (current_char == "("):
              if len(current_node_name) > 0:
                return generate_error(char_index);
              else:
                add_new_tree_level()
              
            elif (current_char in name_quotes):
              if automaton_state == 1 and len(current_node_name) == 0 and len (current_node_attribute) == 0 and len (current_node_annotation) == 0:
                automaton_state = 2
                quote_delimiter = current_char
                continue
              return generate_error(char_index)
            else:
              if (current_char == "["):
                if len (current_node_annotation):
                  return generate_error(char_index)
                else:
                  automaton_state = 4
              else:
                if (automaton_state == 3):
                  current_node_attribute += current_char;
                else:
                  if (space.search(current_char)):
                    continue;
                  if (current_char == ";"):
                    char_index = len(nwk_str)
                    break
                  current_node_name += current_char;
        elif automaton_state == 2: 
            # inside a quoted expression
            if (current_char == quote_delimiter):
              if (char_index < len (nwk_str - 1)):
                if (nwk_str[char_index + 1] == quote_delimiter):
                  char_index+=1
                  current_node_name += quote_delimiter;
                  continue;

              quote_delimiter = 0
              automaton_state = 1
              continue
            else:
              current_node_name += current_char;
        elif automaton_state == 4:
           ##inside a comment / attribute
            if (current_char == "]"):
              automaton_state = 3
            else:
              if (current_char == "["):
                return generate_error(char_index);
              current_node_annotation += current_char;
      except Exception as e:
        return generate_error(char_index);

    if (len (clade_stack) != 1):
      return generate_error(len (nwk_str) - 1);

    if (len (current_node_name)):
        tree_json['name'] = current_node_name;

    return {
      'json': tree_json,
      'error': None
    }
    

node_parents = {}
root_node_name = None

def recurse_tree (node, last = None):  
    node_parents[node['name']] = last
    if 'children' in node:
        for n in node['children']:
            recurse_tree (n, node['name'])
            

recurse_tree (newick_parser (tree, False)['json'])

if annotation_json and fubar is not None:
    hpd = 0.95
    grid = fubar["grid"]
    for site, p in fubar['posterior']["0"].items():
        posteriors = [(i,v) for i,v in enumerate (p[0])]
        hdpc = sorted(posteriors, key=operator.itemgetter(1), reverse = True)
        psum = 0.
        i = 0
                
        alpha = sum ([k[0] * posteriors[i][1] for i, k in enumerate (grid)])
        beta  = sum ([k[1] * posteriors[i][1] for i, k in enumerate (grid)])
        ppp = sum ([posteriors[i][1] for i, k in enumerate (grid) if k[0]<k[1]])
        ppn = sum ([posteriors[i][1] for i, k in enumerate (grid) if k[0]>k[1]])
        non0 = sum ([posteriors[i][1] for i, k in enumerate (grid) if k[0] > 0])
        hpd0 = hpd * non0
        mean_omega  = sum ([k[1]/k[0] * posteriors[i][1] for i, k in enumerate (grid) if k[0]>0]) / non0
        
        omega = []
        
        while (psum <= hpd0):
            idx = hdpc[i][0]
            if grid[idx][0] > 0:
                psum += hdpc[i][1]
                omega.append (grid[idx][1]/grid[idx][0])
            i+=1
            
            
        add_annotation_to_site (int(site), {"fubar" : {
            "a" : alpha,
            "b" : beta,
            "p+" : ppp,
            "p-" : ppn,
            "w" : mean_omega,
            "wl" : min (omega),
            "wu" : max (omega)
        }})


def compute_JH (timing, min_date, max_date):
    #print (timing, file = sys.stderr)
    residue_counts = {}
    mafs_by_date   = {}
    date_cutoff    = min_date + datetime.timedelta(days = 45)
    
    for residue, dates in timing.items():
        residue_counts [residue ] = 0
        for key, value in dates.items():
           #this_date =  datetime.datetime.strptime (key[0], "%Y%m%d")
           this_date = key[0]
           if this_date not in mafs_by_date:
              mafs_by_date [this_date]= {}
           if residue not in mafs_by_date[this_date]:
                mafs_by_date [this_date][residue] = 0
           mafs_by_date[this_date][residue] += value
           if datetime.datetime.strptime (key[0], "%Y%m%d") <= date_cutoff:
                residue_counts[residue] += value
            
    consensus = max(residue_counts.items(), key=operator.itemgetter(1))[0]
    mafs = []
    
    for date, counts in mafs_by_date.items():
        all = sum (counts.values())
        minority = sum ([v for k, v in counts.items() if k != consensus])
        #print (all, minority, counts, [v for k, v in counts.items() if k != consensus], file = sys.stderr)
        mafs.append ([date, minority/all])
  
    bin_count = math.ceil ((max_date - min_date).days/10)
    values_by_bins = [[] for k in range (bin_count)]
    unique_values = set ()
    for v in mafs:
        bin = (datetime.datetime.strptime (v[0], "%Y%m%d")  - min_date).days // 10
        values_by_bins[bin].append (v[1])
        unique_values.add (v[1])
    
    contingency_table = [[0 for k in range (bin_count)] for v in unique_values]
    value_to_index = {}
    for i,v in enumerate (sorted (list(unique_values))):
        value_to_index[v] = i
    
    for i,bin in enumerate (values_by_bins):
        for v in bin:
            contingency_table[value_to_index[v]][i] += 1
           
    value_count = len (unique_values) 
    if value_count == 1:
        return 0.
    row_sums = [sum (row) for row in contingency_table]
    column_sums = [sum ([row[j] for row in contingency_table]) for j in range (bin_count)]
    N = sum (row_sums)
    P = 0
    Q = 0
    
    r3 = sum ([k*k*k for k in row_sums])
    r2 = sum ([k*k for k in row_sums])
    c3 = sum ([k*k*k for k in column_sums])
    c2 = sum ([k*k for k in column_sums])
    
    
    varS = (2.*(N*N*N-r3-c3)+3.*(N*N-r2-c2)+5*N) / 18. + (r3-3*r2+2*N)*(c3-3*c2+2*N)/(9.*N*(N-1)*(N-2)) + (r2-N)*(c2-N)/(2.*N*(N-1))
 
   
    for i,bin in enumerate (values_by_bins):
        for v in bin:
            for j in range (i+1,bin_count):
                for e in values_by_bins[j]:
                    if e < v:
                        Q += 1
                    elif e > v:
                        P += 1
            
    Z = (P-Q) / math.sqrt (varS)            
    return Z
    
site_to_epitope = {}
    
if epitopes and annotation_json:
    this_epitope = collections.deque ()
    site_range = collections.deque ()
    for node,value in slac["branch attributes"]["0"].items():
        if node in slac["tested"]["0"] and slac["tested"]["0"][node] != "test":
            for site in range(sites):
               if "amino-acid" in value:
                    aa_value    = value["amino-acid"][0][site]
                    if len (aa_value) == 1:
                        if aa_value != '-':
                            if len (this_epitope) < 9:
                                this_epitope.append (aa_value)
                                site_range.append (site)
                            else:
                                this_epitope.popleft()
                                this_epitope.append (aa_value)
                                site_range.popleft()
                                site_range.append (site)
                    else:
                        this_epitope.clear()
                        site_range.clear()
                    
               if len (this_epitope) == 9:
                    str_epitope = ''.join (this_epitope)
                    if str_epitope in epitopes:
                        for si,s in enumerate (site_range):
                            if not s in site_to_epitope:
                                site_to_epitope[s] = {}
                                
                            if not str_epitope in site_to_epitope[s]:
                                site_to_epitope[s][str_epitope] = {}
                                
                            if not si in site_to_epitope[s][str_epitope]:
                                site_to_epitope[s][str_epitope][si] = 0
                            
                            site_to_epitope[s][str_epitope][si] += len (dups[node])

            
        

for site in range(sites) if annotation_json else site_list:
    if site in site_list:
        site_list[site]['meme-branches'] = meme["MLE"]["content"]["0"][site][7]
        site_list[site]['substitutions'] = [slac["MLE"]["content"]["0"]['by-site']['RESOLVED'][site][2],slac["MLE"]["content"]["0"]['by-site']['RESOLVED'][site][3]]
        labels      = {}
    evo_composition = {}
    timing      = {}
    ''' 
        for each amino acid, this will record "date" : count for when they were sampled
        timing -> 
            "residue" ->
                "date" -> count
    '''
    check_key           = "%d" % ref_seq_map[site]
    evo_site_annotation = evo_annotation[check_key] if evo_annotation and check_key in evo_annotation else None
    
    substitutions       = None
        
    for node,value in slac["branch attributes"]["0"].items():
        if root_node_name is None and node not in node_parents:
            root_node_name = node
            for k in node_parents:
                if node_parents[k] == 'root':
                    node_parents[k] = node
            
        if "amino-acid" in value:
            aa_value    = value["amino-acid"][0][site]
            codon_value = value["codon"][0][site]
            
            if len (aa_value) == 1 and evo_site_annotation:
                if codon_value not in evo_composition:
                    try:
                        evo_composition[codon_value] = {
                            "support" : evo_annotation[check_key][codon_value], "aa" : aa_value
                            }
                    except Exception as e:
                        evo_composition[codon_value] = {
                            "support" : 0.0, "aa" : aa_value
                            }
                            
            
            if node in date_dups:
                if aa_value not in timing:
                    timing [aa_value] = {}
                for dt, cnt in date_dups[node].items():
                    if not dt in timing [aa_value]:
                        timing [aa_value][dt] = cnt
                    else:
                        timing [aa_value][dt] += cnt
             
            if  value["nonsynonymous substitution count"][0][site] + value["synonymous substitution count"][0][site]:
                if substitutions is None:
                    substitutions = {
                        'cdn' : {},
                        'aa'  : {},
                        'lcdn' : {},
                        'laa' :  {}
                    }
                cdn_pair = "%s|%s" % (value["codon"][0][site], slac["branch attributes"]["0"][node_parents[node]]["codon"][0][site])
                aa_pair = "%s|%s" % (value["amino-acid"][0][site], slac["branch attributes"]["0"][node_parents[node]]["amino-acid"][0][site])
                if slac["tested"]["0"][node] == "test":
                    ks = ('cdn', 'aa')
                else:
                    ks = ('lcdn','laa')
                if cdn_pair in substitutions[ks[0]]:
                    substitutions[ks[0]][cdn_pair] += 1
                else:
                    substitutions[ks[0]][cdn_pair] = 1
                if aa_pair in substitutions[ks[1]]:
                    substitutions[ks[1]][aa_pair] += 1
                else:
                    substitutions[ks[1]][aa_pair] = 1

             
            if site in site_list:
                labels[node] = [aa_value,value["codon"][0][site],value["nonsynonymous substitution count"][0][site],value["synonymous substitution count"][0][site]]
    
    if site in site_list:
        site_list[site]['composition'] = aa_counts_by_site[site]
        site_list[site]['labels'] = labels
        if len (evo_composition):
            site_list[site]['evolutionary_support']     = evo_composition
            site_list[site]['evolutionary_predictions'] = evo_site_annotation
            site_list[site]['counts'] = counts_by_site 
            print ("Site %d" % site, file = sys.stderr)
            print ("Codon\tAmino-Acid\tFrequency\tSupport", file = sys.stderr)
            all_count = sum (counts_by_site[site].values())
            for codon, support in evo_composition.items():
                codon_freq = [v for k,v in counts_by_site[site].items() if k == codon]
                if len (codon_freq):
                    print ("%s\t%s\t%.3g\t%.3g" % (codon,support["aa"],[v for k,v in counts_by_site[site].items()if k == codon][0]/all_count,support["support"]), file = sys.stderr)
            #print (site, evo_composition,aa_counts_by_site[site], file = sys.stderr)
    
    timing_as_array = {}
    
    for aa, t in timing.items():
        timing_as_array [aa] = [[k[0],k[1],country_to_sub[k[1]],v,k[2],k[3]] for k, v in timing[aa].items()]
    
    if site in site_list:
        site_list[site]['timing'] = timing_as_array
              
    jh_z = compute_JH (timing, min_date, max_date)

    if site in site_list:
        site_list[site]['trend'] = jh_z
        
    if annotation_json:
        add_annotation_to_site (site, {"trend" : jh_z})
        if substitutions is not None:
            add_annotation_to_site (site, {"subs" : substitutions})
        if site_to_epitope and site in site_to_epitope:
            add_annotation_to_site (site, {"hla" : site_to_epitope[site]} )
        
    if prime and site in site_list:
        site_list[site]['prime'] = []
        prime_row = prime["MLE"]["content"]["0"][site]
        prime_headers = prime["MLE"]["headers"]
        for idx in [5,7,10,13,16,19]:
            if prime_row[idx] <= import_settings.pvalue :
               if idx == 5:
                    site_list[site]['prime'].append (['Overall', prime_row[idx], 0])
               else:
                    site_list[site]['prime'].append ([prime_headers[idx][1].replace ('p-value for non-zero effect of ',''), prime_row[idx], prime_row[idx-1]])
                    
        
json_out = {
    'sequences' : sequences,
    'bl' : branch_lengths,
    'total sequences' : sum ([len (k) for k in dups.values()]),
    'aminoacid variant sites' : [v for v in aa_variants_by_site if len (v) > 1 and len ([c for c in v.values() if c>1]) > 1],
    'all variant sites' : [v for v in variants_by_site if len (v) > 1 and len ([c for c in v.values() if c>1]) > 1],
    'any variation' : len ([v for v in variants_by_site if len (v) > 1]),
    'sites' : sites,
    'tree' : tree,
    'MAF' : import_settings.MAF,
    'L' : L,
    'p' : import_settings.pvalue,
    'selection' : site_list,
    'map' : ref_seq_map,
    'dN/dS' : {
        'internal' :  meme["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *test*"],
        'leaves'   :  meme["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *background*"]
    }
}

json.dump (json_out, import_settings.output, sort_keys = True, indent = 1)
    
if annotation_json is not None:
    with open (import_settings.overall, "w") as ann:
        json.dump (annotation_json, ann, indent = 1)

