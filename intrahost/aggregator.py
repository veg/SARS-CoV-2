import urllib.request
import csv
from io import StringIO

import argparse
import operator
import sys
import json
import datetime
import os
from math import exp

arguments = argparse.ArgumentParser(description='Aggregate Galaxy SARS-CoV-2 AV data from a master JSON file into variant and metadata annotation for consumption by Observable Dashboards')

arguments.add_argument('-c', '--columns',   help = 'Columns to extract', required = True, type = str, nargs = '*')
arguments.add_argument('-u', '--url', help ='Master URL', required = True, type = str)
arguments.add_argument('-p', '--project', help ='If provided, include only runs that match to the Bioproject[s]; comma separated list; if the first character is !, include everything NOT in the list', required = False, type = str, default = "")
arguments.add_argument('-a', '--annotation', help ='Output metadata annotation to this TSV file', required = True, type=argparse.FileType('w'))
arguments.add_argument('-d', '--diagnostics', help ='Output batch-level diagnostics to this TSV file', required = True, type=argparse.FileType('w'))
arguments.add_argument('-P', '--position_field', help ='The name of variant table column which contains genomic variant position', required = False, type=str, default = 'POS')
arguments.add_argument('-A', '--af_field', help ='The name of variant table column which contains AF', required = False, type=str, default = 'AF')
arguments.add_argument('-f', '--freq', help ='Maximim AF to consider for Poisson filtering', required = False, type=float, default = 0.25)

import_settings = arguments.parse_args()

writer_annotation = csv.writer   (import_settings.annotation,  delimiter = '\t')
writer_annotation.writerow (['run_accession','collection_date','completed_date'])

diag_annotation = csv.writer   (import_settings.diagnostics,  delimiter = '\t')
diag_annotation.writerow (['batch','N','error_rate','poisson_cut'])

exclude = False
filter_chain = []

for i, term in enumerate (import_settings.project.split (',')):
    term = term.strip()
    if len (term):
        if i == 0:
            if term[0] == '!':
                term = term[1:]
                exclude = True
        if len (term):
            filter_chain.append (lambda x: x == term)
                                     
                   
if len (filter_chain) == 0:
    if exclude:
        filter_chain = [lambda x : False]
    else:
        filter_chain = [lambda x : True]

    

writer = csv.writer (sys.stdout, delimiter = '\t')
keep_in = {}

headers = []
shift = 0
for i,n in enumerate (import_settings.columns):
    pieces = n.split ('=')
    if len (pieces) == 1:
        keep_in [n] = (i - shift, )
        headers.append (n)
    else:
        keep_in [pieces[0]] = (i, pieces[1])
        shift += 1

writer.writerow (headers)

def poissonSupport (MLL, support, N):
  el = exp (-MLL);
  lkf = 1;
  points =  [k for k in range (len (support)+1)]
  i = 0;
  sp = []
  for p in points:
      sp.append (el*lkf*N)
      i += 1
      lkf *= MLL/i;
      
  return sp


genome_length = 29903
with urllib.request.urlopen(import_settings.url) as response:
    master_list = json.loads (response.read())
    
    for key, data in master_list.items():
        study = data["study_accession"]
        if (exclude and True in [f(study) for f in filter_chain]) or (not exclude and not True in [f(study) for f in filter_chain]):
            print ("Skipping %s because of project filtering settings" % study, file = sys.stderr)
            continue
        else:   
            print ("Using %s because of project filtering settings" % study, file = sys.stderr)
            
        for i, accession in enumerate (data['samples']):
            writer_annotation.writerow ([accession,data['collection_dates'][i],datetime.datetime.strptime(data['time'],"%Y-%m-%dT%H:%M:%S.%f").strftime ("%Y-%m-%d")])
            
        url = data["report"]["datamonkey_link"]
        N = len (data["samples"])
        print ("Fetching %s (%d samples) " % (url, N), file = sys.stderr)
        with urllib.request.urlopen(url) as response:
        
            #count how often a position has mutations
            hit_count = {}
        
            html = response.read()
            reader = csv.reader (StringIO (html.decode("utf-8")), delimiter = '\t')
            names = next (reader)
            pos_idx = names.index (import_settings.position_field)
            af_idx = names.index (import_settings.af_field)
            field_mapper = {}
            field_filter = []
            for i,n in enumerate(names):
                if n in keep_in:
                    if len (keep_in[n]) == 1:
                        field_mapper[i] = keep_in[n][0]
                    else:
                        field_filter.append ([i,keep_in[n][1]])
        
            for row in reader:
                new_row = ['' for k in keep_in]
                for i,o in field_mapper.items():
                    new_row[o] = row[i]
                if len ([k for k in field_filter if row[k[0]] != k[1]]) > 0:
                    continue
                    
                if float (row[af_idx]) <= import_settings.freq:
                    pos = int (row[pos_idx])
                    if pos in hit_count:
                        hit_count[pos] += 1
                    else:
                        hit_count[pos] = 1
                    
                
                writer.writerow (new_row)   
                
            #convert position=>hits into hits=># of positions
            diag_row = [data['batch_id'],str (N), 'N/A', 'N/A']
            if N >= 10:
                max_hits = max(hit_count.items(), key=operator.itemgetter(1))[1]
                hit_array = [0 for k in range (max_hits+1)]
                hit_array [0] = genome_length - len (hit_count)
                for p, v in hit_count.items():
                    hit_array[v] += 1
                ml_lambda = sum ([i*k for i, k in enumerate (hit_array)]) / genome_length
                expected_counts = poissonSupport (ml_lambda, hit_array, genome_length)
                cutoff = 1
                for i,k in enumerate (hit_array):
                    if i > len (expected_counts): break
                    if k > 0:
                        if (expected_counts[i]-k) / expected_counts[i] < -1:
                            break
                        else:
                            cutoff += 1
                            
                diag_row [2] = str (ml_lambda/genome_length)
                diag_row [3] = str (cutoff) 
                
                #print (expected_counts, file = sys.stderr)
                #print (hit_array,  file = sys.stderr)
                
                print ("Batch ID %s: Poisson cutoff %d, error rate %g" % (data['batch_id'], cutoff, ml_lambda/genome_length), file = sys.stderr)
                
            diag_annotation.writerow (diag_row)