"""
Combine individual alignments into a composite FASTA file

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2020-05-02)
    v0.0.2 (2020-10-20)


"""

import argparse
import sys
import json
import re
import datetime
import os
import math, csv
from   os import  path
import operator
from progress.counter import Counter

from   Bio import SeqIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

arguments = argparse.ArgumentParser(description='Reduce sequences to differences from reference (first sequence in the file)')

arguments.add_argument('-i', '--input',      help = 'Read the FASTA from', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output',     help = 'Write the resulting JSON to', required = True, type = argparse.FileType('w'))
arguments.add_argument('-c', '--cutoff',     help = 'Only report variants reaching this count (if >= 1 or frequency if in [0-1])', required = False, type = float, default = 0.02)
arguments.add_argument('-a', '--annotation', help = 'A JSON file with annotation information', required = True, type = argparse.FileType('r'))

import_settings = arguments.parse_args()

summary_json    = {
                        'reference_base' : [],
                        'records' : [], # each record is an array [anonymized id, sample_date (YYYY-MM-DD), [location1,location2,location3,location4], [variants as pairs of coordinate, SNP]]
                        'non-human' : [] # same as in records except non-human
                  }

variant_counts      = {}
sequence_count      = 0
reference_position  = []

nonhuman_header = re.compile ('^NONHUMAN\:(.+)$')

record_info = json.load (import_settings.annotation)

def lookup_date (id):
    id = '_'.join (id.split ('_')[0:3])
    if id in record_info:
        if 'collected' in record_info[id]:
            try:
                return datetime.datetime.strptime(record_info[id]['collected'], "%Y%m%d").strftime ("%Y%m%d")
            except Exception as e:
                #print (e, file = sys.stderr)
                pass
            
    return None
    
    
    
def lookup_location (id):   
    id = '_'.join (id.split ('_')[0:3])
    if id in record_info:
        if 'location' in record_info[id] and record_info[id]['location']:
           return [record_info[id]['location']['subregion'],record_info[id]['location']['country'],record_info[id]['location']['state'],record_info[id]['location']['locality']]
            
    return None      
      
print ("Loading the FASTA", file = sys.stderr)        
        
with Counter('Loading records: ') as bar:
    for seq_record in SeqIO.parse(import_settings.input, "fasta"):
        seq_id   = seq_record.name
        seq = str (seq_record.seq).upper()
        if len (summary_json['reference_base']) == 0:
            seq_position  = 0
            for n in seq:
                if n == '-':
                    reference_position.append (-1)
                else:
                    reference_position.append (seq_position)
                    seq_position+=1
                    summary_json['reference_base'].append (n)
        else:
            if seq_id != 'CONSENSUS':
                m = nonhuman_header.match (seq_id)
                if m:
                    seq_id = m.group(1)
                    summary_json['non-human'].append (
                    [seq_id,
                        [(k,seq[i]) for (i,k) in enumerate (reference_position) if k >= 0 and summary_json['reference_base'][k] != seq[i] and seq[i] != '-']
                    ])
                else:
                    sequence_count += 1
                    summary_json ['records'].append (["%d" % sequence_count, lookup_date (seq_id), lookup_location (seq_id),
                        [(k,seq[i]) for (i,k) in enumerate (reference_position) if k >= 0 and summary_json['reference_base'][k] != seq[i] and seq[i] != '-']
                    ]);
        
                    for variant in summary_json ['records'][-1][3]:
                        if not variant in variant_counts:
                            variant_counts [variant] = 0
                        variant_counts [variant] += 1
                    bar.next()
                    
                    #print (seq_id, len(summary_json ['sequences'] [seq_id]), file = sys.stderr)

if import_settings.cutoff < 1.:
    import_settings.cutoff *= sequence_count
    
accepted_variants = {}
for v,c in variant_counts.items():
    if c >= import_settings.cutoff:
        accepted_variants[v] = c
        

summary_json['variants'] = [[k[0],k[1],c] for k, c in accepted_variants.items()]
        
filtered_variants = []
for i, variants in enumerate (summary_json ['records']):   
    filtered_variants.append ([variants[0], variants[1], variants[2], [i for k in variants[3] if k in accepted_variants  for i in k]])
    
summary_json ['records'] = filtered_variants

json.dump (summary_json, import_settings.output, indent = 1)
                
       