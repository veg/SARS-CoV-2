"""
Combine individual alignments into a composite FASTA file

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2020-05-02)


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

from   Bio import SeqIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

arguments = argparse.ArgumentParser(description='Reduce sequences to differences from reference (first sequence in the file)')

arguments.add_argument('-i', '--input', help = 'Read the FASTA from', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output', help = 'Write the resulting JSON to', required = True, type = argparse.FileType('w'))
arguments.add_argument('-c', '--cutoff', help = 'Only report variants reaching this count (if >= 1 or frequency if in [0-1])', required = False, type = float, default = 0.02)

import_settings = arguments.parse_args()

summary_json    = {
                        'reference_base' : [],
                        'sequences' : {}
                  }

variant_counts  = {}
sequence_count  = 0
reference_position = []
        
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
        sequence_count += 1
        summary_json ['sequences'][seq_id] = [(k,seq[i]) for (i,k) in enumerate (reference_position) if k >= 0 and summary_json['reference_base'][k] != seq[i] and seq[i] != '-']
        
        for variant in summary_json ['sequences'][seq_id]:
            if not variant in variant_counts:
                variant_counts [variant] = 0
            variant_counts [variant] += 1
        
        print (seq_id, len(summary_json ['sequences'] [seq_id]), file = sys.stderr)
        #sys.exit (0)

if import_settings.cutoff < 1.:
    import_settings.cutoff *= sequence_count
    
accepted_variants = {}
for v,c in variant_counts.items():
    if c >= import_settings.cutoff:
        accepted_variants[v] = c

summary_json['variants'] = [[k[0],k[1],c] for k, c in accepted_variants.items()]
        
filtered_variants = {}
for s, variants in summary_json ['sequences'].items():
    filtered_variants[s] = [k for k in variants if k in accepted_variants]
    
summary_json ['sequences'] = filtered_variants

json.dump (summary_json, import_settings.output, indent = 1)
                
       