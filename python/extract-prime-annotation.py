"""
Extract evolutionary annotation data from PRIME JSON and an alignment to REMAP coordinates

    SLAC (required)
    FEL  (required)
    MEME (required)
    PRIME (optional)

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2020-04-10)


"""

import argparse
import sys
import json
import re
import datetime
import os
import math
from   os import  path
from   Bio import SeqIO
import operator


arguments = argparse.ArgumentParser(description='Summarize selection analysis results.')

arguments.add_argument('-o', '--output', help = 'Write results here', type = str,  required = True)
arguments.add_argument('-p', '--prime',  help = 'PRIME results file', required = False, type = argparse.FileType('r'))
arguments.add_argument('-c', '--coordinates',  help = 'An alignment with reference sequence (assumed to start with NC)', required = True, type = argparse.FileType('r'))
arguments.add_argument('-r', '--region',  help = 'Description of the region', required = True, type = str)
arguments.add_argument('-f', '--offset',  help = 'Offset within the region', required = True, type = int)


import_settings = arguments.parse_args()

prime  = json.load (import_settings.prime)
    
ref_seq_map = None
ref_seq_re = re.compile ("^NC")

# site -> list of codons that have evolutionary support for this sequence
ref_seq_p = 0
ref_seq_id = None


json_out = {}
try:
    output_file = open (import_settings.output, 'r+')
except FileNotFoundError as e:
    output_file = open (import_settings.output, 'w')

try:
    json_out = json.load (output_file)
except Exception as e:
    pass


for seq_record in SeqIO.parse(import_settings.coordinates, "fasta"):
    seq_id   = seq_record.description
    if ref_seq_re.search (seq_id):
        ref_seq = str(seq_record.seq)
        ref_seq_id = seq_id.upper().upper()
        i = 0
        c = 0
        ref_seq_map = {}
        while i < len (ref_seq):
            #ref_seq_map.append (c)
            if ref_seq[i:i+3] != '---':
                ref_seq_map[c] = ref_seq_p
                ref_seq_p += 1
            c += 1
            i+=3
        break
        
for k in prime["tested"]["0"]:
    if k.upper() == ref_seq_id:
        ref_seq_id = k
        break
        
if ref_seq_map is None:
    raise Exception ("Misssing reference sequence for coordinate mapping")


if import_settings.region not in json_out:
    json_out [import_settings.region] = {}

current_segment = 0
current_offset  = 0

#print (len(prime["MLE"]["Imputed States"]["0"]))

for alignment_coord, seq_coord in ref_seq_map.items():
    #print (alignment_coord, current_offset)
    if alignment_coord - current_offset >= len (prime["MLE"]["Imputed States"]["%d" % current_segment]):
        current_offset  += len (prime["MLE"]["Imputed States"]["%d" % current_segment])
        current_segment += 1
        
    site_info = prime["MLE"]["Imputed States"]["%d" % current_segment]['%d' % (alignment_coord-current_offset)]
    if site_info:
        json_out [import_settings.region][import_settings.offset+seq_coord] = site_info[ref_seq_id]
    
output_file.seek(0)  
json.dump (json_out, output_file, indent = 1)
    

