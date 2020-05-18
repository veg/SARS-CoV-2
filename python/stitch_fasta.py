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
from   Bio import SeqIO
import operator


arguments = argparse.ArgumentParser(description='Summarize selection analysis results.')

arguments.add_argument('-d', '--dir',        help = 'Directory with files',    required = True, type = str )
#arguments.add_argument('-d', '--duplicates',   help = 'Duplicate JSON files',      required = True, type = argparse.FileType('r'), nargs = '*')
arguments.add_argument('-o', '--output',       help = 'Write the file to',         required = True, type = argparse.FileType('w'))

genes = ['ORF1a','ORF1b','S','ORF3a','M','ORF6','ORF7a','ORF8','N']

import_settings = arguments.parse_args()

combined_fasta  = {}

for i, gene in enumerate (genes):
    local_set = {}
    
    with open (path.join (import_settings.dir, "sequences.%s.duplicates.json" % gene), "r") as dh:
        dups = json.load (dh)
        
    count = 0
   
    with open (path.join (import_settings.dir, "sequences.%s.compressed.fas" % gene), "r") as fasta:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            seq_id   = seq_record.name
            seq = str (seq_record.seq)
            no_count = '_'.join (seq_id.split ('_')[0:-1])
            local_set [no_count] = seq
            count += len (dups[seq_id])
            for i,dup_id in  dups[seq_id].items():
                if dup_id != seq_id:
                    local_set [dup_id] = seq
            
        
    if len (combined_fasta):
        to_delete = set ()
        for i in combined_fasta:
            if i in local_set:
                combined_fasta[i] += local_set[i]
            else:
                to_delete.add (i)
                
        for i in to_delete:
            del combined_fasta[i]
                
    else:
        for i,s in local_set.items():
            combined_fasta[i] = s

    print (len (combined_fasta), count)
    #sys.exit (0)
    
for i, s in combined_fasta.items():
    print (">%s\n%s\n" % (i,s), file = import_settings.output)
