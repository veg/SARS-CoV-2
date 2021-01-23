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


arguments = argparse.ArgumentParser(description='CTL converter.')

arguments.add_argument('-t', '--tsv',        help = 'Import TSV CTL annotation',    required = True, type = argparse.FileType('r') )
arguments.add_argument('-o', '--output',       help = 'Write the JSON file to',         required = True, type = argparse.FileType('w'))

import_settings = arguments.parse_args()

reader = csv.reader (import_settings.tsv, delimiter = '\t')
epitopes = {}
headers = next (reader)
for line in reader:
    seq = line[3] # epitope sequence
    allele = line[1] 
    
    if seq not in epitopes:
        epitopes[seq] = {}
        
    epitopes[seq][allele] = float (line[4])

print ("Imported %d epitopes" % (len (epitopes)))
json.dump (epitopes, import_settings.output, indent = 1)