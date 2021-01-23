"""
Inject date and location attributes into the variants file

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2021-01-12)


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


import itertools

date_parse_format = "%Y%m%d"


arguments = argparse.ArgumentParser(description='Compute pairwise distances between sequences from a variant file')

arguments.add_argument('-i', '--input', help  = 'Read the variants from', required = True, type = argparse.FileType('r'))
arguments.add_argument('-a', '--annotation', help  = 'Read annotation from', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output', help = 'Write annotated JSON to', required = True, type = argparse.FileType('w'))

import_settings = arguments.parse_args()

variants = json.load (import_settings.input)
annotation = json.load (import_settings.annotation)

annotated_json = {}
sequences = []

baseline = datetime.datetime.strptime (annotation['reference'],"%Y%m%d")

for seq, vars in variants["sequences"].items():
    seqID = int (seq)
    date = annotation['values'][0][annotation['meta'][seqID][0]]
    if date:
        date = (baseline + datetime.timedelta(days=date)).strftime ("%Y%m%d")
    location = annotation['values'][2][annotation['meta'][seqID][2]]
    seq_record = {'V' : vars, 'D' : date, 'L' : location}    
    sequences.append (seq_record)
    

annotated_json ['reference_base'] = variants['reference_base']
annotated_json ['variants'] = variants['variants']
annotated_json ['sequences'] = sequences
       
json.dump (annotated_json, import_settings.output)