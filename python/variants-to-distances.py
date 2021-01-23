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

import itertools

date_parse_format = "%Y%m%d"


arguments = argparse.ArgumentParser(description='Compute pairwise distances between sequences from a variant file')

arguments.add_argument('-i', '--input', help  = 'Read the variants from', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output', help = 'Write the resulting TSV to', required = True, type = argparse.FileType('w'))
arguments.add_argument('-c', '--cutoff', help = 'Only report distances that are no greater than this', required = False, type = float, default = 0.02)
arguments.add_argument('-a', '--annotation', help  = 'Read annotation from', required = True, type = argparse.FileType('r'))
arguments.add_argument('-t', '--attributes', help = 'Write the attribute CSV to', required = True, type = argparse.FileType('w'))

import_settings = arguments.parse_args()

variants = json.load (import_settings.input)
annotation = json.load (import_settings.annotation)



N = len(variants["reference_base"])

variants_as_sets = {}

writer_attr = csv.writer (import_settings.attributes)
writer_attr.writerow (["ID","Country","State"])

for id, seq in variants["sequences"].items():
    accession = id.split('_')[2]
    try:
        date_raw = annotation["epi_isl_%s" % accession]['collected']
        date_check = datetime.datetime.strptime (date_raw,date_parse_format)
        date_raw = date_check.strftime ("%m%d%Y")
        if annotation["epi_isl_%s" % accession]['location']['country'] == "USA":
            if date_raw and date_raw != '20200101':
                id = "%s|%s|%s" % (accession, date_raw,annotation["epi_isl_%s" % accession]['location']['state'])
                writer_attr.writerow ([accession, annotation["epi_isl_%s" % accession]['location']['country'], annotation["epi_isl_%s" % accession]['location']['state'] ])
                variants_as_sets [id] = set (["%d%s" % (v[0],v[1]) for v in seq])
    except:
        pass

print ("Filtered %d variants" % len (variants_as_sets), file = sys.stderr)

writer = csv.writer (import_settings.output)
writer.writerow (["ID1","ID2","Distance"])


for s1, s2 in itertools.combinations(variants_as_sets.items(),2):
    v1 = s1[1]
    v2 = s2[1]
    d = len (v1^v2) / N
    
    if d < import_settings.cutoff:
        writer.writerow ([s1[0],s2[0],str (d)])
                
       