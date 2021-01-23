import csv
import json
import sys
import argparse
import itertools
import shutil
import copy
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

arguments = argparse.ArgumentParser(description='Report which dates have full report')
arguments.add_argument('-f', '--fasta-file',   help = 'fasta to overwrite', required = True, type = argparse.FileType('r'))
arguments.add_argument('-m', '--map-file',   help = 'fasta to filter duplicates', required = True, type = argparse.FileType('r'))
args = arguments.parse_args()

# If one fails, then copy the other to the output. If both fail, then throw an error
map_json = json.load(args.map_file)
seqs = list(SeqIO.parse(args.fasta_file, 'fasta'))

# Overwrite files
orig_fn = args.fasta_file.name
tmp_fn = orig_fn + '.tmp'

# Fix FASTA headers
for seq in seqs:
    old_id = seq.id
    seq.id = map_json[old_id]
    seq.description = map_json[old_id]

args.fasta_file.close()

with open(tmp_fn, 'w') as tmp_fp:
    SeqIO.write(seqs, tmp_fp, "fasta")
    shutil.move(tmp_fn, orig_fn)


