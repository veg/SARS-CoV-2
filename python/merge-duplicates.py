import csv
import json
import sys
import argparse
import itertools
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

arguments = argparse.ArgumentParser(description='Report which dates have full report')
arguments.add_argument('-p', '--protein-duplicates',   help = 'fasta to filter duplicates', required = True, type = argparse.FileType('r'))
arguments.add_argument('-n', '--nuc-duplicates',   help = 'duplicates file', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output', help = 'write compressed fasta here', type = argparse.FileType('w'), default = sys.stdout)
args = arguments.parse_args()

# If one fails, then copy the other to the output. If both fail, then throw an error
def merge(protein_json, nuc_json):
    return protein_json

try:
    protein_json = json.load(args.protein_duplicates)
except:
    protein_json = False

try:
    nuc_json = json.load(args.nuc_duplicates)
except:
    nuc_json = False

if not protein_json and not nuc_json:
    output_json = {}
elif not protein_json:
    output_json = nuc_json
elif not protein_json:
    output_json = nuc_json
else:
    output_json = merge(protein_json, nuc_json)

json.dump(output_json, args.output, indent=4, sort_keys=True)

