import csv
import json
import sys
import argparse
import copy
import itertools
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

arguments = argparse.ArgumentParser(description='Report which dates have full report')
arguments.add_argument('-i', '--input',   help = 'fasta to filter duplicates', required = True, type = argparse.FileType('r'))
arguments.add_argument('-m', '--nuc-input',   help = 'nucleotide fasta to filter duplicates', required = True, type = argparse.FileType('r'))
arguments.add_argument('-d', '--duplicates',   help = 'duplicates file', type = argparse.FileType('w'))
arguments.add_argument('-o', '--output', help = 'write compressed fasta here', type = argparse.FileType('w'), default = sys.stdout)
arguments.add_argument('-n', '--nuc-output', help = 'write nucleotide compressed fasta here', required = True, type = argparse.FileType('w'))

args = arguments.parse_args()

seqs = list(SeqIO.parse(args.input, 'fasta'))
seqs = sorted(seqs, key=lambda x: x.seq)

dupes = {}
for rec in seqs:
    if(rec.seq not in dupes):
        dupes[rec.seq] = [rec]
    else:
        dupes[rec.seq].append(rec)

# Write new fasta from first index of each dupe, and duplicates file if provided
vals = list(dupes.values())

# Write fasta
firsts = [val[0] for val in vals]

# Write dupes
dupe_names = {val[0].name : {"{0}".format(i): val[i].name for i in range(len(val))} for val in vals}


for first in firsts:
    first.id = first.name + "_" + str(len(dupe_names[first.name].keys()))
    first.description = first.name + "_" + str(len(dupe_names[first.name].keys()))

# Second pass to find nearly similar
# Get sequences that are a difference of one
SeqIO.write(firsts, args.output, "fasta")
filtered_seq_names = [seq.name for seq in firsts]

# Compress nucleotide file as well
nuc_seqs = list(SeqIO.parse(args.nuc_input, 'fasta'))
filtered_nuc_seqs = [seq for seq in nuc_seqs if seq.name in filtered_seq_names]

for filtered_nuc_seq in filtered_nuc_seqs:
    filtered_nuc_seq.id = filtered_nuc_seq.name + "_" + str(len(dupe_names[filtered_nuc_seq.name].keys()))
    filtered_nuc_seq.description = filtered_nuc_seq.name + "_" + str(len(dupe_names[filtered_nuc_seq.name].keys()))

SeqIO.write(filtered_nuc_seqs, args.nuc_output, "fasta")

json.dump(dupe_names, args.duplicates, indent=4, sort_keys=True)

