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
arguments.add_argument('-i', '--input', help = 'fasta to filter duplicates', required = True, type = argparse.FileType('r'))
arguments.add_argument('-m', '--nuc-input', help='nucleotide fasta to filter duplicates', required = True, type = argparse.FileType('r'))
arguments.add_argument('-d', '--duplicates', help='duplicates file', type = argparse.FileType('w'))
arguments.add_argument('-g', '--nucleotide-duplicates', help='nucleotide-duplicates file', type = argparse.FileType('w'))
arguments.add_argument('-o', '--output', help = 'write compressed fasta here', type = argparse.FileType('w'), default = sys.stdout)
arguments.add_argument('-n', '--nuc-output', help = 'write nucleotide compressed fasta here', required = True, type = argparse.FileType('w'))

args = arguments.parse_args()


def get_dupes(seqs):

    seqs = sorted(seqs, key=lambda x: x.seq)

    dupes = {}
    for rec in seqs:
        if(rec.seq not in dupes):
            dupes[rec.seq] = [rec]
        else:
            dupes[rec.seq].append(rec)

    return dupes

p_seqs = list(SeqIO.parse(args.input, 'fasta'))
dupes = get_dupes(p_seqs)

nuc_seqs = list(SeqIO.parse(args.nuc_input, 'fasta'))
nuc_dupes = get_dupes(nuc_seqs)

# Write new fasta from first index of each dupe, and duplicates file if provided
vals = list(dupes.values())

# Write fasta
firsts = [val[0] for val in vals]

# Write dupes
dupe_names = {val[0].name : {"{0}".format(i): val[i].name for i in range(len(val))} for val in vals}
json.dump(dupe_names, args.duplicates, indent=4, sort_keys=True)

# Write nuc_dups
nuc_vals = list(nuc_dupes.values())
nuc_dupe_names = {val[0].name : {"{0}".format(i): val[i].name for i in range(len(val))} for val in nuc_vals}
json.dump(nuc_dupe_names, args.nucleotide_duplicates, indent=4, sort_keys=True)

for first in firsts:
    first.id = first.name
    first.description = first.name
    # first.id = first.name + "_" + str(len(dupe_names[first.name].keys()))
    # first.description = first.name + "_" + str(len(dupe_names[first.name].keys()))

# Second pass to find nearly similar
# Get sequences that are a difference of one
SeqIO.write(firsts, args.output, "fasta")

nuc_firsts = [val[0] for val in nuc_vals]
filtered_seq_names = [seq.name for seq in nuc_firsts]

for first in nuc_firsts:
    first.id = first.name
    first.description = first.name
    # first.id = first.name + "_" + str(len(nuc_dupe_names[first.name].keys()))
    # first.description = first.name + "_" + str(len(nuc_dupe_names[first.name].keys()))

SeqIO.write(nuc_firsts, args.nuc_output, "fasta")
