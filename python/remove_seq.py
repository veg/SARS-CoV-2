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
arguments.add_argument('-i', '--input',   help = 'fasta to remove sequences from', required = True, type = argparse.FileType('r'))
arguments.add_argument('-r', '--reference',   help = 'msa of sequences to remove', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output', help = 'output with removed sequences', type = argparse.FileType('w'), default = sys.stdout)

args = arguments.parse_args()

seqs = list(SeqIO.parse(args.input, 'fasta'))
to_remove = list(SeqIO.parse(args.reference, 'fasta'))
to_remove_names = [x.name for x in to_remove]

to_write = list(filter(lambda x: x.name not in to_remove_names, seqs))

# Second pass to find nearly similar
SeqIO.write(to_write, args.output, "fasta")
