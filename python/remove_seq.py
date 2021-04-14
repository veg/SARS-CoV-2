import csv
import json
import sys
import argparse
import copy
import itertools
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

def remove_reference_seq(input_fn, reference_fn, output_fn):

    input_fh = open(input_fn, 'r')
    reference_fh = open(reference_fn, 'r')

    seqs = list(SeqIO.parse(input_fn, 'fasta'))
    to_remove = list(SeqIO.parse(reference_fn, 'fasta'))
    to_remove_names = [x.name for x in to_remove]

    to_write = list(filter(lambda x: x.name not in to_remove_names, seqs))

    # Second pass to find nearly similar
    SeqIO.write(to_write, output_fn, "fasta")

def reserve_only_original_input(input_fn, original_fn, output_fn):

    input_fh = open(input_fn, 'r')
    original_fh = open(original_fn, 'r')

    seqs = list(SeqIO.parse(input_fn, 'fasta'))

    to_keep = list(SeqIO.parse(original_fn, 'fasta'))
    to_keep_names = [x.name for x in to_keep]

    to_write = []
    uniq_names = []
    for seq in seqs:
        if seq.name in to_keep_names and seq.name not in uniq_names:
            to_write.append(seq)
            uniq_names.append(seq.name)

    # Second pass to find nearly similar
    SeqIO.write(to_write, output_fn, "fasta")


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-i', '--input',   help = 'fasta to remove sequences from', required = True, type = str)
    arguments.add_argument('-r', '--reference',   help = 'msa of sequences to remove', required = True, type = str)
    arguments.add_argument('-o', '--output', help = 'output with removed sequences', type = str, default = sys.stdout)
    args = arguments.parse_args()
    remove_reference_seq(args.input, args.reference, args.output)
    #reserve_only_original_input(args.input, args.reference, args.output)

