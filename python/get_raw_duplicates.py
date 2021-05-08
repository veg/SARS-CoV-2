import csv
import json
import sys
import argparse
import copy
import itertools
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

def get_dupes(seqs):

    seqs = sorted(seqs, key=lambda x: x.seq)

    dupes = {}
    for rec in seqs:
        if(rec.seq not in dupes):
            dupes[rec.seq] = [rec]
        else:
            dupes[rec.seq].append(rec)

    return dupes

def get_raw_duplicates(input, nuc_input):

    p_seqs = list(SeqIO.parse(input, 'fasta'))
    dupes = get_dupes(p_seqs)

    nuc_seqs = list(SeqIO.parse(nuc_input, 'fasta'))
    nuc_dupes = get_dupes(nuc_seqs)

    # Write new fasta from first index of each dupe, and duplicates file if provided
    vals = list(dupes.values())

    # Write fasta
    firsts = [val[0] for val in vals]

    # Write dupes
    dupe_names = {val[0].name : {"{0}".format(i): val[i].name for i in range(len(val))} for val in vals}

    # Write nuc_dups
    nuc_vals = list(nuc_dupes.values())
    nuc_dupe_names = {val[0].name : {"{0}".format(i): val[i].name for i in range(len(val))} for val in nuc_vals}

    for first in firsts:
        first.id = first.name
        first.description = first.name
        # first.id = first.name + "_" + str(len(dupe_names[first.name].keys()))
        # first.description = first.name + "_" + str(len(dupe_names[first.name].keys()))

    nuc_firsts = [val[0] for val in nuc_vals]
    filtered_seq_names = [seq.name for seq in nuc_firsts]

    for first in nuc_firsts:
        first.id = first.name
        first.description = first.name
        # first.id = first.name + "_" + str(len(nuc_dupe_names[first.name].keys()))
        # first.description = first.name + "_" + str(len(nuc_dupe_names[first.name].keys()))


    # Protein seq map
    prot_seq_map = {record.id : record for record in p_seqs}
    prot_nuc_firsts = [prot_seq_map[nuc_first.id] for nuc_first in nuc_firsts if nuc_first.id in prot_seq_map.keys()]

    return (dupe_names, nuc_dupe_names, prot_nuc_firsts, nuc_firsts)

def get_nuc_raw_duplicates(input):

    nuc_seqs = list(SeqIO.parse(input, 'fasta'))
    dupes = get_dupes(nuc_seqs)

    # Write new fasta from first index of each dupe, and duplicates file if provided
    vals = list(dupes.values())

    # Write fasta
    firsts = [val[0] for val in vals]

    # Write dupes
    dupe_names = {val[0].name : {"{0}".format(i): val[i].name for i in range(len(val))} for val in vals}

    for first in firsts:
        first.id = first.name
        first.description = first.name

    return (dupe_names, firsts)


def write_raw_duplicates(input, nuc_input, duplicates, nucleotide_duplicates):
    input_fh = open(input, 'r')
    nuc_input_fh = open(nuc_input, 'r')

    dupe_names, nuc_dupe_names, prot_nuc_firsts, nuc_firsts = get_raw_duplicates(input_fh, nuc_input_fh)

    duplicates_fh = open(duplicates, 'w')
    nucleotide_duplicates_fh = open(nucleotide_duplicates, 'w')
    json.dump(dupe_names, duplicates_fh, indent=4, sort_keys=True)
    json.dump(nuc_dupe_names, nucleotide_duplicates_fh, indent=4, sort_keys=True)

def write_nuc_raw_duplicates(input, duplicate_output, uniques_output):
    input_fh = open(input, 'r')
    dupe_names, nuc_firsts = get_nuc_raw_duplicates(input_fh)
    duplicates_fh = open(duplicate_output, 'w')
    json.dump(dupe_names, duplicates_fh, indent=4, sort_keys=True)
    SeqIO.write(nuc_firsts, uniques_output, "fasta")

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-i', '--input', help = 'fasta to filter duplicates', required = True, type = argparse.FileType('r'))
    arguments.add_argument('-m', '--nuc-input', help='nucleotide fasta to filter duplicates', required = True, type = argparse.FileType('r'))
    arguments.add_argument('-d', '--duplicates', help='duplicates file', type = argparse.FileType('w'))
    arguments.add_argument('-g', '--nucleotide-duplicates', help='nucleotide-duplicates file', type = argparse.FileType('w'))
    arguments.add_argument('-o', '--output', help = 'write compressed fasta here', type = argparse.FileType('w'), default = sys.stdout)
    arguments.add_argument('-n', '--nuc-output', help = 'write nucleotide compressed fasta here', required = True, type = argparse.FileType('w'))
    args = arguments.parse_args()

    dupe_names, nuc_dupe_names, prot_nuc_firsts, nuc_firsts = get_raw_duplicates(args.input, args.nuc_input)

    json.dump(dupe_names, args.duplicates, indent=4, sort_keys=True)
    json.dump(nuc_dupe_names, args.nucleotide_duplicates, indent=4, sort_keys=True)

    # Write protein sequences based on nucleotide dupes
    # Get p_seqs based on nucleotide dupes
    SeqIO.write(prot_nuc_firsts, args.output, "fasta")
    SeqIO.write(nuc_firsts, args.nuc_output, "fasta")



