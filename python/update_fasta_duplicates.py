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

def update_fasta_duplicates(fasta_file, map_file):
    # If one fails, then copy the other to the output. If both fail, then throw an error
    with open(map_file, 'r') as map_file_fh:
        map_json = json.load(map_file_fh)

    seqs = list(SeqIO.parse(fasta_file, 'fasta'))

    # Overwrite files
    orig_fn = fasta_file
    tmp_fn = orig_fn + '.tmp'

    # Fix FASTA headers
    # Create new map that does not have copy numbers
    just_id_map = {'_'.join(k.split('_')[:3]): v for k,v in map_json.items()}

    for seq in seqs:
        old_id = '_'.join(seq.id.split('_')[:3])
        seq.id = just_id_map[old_id]
        seq.description = just_id_map[old_id]

    with open(tmp_fn, 'w') as tmp_fp:
        SeqIO.write(seqs, tmp_fp, "fasta")
        shutil.move(tmp_fn, orig_fn)

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-f', '--fasta-file',   help = 'fasta to overwrite', required = True, type = str)
    arguments.add_argument('-m', '--map-file',   help = 'fasta to filter duplicates', required = True, type = str)
    args = arguments.parse_args()
    update_fasta_duplicates(args.fasta_file, args.map_file)

