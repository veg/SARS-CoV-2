import csv
import json
import sys
import argparse
import itertools
import shutil
import copy
import os
import multiprocessing
from multiprocessing import Pool
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

fn = './data/to-import/sequences.fasta'
chunk_by = 50000
seqs = list(SeqIO.parse(fn, 'fasta'))
chunked_seqs = grouper(seqs, chunk_by)
index = 1

def write_seqs(seqs):
    # First index is the directory
    # The second is the id
    split = seq.description.split('|')
    sub_fn = './data/to-import/sequences.' + str(index) + '-' + str(index + chunk_by) + '.fasta'
    index = index + chunk_by + 1
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn, "w") as f:
        SeqIO.write(seq, f, "fasta")


# Write each file to directories based on headers
cpus = multiprocessing.cpu_count()
with Pool(cpus) as p:
    p.map(write_seqs, chunked_seqs)

