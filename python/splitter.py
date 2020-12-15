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

def write_seq(seq):
    # First index is the directory
    # The second is the id
    split = seq.description.split('|')
    fn = '/'.join(split[:2]).lower().replace(' ', '_') + '.fas'
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    with open(fn, "w") as f:
        SeqIO.write(seq, f, "fasta")



fn = './2020-12-03.fasta'
seqs = list(SeqIO.parse(fn, 'fasta'))

# Write each file to directories based on headers
cpus = multiprocessing.cpu_count()
with Pool(cpus) as p:
    p.map(write_seq, seqs)

