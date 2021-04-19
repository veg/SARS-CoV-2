import csv
import json
import sys
import argparse
import itertools
import shutil
import copy
import os
import random
import multiprocessing
from multiprocessing import Pool
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

def sampler(seq_fn, output_fn, num):
    seqs = list(SeqIO.parse(seq_fn, 'fasta'))
    reference_seqs = random.choices(seqs,k=num)

    with open(output_fn, "w") as f:
        SeqIO.write(reference_seqs, f, "fasta")

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-i', '--input',   help = 'fasta to update', required = True, type = str)
    arguments.add_argument('-o', '--output',   help = 'gene region', required = True, type = str)
    arguments.add_argument('-n', '--num',   help = 'number sequences', required = True, type = int)
    args = arguments.parse_args()
    sampler(args.input, args.output, args.num)

