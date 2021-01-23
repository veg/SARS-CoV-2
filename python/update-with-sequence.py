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

import pymongo
from pymongo import MongoClient

arguments = argparse.ArgumentParser(description='Report which dates have full report')
arguments.add_argument('-i', '--input',   help = 'fasta to update', required = True, type = argparse.FileType('r'))

args = arguments.parse_args()
db = MongoClient()

def update_record(seq):
    seq_str = str(seq.seq)
    try:
        epi_id = seq.description.split('|')[1].lower()
    except:
        print("could not process " + seq.description)
        return
    print(epi_id)
    db.gisaid.records.update_one({'id': epi_id}, {'$set': {'seq': seq_str}})

seqs = list(SeqIO.parse(args.input, 'fasta'))

cnt = 0
for seq in seqs:
    cnt = cnt + 1
    update_record(seq)

