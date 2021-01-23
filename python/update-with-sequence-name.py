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

def get_missing_seqs():
    return [d['name'] for d in db.gisaid.records.find({'seq': None})]

seqs = list(SeqIO.parse(args.input, 'fasta'))
missing_names = get_missing_seqs()

def update_record(seq):
    seq_str = str(seq.seq)
    try:
        name = seq.description
    except:
        print("could not process " + seq.description)
        return
    if(name in missing_names):
        db.gisaid.records.update_one({'name': name}, {'$set': {'seq': seq_str}})

cnt = 0
for seq in seqs:
    cnt = cnt + 1
    update_record(seq)

