import csv
import json
import sys
import argparse
import itertools
import shutil
import copy
import os
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO

from pymongo import MongoClient

arguments = argparse.ArgumentParser(description='Report which dates have full report')
arguments.add_argument('-i', '--input',   help = 'fasta to update', required = True, type = argparse.FileType('r'))

args = arguments.parse_args()
db = MongoClient(host='192.168.0.4')

def get_missing_seqs():
    return [d['name'] for d in db.gisaid.records.find({'seq': None})]

# Add 100,000 limit
seqs = list(SeqIO.parse(args.input, 'fasta'))
missing_names = get_missing_seqs()
print('got ' + str(len(missing_names)) + ' missing seqs in the database')

def update_record(seq):
    seq_str = str(seq.seq)
    try:
        name = '/'.join(seq.description.split('/')[1:3]) + '/' + ''.join(list(seq.description.split('/')[3])[:4])
    except:
        print("could not process " + seq.description)
        return
    # Get description
    if(name in missing_names):
        result = db.gisaid.records.update_one({'name': name}, {'$set': {'seq': seq_str}})
        print('updated ' + name)
        print(result.raw_result)

cnt = 0
for seq in seqs:
    cnt = cnt + 1
    update_record(seq)

