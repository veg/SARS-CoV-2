# -*- coding: utf-8 -*-

import csv
import json
import sys
import argparse
import itertools
import shutil
import copy
import os
import multiprocessing
import unicodedata
from multiprocessing import Pool
from datetime import date, timedelta
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pymongo
from pymongo import MongoClient


HOST= "Human"
MINLENGTH=28000

arguments = argparse.ArgumentParser(description='Report which dates have full report')
arguments.add_argument('-o', '--output',   help = 'fasta output', type = argparse.FileType('w', encoding='utf-8'))

args = arguments.parse_args()
db = MongoClient()

def sequence_name(record):
    #hCoV-19/England/20102000906/2020
    def value_or_null (v):
        if v is not None:
            return v.split (" ")[0]
        return "null"
    def location (v):
        for k in ['state','country','subregion']:
            if v['location'][k]:
                return v['location'][k].replace (' ', '_')

    fields = [record['id'],location(record),value_or_null (record['collected']),value_or_null (record['technology'])]
    return unicodedata.normalize('NFKD', "/".join (fields))

# Query for human host and sequence length greater than 28000, and sequence populated

# LIMIT=100000
# LIMIT=10
records = list(db.gisaid.records.find({"host":HOST, "length": {"$gt": MINLENGTH }, "seq":{"$exists":True}}))

seq_records = [SeqRecord(Seq(rec["seq"]),id=sequence_name(rec),name='',description='') for rec in records]

# Write to fasta
SeqIO.write(seq_records, args.output, "fasta")
