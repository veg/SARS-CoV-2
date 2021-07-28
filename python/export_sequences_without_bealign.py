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

def sequence_name(record):
    #hCoV-19/England/20102000906/2020
    def value_or_null (v):
        if v is not None:
            if(isinstance(v, date)):
                return v.strftime('%Y-%m-%d')
            else:
                return v.split (" ")[0]
        return "null"
    def location (v):
        for k in ['state','country','subregion']:
            if v['location'][k]:
                return v['location'][k].replace (' ', '_')

    if('originalCollected' in record.keys()):
        fields = [record['id'], location(record), value_or_null(record['originalCollected']), value_or_null(record['technology'])]
    else:
        fields = [record['id'], location(record), value_or_null(record['collected']), value_or_null(record['technology'])]
    return unicodedata.normalize('NFKD', "/".join (fields))

def export_sequences_without_bealign(gene, output_fn):

    db = MongoClient(host='192.168.0.4')

    get_bealign_key = lambda gene: '.'.join(['bealign', gene])
    get_qc_passed_key = lambda gene: '.'.join(['qc', gene, 'passed'])

    bealign_key = get_bealign_key(gene)
    qc_passed_key = get_qc_passed_key(gene)

    acceptable = ['collected', 'originalCollected', 'host', 'id', 'location', 'name', 'technology', 'type', 'nextstrainClade', 'pangolinLineage', 'gisaidClade', 'seq']

    HOST= "Human"
    MINLENGTH=28000
    mongo_query = { "host" : HOST,  "length": {"$gt": MINLENGTH }, "seq": {"$exists":True} }
    mongo_query[bealign_key] = { "$exists": False }

    # Ensure that the gene has passed quality control. We may use pre-msa output at a later date.
    mongo_query[qc_passed_key] = True

    # Query for human host and sequence length greater than 28000, and sequence populated
    records = list(db.gisaid.records.find(mongo_query, limit=75000))

    # Need to write prot_key
    seq_records = [SeqRecord(Seq(rec['seq']),id=sequence_name(rec),name='',description='') for rec in records]

    # Write to fasta
    with open(output_fn, 'w', encoding='utf-8') as output_fh:
        SeqIO.write(seq_records, output_fh, "fasta")

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-o', '--output', help = 'fasta output', type = str)
    arguments.add_argument('-g', '--gene', help = 'gene output', type = str)
    args = arguments.parse_args()
    export_sequences_without_bealign(args.gene, args.output)

