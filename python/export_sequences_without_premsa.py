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
from dateutil.relativedelta import *
from multiprocessing import Pool
import datetime
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

def export_sequences(gene, output_fn):

    db = MongoClient(host='129.32.209.134')

    get_nuc_key = lambda gene: gene + '_premsa_nuc_seq'
    get_prot_key = lambda gene: gene + '_premsa_protein_seq'
    get_qc_key = lambda gene: 'qc.' + gene + '.passed'

    nuc_key = get_nuc_key(gene)
    prot_key = get_prot_key(gene)
    qc_key = get_qc_key(gene)

    acceptable = ['collected', 'originalCollected', 'host', 'id', 'location', 'name', 'technology', 'type', 'nextstrainClade', 'pangolinLineage', 'gisaidClade', 'seq']
    HOST= "Human"
    MINLENGTH=28000
    TODAY = datetime.datetime.today()
    LAST_MONTH = TODAY-relativedelta(months=+1, day=1)


    mongo_query = { "host" : HOST,  "length": {"$gt": MINLENGTH }, "seq": {"$exists":True} }
    mongo_query[nuc_key] = { "$exists": False }
    mongo_query[prot_key] = { "$exists": False }
    mongo_query[qc_key] = { "$exists": False }

    # Speed up query by only looking at submitted from last month
    mongo_query["submitted"] = { "$gte": LAST_MONTH }

    # Query for human host and sequence length greater than 28000, and sequence populated

    # LIMIT=100000
    # LIMIT=10
    records = list(db.gisaid.records.find(mongo_query, limit=10000).sort([( '$natural', -1 )]))
    seq_records = [SeqRecord(Seq(rec["seq"]),id=sequence_name(rec),name='',description='') for rec in records]

    # Write to fasta
    with open(output_fn, 'w', encoding='utf-8') as output_fh:
        SeqIO.write(seq_records, output_fh, "fasta")

def export_sequences_without_reference(gene, output_fn, nuc_output_fn):

    db = MongoClient(host='129.32.209.134')

    get_nuc_key = lambda gene: gene + '_premsa_nuc_seq'
    get_prot_key = lambda gene: gene + '_premsa_protein_seq'
    get_qc_key = lambda gene: 'qc.' + gene + '.passed'

    nuc_key = get_nuc_key(gene)
    prot_key = get_prot_key(gene)
    qc_key = get_qc_key(gene)
    duplicate_key = 'duplicate_of_by_gene.' + gene

    key_to_export = ".".join(["reference_alignment", gene])
    acceptable = ['collected', 'originalCollected', 'host', 'id', 'location', 'name', 'technology', 'type', 'nextstrainClade', 'pangolinLineage', 'gisaidClade', 'seq']

    HOST= "Human"
    MINLENGTH=28000

    mongo_query = { "host" : HOST,  "length": {"$gt": MINLENGTH }, "seq": {"$exists":True} }

    mongo_query[nuc_key] = { "$exists": True }
    mongo_query[prot_key] = { "$exists": True }
    mongo_query[qc_key] = { "$exists": True }
    mongo_query[key_to_export] = { "$exists": False }
    mongo_query[duplicate_key] = { "$exists": False }

    # Not a duplicate

    # Query for human host and sequence length greater than 28000, and sequence populated
    records = list(db.gisaid.records.find(mongo_query, limit=99000))

    # Need to write prot_key
    seq_records = [SeqRecord(Seq(rec[prot_key]),id=sequence_name(rec),name='',description='') for rec in records]
    nuc_seq_records = [SeqRecord(Seq(rec[nuc_key]),id=sequence_name(rec),name='',description='') for rec in records]

    # Write to fasta
    with open(output_fn, 'w', encoding='utf-8') as output_fh:
        SeqIO.write(seq_records, output_fh, "fasta")

    with open(nuc_output_fn, 'w', encoding='utf-8') as nuc_output_fh:
        SeqIO.write(nuc_seq_records, nuc_output_fh, "fasta")


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-o', '--output',   help = 'fasta output', type = str)
    arguments.add_argument('-g', '--gene',   help = 'gene output', type = str)
    args = arguments.parse_args()
    export_sequences(args.gene, args.output)

