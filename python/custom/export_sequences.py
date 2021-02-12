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

    fields = [record['id'], location(record), value_or_null(record['originalCollected']), value_or_null(record['technology'])]
    return unicodedata.normalize('NFKD', "/".join (fields))

def export_sequences(config):

    db = MongoClient()

    acceptable = ['collected', 'originalCollected', 'host', 'id', 'location', 'name', 'technology', 'type', 'nextstrainClade', 'pangolinLineage', 'gisaidClade', 'seq']
    HOST= "Human"
    MINLENGTH=28000

    mongo_query = { "host" : HOST,  "length": {"$gt": MINLENGTH }, "seq": {"$exists":True} }

    output_fn = config["sequence-output"]

    if("clade-type" in config.keys()):
        clade_type = config["clade-type"]
    else:
        clade_type = "pangolinLineage"

    if("clades" in config.keys()):
        # db.inventory.find ( { quantity: { $in: [20, 50] } } )
        mongo_query[clade_type] = { "$in": config["clades"] }

    # Query for human host and sequence length greater than 28000, and sequence populated

    # LIMIT=100000
    # LIMIT=10
    records = list(db.gisaid.records.find(mongo_query))
    seq_records = [SeqRecord(Seq(rec["seq"]),id=sequence_name(rec),name='',description='') for rec in records]

    # Write to fasta
    with open(output_fn, 'w', encoding='utf-8') as output_fh:
        SeqIO.write(seq_records, output_fh, "fasta")

if __name__ == "__main__":

    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-o', '--output',   help = 'fasta output', type = str)
    args = arguments.parse_args()
    config = {"sequence-output" : args.output }
    config["clades"] = ["B.1.351", "P.1"]
    config["clade-type"] = "pangolinLineage"
    export_sequences(config)

