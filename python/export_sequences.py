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
from datetime import date, timedelta, datetime
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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

    # if('originalCollected' in record.keys()):
    #     fields = [record['id'], location(record), value_or_null(record['originalCollected']), value_or_null(record['technology'])]
    # else:
    #     fields = [record['id'], location(record), value_or_null(record['collected']), value_or_null(record['technology'])]

    # Just export ID
    return unicodedata.normalize('NFKD', record['id'])

def export_sequences(config):

    db = MongoClient(host='129.32.209.134')

    acceptables = ['collected', 'originalCollected', 'host', 'id', 'location', 'name', 'technology', 'type', 'nextstrainClade', 'pangolinLineage', 'gisaidClade', 'seq']
    HOST= "Human"
    MINLENGTH=28000

    mongo_query = { "host" : HOST,  "length": {"$gt": MINLENGTH }, "seq": {"$exists":True} }

    output_fn = config["sequence-output"]

    if("clade-type" in config.keys()):
        clade_type = config["clade-type"]
    else:
        clade_type = "pangolinLineage"

    if("clades" in config.keys()):
        mongo_query[clade_type] = { "$in": config["clades"] }
    elif("ignore-clades" in config.keys()):
        mongo_query[clade_type] = { "$nin": config["ignore-clades"] }

    if("collection-date-range" in config.keys()):
        start_date = config["collection-date-range"][0]
        end_date = config["collection-date-range"][1]
        mongo_query["collected"] = { "$gt": datetime.strptime(start_date, "%Y-%m-%d"), "$lt": datetime.strptime(end_date, "%Y-%m-%d") }
        mongo_query["originalCollected"] = { "$regex": "[0-9]{4}-[0-9]{2}" }

    db_mongo_query = db.gisaid.records.find(mongo_query, acceptables);

    if("get-latest-by-collection-date") in config.keys():
        # Add sort and limit
        limit = config["get-latest-by-collection-date"]
        db_mongo_query = db_mongo_query.sort([("collected", -1)]).limit(limit)

    # Query for human host and sequence length greater than 28000, and sequence populated
    records = list(db_mongo_query)
    seq_records = [SeqRecord(Seq(rec["seq"]),id=sequence_name(rec),name='',description='') for rec in records]

    if len(seq_records) > 0:
        # Write to fasta
        with open(output_fn, 'w', encoding='utf-8') as output_fh:
            SeqIO.write(seq_records, output_fh, "fasta")

def export_premsa_sequences(config, nuc_output_fn, prot_output_fn, gene):
    '''
    config
    nuc_output_fn -- nucleotide output filename
    protein_output_fn -- protein output filename
    gene -- 'region of SARS-CoV-2
    '''

    db = MongoClient(host='129.32.209.134')

    acceptables = ['collected', 'originalCollected', 'host', 'id', 'location', 'name', 'technology', 'type', 'nextstrainClade', 'pangolinLineage', 'gisaidClade']
    HOST= "Human"
    MINLENGTH=28000
    only_uniques = True
    key_to_export = ""

    if("only-uniques" in config.keys()):
        only_uniques = config["only-uniques"]

    # Get QC key
    validation_key = 'qc.' + gene + '.passed'
    not_duplicate_key = 'qc.' + gene + '.duplicate_of'
    not_duplicate_val = {"$exists":False}

    second_duplicate_key = 'duplicate_of_by_gene.' + gene
    second_duplicate_val = [{second_duplicate_key:{"$exists":False}},{second_duplicate_key:'reference'}]

    nuc_key_to_export = gene + '_premsa_nuc_seq'
    prot_key_to_export = gene + '_premsa_protein_seq'

    mongo_query = { "host" : HOST, validation_key: True}
    acceptables.extend([nuc_key_to_export, prot_key_to_export])

    if(only_uniques):
        mongo_query[not_duplicate_key] = not_duplicate_val
        mongo_query["$or"] = second_duplicate_val

    if("clade-type" in config.keys()):
        clade_type = config["clade-type"]
    else:
        clade_type = "pangolinLineage"

    if("clades" in config.keys()):
        # db.inventory.find ( { quantity: { $in: [20, 50] } } )
        mongo_query[clade_type] = { "$in": config["clades"] }

    elif("ignore-clades" in config.keys()):
        mongo_query[clade_type] = { "$nin": config["ignore-clades"] }

    db_mongo_query = db.gisaid.records.find(mongo_query, acceptables)

    if("get-latest-by-collection-date") in config.keys():
        # Add sort and limit
        limit = config["get-latest-by-collection-date"]
        db_mongo_query = db_mongo_query.sort([("collected", -1)]).limit(limit)

    # Query for human host and sequence length greater than 28000, and sequence populated
    records = list(db_mongo_query)

    # Filter sequences down to those that have been processed
    recs_with_nucs = filter(lambda x: nuc_key_to_export in x.keys(), records)
    nuc_seq_records = [SeqRecord(Seq(rec[nuc_key_to_export]),id=sequence_name(rec),name='',description='') for rec in recs_with_nucs]

    recs_with_prot = filter(lambda x: prot_key_to_export in x.keys(), records)
    prot_seq_records = [SeqRecord(Seq(rec[prot_key_to_export]),id=sequence_name(rec),name='',description='') for rec in recs_with_prot]

    if len(nuc_seq_records) <= 0 or len(prot_seq_records) <= 0:
        return False

    # Write to fasta
    with open(nuc_output_fn, 'w', encoding='utf-8') as nuc_output_fh:
        SeqIO.write(nuc_seq_records, nuc_output_fh, "fasta")

    # Write to fasta
    with open(prot_output_fn, 'w', encoding='utf-8') as prot_output_fh:
        SeqIO.write(prot_seq_records, prot_output_fh, "fasta")

    return True

def export_bealign_sequences(config, nuc_output_fn, gene):
    '''
    config -- TODO: Outline what it supports.
    nuc_output_fn -- nucleotide output filename
    gene -- region of SARS-CoV-2
    '''

    db = MongoClient(host='129.32.209.134')

    acceptables = ['collected', 'originalCollected', 'host', 'id', 'location', 'name', 'technology', 'type', 'nextstrainClade', 'pangolinLineage', 'gisaidClade']
    HOST= "Human"
    MINLENGTH=28000
    validation_key = 'qc.' + gene + '.passed'

    # Get QC key
    bealign_key_to_export = 'bealign.' + gene

    mongo_query = { "host" : HOST, validation_key: True }
    acceptables.extend([bealign_key_to_export])

    if("clade-type" in config.keys()):
        clade_type = config["clade-type"]
    else:
        clade_type = "pangolinLineage"

    if("clades" in config.keys()):
        mongo_query[clade_type] = { "$in": config["clades"] }

    elif("ignore-clades" in config.keys()):
        mongo_query[clade_type] = { "$nin": config["ignore-clades"] }

    if("collection-date-range" in config.keys()):
        start_date = config["collection-date-range"][0]
        end_date = config["collection-date-range"][1]
        mongo_query["collected"] = { "$gt": datetime.strptime(start_date, "%Y-%m-%d"), "$lt": datetime.strptime(end_date, "%Y-%m-%d") }
        # mongo_query["originalCollected"] = { "$regex": "[0-9]{4}-[0-9]{2}" }

    db_mongo_query = db.gisaid.records.find(mongo_query, acceptables)

    if("get-latest-by-collection-date") in config.keys():
        # Add sort and limit
        limit = config["get-latest-by-collection-date"]
        db_mongo_query = db_mongo_query.sort([("collected", -1)]).limit(limit)

    # if(True):
    #     db_mongo_query.limit(5000)

    # Query for human host and sequence length greater than 28000, and sequence populated
    records = list(db_mongo_query)


    # Filter records that have bealign
    # Filter sequences down to those that have been processed
    nuc_seq_records = []
    for rec in records:
        if 'bealign' in rec.keys():
            if gene in rec['bealign'].keys():
                nuc_seq_records.append(SeqRecord(Seq(rec['bealign'][gene]),id=sequence_name(rec),name='',description=''))

    if len(nuc_seq_records) > 0:
        # Write to fasta
        with open(nuc_output_fn, 'w', encoding='utf-8') as nuc_output_fh:
            SeqIO.write(nuc_seq_records, nuc_output_fh, "fasta")
        return True

    return False

if __name__ == "__main__":

    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-o', '--output',   help = 'fasta output', type = str)
    args = arguments.parse_args()

    config = {}
    config["sequence-output"] = args.output
    # config['get-latest-by-collection-date'] = 100000
    config['only-uniques'] = False
    config["collection-date-range"] = ("2020-03-01", "2020-04-01")
    # config["clades"] = ["B.1.351"]
    # config["clades"] = ["B.1.427", "B.1.429"]
    # config["clades"] = ["B.1.1.7"]
    # config["ignore-clades"] = ["B.1.351", "P.1", "B.1.1.7"]
    # config["clade-type"] = "pangolinLineage"

    # export_sequences(config)
    # export_premsa_sequences(config, 'nuc.fas', 'prot.fas', '3C')
    # export_postmsa_sequences(config, 'aligned.fas', '3C')
    export_bealign_sequences(config, args.output, 'leader')


