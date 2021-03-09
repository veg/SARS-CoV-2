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

def get_duplicates(records, gene):
    references = {}
    duplicates = []
    for rec in records:
        if 'duplicate_of_by_gene' not in rec.keys():
            references[rec['id']] = []
        elif gene not in rec['duplicate_of_by_gene'].keys():
            references[rec['id']] = []
        elif rec['duplicate_of_by_gene'][gene] == 'reference':
            references[rec['id']] = []
        else:
            duplicates.append((rec['id'], rec['duplicate_of_by_gene'][gene]))

    for dupe in duplicates:
        try:
            references[dupe[1]].append(dupe[0])
        except:
            print(dupe[1] + " not in references, which means there was a circular or non-root reference somewhere. Adding it to list of references, but please take care.")
            references[dupe[1]] = [dupe[0]]

    return references

def export_duplicates(output_fn, gene):

    db = MongoClient(host='192.168.0.4')
    acceptable = ['id', 'duplicate_of_by_gene']

    # transform acceptable into mongo query
    acceptable_dict = { k: 1 for k in acceptable}

    mongo_query = {}

    records = list(db.gisaid.records.find(mongo_query, acceptable_dict))
    records = [{k: v for k, v in rec.items() if k in acceptable} for rec in records]

    # split items with duplicate and those that have reference/missing
    duplicates = get_duplicates(records, gene)

    # Write attributes.csv file and MASTER-NO-JSON file
    with open(output_fn, 'w') as output_fh:
        json.dump(duplicates, output_fh, indent=4, sort_keys=True)

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-o', '--output', help = 'json output', type = str)
    arguments.add_argument('-g', '--gene', help = 'gene to output', type = str)
    args = arguments.parse_args()
    export_duplicates(args.output, args.gene)

