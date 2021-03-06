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

from pymongo import MongoClient, InsertOne, DeleteMany, ReplaceOne, UpdateOne, UpdateMany

db = MongoClient(host='192.168.0.4')

def update_dupe_records(doc, gene):
    dupe_key = ".".join(["duplicate_of_by_gene", gene])
    old_key = 'qc.' + gene + '.duplicate_of'
    reference = doc['qc'][gene]['duplicate_of']
    dupe_val = { dupe_key : reference }
    id = doc['id']
    return UpdateOne({'id': id}, {'$set': dupe_val,  "$unset" : { old_key : 1 } })

def get_items_with_old_key(gene):
    old_key = 'qc.' + gene + '.duplicate_of'
    return list(db.gisaid.records.find({old_key : { "$exists" : True }}))

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Mark duplicates in MongoDB')
    arguments.add_argument('-g', '--gene',   help = 'gene region', required = True, type = str)
    args = arguments.parse_args()
    docs_to_update = get_items_with_old_key(args.gene)
    items = [update_dupe_records(doc, args.gene) for doc in docs_to_update]
    results = db.gisaid.records.bulk_write(items)
    print("Updated " +  str(results.modified_count) + " of " + str(len(items)) + " items to update")
