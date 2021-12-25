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
from pymongo import InsertOne, DeleteMany, ReplaceOne, UpdateOne

def parse_premsa_log(log_file, gene):
    # Open file and parse
    with open(log_file, 'r') as log_file_fh:
        lines = log_file_fh.readlines()

    # Get all lines that start with warning
    warnings = [l.split('Sequence')[1].strip() for l in lines if "WARNING" in l]

    # Parse sequence ids and the message
    # Return list of errors
    return [("_".join(warning.split()[0].split("_")[:3]), warning) for warning in warnings]

def update_mongo(items, gene):
    '''
    For each item of the form
    ('epi_isl_860592_KwaZulu_Natal_2021_01_02_null', 'epi_isl_860592_KwaZulu_Natal_2021_01_02_null has too many ambiguous nucleotides; try setting --N-fraction flag to a higher value')
    '''

    db = MongoClient(host='129.32.209.134')

    bulk_updates = []
    for item in items:
        epi_id = { 'id' : item[0] }
        key_to_update = ".".join(["qc", gene])

        setter = { "$set": {
                                key_to_update : { "passed": False, "msg": item[1] }
                           }
                    }

        update_one = UpdateOne(epi_id, setter)
        bulk_updates.append(update_one)

    results = db.gisaid.records.bulk_write(bulk_updates)

    return results

def mark_troubled(log_file, gene):
    items = parse_premsa_log(log_file, gene)
    results = update_mongo(items, gene)
    print("Updated " +  str(results.modified_count) + " of " + str(len(items)) + " items found in the log")

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-l', '--log-file', help = 'premsa log file', type = str)
    arguments.add_argument('-g', '--gene', help = 'gene', type = str)
    args = arguments.parse_args()
    mark_troubled(args.log_file, args.gene)

