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


def update_record(seq, gene):

    seq_str = str(seq.seq)

    try:
        epi_id = '_'.join(seq.description.split('_')[:3]).lower()
    except:
        print("could not process " + seq.description)
        return

    # print(epi_id)
    key_to_update = ".".join(["reference_alignment", gene])
    return UpdateOne({'id': epi_id}, {'$set': {key_to_update: seq_str}})

def store_postmsa_file(input, gene):

    db = MongoClient(host='192.168.0.4')

    input_fh = open(input, 'r')
    seqs = list(SeqIO.parse(input_fh, 'fasta'))
    to_update = [update_record(seq, gene) for seq in seqs]

    if(len(to_update)):
        results = db.gisaid.records.bulk_write(to_update)
    print("Updated " +  str(results.modified_count) + " of " + str(len(to_update)) + " items to update")


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-i', '--input',   help = 'fasta to update', required = True, type = str)
    arguments.add_argument('-g', '--gene',   help = 'gene region', required = True, type = str)
    args = arguments.parse_args()
    store_postmsa_file(args.input, args.gene)

