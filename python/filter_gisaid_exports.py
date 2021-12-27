import csv
import json
import sys
import argparse
import itertools
import shutil
import copy
import os
import glob
import multiprocessing
from multiprocessing import Pool
from datetime import date, timedelta, datetime
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser

import pymongo
from pymongo import MongoClient

def get_most_recent_submission_date():
    db = MongoClient(host='129.32.209.134')
    acceptable = ['id', 'submitted']
    acceptable_dict = { k: 1 for k in acceptable}
    mongo_query = {}
    most_recent_submission_date = list(db.gisaid.records.find(mongo_query, acceptable_dict).sort([("submitted", -1)]).limit(1))
    return most_recent_submission_date[0]['submitted']

def split_metadata_based_on_date(input, output, date):
    # read csv
    with open(input) as input_fh:
       # Unfortunately, not sorted by submission date
       metadata = list(csv.reader(input_fh, delimiter='\t', quotechar='|'))
       i = metadata[0].index('date_submitted')
       sorted_meta = sorted(metadata[1:], key=lambda x:  datetime.strptime(x[i], "%Y-%m-%d"))
       # Find first index of date and split it out
       new_meta = list(itertools.dropwhile(lambda x: datetime.strptime(x[i],"%Y-%m-%d") < date, sorted_meta))

       # Write new meta to output file
       with open(output, 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter='\t',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            # Header file
            spamwriter.writerow(metadata[0])
            for meta in new_meta:
                spamwriter.writerow(meta)

       # Return new_meta for further use
       return metadata[0] + new_meta

def split_fasta_based_on_names(input, output, meta):
    # Get index of name
    i = meta[0].index('strain')
    names = [x[i] for x in meta[1:]]
    seq_records = []

    # *DO NOT* load entire FASTA into memory, just pick off based on names
    with open(input) as handle:
        for rec in SimpleFastaParser(handle):
            strain_name = '/'.join(rec[0].split('/')[1:]).split('|')[0]
            if rec[0] in names:
                seq_records.append(SeqRecord(Seq(rec[1]),id=strain_name,name='',description=''))

    # Write new FASTA with just new records
    with open(output, 'w', encoding='utf-8') as output_fh:
        SeqIO.write(seq_records, output_fh, "fasta")

def filter_gisaid_exports(meta_input, meta_output, fasta_input, fasta_output):
    most_recent_submission_date = get_most_recent_submission_date()

    # Return names of new metadata as well
    new_meta = split_metadata_based_on_date(meta_input, meta_output, most_recent_submission_date)
    split_fasta_based_on_names(fasta_input, fasta_output, new_meta)

def filter_gisaid_exports_by_dir(dir, meta_output, fasta_output):
    meta_input = max(glob.glob(dir + '*.tsv'), key=os.path.getctime)
    fasta_input = max(glob.glob(dir + '*.fasta'), key=os.path.getctime)
    most_recent_submission_date = get_most_recent_submission_date()
    # Return names of new metadata as well
    new_meta = split_metadata_based_on_date(meta_input, meta_output, most_recent_submission_date)
    split_fasta_based_on_names(fasta_input, fasta_output, new_meta)


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-m', '--meta-input', help = 'tsv metadata', type = str)
    arguments.add_argument('-n', '--meta-output', help = 'output tsv metadata', type = str)
    arguments.add_argument('-f', '--fasta-input', help = 'fasta input', type = str)
    arguments.add_argument('-g', '--fasta-output', help = 'fasta output', type = str)
    args = arguments.parse_args()
    filter_gisaid_exports(args.meta_input, args.meta_output, args.fasta_input, args.fasta_output)
    # filter_gisaid_exports_by_dir('./data/to-import/', args.meta_output, args.fasta_output)


