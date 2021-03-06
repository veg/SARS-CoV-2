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


def update_record(nuc_seq_key, nuc_seq, prot_seq_key, prot_seq, gene):
    prot_seq_str = ''

    if(nuc_seq):
        nuc_seq_str = str(nuc_seq.seq)
    if(prot_seq):
        prot_seq_str = str(prot_seq.seq)

    try:
        epi_id = '_'.join(nuc_seq.description.split('_')[:3]).lower()
    except:
        print("could not process " + nuc_seq.description)
        return

    # print(epi_id)
    key_to_update = ".".join(["qc", gene])
    return UpdateOne({'id': epi_id}, {'$set': {nuc_seq_key: nuc_seq_str, prot_seq_key: prot_seq_str, key_to_update: { "passed" : True }}})

def store_premsa_file(nuc_input, prot_input, gene):

    db = MongoClient(host='192.168.0.4')

    nuc_input_fh = open(nuc_input, 'r')
    prot_input_fh = open(prot_input, 'r')

    seqs = list(SeqIO.parse(nuc_input_fh, 'fasta'))

    items_to_update = {}
    for seq in seqs:
        items_to_update[seq.id] = {}
        items_to_update[seq.id][gene + '_premsa_nuc_seq'] = seq
        items_to_update[seq.id][gene + '_premsa_protein_seq'] = ''


    seqs = list(SeqIO.parse(prot_input_fh, 'fasta'))

    for seq in seqs:
        items_to_update[seq.id][gene + '_premsa_protein_seq'] = seq

    to_update = [update_record(gene + '_premsa_nuc_seq', v[gene + '_premsa_nuc_seq'], gene + '_premsa_protein_seq', v[gene + '_premsa_protein_seq'], gene) for k,v in items_to_update.items()]
    if(len(to_update)):
        results = db.gisaid.records.bulk_write(to_update)
    print("Updated " +  str(results.modified_count) + " of " + str(len(items_to_update)) + " items to update")


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-n', '--nuc-input',   help = 'fasta to update', required = True, type = str)
    arguments.add_argument('-p', '--prot-input',   help = 'fasta to update', required = True, type = str)
    arguments.add_argument('-t', '--type',   help = 'gene region', required = True, type = str)
    args = arguments.parse_args()
    store_premsa_file(args.nuc_input, args.prot_input, args.type)

