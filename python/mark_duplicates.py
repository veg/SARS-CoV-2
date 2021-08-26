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

from datetime import date, timedelta, datetime
from dateutil.relativedelta import *
from dateutil.rrule import *
from dateutil.parser import *

from operator import itemgetter
from Bio import SeqIO
from pymongo import MongoClient, InsertOne, DeleteMany, ReplaceOne, UpdateOne, UpdateMany

def update_dupe_records(dupes, reference, gene):

    # Duplicates are computed by gene
    dupe_key = ".".join(["duplicate_of_by_gene", gene])
    dupe_val = { dupe_key : reference }
    return UpdateMany({'id': {"$in" : dupes}}, {'$set': dupe_val })

def update_reference_records(reference, gene):

    # Duplicates are computed by gene
    dupe_key = ".".join(["duplicate_of_by_gene", gene])
    dupe_val = { dupe_key : 'reference'}
    return UpdateOne({'id': reference}, {'$set': dupe_val })


def mark_duplicates(dupe_input, gene):

    db = MongoClient(host='192.168.0.4')
    dupes = json.loads(open(dupe_input, 'r').read())

    # Shave meta
    shave_meta = lambda x: '_'.join(x.split('_')[:3])

    # Transform ids to only include epi_isl_xxxxxxx
    fmt_dupes = { shave_meta(ref): list(map(shave_meta, dupe.values()))[1:]  for ref, dupe in dupes.items() if len(dupe) > 1}

    # Get documents for all references
    ref_keys = list(fmt_dupes.keys())

    # key references by id
    get_update_queries = [update_dupe_records(v, k, gene) for k,v in fmt_dupes.items()]
    get_update_reference_queries = [update_reference_records(k, gene) for k in ref_keys]

    results = db.gisaid.records.bulk_write(get_update_queries + get_update_reference_queries)
    print("Updated " +  str(results.modified_count) + " of " + str(sum([len(v) for k,v in fmt_dupes.items()]) + len(ref_keys)) + " items to update")

def aggregate_duplicates(gene, start_date, end_date):

    db = MongoClient(host='192.168.0.4')

    gene_id = gene + '_premsa_nuc_seq'

    DEV_DB = "dev_duplicate_agg_"
    PROD_DB = "record_duplicate_agg_"

    #db_prefix = PROD_DB
    db_prefix = DEV_DB

    aggregate_query = [
        { "$match": { "collected": { "$gte": datetime.strptime(start_date, "%Y-%m-%d"), "$lte": datetime.strptime(end_date, "%Y-%m-%d") } } },
        { "$group" : { "_id": {gene_id: "$" + gene_id}, "uniqueIds": {"$addToSet": "$id"}, "count": {"$sum": 1 } } },
        { "$match": { "count": {"$gt": 1} } },
        { "$sort": { "count" : -1 } },
        { "$merge": {
           "into": db_prefix + gene,
           "on" : "_id",
           "whenMatched": [ { "$set": {
              "uniqueIds": { "$setUnion": ["$uniqueIds", "$$new.uniqueIds"] },
              "count": { "$size": { "$setUnion": ["$uniqueIds", "$$new.uniqueIds"] } },
           } } ],
           "whenNotMatched": "insert"
       }}
    ]

    records = db.gisaid.dev.aggregate(aggregate_query, allowDiskUse=True)
    #records = db.gisaid.records.aggregate(aggregate_query, allowDiskUse=True)


if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Mark duplicates in MongoDB')
    # arguments.add_argument('-d', '--dupe-input',   help = 'fasta to update', required = True, type = str)
    arguments.add_argument('-t', '--type',   help = 'gene region', required = True, type = str)
    args = arguments.parse_args()
    # mark_duplicates(args.dupe_input, args.type)

    TODAY = date.today()
    # THISMONTH = TODAY-relativedelta(day=1)
    # ONEMONTHAGO = TODAY-relativedelta(day=1)

    THISMONTH = TODAY
    ONEMONTHAGO = TODAY

    print(THISMONTH)
    print(ONEMONTHAGO)

    starts = [dt.strftime('%Y-%m-%d') for dt in rrule(DAILY,dtstart=parse("20191101T000000"), until=ONEMONTHAGO)]
    ends = [dt.strftime('%Y-%m-%d') for dt in rrule(DAILY,dtstart=parse("20191101T000000"), until=THISMONTH)]
    months = set(list(zip(starts,ends)))

    for x in months:
        print(x)
        # try:
        aggregate_duplicates(args.type, x[0], x[1])
        # except:
            # print('could not aggregate counts for ' + x[0])


