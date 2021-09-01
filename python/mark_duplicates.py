import json
import argparse
import itertools
import math

from datetime import date, timedelta, datetime
from dateutil.relativedelta import *
from dateutil.rrule import *
from dateutil.parser import *

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


def verify_duplicates(gene):
    db = MongoClient(host='192.168.0.4')

    original_records = set([rec['id'] for rec in list(db.gisaid.dev.find({}, { 'id' : 1, '_id': False }))])
    records = set(itertools.chain(*[rec['uniqueIds'] for rec in list(db.gisaid.dev_duplicate_agg_S.find({}, { 'uniqueIds' : 1, '_id': False }))]))

    missing_ids = list(original_records.difference(records))
    seq = list(db.gisaid.dev.find({'id': missing_ids[0]}))[0]['S_premsa_nuc_seq']
    hi = list(db.gisaid.dev_duplicate_agg_S.find({'_id.S_premsa_nuc_seq' : seq}, { 'uniqueIds' : 1, '_id': False }))

    import pdb; pdb.set_trace()
    print(records)

def get_all_count():
    db = MongoClient(host='192.168.0.4')
    return db.gisaid.records.find({}, {'_id': 1}).count()

def aggregate_duplicates(gene, month, count_index, starting_id, month_index):

    db = MongoClient(host='192.168.0.4')

    gene_id = gene + '_premsa_nuc_seq'

    DEV_DB = "dev_duplicate_agg_"
    PROD_DB = "record_duplicate_agg_"

    start_date = month[0]
    end_date = month[1]

    # Create new documents for every 500k sequences
    COLLECTION_INDEX = 500000

    db_string = PROD_DB + gene
    # db_string = DEV_DB + gene

    didx = '_'.join([str(month_index), str(math.floor(count_index/COLLECTION_INDEX))])

    # Use date instead of skip
    aggregate_query = [
        { "$match": { "collected": { "$gte": datetime.strptime(start_date, "%Y-%m-%d"), "$lte": datetime.strptime(end_date, "%Y-%m-%d") } } },
        { "$skip": count_index },
        { "$limit" : 50000 },
        { "$group" : { "_id": {didx + '_' + gene_id: "$" + gene_id}, "uniqueIds": {"$addToSet": "$id"}, "count": {"$sum": 1 } } },
        { "$match": { "count": {"$gt": 0} } },
        { "$sort": { "count" : -1 } },
        { "$merge": {
           "into": db_string,
           "on" : "_id",
           "whenMatched": [ { "$set": {
              "uniqueIds": { "$setUnion": ["$uniqueIds", "$$new.uniqueIds"] },
              "count": { "$size": { "$setUnion": ["$uniqueIds", "$$new.uniqueIds"] } },
              "month": month_index,
              "index": math.floor(count_index/COLLECTION_INDEX),
              "seq" : "$" + gene_id
           } } ],
           "whenNotMatched": "insert"
       }}
    ]

    # records = db.gisaid.dev.aggregate(aggregate_query, allowDiskUse=True)
    records = db.gisaid.records.aggregate(aggregate_query, allowDiskUse=True)

def get_month_id(month):
    try:
        db = MongoClient(host='192.168.0.4')
        records = list(db.gisaid.records.find({ "collected": { "$gte": datetime.strptime(month[0], "%Y-%m-%d"), "$lte": datetime.strptime(month[1], "%Y-%m-%d") } }, { 'id' : 1 }).limit(1))
        earliest_id = records[0]['_id']
        count = db.gisaid.records.find({ "collected": { "$gte": datetime.strptime(month[0], "%Y-%m-%d"), "$lte": datetime.strptime(month[1], "%Y-%m-%d") } }, { 'id' : 1 }).count()
        return earliest_id, count
    except:
        return False, 0



if __name__ == "__main__":

    arguments = argparse.ArgumentParser(description='Mark duplicates in MongoDB')
    # arguments.add_argument('-d', '--dupe-input',   help = 'fasta to update', required = True, type = str)
    arguments.add_argument('-t', '--type',   help = 'gene region', required = True, type = str)
    args = arguments.parse_args()

    # verify_duplicates(args.type)

    TODAY = date.today()
    THISMONTH = TODAY-relativedelta(day=31)
    ONEMONTHAGO = TODAY-relativedelta(months=+1, day=1)


    starts = [dt.strftime('%Y-%m-%d') for dt in rrule(MONTHLY,dtstart=parse("20191101T000000"), until=ONEMONTHAGO)]
    ends = [dt.strftime('%Y-%m-%d') for dt in rrule(MONTHLY,dtstart=parse("20191201T000000"), until=THISMONTH)]
    months = sorted(set(list(zip(starts,ends))), reverse=False)

    for idx in range(len(months)):
        month = months[idx]
        earliest_id, count = get_month_id(month)
        # split by 30000 each
        print(earliest_id)
        print(count)
        lst = list(range(count))[0::50000]
        # Need to get each month, then first id with that month and count, then loop
        for x in lst:
            print('processing ' + str(x))
            try:
                aggregate_duplicates(args.type, month, x, earliest_id, idx)
            except:
                print('could not aggregate counts for ' + x)


