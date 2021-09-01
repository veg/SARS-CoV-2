import json
import argparse
from datetime import date, timedelta, datetime
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

    start_date='2020-06-01'
    end_date='2021-09-01'

    aggregate_query = [
        { "$match": { "submitted": { "$gte": datetime.strptime(start_date, "%Y-%m-%d"), "$lte": datetime.strptime(end_date, "%Y-%m-%d") } } },
        { "$group" : { "_id": {"S_premsa_nuc_seq": "$S_premsa_nuc_seq"}, "uniqueIds": {"$addToSet": "$id"}, "count": {"$sum": 1 } } },
        { "$match": { "count": {"$gt": 1} } },
        { "$sort": { "count" : -1 } },
        { "$merge": {
           "into": "dev_duplicate_agg",
           "on" : "_id",
           "whenMatched": [ { "$set": {
              "uniqueIds": { "$setUnion": ["$uniqueIds", "$$new.uniqueIds"] },
              "count": { "$size": { "$setUnion": ["$uniqueIds", "$$new.uniqueIds"] } },
           } } ],
           "whenNotMatched": "insert"
       }}
    ]

    # records = list(db.gisaid.dev.aggregate(aggregate_query, aggregate_opts))
    records = db.gisaid.dev.aggregate(aggregate_query, allowDiskUse=True)
    # import pdb;pdb.set_trace()
    records = [{k: v for k, v in rec.items()} for rec in records]
    print(records)

    # Now, instead of getting the entire collection, let's loop through submit dates and merge into dev_duplicates


    with open(output_fn, 'w') as output_fh:
        json.dump(records, output_fh, indent=4, sort_keys=True)

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-o', '--output', help = 'json output', type = str)
    arguments.add_argument('-g', '--gene', help = 'gene to output', type = str)
    args = arguments.parse_args()
    export_duplicates(args.output, args.gene)

