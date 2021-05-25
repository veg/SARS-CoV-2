import csv
import unicodedata
import json
import argparse
import os
from datetime import date, datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pymongo import MongoClient

genes = ["S","M","N","E","ORF3a","ORF6","ORF7a","ORF8","leader","nsp2","nsp3","nsp4","3C","nsp6","nsp7","nsp8","nsp9","nsp10","RdRp","helicase","exonuclease","endornase","methyltransferase"]

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


    # Just export ID
    return unicodedata.normalize('NFKD', record['id'])

def export_random(config):

    db = MongoClient(host='192.168.0.4')
    acceptable = ['address', 'genbank_accession', 'age', 'assembly', 'authors', 'originalCollected', 'coverage', 'length', 'gender', 'sex', 'host', 'id', 'lab', 'originating_lab', 'location', 'name', 'passage', 'seqLength', 'originalSubmitted', 'submitter', 'submitting_lab', 'technology', 'type', 'nextstrainClade', 'pangolinLineage', 'gisaidClade', 'bealign']
    acceptable_meta_output = ['address', 'genbank_accession', 'age', 'assembly', 'authors', 'originalCollected', 'coverage', 'length', 'gender', 'sex', 'host', 'id', 'lab', 'originating_lab', 'location', 'name', 'passage', 'seqLength', 'originalSubmitted', 'submitter', 'submitting_lab', 'technology', 'type', 'nextstrainClade', 'pangolinLineage', 'gisaidClade']

    # transform acceptable into mongo query
    acceptable_dict = { k: 1 for k in acceptable}
    HOST= "Human"

    output_dir = config["output-dir"]
    meta_output_fn = os.path.join(output_dir, 'meta.json')

    # All genes must have passed quality control
    mongo_query = {'qc.' + gene + '.passed': True for gene in genes}
    mongo_query['host'] = HOST

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
        mongo_query["collected"] = { "$gte": datetime.strptime(start_date, "%Y-%m-%d"), "$lte": datetime.strptime(end_date, "%Y-%m-%d") }
        mongo_query["originalCollected"] = { "$regex": "[0-9]{4}-[0-9]{2}" }

    # Construct random sample query
    db_mongo_query = db.gisaid.records.aggregate([
       { "$match": mongo_query},
       { "$project": acceptable_dict },
       { "$sample": { "size": config["size"] } }
    ])

    records = list(db_mongo_query)
    print(len(records))
    records = [{k: v for k, v in rec.items() if k in acceptable} for rec in records]
    meta_records = [{k: v for k, v in rec.items() if k in acceptable_meta_output} for rec in records]
    meta_output_json = { row["id"] : row for row in meta_records }

    # Write attributes.csv file and MASTER-NO-JSON file
    with open(meta_output_fn, 'w') as output_fh:
        json.dump(meta_output_json, output_fh, indent=4, sort_keys=True)

    # Write sequences with bealign
    for gene in genes:
        nuc_output_fn = os.path.join(output_dir, 'sequences.' + gene + '.fas')
        nuc_seq_records = []
        for rec in records:
            if 'bealign' in rec.keys():
                if gene in rec['bealign'].keys():
                    nuc_seq_records.append(SeqRecord(Seq(rec['bealign'][gene]),id=sequence_name(rec),name='',description=''))

        # Write to fasta
        with open(nuc_output_fn, 'w', encoding='utf-8') as nuc_output_fh:
            SeqIO.write(nuc_seq_records, nuc_output_fh, "fasta")

if __name__ == "__main__":
    arguments = argparse.ArgumentParser(description='Report which dates have full report')
    arguments.add_argument('-o', '--output', help = 'output directory', type = str)
    args = arguments.parse_args()
    config = {"output-dir" : args.output }
    # config["clades"] = ["B.1.351"]
    # config["clades"] = ["B.1.427", "B.1.429"]
    # config["clades"] = ["P.1"]
    # config["clades"] = ["B.1.1.7"]
    # config["ignore-clades"] = ["B.1.351", "P.1", "B.1.1.7"]
    # config["clade-type"] = "pangolinLineage"
    config["collection-date-range"] = ("2019-12-01", "2020-02-28")
    config["size"] = 10
    export_random(config)


