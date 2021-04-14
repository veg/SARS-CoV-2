import sys
import json
import argparse
import operator
import datetime
import csv

arguments = argparse.ArgumentParser(description='Compress and obfuscate GISAID records.')

arguments.add_argument('-i', '--input',    help = 'Input database', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output',   help = 'Output file', required = True, type = argparse.FileType('w'))
arguments.add_argument('-m', '--map',      help = 'ID-accession map', required = True, type = argparse.FileType('a+'))
arguments.add_argument('-d', '--date',     help = 'Reference date (YYYYMMDD)', required = False, type = str, default = "20191201")

settings = arguments.parse_args()

date_parse_format = "%Y%m%d"
date_parse_format_db = "%Y-%m-%d"

db = json.load (settings.input)
obfuscated = {'reference' : settings.date,
              'values' : []}

ref_date = datetime.datetime.strptime (settings.date, date_parse_format)

record_data = {
    'date'     : {}, 
    'subregion' : {},
    'country': {},
    'state': {},
    'locality' : {},
    'pangolinLineage' : {},
    'nextstrainClade' : {}
}

id_map = []

try:
    settings.map.seek (0)
    id_map = json.load (settings.map)
except Exception as e:
    pass
    
existing_ids = {}
for i,id in enumerate (id_map):
    existing_ids[id] = i
    if id not in db:
        raise Exception ("Accession %s present in the mapping file has been deleted from the database. It is not possible to maintain a consistent map of accessions to IDs in this setting. Run with a new mapping file" % id)
    
settings.map.truncate (0)

loc_tags = ['subregion', 'country', 'state', 'locality']
field_order = ['date', 'subregion', 'country', 'state', 'locality','pangolinLineage','nextstrainClade']
old_clades = set (['19A','19B','20A'])

now = datetime.datetime.now()
wk = datetime.timedelta(weeks=1)

print ("Loaded %d records" % len (db), file = sys.stdout)

for id, record in db.items():
    date_check = None
    try:
        if 'collected' in record:
            date_check = datetime.datetime.strptime (record['collected'],date_parse_format_db)
        if 'originalCollected' in record:
            date_check = datetime.datetime.strptime (record['originalCollected'],date_parse_format_db)
        if date_check.year < 2019 or date_check.year == 2019 and date_check.month < 10 or date_check >= now or (date_check.year == 2020 and date_check.month == 1 and date_check.day == 1): 
            date_check = None
        else: 
            if 'originalSubmitted' in record:
                date_submitted = datetime.datetime.strptime (record['originalSubmitted'],date_parse_format_db)
                if (date_submitted-date_check) / wk > 52 and record["nextstrainClade"] not in old_clades:
                    print ("Over a year between sampling and submission", record["id"], record["pangolinLineage"], record["nextstrainClade"], record["location"]["country"])
                    date_check = None
            if date_check:     
                date_check = (date_check-ref_date).days
                if date_check < 0:
                    date_check = None
                    
    except Exception as e:
        #print (e)
        date_check = None
        
    if date_check in record_data['date']:
        record_data['date'][date_check]  += 1
    else:
        record_data['date'][date_check] = 1
        
    record['date'] = date_check
        
    loc_data = [None,None,None,None]
        
    if 'location' in record:
        for i,k in enumerate (loc_tags):
            loc_data [i] = record['location'][k]
            
            
    for i,k in enumerate (loc_tags):
        record[k] = loc_data[i]
        if loc_data[i] in record_data[k]:
            record_data[k][loc_data[i]] += 1
        else:
            record_data[k][loc_data[i]]  = 1

    for key in ['pangolinLineage','nextstrainClade']:
        lineage = None
        if key in record:
            lineage = record[key]
        else:
            record[key] = None
        
        if lineage in record_data[key]:
            record_data[key][lineage]  += 1
        else:
            record_data[key][lineage] = 1
            
    
for field in field_order:
    sorted_dict = (sorted(list(record_data[field].items()), reverse = True, key=operator.itemgetter(1)))
    obfuscated["values"].append ([k[0] for k in sorted_dict])
    for i,k in enumerate (sorted_dict):
        record_data [field][k[0]] = i
        
    
current_id = len (existing_ids)
id_map = [None for k in range (len (db))]
obfuscated['meta'] = [[] for  k in range (len (db))]

for id, record in db.items():
    if id in existing_ids:
        i = existing_ids[id]
    else:
        i = current_id
        current_id += 1
        
    id_map [i] = id    
    obfuscated['meta'][i] = [
        record_data [field][record[field]] for field in field_order
    ]
    
json.dump (obfuscated, settings.output, separators=(',', ':'), indent = None)
json.dump (id_map, settings.map, separators=(',', ':'), indent = None)
   
        

