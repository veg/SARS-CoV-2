"""
Convert CSV data into a JSON annotation file

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2021-01-12)


"""

import argparse
import sys
import json
import datetime
import math, csv, re
from   os import  path


date_output_format = "%Y-%m-%d"


arguments = argparse.ArgumentParser(description='Convert CSV data to JSON annotation')

arguments.add_argument('-i', '--input', help  = 'Read the metadata from', required = True, type = argparse.FileType('r'))
arguments.add_argument('-d', '--date', help  = 'Read dates in the following format', required = False, type = str, default = "%Y-%m-%d")
arguments.add_argument('-o', '--output', help = 'Write JSON to', required = True, type = argparse.FileType('w'))

import_settings = arguments.parse_args()

meta_data = csv.reader(import_settings.input, delimiter = ',')

annotated_json = {}
valid_accession = re.compile ("^[0-9A-Z]+$")

next (meta_data)
for l in meta_data:
    id = l[0]
    if valid_accession.match (id):
        location = l[1].split (':')[0].strip()
        try:
            sampled = datetime.datetime.strptime (l[2], import_settings.date)
        except:
            sampled = None
        
        annotated_json['epi_isl_%s' % id] = {
            'collected' : sampled.strftime (date_output_format) if sampled else None,
            'location' : {
                'country' : location if len (location) else None,
                'locality': None,
                'state' : None,
                'subregion' : None
            }
         }
    else:
        print ("%s is not a valid accession number" % id, file = sys.stderr)
        
    
       
json.dump (annotated_json, import_settings.output, indent = 0)