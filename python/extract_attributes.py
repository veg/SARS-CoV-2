"""
Extract sequence attributes as a CSV

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2020-03-22)


"""

import argparse
import sys
import json
import re
import datetime
import os
import csv




#-------------------------------------------------------------------------------


arguments = argparse.ArgumentParser(description='Extract attributes as CSV.')

arguments.add_argument ('-j', '--json', help = 'The master JSON cache file where sequence records live', required = True, type = str)
arguments.add_argument ('-f', '--field', help = 'Fields to extract in addition to the basics', action = 'append', default = ['id','location','submitted','collected'])

import_settings = arguments.parse_args()

current_sequence_db = {}
try:
    with open (import_settings.json, 'r') as cache:
        try:
            current_sequence_db = json.load (cache)
        except ValueError:
               print ('Error parsing JSON; please fix the file format and try again', file = sys.stderr)
               sys.exit (1)
except:
    pass

columns = []

csv_writer = csv.writer (sys.stdout)
first = True

for id, record in current_sequence_db.items():
    if first:
        for f in import_settings.field:
            if type (record[f]) == dict:
                for fk in sorted (record[f].keys()):
                    columns.append (fk)
            else:
                columns.append (f)
        csv_writer.writerow (columns)
                
        first = False
    
    row = []
    for f in import_settings.field:
        if type (record[f]) == dict:
            for fk in sorted (record[f].keys()):
                row.append (record[f][fk])
        else:
            row.append (record[f])
    csv_writer.writerow (row)
    
