import sys
import json
import argparse
import operator
import csv
import re
import os
from datetime import datetime

arguments = argparse.ArgumentParser(description='Summarize selection analysis over time.')

arguments.add_argument('-d', '--directory',   help = 'Directory to scan', required = True, type = str)
settings = arguments.parse_args()

time_format = "%Y-%m-%d"
report_name = "report.json"
dates = []

for root, dirs, files in os.walk(settings.directory):
    for analysis_file  in files:
        if analysis_file == report_name:
            try:
                date = datetime.strptime(os.path.split (root)[1], "%Y-%m-%d")
                dates.append ([date, os.path.join (root, analysis_file)])
            except ValueError as ve:
                pass
        
 
all_dates = sorted (dates, key = lambda x: x[0])

with open (all_dates[-1][1], "r") as fh:
    reference_json = json.load (fh)
    
reference_json['dates'] = [k[0].strftime (time_format) for k in all_dates]


for gene,data in reference_json.items():
    if gene != 'dates':
        #['sequences', 'total sequences', 'aminoacid variant sites', 'all variant sites', 'duplicate_map', 'any variation', 'sites', 'tree', 'MAF', 'L', 'p', 'selection', 'map', 'dN/dS']
        
        
        data ['temporal'] = {
            'sequences' :           [data['sequences']],
            'total sequences':      [data['total sequences']],
            'variants' :            [],
            'variants2x' :          [],
            'variantsAA' :          [],
            'variantsAA2x' :        [],
            'L' : [data['L']],
            'dN/dS internal' : [data['dN/dS']['internal'][0][0]],
            'dN/dS leaves' : [data['dN/dS']['leaves'][0][0]]
        }
        
        print (data ['temporal'])
        
        sys.exit (0)