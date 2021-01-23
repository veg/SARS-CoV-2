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

writer = csv.writer (sys.stdout)
writer.writerow (['Gene','Date','Value','Type'])

for date in all_dates:
    with open (date[1], "r") as fh:
        summary_json = json.load (fh)
        for gene,data in summary_json.items():
            if gene != 'dates':
                #['sequences', 'total sequences', 'aminoacid variant sites', 'all variant sites', 'duplicate_map', 'any variation', 'sites', 'tree', 'MAF', 'L', 'p', 'selection', 'map', 'dN/dS']
    
                df =  date[0].strftime("%Y-%m-%d")
                if data['dN/dS']['internal']:
                   writer.writerow  ([gene, df, "%g" % data['dN/dS']['internal'][0][0], "dN/dS Internal"])
                if data['dN/dS']['leaves']:
                    writer.writerow  ([gene, df, "%g" % data['dN/dS']['leaves'][0][0], "dN/dS Leaves"])
                writer.writerow  ([gene, df, "%g" % data['L'], "Internal Branch Length"])
                writer.writerow  ([gene, df, "%d" % data['sequences'], "Unique Haplotypes"])
                writer.writerow  ([gene, df, "%d" % data['total sequences'], "Sequences"])
                writer.writerow  ([gene, df, "%g" % (data['any variation'] / data['sites']), "Variable sites"])
                writer.writerow  ([gene, df, "%g" % (len (data['aminoacid variant sites']) / data['sites']), "Variable amino-acids"])
                writer.writerow  ([gene, df, "%g" % (len ([i for i,s in data['selection'].items() if 'kind' in s and s['kind'] == "negative"]) / data['sites']), "Negatively Selected"])
                writer.writerow  ([gene, df, "%g" % (len ([i for i,s in data['selection'].items() if 'kind' in s and s['kind'] == "positive" or s['meme'] <= 0.05]) / data['sites']), "Positively Selected"])
    
    
                '''
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
                '''
    
