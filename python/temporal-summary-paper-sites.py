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
writer.writerow (['coordinate','date','p', 'gene', 'branches', 'fel', 'kind', 'maf'])

for date in all_dates:
    with open (date[1], "r") as fh:
        summary_json = json.load (fh)
        for gene,data in summary_json.items():
            if gene != 'dates':
                #['sequences', 'total sequences', 'aminoacid variant sites', 'all variant sites', 'duplicate_map', 'any variation', 'sites', 'tree', 'MAF', 'L', 'p', 'selection', 'map', 'dN/dS']
    
                df =  date[0].strftime("%Y-%m-%d")
                
                for s, site_data in data["selection"].items():
                    si = int (s)
                    if si >= 0:
                        coordinate = "%s %d" % (gene, data["map"][int(s)]+1)
                        if site_data['meme'] <= 0.05 or site_data['fel'] <= 0.5:
                            writer.writerow ([coordinate, df, "%g" % site_data['meme'], gene, "%g" % site_data['meme-branches'], "%g" % site_data['fel'],  site_data['kind'] if 'kind' in site_data else None, site_data['MAF']])
    
    
                
    
