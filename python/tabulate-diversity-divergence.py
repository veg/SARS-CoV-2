"""
Read in sequence information, a TN93 distance file, and report mean diversity and divergence
over a time period, using a sliding window

Date, Region, Diversity, Divergence

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2020-03-23)


"""

import argparse
import sys
import json
import re
import datetime
import os
import operator
import csv
import math


#-------------------------------------------------------------------------------


arguments = argparse.ArgumentParser(description='Tabulate diversity and divergence.')

arguments.add_argument ('-j', '--json', help = 'The master JSON cache file where sequence records live', required = True, type = str)
arguments.add_argument ('-t', '--tn93', help = 'The TN93 distance file',  required = True, type = str)
arguments.add_argument ('-d', '--duration', help = 'The binning duration in days',  required = False, type = int, default = 7)

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
    
sequences_with_dates = {}

for id, record in current_sequence_db.items():
    try:
        sequences_with_dates[id] =  datetime.datetime.strptime (record['collected'], "%Y%m%d")
        if sequences_with_dates[id].year < 2019 or sequences_with_dates[id].year == 2019 and sequences_with_dates[id].month < 10: 
            del sequences_with_dates[id]
    except:
        pass
        
date_range = [min (sequences_with_dates.values()),max (sequences_with_dates.values())]
step = datetime.timedelta (days = import_settings.duration)
bins =  int (math.ceil ((date_range[1]-date_range[0]).days / import_settings.duration))


diversity_by_bin  = [[0.,0] for k in range (bins+1)]
divergence_by_bin = [[0.,0] for k in range (bins+1)]

diversity_by_bin_country  = {}
divergence_by_bin_country = {}

def get_location (id):
    record = current_sequence_db[id]
    if 'country' in record['location'] :
        return record['location']['country']
    return None

with open (import_settings.tn93, 'r') as distances:
    csv_reader = csv.reader (distances)
    next (csv_reader)
    for line in csv_reader:
        seq1 = '_'.join (line[0].split ('_')[0:3])
        seq2 = '_'.join (line[1].split ('_')[0:3])
        if seq1 in current_sequence_db and seq2 in current_sequence_db:
            if seq1 in sequences_with_dates and seq2 in sequences_with_dates:
                bin1 = int (((sequences_with_dates[seq1] - date_range[0]).days / import_settings.duration) // 1)
                bin2 = int (((sequences_with_dates[seq2] - date_range[0]).days / import_settings.duration) // 1)
                if bin1 == bin2:
                    diversity_by_bin [bin1][1] += 1
                    diversity_by_bin [bin1][0] += float (line[2])
                    c = get_location (seq1)
                    c2 = get_location (seq2)
                    if c and c == c2:
                        if not c in diversity_by_bin_country:
                            diversity_by_bin_country [c] = [[0.,0] for k in range (bins+1)]
                        diversity_by_bin_country [c][bin1][1] += 1
                        diversity_by_bin_country [c][bin1][0] += float (line[2])
        else:
            db_seq = seq1 if seq1 in current_sequence_db else seq2
            if db_seq in sequences_with_dates:
                which_bin = int (((sequences_with_dates[db_seq] - date_range[0]).days / import_settings.duration) // 1)
                divergence_by_bin [which_bin][1] += 1
                divergence_by_bin [which_bin][0] += float (line[2])
                c = get_location (db_seq)
                if (c):
                    if not c in divergence_by_bin_country:
                        divergence_by_bin_country [c] = [[0.,0] for k in range (bins+1)]
                    divergence_by_bin_country [c][which_bin][1] += 1
                    divergence_by_bin_country [c][which_bin][0] += float (line[2])
        

current_date = date_range[0]
output       = csv.writer (sys.stdout)
output.writerow (["Date","Location","N","Diversity","Divergence"])

for i,k in enumerate (divergence_by_bin):
    if k[1] > 1:
        row = [current_date.strftime ("%Y%m%d"), "World", str (k[1]), str (diversity_by_bin[i][0]/diversity_by_bin[i][1]) , str (k[0]/k[1])]
        output.writerow (row)
        
        for country, divergences in divergence_by_bin_country.items():
            if divergences[i][1] > 1:
                row = [current_date.strftime ("%Y%m%d"), country, str (divergences[i][1]), str (diversity_by_bin_country[country][i][0]/diversity_by_bin_country[country][i][1]) , str (divergences[i][0]/divergences[i][1])]
                output.writerow (row)
   
    current_date += step
    
