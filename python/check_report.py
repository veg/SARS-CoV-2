import csv
import sys
import argparse
import itertools
from datetime import date, timedelta
from operator import itemgetter

## Print dates with complete reports
genes=['leader','nsp2','nsp3','nsp4','3C','nsp6','nsp7','nsp8','nsp9','nsp10','helicase','exonuclease','endornase','S','E','M','N','ORF3a','ORF6','ORF7a','ORF8','RdRp','methyltransferase']

arguments = argparse.ArgumentParser(description='Report which dates have full report')
arguments.add_argument('-i', '--input',   help = 'summary results file', required = True, type = argparse.FileType('r'))
arguments.add_argument('-o', '--output', help = 'Write results here', type = argparse.FileType('w'), default = sys.stdout)

args = arguments.parse_args()

results = [row for row in csv.DictReader(args.input, delimiter=',', quotechar='|')]
results = sorted(results, key=itemgetter('analysis_date'))

# group by dates
grouped = itertools.groupby(results, lambda x: x["analysis_date"])

valid_dates = []

for key, group in grouped:

    group = list(group)

    # Check for at least the number of genes
    if(len(group) < len(genes)):
        break

    # Check that each gene has num_seqs field
    valids = len(list(filter(lambda x: len(x['num_seqs']) > 0, group)))
    if valids >= len(genes):
        valid_dates.append(key)

## Print out good dates
csvwrite = csv.writer(args.output, delimiter=',')
for valid_date in valid_dates:
    csvwrite.writerow ([valid_date])


