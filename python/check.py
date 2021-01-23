import sys
import json
import argparse
import operator

arguments = argparse.ArgumentParser(description='Summarize selection analysis results.')

arguments.add_argument('-f', '--file',   help = 'File to process', required = True, type = str)
import_settings = arguments.parse_args()

with open (import_settings.file, "r") as fh:
    data = json.load (fh)
    for site, info in data.items():
        if ('bSC2' in info and 'cdn' in info):
            mx = max(info['cdn'].items(), key=operator.itemgetter(1))[0]
            if info['bSC2'] != mx:
                if len (set (mx) - set ('ACGT')) == 0:
                    print (json.dumps(info, indent = 1), "\n#####\n")
 