import csv
import argparse
import operator
import sys
import json

arguments = argparse.ArgumentParser(description='Read a TSV file and compress in into a JSON with indexed key storage. This assumes a sufficient degree of field value repetitions')

arguments.add_argument('-i', '--input', help ='The TSV file to compress', required = True, type = str)
arguments.add_argument('-o', '--output', help ='The JSON file to compress the input file to', required = True, type = str)
arguments.add_argument('-s', '--switch', help ='If compression is below this cutoff, store values directly', required = False, type = float, default = 10.)

args = arguments.parse_args()

print ("Analyzing the TSV file", file = sys.stderr)
with open (args.input, "r") as fh:
    reader = csv.reader (fh, delimiter = '\t')
    headers = next (reader)
    labeled_columns = len (headers)
    header_uniques = {}
    for i,h in enumerate (headers):
        header_uniques[i] = {}
    l = 0
    for line in reader:
         for i,field in enumerate (line[:labeled_columns]):
            field_value = field
            try:
                field_value = int (field)
            except: 
                try:
                    field_value = float (field)
                except: 
                    pass
            
         
            if not field_value in header_uniques[i]:
                header_uniques[i][field_value] = 0
                
            header_uniques[i][field_value] += 1
            
         l += 1
         if l % 50000 == 0:
            print ("...Read %d lines" % l)
            
    print ("Read %d lines with %d columns" % (l, len (headers)), file = sys.stderr)
    for i,h in enumerate(headers):
        ratio = l / (0. + len (header_uniques[i]))
        print ("\tColumn %s has %d unique values (%g compression)" % (h, len (header_uniques[i]), ratio), file = sys.stderr)
        if ratio <= args.switch:
            header_uniques[i] = None
            print ("\t\tColumn %s has will be stored using direct values" % h)
        else:
            sort_headers = sorted(header_uniques[i].items(), key=lambda item: -item[1])
            header_uniques[i] = {}
            for i2, pair in enumerate (sort_headers):
                header_uniques[i][pair[0]] = i2
            
        
print ("Compressing to JSON", file = sys.stderr)

with open (args.input, "r") as fh:
    reader = csv.reader (fh, delimiter = '\t')
    next (reader)
    outputJSON = {
        'cols' : headers,
        'rows' :l,
        'keys' : []
    }
    
    flat_values = []
    
    for i,h in enumerate(headers):
        indexed_values = {}
        col_keys = []
        if header_uniques [i]:
            for v,i2 in header_uniques [i].items():
                indexed_values[i2] = v
            for k in range (len (header_uniques [i])):
                col_keys.append (indexed_values[k])
            outputJSON['keys'].append (col_keys)
        else:
            outputJSON['keys'].append (None)
        
    for line in reader:
        for i,field in enumerate (line[:labeled_columns]):
        
            field_value = field
            try:
                field_value = int (field)
            except:             
                try:
                    field_value = float (field)
                except: 
                    pass
                
            if header_uniques[i]:
                flat_values.append (header_uniques[i][field_value])
            else:   
                flat_values.append (field_value)
                 
    outputJSON ['values'] = flat_values
    print ("Writing %d flat attribute values to JSON" % (len (flat_values)), file = sys.stderr)
    
    with open (args.output, "w") as jh:
        json.dump (outputJSON, jh, separators=(',', ':'))
