"""
Combine individual alignments into a composite FASTA file

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2020-05-02)


"""

import argparse
import sys
import json
import re
import datetime
import os
import math, csv
from   os import  path
from   Bio import SeqIO
import operator


arguments = argparse.ArgumentParser(description='Summarize selection analysis results.')

arguments.add_argument('-f', '--first',        help = 'Earlier file',    required = True, type = argparse.FileType('r'))
arguments.add_argument('-s', '--second',   help = 'Later file',      required = True, type = argparse.FileType('r'))

import_settings = arguments.parse_args()

combined_fasta  = {}

first = json.load (import_settings.first)
second = json.load (import_settings.second)

for id, data in first.items():
    if second[id]['sequence'].upper () != data['sequence'].upper():
        newer = second[id]['sequence'].upper ()
        #print (">%s_old\n%s\n>%s_new\n%s\n\n" % (id, data['sequence'].upper(), id, newer))
        for i, char in enumerate (data['sequence'].upper()):
            if i < len (newer):
                if char != newer[i] and char != 'N' and newer[i] != 'N' :
                    print ("%s %d %s %s" %(id, i, char, newer[i]))