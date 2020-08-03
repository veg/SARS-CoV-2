from Bio import SeqIO
import sys
import argparse

parser = argparse.ArgumentParser(description='Merge two gisaid FASTA files')
parser.add_argument('-l', '--last', type=argparse.FileType('r'))
parser.add_argument('-c', '--current', type=argparse.FileType('r'))
parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout)

args = parser.parse_args()

last = list(SeqIO.parse(args.last, 'fasta'))

# collect descriptions
last_descs = [l.description for l in last]
current = list(SeqIO.parse(args.current, 'fasta'))

cnt=0

for c in current:
    if c.description not in last_descs:
        last.append(c)
        cnt+=1

print("Added {} sequences of {} current sequences for a new merged total of {}".format(cnt, len(current), len(last)))

# Get difference between two sets and merge.
SeqIO.write(last, args.outfile, "fasta")
