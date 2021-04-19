"""
Perform a comparative analysis of SARS-CoV-2 clades

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2021-01-02)
    v0.0.2 (2021-04-15) -M mode, reordering execution modules


"""

import argparse
import sys
import json
import re, time, shutil
import datetime, random
import os
import math, csv
from   os import  path
import operator
import compress_json
import subprocess
from progress.bar import Bar
from termcolor import colored, cprint

from   Bio import SeqIO

task_runners = {
    'hyphy' : '/Users/sergei/Development/hyphy/hyphy LIBPATH=/Users/sergei/Development/hyphy/res',
    'hyphy-mpi' : '/Users/sergei/Development/hyphy/hyphy LIBPATH=/Users/sergei/Development/hyphy/res',
    'bealign' : 'bealign' ,
    'bam2msa' : 'bam2msa' ,
    'tn93-cluster' : '/Users/sergei/Development/tn93/tn93-cluster',
    'tn93' : '/Users/sergei/Development/tn93/tn93',
    'raxml-ng' : 'raxml-ng --force '
}
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord




arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')

arguments.add_argument('-i', '--input',            help = 'FASTA file to process',                         required = True, type = str )
arguments.add_argument('-o', '--output',           help = 'Directory for output and working files',        required = True, type = str)
arguments.add_argument('-g', '--gene',             help = 'Gene/ORF name',                                 required = True, type = str)
arguments.add_argument('-r', '--reference',        help = 'Reference sequence for gene',                   required = True, type = str)
arguments.add_argument('-O', '--other',            help = 'FASTA file with other sequences to analyze',    required = True, type = str )
arguments.add_argument('-m', '--max_reference',    help = 'The maximum number of reference sequences to retain',    required = False, type = int, default = 500)
arguments.add_argument('-q', '--max_query',        help = 'The maximum number of query sequences to retain',    required = False, type = int, default = 200)
arguments.add_argument('--threshold_query',        help = 'Distance threshold for clustering query sequences',    required = False, type = float, default = 0.0005)
arguments.add_argument('--threshold_ref',          help = 'Initial threshold for clustering reference sequences',    required = False, type = float, default = 0.001)
arguments.add_argument('-P', '--paths',            help = 'Read command settings from this file',    required = False, type = str)
arguments.add_argument('-L', '--label',            help = 'Label query sequences in tree as follows', required = True, type = str )
arguments.add_argument('-A', '--ambigs',           help = 'Remove ambigs from sequence files', action='store_true' )
arguments.add_argument('-M', '--minimal',           help = 'Run only the minimal analyses', action='store_true' )

settings = arguments.parse_args()

strike_ambigs = settings.ambigs
print ("Strike ambigs : ", strike_ambigs, file = sys.stderr)
run_minimal = settings.minimal
print ("Run minimal : ", run_minimal, file = sys.stderr)

if settings.paths:
    with open (settings.paths, "r") as fh:
        task_runners = json.load (fh)

input_stamp = os.path.getmtime(settings.input)

print ("Input file %s modified on " % settings.input, time.ctime(input_stamp))
with open (settings.reference) as fh:
    for l in fh:
        if l[0] == '>':
            _ref_seq_name = l[1:].split (' ')[0].strip()
            break
            
print ("Reference seq_name %s" % _ref_seq_name)
        

def check_result_file (file_name, tag):
    full_file_name = os.path.join (settings.output, "%s.%s" % (settings.gene, file_name))
    if os.path.isfile(full_file_name) and os.path.getsize(full_file_name) > 0:
        file_stamp =  os.path.getmtime(full_file_name)
        if file_stamp >= input_stamp:
            print (colored('%s already done in %s' % (tag, full_file_name), 'green'))
            file_stamp = input_stamp
            return (full_file_name, True) # already done
        print (colored('%s exists in %s but is OLDER than the input/parent file' % (tag, full_file_name), 'red'))            
        
    if not os.path.isdir (settings.output):
        os.mkdir (settings.output)
    return (full_file_name, False)
    
    
def run_command (exec, arguments, filename, tag):    
    print (colored('Running ... %s\n' % (tag), 'cyan'))
    print ("\t", colored('Command ... %s\n' % (" ".join ([exec] + arguments)), 'yellow'))
    result = os.system (" ".join ([exec] + arguments))
    if result != 0:
        #raise Exception ('Command exection failed code %s' % result)
        print ('Command exection failed code %s' % result)
        return input_stamp
    return os.path.getmtime(filename)
    
def cluster_to_fasta (in_file, out_file, ref_seq = None):
    with open (in_file, "r") as fh:
        cluster_json = json.load (fh)
        print (colored('Running ... converting representative clusters to .FASTA\n', 'cyan'))
        check_uniq = set ()
        with open (out_file, "w") as fh2:
            for c in cluster_json:
                cc = c['centroid'].split ('\n')
                if ref_seq:
                    if ref_seq in c['members']:
                        cc[0] = ">" + ref_seq
                        print ("\n".join (cc), file = fh2)
                        continue
                seq_id = cc[0]
                while seq_id in check_uniq:
                        seq_id = cc[0] + '_' + ''.join(random.choices ('0123456789abcdef', k = 10))
                check_uniq.add (seq_id) 
                print (seq_id,"\n",cc[1], file = fh2)
           
    return (os.path.getmtime(out_file), len (cluster_json))
        


query_bam, done = check_result_file ("query.bam", "Aligning query sequences to reference")
if not done:
    input_stamp = run_command (task_runners['bealign'], ['-r', settings.reference, '-m', 'HIV_BETWEEN_F', settings.input, query_bam], query_bam, "create .BAM alignment to reference")
    
query_msa, done = check_result_file ("query.msa", "Mapping BAM to MSA")
if not done:
    input_stamp = run_command (task_runners['bam2msa'], [query_bam, query_msa], query_msa, "map BAM to MSA")
    
if strike_ambigs:
    query_msa_cpy = query_msa + ".bk" 
    shutil.copy (query_msa, query_msa_cpy)
    run_command (task_runners['hyphy'], [os.path.abspath ("lib/strike-ambigs.bf"), '--alignment', os.path.abspath (query_msa_cpy), '--output', os.path.abspath (query_msa)  ], query_msa, "remove ambigs")


query_cluster, done = check_result_file ("query.json", "Extracting representative clusters")

#if not done:
#    input_stamp = run_command (task_runners['tn93-cluster'], ['-f', '-o', query_cluster, '-t', '%g' % settings.threshold_query, query_msa], query_cluster, "extract representative clusters")
    
ref_threshold = settings.threshold_query
ref_step      = ref_threshold*0.25
    
query_compressed, done = check_result_file ("query_compressed.fas", "Converting representative clusters to .FASTA")
if not done:
    query_msa_cpy = query_msa + ".bk" 
    shutil.copy (query_msa, query_msa_cpy)
    while True:
        input_stamp = run_command (task_runners['tn93-cluster'], ['-f', '-o', query_cluster, '-t', "%g" % ref_threshold, query_msa_cpy], query_cluster, "extract representative clusters at threshold %g" % ref_threshold)
        input_stamp, cluster_count = cluster_to_fasta (query_cluster, query_msa_cpy)
        print (colored('Found %d clusters' % cluster_count, 'yellow'))
        if cluster_count <= settings.max_query:
            shutil.copy (query_msa_cpy, query_msa)
            break
        else:
            ref_threshold += ref_step
    os.remove (query_msa_cpy)

    input_stamp, cluster_count = cluster_to_fasta (query_cluster, query_compressed)
    

ref_bam, done = check_result_file ("reference.bam", "Aligning referemce sequences to reference")
if not done:
    input_stamp = run_command (task_runners['bealign'], ['-r', settings.reference, '-m', 'HIV_BETWEEN_F', '-K', settings.other, ref_bam], ref_bam, "create .BAM alignment to reference")
    
ref_msa, done = check_result_file ("reference.msa", "Mapping reference BAM to MSA")
if not done:
    input_stamp = run_command (task_runners['bam2msa'], [ref_bam, ref_msa], ref_msa, "map reference BAM to MSA")
      
ref_cluster, done = check_result_file ("ref.json", "Extracting representative clusters from the reference")
ref_threshold = settings.threshold_ref
ref_step      = ref_threshold*0.25

if not done:
    ref_msa_cpy = ref_msa + ".bk" 
    shutil.copy (ref_msa, ref_msa_cpy)
    while True:
        input_stamp = run_command (task_runners['tn93-cluster'], ['-f', '-o', ref_cluster, '-t', "%g" % ref_threshold, ref_msa_cpy], ref_cluster, "extract representative clusters at threshold %g" % ref_threshold)
        input_stamp, cluster_count = cluster_to_fasta (ref_cluster, ref_msa_cpy, _ref_seq_name)
        print (colored('Found %d clusters' % cluster_count, 'yellow'))
        if cluster_count <= settings.max_reference:
            shutil.copy (ref_msa_cpy, ref_msa)
            break
        else:
            ref_threshold += ref_step
    os.remove (ref_msa_cpy)
        
if strike_ambigs:
    ref_msa_cpy = ref_msa + ".bk" 
    shutil.copy (ref_msa, ref_msa_cpy)
    run_command (task_runners['hyphy'], [os.path.abspath ("lib/strike-ambigs.bf"), '--alignment', os.path.abspath (ref_msa_cpy), '--output', os.path.abspath (ref_msa)  ], ref_msa, "remove ambigs")


combined_msa, done = check_result_file ("combined.fas", "Combining and filtering reference")
if not done:
    pairwise = combined_msa + ".csv"
    input_stamp = run_command (task_runners['tn93'], ['-o', pairwise, '-s', ref_msa, '-t', "%g" % (settings.threshold_query*2.), query_compressed], pairwise, "filtering reference sequuences that are closer than %g to any query cluster" % (settings.threshold_query*2.))
    with open (pairwise) as fh:
        reader = csv.reader (fh, delimiter = ',')
        next (reader)
        seqs_to_filter = set ()
        for l in reader:
            seqs_to_filter.add (l[1])
        if _ref_seq_name in seqs_to_filter:
            seqs_to_filter.remove (_ref_seq_name)
    shutil.copy (query_compressed, combined_msa)
    with open (combined_msa, "a+") as fh:
        check_uniq = set () 
        for seq_record in SeqIO.parse(ref_msa, "fasta"):
            if not seq_record.name in seqs_to_filter:
                if seq_record.name == _ref_seq_name:
                    print ("\n>%s\n%s" % ("REFERENCE", str(seq_record.seq)), file = fh)   
                    check_uniq.add ('REFERENCE')             
                else:
                    seq_id = seq_record.name
                    while seq_id in check_uniq:
                        seq_id = seq_record.name + '_' + ''.join(random.choices ('0123456789abcdef', k = 10))
                    check_uniq.add (seq_id) 
                    print ("\n>%s\n%s" % (seq_id, str(seq_record.seq)), file = fh)
           
    os.remove (pairwise)        
    input_stamp = os.path.getmtime(combined_msa)
    
    
combined_tree, done = check_result_file ("combined.fas.raxml.bestTree", "Inferring tree topology")
if not done:
    input_stamp = run_command (task_runners['raxml-ng'], ["--model GTR", "--tree pars{3}", "--msa", "%s" % combined_msa], combined_tree, "infer tree topology")


labeled_tree, done = check_result_file ("int.nwk", "Labeling the tree")
labeled_tree_clade, done_clade = check_result_file ("clade.nwk", "Labeling the tree")
labeled_tree_full, done_clade = check_result_file ("full.nwk", "Labeling the tree")

if not done or not done_clade:
    input_stamp = run_command (task_runners['hyphy'], [os.path.abspath ("lib/annotator.bf"), combined_tree, 'REFERENCE', query_compressed, settings.label, os.path.abspath(check_result_file ("", "")[0])], labeled_tree, "label tree topology")




slac, done = check_result_file ("SLAC.json", "SLAC analysis")
if not done:
    input_stamp = run_command (task_runners['hyphy'], ['slac', '--alignment', os.path.abspath (combined_msa), "--tree", os.path.abspath (labeled_tree), "--output", os.path.abspath (slac), "--samples", "0"], slac, "SLAC analysis")



prot_msa, done = check_result_file ("aa.fas", "Translate to amino-acids")

if not done or not done_clade:
    input_stamp = run_command (task_runners['hyphy'], ["conv", 'ENV="NORMALIZE_SEQUENCE_NAMES=0"' , "Universal", '"Keep Deletions"', os.path.abspath(combined_msa), os.path.abspath(prot_msa)], prot_msa, "translate to amino-acids")


bgm, done = check_result_file ("combined.fas.BGM.json", "BGM analysis")
if not done:
    input_stamp = run_command (task_runners['hyphy'], ['bgm', '--alignment', os.path.abspath (combined_msa), "--tree", os.path.abspath (labeled_tree), "--output", os.path.abspath (bgm), "--branches", settings.label], bgm, "BGM analysis")

fel_cache, done = check_result_file ("fit-cache", "")
fel, done = check_result_file ("FEL.json", "FEL analysis on the clade of interest")
if not done:
    input_stamp = run_command (task_runners['hyphy-mpi'], ['fel', '--alignment', os.path.abspath (combined_msa), "--tree", os.path.abspath (labeled_tree), "--output", os.path.abspath (fel), "--branches", settings.label], fel, "FEL analysis on the clade of interest")

meme, done = check_result_file ("MEME.json", "MEME analysis on the clade of interest")
if not done:
    input_stamp = run_command (task_runners['hyphy-mpi'], ['meme', '--alignment', os.path.abspath (combined_msa), "--tree", os.path.abspath (labeled_tree), "--output", os.path.abspath (meme), "--branches", settings.label], meme, "MEME analysis on the clade of interest")

meme_full, done = check_result_file ("MEME-full.json", "MEME analysis on the clade of interest")
if not done:
    input_stamp = run_command (task_runners['hyphy-mpi'], ['meme', '--alignment', os.path.abspath (combined_msa), "--tree", os.path.abspath (labeled_tree_full), "--output", os.path.abspath (meme_full), "--branches", settings.label], meme_full, "MEME analysis on the clade of interest")


prime, done = check_result_file ("PRIME.json", "PRIME analysis on the clade of interest")
if not done:
    input_stamp = run_command (task_runners['hyphy-mpi'], ['prime', '--alignment', os.path.abspath (combined_msa), "--tree", os.path.abspath (labeled_tree), "--output", os.path.abspath (prime), "--branches", settings.label], prime, "PRIME analysis on the clade of interest")

fade, done = check_result_file ("FADE.json", "FADE analysis")
if not done:
    input_stamp = run_command (task_runners['hyphy'], ['fade', '--alignment', os.path.abspath (prot_msa), "--tree", os.path.abspath (labeled_tree_clade), "--output", os.path.abspath (fade), "--branches", settings.label], fade, "FADE")

if not run_minimal:
    absrel, done = check_result_file ("ABSREL.json", "ABSREL analysis")
    if not done:
        input_stamp = run_command (task_runners['hyphy'], ['absrel', '--alignment',  os.path.abspath (combined_msa), "--tree", os.path.abspath (labeled_tree), "--output", os.path.abspath (absrel), "--branches", settings.label], absrel, "ABSREL")

    busted, done = check_result_file ("BUSTED.json", "BUSTED analysis")
    if not done:
        input_stamp = run_command (task_runners['hyphy'], ['busted', '--alignment', os.path.abspath (combined_msa), "--tree", os.path.abspath (labeled_tree_clade), "--output", os.path.abspath (busted), "--branches", settings.label, "--starting-points", "10"], busted, "BUSTED")

    relax, done = check_result_file ("RELAX.json", "RELAX analysis")
    if not done:
        input_stamp = run_command (task_runners['hyphy'], ['relax', '--models', 'Minimal', '--alignment',  os.path.abspath (combined_msa), "--tree", os.path.abspath (labeled_tree_clade), "--output", os.path.abspath (relax), "--test", settings.label, "--reference", "Reference", "--starting-points", "10", "--srv", "Yes"], relax, "RELAX")
    
cfel, done = check_result_file ("CFEL.json", "CFEL analysis on the clade of interest")
if not done:
    input_stamp = run_command (task_runners['hyphy-mpi'], ['contrast-fel', '--alignment', os.path.abspath (combined_msa), "--tree", os.path.abspath (labeled_tree_clade), "--output", os.path.abspath (cfel), "--branch-set", settings.label,  "--branch-set", "Reference"], labeled_tree, "CFEL analysis on the clade of interest")


