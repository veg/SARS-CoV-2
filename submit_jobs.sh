#!/bin/bash
fdate=$(date +"%Y-%m-%d")
#fdate=$1
#fdate='2020-09-01'

FQUEUE='epyc2'
QUEUE='epyc'
OUTPUT_DIR=`pwd`/logs/$fdate/
mkdir $OUTPUT_DIR

# Gene
qsub -l nodes=1:ppn=32 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate S 32" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate M 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate N 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate E 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF3a 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF6 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF7a 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF7b 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF8 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh 
#qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF1a 128" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
#qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF1b 128" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh

# Products
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate leader 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp2 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp3 128" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp4 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate 3C 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp6 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp7 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp8 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp9 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp10 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate RdRp 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate helicase 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate exonuclease 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate endornase 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate methyltransferase 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF10 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
