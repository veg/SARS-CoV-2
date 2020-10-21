#!/bin/bash
fdate=$(date +"%Y-%m-%d")
#fdate=$1
#fdate='2020-09-22'
#fdate='2020-10-12'
#fdate='2020-10-05'

FQUEUE='epyc2'
QUEUE='epyc'
OUTPUT_DIR=`pwd`/logs/$fdate/
mkdir $OUTPUT_DIR

# Gene
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate S 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate M 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate N 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate E 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF3a 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF6 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF7a 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF7b 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF8 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh 
#qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF1a 128" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
#qsub -l nodes=2:ppn=64 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF1b 128" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh

# Products
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate leader 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp2 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=16 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $FQUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp3 16" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp4 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate 3C 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp6 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp7 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp8 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp9 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate nsp10 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate RdRp 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate helicase 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate exonuclease 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate endornase 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate methyltransferase 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -l nodes=1:ppn=8 -d `pwd` -o $OUTPUT_DIR -e $OUTPUT_DIR -q $QUEUE -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF10 8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh

