#!/bin/bash
fdate=$(date +"%Y-%m-%d")
#fdate="2020-04-05"

qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate S" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate M" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate N" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF3a" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF7a" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh 
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF1a" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF1b" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
