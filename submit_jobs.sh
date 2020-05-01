#!/bin/bash
#fdate=$(date +"%Y-%m-%d")
fdate="2020-04-29"

qsub -o `pwd`/ORF1a.$fdate.stdout -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF1a" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -o `pwd`/ORF1b.$fdate.stdout -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF1b" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -o `pwd`/S.$fdate.stdout -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate S" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -o `pwd`/M.$fdate.stdout -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate M" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -o `pwd`/N.$fdate.stdout -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate N" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -o `pwd`/ORF3a.$fdate.stdout -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF3a" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -o `pwd`/ORF7a.$fdate.stdout -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF7a" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
qsub -o `pwd`/ORF8.$fdate.stdout -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh 
qsub -o `pwd`/ORF6.$fdate.stdout -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF6" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh
