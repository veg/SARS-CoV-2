#!/bin/bash
fdate=$(date +"%Y-%m-%d")

#qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate S" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh > S.log 2>S.log
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate M" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh > M.log 2>M.log
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate N" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh > N.log 2>N.log
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF3a" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh > ORF3a.log 2>ORF3a.log
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF7a" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh > ORF7a.log 2>ORF7a.log
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF8" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh > ORF8.log 2>ORF8.log
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF1a" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh > ORF1a.log 2>ORF1a.log
qsub -d `pwd` -q epyc -F "/data/shares/veg/SARS-CoV-2/SARS-CoV-2/data/fasta/$fdate ORF1b" /data/shares/veg/SARS-CoV-2/SARS-CoV-2/scripts/extract_genes.sh > ORF1b.log 2>ORF1b.log
