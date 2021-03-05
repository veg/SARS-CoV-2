## SARS-CoV-2 Clade analysis pipeline

### Overview.

The pipeline uses several open-source tools to prepare SARS-CoV-2 full length genome data for clade selection analysis using [HyPhy](hyphy.org), assemble the results of these analysis using a Python3 script into a monolithic JSON file which can then be visualized with an [Observable notebook](https://observablehq.com/@spond/revised-sars-cov-2-analytics-page). The results are then visualized to pick out patterns at each site, such as positive and negative selection, the history of selection at that site (when that site started becoming selected for), haplotype information, global distribution at that site, etc.  

<!-- these are commented out -->
<!--1. A (large) .FASTA file of unaligned genomes with sequence names formatted as `epi_isl_XXXX|...` where `XXXX` is a unique accession number and the values after the pipe are ignored.-->
<!--2.--> 
<!--3. A (large) .JSON file in a format described below which contains **metadata** about the sequences; the two are linked via the unique accession number. -->
<!--The metadata file is not necessary for sequence analysis, but it is "baked" into several processing scripts (will be unwound later), and may break things downsteam.-->

### Installation and dependancies

#### Environment 
Install the env using: ```pip3 install -r python/requirements.txt```. This will install:
```
numpy==1.19.4
compress_json==1.0.4
progress==1.5
Bio>=0.3.0
BioExt==0.19.8
mappy==2.17
seqmagick==0.8.4
termcolor==1.1.0
```

##### Binary packages

1. HyPhy 2.5.26+ (use latest version): github.com/veg/hyphy.git
2. MAFFT: http://mafft.cbrc.jp/alignment/software/ (can be swapped out for others with trival script changes; alignment is done on amino-acid sequences)
3. RapidNJ: https://birc.au.dk/software/rapidnj/ (can be swapped out for others with trival script changes)
4. jq: https://stedolan.github.io/jq/ JSON processing

##### HyPhy libraries

Assumes the availability of [HyPhy analyses](http://github.com/veg/hyphy-analyses) on the file system

### Pipeline 

1. FIRST change the ```paths``` in ```config.json```
2. Run ```submit_jobs.sh``` which sends jobs to ```run_gene.sh```
    * ```submit_jobs.sh```:
        1. this loops over each gene in SARS-CoV-2 and ```qsub```s that gene to ```run_gene.sh```
3. Then in ```run_gene.sh``` passes each gene as an argument to the ```run.py``` script. That python script runs on EACH gene iteratively:
    * Required Parameters:
        1. ```-i```, ```--input```: FASTA file to process (this is the user specific clade file)
        2. ```-o```, ```--output```: Directory for output and working files
        3. ```-g```, ```--gene```: Gene/ORF name
        4. ```-r``` ,```--reference```
        5. ```-O```, ```--other```: This is the codon-alignment for each gene from the supplied reference dataset (**might want to address this and keep it on our end, have a bubble list of references to select from?**) 
        6. ```-L```, ```--label```: Label query sequences in tree as follows

    * Optional Parameters:
        1. ```-m```, ```--max_reference```: the maximum number of reference sequences to retain | ```200```
        2. ```-q```, ```--max_query```: The maximum number of query sequences to retain | ```500```
        3. ```--threshold_query```: Distance threshold for clustering query sequences | ```0.0001```
        4. ```-A```, ```--ambigs```: Removes ambigs from sequence files, and stores them

4. Now: in ```run.py```, for each gene (in this order) the script will run:
    1. ```bealign``` to align query sequences to reference
    2. ```bam2msa``` to make a bam file
    3. ```tn93-cluster``` to extract representative clusters, this helps to compress the alignment
    4. ```tn93``` to pairwise cluster sequences in the combined MSA
    3. ```raxml``` to infer tree topology
    4. ```lib/annotator.bf``` to label the branches in the clade to analyze
    5. ```SLAC``` 
    6. run ```hyphy``` command to translate codons to amino acids
    7. ```BGM```
    8. ```BUSTED```
    9. ```FEL```
    10. ```MEME``` and ```MEME-full```
    11. ```PRIME```
    12. ```CFEL```
    13. ```FADE```
    14. ```ABSREL```
    15. ```RELAX```

5. Visualize results: ```lib/generate-report.py```:
    1. asdfasf

6. Go to this [Observable notebook](https://observablehq.com/@spond/revised-sars-cov-2-analytics-page). Here you will find information about YY
        

