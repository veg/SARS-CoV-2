
## SARS-CoV-2 selection pipeline

> **Disclaimer** This pipeline was a quick and dirty effort that continued evolving but was never properly productized (not even close). Expect breakage, tedious installation issues and poor documentation. We will work on improving it, but our efforts are geared towards speeding up core components and improving visualization.

### Overview.

The pipeline is uses several open-source tools to prepare SARS-CoV-2 full length genome data for selection analysis using [HyPhy](hyphy.org), assemble the results of these analysis using a Python3 script into a monolithic JSON file which can then be visualized with an [Observable notebook](https://observablehq.com/@spond/revised-sars-cov-2-analytics-page). The pipeline is used on [GISAID](gisaid.org) data which consists of 

1. A (large) .FASTA file of unaligned genomes with sequence names formatted as `epi_isl_XXXX|...` where `XXXX` is a unique accession number and the values after the pipe are ignored.
3. A (large) .JSON file in a format described below which contains **metadata** about the sequences; the two are linked via the unique accession number. 

The metadata file is not necessary for sequence analysis, but it is "baked" into several processing scripts (will be unwound later), and may break things downsteam.

### Installation and dependancies

##### Binary packages

1. HyPhy 2.5.26+ (use latest version): github.com/veg/hyphy.git
2. MAFFT: http://mafft.cbrc.jp/alignment/software/ (can be swapped out for others with trival script changes; alignment is done on amino-acid sequences)
3. RapidNJ: https://birc.au.dk/software/rapidnj/ (can be swapped out for others with trival script changes)
4. jq: https://stedolan.github.io/jq/ JSON processing

These can be installed anywhere on the system with paths specified in `extract_gene.sh`

##### Python modules and packages

```
numpy==1.19.4
compress_json==1.0.4
progress==1.5
Bio==0.3.0
BioExt==0.19.8
mappy==2.17
```

Install using `pip3 install -r python/requirements.txt` (see below for tree structure)

##### HyPhy libraries

Assumes the availability of [HyPhy analyses](http://
github.com/veg/hyphy-analyses) on the file system

### Configuration

This is the ugly part, sorry.
Once the pipeline is installed, two shell scripts need to be edited to point to the right paths (see below): `submit_jobs.sh` , `extract_gene.sh`


### Directory structure

The pipeline assumes a specific simple directory structure

```
|-- data
|   |-- evo-annotation.json
|   |-- ctl
|   |-- fasta
|   |   |-- YYYY-MM-DD
...
|-- logs
|-- python
|-- reference_genes
|-- scripts
```

The input file FASTA, named `sequences`, for an analysis run on a specific date (YYYY-MM-DD) is placed in the `fasta/YYYY-MM-DD`  directory as described above. Temporary and individual result files are placed in the same directory by the pipeline. 


#### Read-only data

`data/evo-annotation.json` contains precomputed annotation of SARS-CoV-2 positions based on the evolutionary history of related beta-CoV

`reference_genes` contains reference genes and proteins of SARS-CoV-2

`data/ctl` contains precomputed data on CTL epitope coverage


### Workflow for sequence analysis

1. Create a data analysis directory in the format YYYY-MM-DD in `data/fasta/`. Copy unaliged whole length genomes to be analyzed, as a FASTA file named `sequences` to that directory.

> **Important.** format your sequence names as `epi_isl_XXXXX|anything` where `XXXXX` is a unique numerical ID (could be any lengths). These are also the IDs that will be used to match sequences with their metadata.

##### To test, copy the file `sequences.fasta` from `test/` to `data/fasta/test/2021-01-01/sequences` (no extension)

2. Each gene (or ORF1a/b product like `nsp2`) is analyzed separately. The [reference genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512) partitoned into these segments is found in the [reference_genes](https://github.com/veg/SARS-CoV-2/tree/compact/reference_genes) directory.

3. The top-level dispatcher `bash` script is `scripts/submit_jobs.sh` is invoked from the root directory of the pipeline (e.g. to use OpenMPI for indivdual jobs but not to use a scheduler, and to use 4 processors for smaller jobs and 8 processors for larger jobs)

```bash scripts/submit_jobs.sh -s4 -L8 -p data/fasta/YYYY-MM-DD```

Optional arguments can be viewed by calling

```
$bash scripts/submit_jobs.sh -h

Top level job dispatcher for the SARS-CoV-2 selection analysis pipeline

Syntax: ./submit_jobs [-hm] [-b dir] [-l dir] [-s number] [-L number] [-q queue] analysis_directory 
options:
h     Print this Help.
m     Enable MPI mode and use scheduler (requires mpirun and qsub; otherwise run serially).
p     Enable MPI mode but do not use scheduler (requires mpirun)
b     path to the root directory of the pipeline (default = /Users/sergei/Dropbox/COVID-19)
l     path to the directory for MPI job logs (default /Users/sergei/Dropbox/COVID-19/logs/analysis_dir_name)
s     allocate this many CPUs to 'small' jobs (default 16)
L     allocate this many CPUs to 'large' jobs (default 48)
q     MPI-queue send jobs to this MPI queue (default : no queue specification)
positional arguments:
[req] analysis_directory path to the directory with the sequence file

```

For each individual gene/product (e.g. `S`, `3C`, `nsp3`), `submit_jobs.sh` invokes a wrapper `bash` script `scripts/extract_genes.sh`.


> **Important** When run with the `-m` mode, this script will submit jobs to an MPI queue with `qsub`. Edit [the submitter function](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/submit_jobs.sh#L91) to change this behavior. The default behavior is to execute the jobs serially (with or without MPI, this is controlled by the `p` flag), so you may want to start the process from inside a `screen` session.

4. QA, align, compress, build trees, perform selection analyses on segments: `scripts/extract_genes.sh` 

`scripts/extract_genes.sh` is called by `scripts/extract_genes.sh` but call also be called directly 

```
QA, align, compress, build trees, perform selection analyses on segments

Syntax: ./extract_genes.sh [-h] data_directory working_directory gene number_cpu MP_flag
options:
h     Print this Help.
positional arguments:
data_directory     path to the directory with the data
base_driectory     path to the root directory of the pipeline 
gene               a SARS-CoV-2 gene/ORF/product (one of S,M,N,E,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF10,ORF1a,ORF1b,leader,nsp2,nsp3,nsp4,3C,nsp6,nsp7,nsp8,nsp9,nsp10,RdRp,helicase,exonuclease,endornase,methyltransferase)
number_cpu         allocate this many CPUs to a jobs
[opt] MP_flag      if set to MP, run the jobs WITHOUT MPI (otherwise run them with MPI)
```

> **Important** please ensure that the paths in the [appropriate block of the script](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L92) point to the right exectuables on your system.

##### To test, call 
```
bash scripts/extract_gene.sh `pwd`/data/fasta/tests/2021-01-18/ `pwd` S 4 MP
```

or (to use MPI via `mpirun`)

```
bash scripts/extract_gene.sh `pwd`/data/fasta/tests/2021-01-18/ `pwd` S 4 
```

The following steps are executed. The script will check for existence of non-empty result files and skip steps where these already exist. These files will be named `sequence.gene.xxx` and look like the following (e.g. for gene `S`)

```
data/fasta/test/2021-01-18
├── sequences
├── sequences.S.FEL.json
├── sequences.S.MEME.json
├── sequences.S.SLAC.json
├── sequences.S.cache
├── sequences.S.compressed.fas
├── sequences.S.compressed.filtered.fas
├── sequences.S.compressed.filtered.fas.rapidnj.bestTree
├── sequences.S.compressed.filtered.sto
├── sequences.S.duplicates-merged.json
├── sequences.S.duplicates.json
├── sequences.S.edits.json
├── sequences.S.filtered.json
├── sequences.S.map.json
├── sequences.S.msa
├── sequences.S.variants.csv
├── sequences.S.variants.json
├── sequences.S_nuc.compressed.fas
├── sequences.S_nuc.fas
├── sequences.S_nucleotide.duplicates.json
├── sequences.S_protein.compressed.fas
├── sequences.S_protein.duplicates.json
├── sequences.S_protein.fas
└── sequences.S_raw_nucleotide.duplicates.json
```

4.1 [Extract gene segment](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L159). 

Extract the appropriate gene segment from input genomes using codon-aware alignment, filter based on `N` content, translate to amino-acids. 

Inputs:

* genome file
* reference sequence (nucleotide)

Settings:

* padded coordinate-range where the gene is expected to be
* fraction of `N` to tolerate
* whether or not to remove sequences with premature stop codons


Outputs:

* `sequences.[gene].protein.fas` : sequences passing filter translated to amino-acids
* `sequences.[gene].nuc.fas` : sequences passing filter 

4.2 [Compress duplicates in protein and nucleotide space](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L174). 


Inputs:

* `sequences.[gene].protein.fas` : sequences passing filter translated to amino-acids
* `sequences.[gene].nuc.fas` : sequences passing filter 


Outputs:

* `sequences.[gene]_protein.compressed.fas`  : unique protein haplotypes
* `sequences.[gene].nuc.compressed.fas` :  unique nucleotide haplotypes
* `sequences.[gene]_protein.duplicates.json`  : keeps track of what other sequences are being represented by the retained unique sequence
* `sequences.[gene].nuc.duplicates.json` :  unique nucleotide haplotypes

All `duplicates.json` file have the following format

`retained sequence` : {"0" : copy 1`, "1": copy 2, ...}``

```
{
    "epi_isl_LR877727": {
        "0": "epi_isl_LR877727"
    },
    "epi_isl_LR877733": {
        "0": "epi_isl_LR877733",
        "1": "epi_isl_LR877737",
        "2": "epi_isl_LR877757",
        "3": "epi_isl_LR877776"
    },
...
```

4.3 [Align compressed protein sequences using MAFFT](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L185). 

Inputs:

* `sequences.[gene]_protein.compressed.fas` : unique protein haplotypes

Outputs:

* `sequences.[gene].msa` : a multiple sequence alignment of protein sequences

4.4 [Map codon sequences onto protein alignments](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L199). 

Also, replace all sequence letters that are not `ACGTN-` with `N`.

Inputs:

* `sequences.[gene].msa` : a multiple sequence alignment of protein sequences
* `sequences.[gene].nuc.fas` : sequences passing filter 


Outputs:

* `sequences.[gene].compressed.fas ` : a multiple sequence alignment of unique codon sequences
* `sequences.[gene].nucleotide.duplicates.json` : keeps track of what other sequences are being represented by the retained unique sequence

4.5 [Merge duplicate data](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L209). 

Inputs:

* `sequences.[gene]_raw_nucleotide.duplicates.json` 
* `sequences.[gene]_nucleotide.duplicates.json` 

Inputs/modified:

* `sequences.[gene].compressed.fas ` : a multiple sequence alignment of unique codon sequences
* `sequences.[gene].duplicates.json` : the master duplicate tracker file


4.6 [Perform additional error filtering](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L224). 

This step examines the nucelotide MSA, estimates the probability that singleton low frequency variants are errors, and replaces those with `-`. Low frequency variants that co-occur are not removed. Sequences that have too many low frequency variants are removed. 

Input:

* `sequences.[gene].compressed.fas ` : a multiple sequence alignment of unique codon sequences
* `sequences.[gene].duplicates.json` : the master duplicate tracker file

Output:

* `sequences.[gene].variants.csv` : a description of position by position nucleotide composition of the alignment

```
Site,Consensus,A,G,C,T,ambig,N,gap
...
2171,C,0,1250,0,42,0,0,0
...
```

* `sequences.[gene].variants.json` : a description of which variants occur in which sequence

```
...
epi_isl_MW451205_1_null_1":{
   "1500":5,
   "1708":5,
   "1760":1,
   "203":7,
   "2044":5,
   "2149":5,
   "2946":5,
   "3354":5,
   "copies":1
  },
...
```

* `sequences.[gene].edits.json` : masking operations by sequence

```
...
"epi_isl_MW451205_1_null_1":{
   "1760":"C"
  },
 "epi_isl_MW480858_1_null_1":"removed",
 "epi_isl_MW480894_1_null_3":{
  },
...
```

* `sequences.[gene].compressed.filtered.fas ` : a multiple sequence alignment of unique codon sequences with putative errors masked and sequences that have too many rare variants removed

This stage also compresses `sequences.[gene].duplicates.json``

4.7 [Perform neighbor joining inference](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L240). 

Inputs:

* `sequences.[gene].compressed.filtered.fas ` : a multiple sequence alignment of unique codon sequences with putative errors masked and sequences that have too many rare variants removed

Outputs:

* `sequences.[gene].compressed.filtered.rapidnj.bestTree` : a Newick tree string

4.8 [Run the SLAC analysis (all branches)](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L256). 

Substitution mapping, ancestral state inference, dN and dS counting.

Inputs:

* `sequences.[gene].compressed.filtered.fas` : a multiple sequence alignment of unique codon sequences with putative errors masked and sequences that have too many rare variants removed

* `sequences.[gene].compressed.filtered.rapidnj.bestTree` : a Newick tree string


Outputs:

* `sequences.[gene].SLAC.json.gz` : compressed detailed SLAC analysis results

4.9 [Run the FEL analysis (internal branches)](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L267). 

Find sites subject to pervasive positive/negative selection along internal tree branches.

Inputs:

* `sequences.[gene].compressed.filtered.fas` : a multiple sequence alignment of unique codon sequences with putative errors masked and sequences that have too many rare variants removed

* `sequences.[gene].compressed.filtered.rapidnj.bestTree` : a Newick tree string


Outputs:

* `sequences.[gene].FEL.json.gz` : compressed detailed FEL analysis results

4.10 [Run the MEME analysis (internal branches)](https://github.com/veg/SARS-CoV-2/blob/compact/scripts/extract_genes.sh#L278). 

Find sites subject to pervasive positive/negative selection along internal tree branches.

Inputs:

* `sequences.[gene].compressed.filtered.fas` : a multiple sequence alignment of unique codon sequences with putative errors masked and sequences that have too many rare variants removed

* `sequences.[gene].compressed.filtered.rapidnj.bestTree` : a Newick tree string


Outputs:

* `sequences.[gene].MEME.json.gz` : compressed detailed MEME analysis results


### Workflow for processing and annotation

To create a full report, the pipeline requires to additional pieces of information about each sequence, its collection date and collection location. Some, but not all can be set to null. The format of the file is 

```
...
{
"epi_isl_MW505982": {
    "collected": "2020-09-24",
    "location": {
        "country": "France",
        "locality": null,
        "state": null,
        "subregion": null
    }
...
```

See an example in `test/annotation.json`. There is also a simple convenience script in `python/meta-to-json.py` to convert a CSV file like `test/sequences.csv` to this `JSON format`

To process a directory of results (like `data/fasta/test/2021-01-18`) call a shell script (from the base directory of the pipeline) you need to place the `annotation.json` (named exactly like that) in the data directory (e.g. `data/fasta/test/2021-01-18/`) and call

```
bash scripts/summarize_genes.sh `pwd`/data/fasta/test/2021-01-18/sequences
```

## Results

Depending on the size of the data, this may take a while to process. The end result will be several files that can be either examined directly, or via attendant visualizations. 

### Site by site report (comparative-annotation.json)

The `comparative-annotation.json` file is placed in the data directory (e.g., `data/fasta/test/2021-01-18/`), and contains, for each analyzed site in the alignment, the following information (comments in-line).

Precomputed results below are copied from `data/evo-annotation.json` and **do not depend on the data analyzed by the scripts**

```
"23401": {  | genomic coordinate for the start of the codon (nuc)
  "G": "S", | gene or segment 
  "S": 614, | site in gene or segment
  "bCFEL": { | precomputed results for differences in selective pressure between nCoV lineages and bats (contrast-FEL)
   "p": 0.9996469714423809, | p-value
   "a": 0.2551011391587314, | syn rate
   "b-nCOV": 0, | non-syn rate for nCOV
   "b": 0 | non-syn rate for bats/pangolins
  },
  "bFEL": {  precomputed results for selection in the nCOV clade (FEL)
   "a": 0.255037454036329, | syn rate
   "b": 0, | non-syn rate
   "p": 0.1433611348116214 | p-value
  },
  "bMEME": { precomputed results for epdosic selection in the nCOV clade (MEME)
   "p": 0.6666666666666666, | p-value
   "a": 0.2550391508973095, | syn rate
   "b+": 0.01635397993808615, | non-syn rate unrestricted
   "w+": 1.110223024625157e-16, | weight
   "b-": 0, | non-syn rate restricted to be ≤ syn.rate
   "w-": 0.9999999999999999, | weight
   "br": 0 | branches under episodic selection
  },
  "bSC2": "GAT", | precomputedc odon in the reference SARS-CoV-2 sequence
  "bSC2-aa": "D", | precomputed aa in the reference SARS-CoV-2 sequence
  "bcdn": { | precomputed codon composition of the site
   "nCOV": { | nCoV clade
    "GAT": 11
   },
   "others": { | others
    "GAT": 57,
    "GAC": 1
   }
  },
  "baa": {  | precomputed a.a composition of the site
   "nCOV": {
    "D": 11
   },
   "others": {
    "D": 58
   }
  },
  "evo": { | precomputed evolutionary predictions for what will be obsevred in SARS-CoV-2 at this site
   "GAT": 0.9873479741958492,
   "GAC": 0.01265202580415079
  },
  "cdn": { | codon composition of this site in the analyzed data
   "GGT": 914,
   "GAT": 298
  },
  "aa": { | a.a composition of this site in the analyzed data
   "G": 914,
   "D": 298
  },
  "SLAC": { | SLAC report
   "N": 11, | # of non-syn changes
   "S": 0,  | # of syn changes
   "EN": 2.053278866216693, | non-syn sites
   "ES": 0.9467211337833055, | syn sites
   "p": 0.315573711261102, | p-value for dN/dS ≠ 1
   "NmS": 260.2139993512735 | scaled dN - dS
  },
  "FEL": { | FEL report
   "a": 0.8958182295483623, | syn rate
   "b": 114.4993776292876,  | non-syn rate
   "p": 0.1051838515962932  | p-value
  },
  "MEME": { | MEME report
   "p": 0.1224233017791853,  | p-value
   "a": 0.2505266788143614,  | syn rate
   "b+": 259.5369555266688,  | non-syn rate unrestricted
   "w+": 0.4947082115422307, | weight
   "b-": 0.1290370800966509, | non-syn rate restricted to be ≤ syn.rate
   "w-": 0.5052917884577693, | weight
   "br": 2 | branches under episodic selection
  },
  "trend": 7.38492838742239, | temportal trend Z-score
  "subs": { | inferred substitutions 
   "cdn": { | internal branches, codons, from|to : count
    "GAT|GGT": 2
   },
   "aa": { | internal branches, a.a., from|to : count
    "D|G": 2
   },
   "lcdn": { | terminal branches, codons, from|to : count
    "GAT|GGT": 8,
    "GGT|GAT": 1
   },
   "laa": { | terminal branches, a.a., from|to : count
    "D|G": 8,
    "G|D": 1
   }
  }
 }
```
### Files for visualization 

Placed in the analysis directory (e.g. `data/fasta/test/2021-01-18/`)

Compressed report for the Observable notebook
`report.json`

Placed in the parent of the analysis directory (e.g. `data/fasta/test/`)

`temporal-gene-properties.csv`
`temporal-gene-sites.csv`

#### How to use visualization

Make the root directory for the analyses (e.g. `data/fasta/test/`) accessible via `http` (must allow CORS [cross-origin requests]). Visualization consumes: 


From `data/fasta/test/YYYY-MM-DD`

    `report.json`, `sequence-annotation.json`, `sequence-information.json`

From `data/fasta/test/`

    `temporal-gene-properties.csv`, `temporal-gene-sites.csv`
    
    
To view results, view 

`https://observablehq.com/@spond/revised-sars-cov-2-analytics-page?base=baseURL&dir=dataURL&time=&title=test`

Where `baseURL` points at `data/fasta/test/` and `dir` is the `YYYY-MM-DD` date

For example, to view the results of the test data set run, use 

https://observablehq.com/@spond/revised-sars-cov-2-analytics-page?base=https%3A%2F%2Fraw.githubusercontent.com%2Fveg%2FSARS-CoV-2%2Fcompact%2Fdata%2Ffasta%2Ftest&dir=2021-01-18&title=Test&time=









