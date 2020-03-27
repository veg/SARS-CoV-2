### Daily analyses of SARS-CoV-2 genomic data

> This project is a part of a larger effort with the Galaxy team: [covid19.galaxyproject.org](http://covid19.galaxyproject.org])

#### TL;DR

1. [Analysis of all current SARS-CoV-2 genomes for evidence of natural selection](https://observablehq.com/@spond/natural-selection-analysis-of-sars-cov-2-covid-19)
2. [Divergence and diversity of SARS-CoV-2 genomes over time overall and by region](https://observablehq.com/@spond/current-state-of-sars-cov-2-evolution)


#### Analysis pipeline

1. We collect data from the <img src="https://www.gisaid.org/fileadmin/gisaid/img/schild.png" alt="gisaid-logo" width="65"> database daily. These are mostly full genome sequences collected from different platforms and different regions. See here for a [summary of the sequence data](https://observablehq.com/@stevenweaver/case-vs-sequence-count). The metadata on the sequences can be found in `data/db/master-no-fasta.json`; in accordance with GISAID data usage policies, we do not distribute sequence data here.

2. We extract full genome human sequences and map them to the [reference genes](https://www.ncbi.nlm.nih.gov/nuccore/?term=COVID) using a simple [codon-aware pipeline](https://github.com/veg/hyphy-analyses/tree/master/codon-msa). 
At this step we also compress the data to retain a single copy of each unique haplotype in the gene, and filter out sequences that have too many (>0.5%) uncalled/unresolved (`N`) bases.

3. We reconstruct ML phylogenies on compressed data using [raxml-ng](https://github.com/amkozlov/raxml-ng)

4. We estimate gene-by-gene distances to compute diversity and divergence using [TN93](https://github.com/veg/tn93), summarized [here](https://observablehq.com/@spond/current-state-of-sars-cov-2-evolution)

5. We run several [HyPhy](https://github.com/veg/hyphy) dN/dS based selection analyses on each gene. We restrict these analyses to internal branches of the tree [filter within-host evolution](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0020062).	
>When analyzing intra-species or intra-host data, dN/dS estimates may be inflated due to the fact that not all observed sequence variation is due to substitutions, but some are simply mutations that have not yet been filtered by selection. In other words, dN/dS may be elevated by intra-species / intra- host polymorphism that need not be attributable by positive selection. One simple approach to mitigating this undesirable effect is to restrict site-specific analyses to Internal branches only. This is because internal branches encompass at least one step that is visible to selection (transmission and/or multiple rounds of replication), and are less likely to contain spurious polymorphic variants.

6. These analyses include [SLAC and FEL](https://www.ncbi.nlm.nih.gov/pubmed/15703242), [MEME](https://www.ncbi.nlm.nih.gov/pubmed/22807683), and PRIME (the latter allows to test for conservation/change in specific biochemical properties at site) to identify which sites may be experiencing positve selection, and what properties may be important to preserve/change during these changes. The up-to-date summary is hosted [here](https://observablehq.com/@spond/natural-selection-analysis-of-sars-cov-2-covid-19)
