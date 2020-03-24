Workflow

1). Scrape sequences from GISAID into `json` files

2). Import into `data/db/master.json` using `python/import-batch.py`. Command used 

``` 
python3 python/import-batch.py -j data/db/master.json -f \
~/Downloads/gisaid.json  --update -L data/db/README.md
```

Extract only attribute data

```
jq 'del(.[].sequence)' data/db/master.json > data/db/master-no-fasta.json
```

Extract attributes

``` 
python3 python/extract-attributes.py -j data/db/master.json > data/attributes.csv
```

3). Extract daily snapshot all human sequeneces that are full genomes (`YYYYMMDD`)

```
python3 python/extract-sequences.py -j data/db/master.json -f \
"host" "re" "[hH]uman" -f "sequence" ">" 28000 > data/fasta/YYYYMMDD/sequencs
```

4). Run the processor shell script 

```
sh scripts/extract_genes.sh data/fasta/YYYYMMDD/sequences
```