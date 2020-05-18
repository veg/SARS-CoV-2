IMPORT=$1
DATE=$2

python3 python/import-batch.py -j data/db/master.json -f $IMPORT --update -L data/db/README.md
jq 'del(.[].sequence)' data/db/master.json > data/db/master-no-fasta.json
mkdir data/fasta/${DATE}
python3 python/extract-sequences.py -j data/db/master.json -f "host" "re" "[hH]uman" -f "sequence" ">" 28000 > data/fasta/${DATE}/sequences

