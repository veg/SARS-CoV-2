FILE=$1
MAFFT=/usr/local/bin/mafft
rm data/mafs.csv

for GENE in {S,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF8}; do
    echo $GENE
    if [ -s ${FILE}.${GENE}.withref.fas ]
    then 
        echo "Already has alignment with reference"
    else
        echo ${FILE}.${GENE}.compressed.fas
        $MAFFT --add data/reference_genes/${GENE}.fas --reorder ${FILE}.${GENE}.compressed.fas > ${FILE}.${GENE}.withref.fas
        cp ${FILE}.${GENE}.withref.fas ${FILE}.${GENE}.bkup.withref.fas
    fi 
    python3 python/summarize-gene.py -D data/db/master-no-fasta.json -d ${FILE}.${GENE}.duplicates.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.withref.fas -E data/evo_annotation.json -F $GENE -A data/mafs.csv > ${FILE}.${GENE}.json
done;