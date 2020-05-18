FILE=$1
MAFFT=/usr/local/bin/mafft
rm data/mafs.csv
rm data/evo_freqs.csv

#gene_coordinates = [[265,13482, 'ORF1a', 0],
#  [13467,21554, 'ORF1b',-1],
#  [21562,25383,'S',0],
#  [25392,26219,'ORF3a',0],
#  [26244,26471,'E',0],
#  [26522,27190,'M',0],
#  [27201,27386,'ORF6',0],
#  [27393,27758,'ORF7a',0],
#  [27893,28258,'ORF8',0],
#  [28273,29532,'N',0],
#  [29557,29673,'ORF10',0]
#]

genes=(S M N ORF3a ORF6 ORF7a ORF8 ORF1a ORF1b)
# ORF1a ORF1b ORF3a ORF6 ORF7a ORF8)
offsets=(21562 26522 28273 25392 27201 27393 27893 265 13471)
#genes=(ORF3a)
#offsets=(28273) 

#for GENE in {S,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF8}; do

ANNOTATION=${FILE}.annotation.json

cp data/comparative-annotation-between.json ${FILE}.annotation.json

for i in ${!genes[@]}; do
    GENE=${genes[i]}
    OFFSET=${offsets[i]}
    echo $GENE
    if [ -s ${FILE}.${GENE}.withref.fas ]
    then 
        echo "Already has alignment with reference"
    else
        echo ${FILE}.${GENE}.compressed.fas
        $MAFFT --add data/reference_genes/${GENE}.fas --reorder ${FILE}.${GENE}.compressed.fas > ${FILE}.${GENE}.withref.fas
        cp ${FILE}.${GENE}.withref.fas ${FILE}.${GENE}.bkup.withref.fas
    fi 
    python3 python/summarize-gene.py -T data/ctl/epitopes.json -D data/db/master-no-fasta.json -d ${FILE}.${GENE}.duplicates.json -u ${FILE}.${GENE}.compressed.fas.FUBAR.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.withref.fas -E data/evo_annotation.json -F $GENE -A data/mafs.csv -V data/evo_freqs.csv -S $OFFSET -O $ANNOTATION > ${FILE}.${GENE}.json
done;

cp $ANNOTATION data/comparative-annotation.json