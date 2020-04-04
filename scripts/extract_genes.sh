FILE=$1
NP=10
HYPHY=/usr/local/bin/hyphy
HYPHYMPI=/usr/local/bin/HYPHYMPI
MAFFT=/usr/local/bin/mafft
RAXML=/usr/local/bin/raxml-ng
TN93=/usr/local/bin/tn93

function RunAGene {



GENE=$1
REFERENCE_SEQUENCE=$2
TRIM_FROM=$3
TRIM_TO=$4
N_FRAC=$5

if [ -s ${FILE}.${GENE}.1.json ]
then 
   echo "$GENE alignment already processed"
else
    echo "EXTRACTING $GENE"
    if [ -s ${FILE}.${GENE}_protein.fas ] 
    then
        echo "Already extracted"
    else
        mpirun -np $NP --use-hwthread-cpus $HYPHYMPI  /Users/sergei/Development/hyphy-analyses/codon-msa/pre-msa.bf --input $FILE --reference $REFERENCE_SEQUENCE --trim-from $TRIM_FROM --trim-to $TRIM_TO --N-fraction $N_FRAC
        mv ${FILE}_protein.fas ${FILE}.${GENE}_protein.fas
        mv ${FILE}_nuc.fas ${FILE}.${GENE}_nuc.fas
    fi
    
    if [ -s ${FILE}.${GENE}.msa ] 
    then
        echo "Already aligned"
    else
        $MAFFT ${FILE}.${GENE}_protein.fas > ${FILE}.${GENE}.msa
    fi
        
    if [ -s ${FILE}.${GENE}.duplicates.json ] 
    then
        echo "Already reverse translated"
    else
        $HYPHY /Users/sergei/Development/hyphy-analyses/codon-msa/post-msa.bf --protein-msa ${FILE}.${GENE}.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.fas --output ${FILE}.${GENE}.compressed.fas --duplicates ${FILE}.${GENE}.duplicates.json
        $HYPHY /Users/sergei/Development/hyphy-analyses/codon-msa/post-msa.bf --protein-msa ${FILE}.${GENE}.msa --nucleotide-sequences ${FILE}.${GENE}_nuc.fas --compress No --output ${FILE}.${GENE}.all.fas    
    fi
    
     
    if [ -s ${FILE}.${GENE}.withref.fas ]
    then 
        echo "Already has alignment with reference"
    else
        $MAFFT --add $REFERENCE_SEQUENCE --reorder ${FILE}.${GENE}.all.fas > ${FILE}.${GENE}.withref.fas
    fi 
    
    if [ -s ${FILE}.${GENE}.tn93 ] 
    then
        echo "Already computed TN93"
    else
        $TN93 -q -t 0.05 ${FILE}.${GENE}.withref.fas > ${FILE}.${GENE}.tn93 2> ${FILE}.${GENE}.tn93.json
        python3 python/tabulate-diversity-divergence.py -j data/db/master.json -t ${FILE}.${GENE}.tn93 > data/evolution.${GENE}.csv
    fi
   
    if [ -s ${FILE}.${GENE}.compressed.fas.raxml.bestTree ] 
    then
        echo "Already has tree"
    else
        $RAXML --msa ${FILE}.${GENE}.compressed.fas --model GTR+G --force
    fi
    
    if [ -s ${FILE}.${GENE}.SLAC.json ] 
    then
        echo "Already has SLAC results"
    else
        mpirun --use-hwthread-cpus -np $NP $HYPHYMPI slac --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.SLAC.json
    fi

    if [ -s ${FILE}.${GENE}.FEL.json ] 
    then
        echo "Already has FEL results"
    else
        mpirun --use-hwthread-cpus -np $NP $HYPHYMPI fel --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.FEL.json
    fi

    if [ -s ${FILE}.${GENE}.MEME.json ] 
    then
        echo "Already has MEME results"
    else
        mpirun --use-hwthread-cpus -np $NP $HYPHYMPI meme --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.MEME.json
    fi
    
    if [ -s ${FILE}.${GENE}.PRIME.json ] 
    then
        echo "Already has PRIME results"
    else
        mpirun --use-hwthread-cpus -np $NP $HYPHYMPI prime --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.PRIME.json
    fi
    
    python3 python/summarize-gene.py -D data/db/master-no-fasta.json -d ${FILE}.${GENE}.duplicates.json -s ${FILE}.${GENE}.SLAC.json -f ${FILE}.${GENE}.FEL.json -m ${FILE}.${GENE}.MEME.json -P 0.1 -p ${FILE}.${GENE}.PRIME.json --output  ${FILE}.${GENE}.json -c ${FILE}.${GENE}.withref.fas

    if [ -s ${FILE}.${GENE}.BGM.json ] 
    then
        echo "Already has BGM results"
    else
       $HYPHY bgm --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --min-subs 2 --type codon
       mv ${FILE}.${GENE}.compressed.fas.BGM.json ${FILE}.${GENE}.BGM.json
    fi

    if [ -s ${FILE}.${GENE}.FADE.json ] 
    then
        echo "Already has FADE results"
    else
        $HYPHY scripts/reroot-on-oldest.bf --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --csv data/attributes.csv --output ${FILE}.${GENE}.compressed.fas.rooted
        $HYPHY conv Universal "Keep Deletions" ${FILE}.${GENE}.compressed.fas  ${FILE}.${GENE}.compressed.fas.prot
        $HYPHY fade --alignment ${FILE}.${GENE}.compressed.fas.prot --tree ${FILE}.${GENE}.compressed.fas.rooted 
    fi


    #if [ -s ${FILE}.${GENE}.BUSTED.json ] 
    #then
    #    echo "Already has BUSTED results"
    #else
    #    $HYPHY busted --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.BUSTED.json --rates 2 --syn-rates 2 --starting-points 10
    #fi


    #if [ -s ${FILE}.${GENE}.ABSREL.json ] 
    #then
    #    echo "Already has ABSREL results"
    #else
    #    mpirun -np $NP $HYPHYMPI absrel --alignment ${FILE}.${GENE}.compressed.fas --tree ${FILE}.${GENE}.compressed.fas.raxml.bestTree --branches Internal --output ${FILE}.${GENE}.ABSREL.json
    #fi
    
    
fi

}


RunAGene "S" "data/reference_genes/S.fas" "20000" "27000" 0.005
RunAGene "ORF1a" "data/reference_genes/ORF1a.fas" "1" "15000" 0.001
RunAGene "M" "data/reference_genes/M.fas" "25000" "30000" 0.01
RunAGene "N" "data/reference_genes/N.fas" "26000" "35000" 0.01
RunAGene "ORF3a" "data/reference_genes/ORF3a.fas" "24000" "27000" 0.05
RunAGene "ORF6" "data/reference_genes/ORF6.fas" "26000" "30000" 0.05
RunAGene "ORF7a" "data/reference_genes/ORF7a.fas" "26000" "35000" 0.05
RunAGene "ORF8" "data/reference_genes/ORF8.fas" "26000" "35000" 0.05
RunAGene "ORF1b" "data/reference_genes/ORF1b.fas" "12000" "24000" 0.001


