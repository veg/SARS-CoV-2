rm data/evo_annotation.json

python3 python/extract-prime-annotation.py  -c data/betas/sequences.nsp2.compressed.fas -p data/betas/sequences.nsp2.compressed.fas.PRIME.json -r ORF1a  -f 180 -o data/evo_annotation.json
python3 python/extract-prime-annotation.py  -c data/betas/sequences.3C.compressed.fas -p data/betas/sequences.3C.compressed.fas.PRIME.json -r ORF1a  -f 3263 -o data/evo_annotation.json 
python3 python/extract-prime-annotation.py  -c data/betas/sequences.leader.compressed.fas -p data/betas/sequences.leader.compressed.fas.PRIME.json -r ORF1a  -f 0 -o data/evo_annotation.json
python3 python/extract-prime-annotation.py  -c data/betas/sequences.nsp4.compressed.fas -p data/betas/sequences.nsp4.compressed.fas.PRIME.json -r ORF1a  -f 2763 -o data/evo_annotation.json 
python3 python/extract-prime-annotation.py  -c data/betas/sequences.nsp3.compressed.fas -p data/betas/sequences.nsp3.compressed.fas.PRIME.json -r ORF1a  -f 818 -o data/evo_annotation.json 
python3 python/extract-prime-annotation.py  -c data/betas/nsp6_nuc.fas -p data/betas/nsp6.PRIME.json -r ORF1a  -f 3569 -o data/evo_annotation.json 
python3 python/extract-prime-annotation.py  -c data/betas/sequences.nsp7.compressed.fas -p  data/betas/sequences.nsp7.compressed.fas.PRIME.json -r ORF1a  -f 3859 -o data/evo_annotation.json 
python3 python/extract-prime-annotation.py  -c data/betas/sequences.nsp8.compressed.fas -p data/betas/sequences.nsp8.compressed.fas.PRIME.json -r ORF1a  -f 3942 -o data/evo_annotation.json 
python3 python/extract-prime-annotation.py  -c data/betas/sequences.nsp9.compressed.fas -p data/betas/sequences.nsp9.compressed.fas.PRIME.json -r ORF1a  -f 4140 -o data/evo_annotation.json 
python3 python/extract-prime-annotation.py  -c data/betas/sequences.nsp10.compressed.fas -p data/betas/sequences.nsp10.compressed.fas.PRIME.json -r ORF1a  -f 4253 -o data/evo_annotation.json 

python3 python/extract-prime-annotation.py  -c data/betas/sequences.RdRp.compressed.fas -p data/betas/sequences.RdRp.compressed.fas.PRIME.json -r ORF1b  -f -10 -o data/evo_annotation.json 
python3 python/extract-prime-annotation.py  -c data/betas/sequences.helicase.compressed.fas -p data/betas/sequences.helicase.compressed.fas.PRIME.json -r ORF1b  -f 922 -o data/evo_annotation.json 
python3 python/extract-prime-annotation.py  -c data/betas/sequences.endornase.compressed.fas -p data/betas/sequences.endornase.compressed.fas.PRIME.json -r ORF1b  -f 2050 -o data/evo_annotation.json  
python3 python/extract-prime-annotation.py  -c data/betas/sequences.methyltransferase.compressed.fas -p data/betas/sequences.methyltransferase.compressed.fas.PRIME.json -r ORF1b  -f 2396 -o data/evo_annotation.json 

python3 python/extract-prime-annotation.py  -c data/betas/S.aligned.fas -p data/betas/S.aligned.fas.best-gard.PRIME.json -r S  -f 0 -o data/evo_annotation.json
#  
python3 python/extract-prime-annotation.py  -c data/betas/sequences.ORF6.compressed.fas -p data/betas/sequences.ORF6.compressed.fas.PRIME.json -r ORF6 -f 0 -o data/evo_annotation.json
python3 python/extract-prime-annotation.py  -c data/betas/sequences.ORF8.compressed.fas -p  data/betas/sequences.ORF8.compressed.fas.PRIME.json -r ORF8 -f 0 -o data/evo_annotation.json
python3 python/extract-prime-annotation.py  -c data/betas/sequences.ORF7a.compressed.fas -p data/betas/sequences.ORF7a.compressed.fas.PRIME.json -r ORF7a -f 0 -o data/evo_annotation.json
python3 python/extract-prime-annotation.py  -c data/betas/sequences.ORF3a.compressed.fas -p data/betas/sequences.ORF3a.compressed.fas.PRIME.json -r ORF3a -f 0 -o data/evo_annotation.json

python3 python/extract-prime-annotation.py  -c data/betas/sequences.N.compressed.fas -p data/betas/sequences.N.compressed.fas.PRIME.json -r N -f 0 -o data/evo_annotation.json
python3 python/extract-prime-annotation.py  -c data/betas/sequences.M.compressed.fas -p data/betas/sequences.M.compressed.fas.PRIME.json -r M -f 0 -o data/evo_annotation.json