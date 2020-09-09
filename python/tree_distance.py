import dendropy
from dendropy.calculate import treecompare

file1 = './data/fasta/2020-05-09/sequences.leader.compressed.fas.rapidnj.bestTree'
file2 = './data/fasta/2020-05-10/sequences.leader.compressed.fas.rapidnj.bestTree'
file3 = './data/fasta/2020-05-11/sequences.leader.compressed.fas.rapidnj.bestTree'
file4 = './data/fasta/2020-05-12/sequences.leader.compressed.fas.rapidnj.bestTree'

s1 = open(file1).read()
s2 = open(file2).read()
s3 = open(file3).read()
s4 = open(file4).read()

tns = dendropy.TaxonNamespace()

tree1 = dendropy.Tree.get(
        data=s1,
        schema='newick',
        taxon_namespace=tns)
tree2 = dendropy.Tree.get(
        data=s2,
        schema='newick',
        taxon_namespace=tns)

tree3 = dendropy.Tree.get(
        data=s3,
        schema='newick',
        taxon_namespace=tns)

tree4 = dendropy.Tree.get(
        data=s4,
        schema='newick',
        taxon_namespace=tns)


## Euclidean distance = 2.22326363775
print(treecompare.euclidean_distance(tree1, tree2))
print(treecompare.euclidean_distance(tree1, tree3))
print(treecompare.euclidean_distance(tree2, tree3))
print(treecompare.euclidean_distance(tree3, tree4))


print(treecompare.false_positives_and_negatives(tree1, tree2))
