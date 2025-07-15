# Utility python script for the <godon_stats.vHive.sh> bash script
# Mark branches under positive selection on the tree file
# Threshold to consider positive selection defined in qval_threshold

from ete3 import Tree
import os
import csv
import sys


R_results_csv = sys.argv[1]
godon_out = sys.argv[2]
raxmltree = sys.argv[3]
outtree = sys.argv[4]
treename = outtree.split("/")[-1]
ps = 0
sf = str(godon_out).split('.')[1]


with open(raxmltree, 'r') as raxtree:
    raxml_tree = raxtree.read()
out_tree = Tree(raxml_tree, format=1)

dico_branch_num = {}

with open(godon_out, 'r') as godon_file:
    for line_godon in godon_file:
        if line_godon.startswith("Testing branch "):
            branch_num = line_godon.split(" ")[2]
            branch_num = branch_num.split("/")[0]
            gtree = next(godon_file)[19:]
            godon_tree = Tree(gtree, format=1)
            
            count = 1
            for branch in godon_tree.traverse():
                if "#1" in branch.name:
                    dico_branch_num[count] = branch_num
                count += 1

count = 1
for branch in out_tree.traverse():
    if count in dico_branch_num:
        branch.name = branch.name + "_" + sf + "_branch" + dico_branch_num[count]
    count += 1
out_tree.write(format=1, outfile=outtree)
