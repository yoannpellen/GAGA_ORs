# Utility python script for the <generax.sh> bash script
# Various functions to prepare files for Generax and GSEA (actually used for hypergeometric test now)

from Bio import SeqIO
from ete3 import Tree
import argparse
import csv
import os

# Description for script whith '-h'
parser = argparse.ArgumentParser(description = "Different functions needed to deal with Guidance2")
optional_args = parser.add_argument_group(title = "Arguments (depending on the function called)")

# Define arguments with help message
optional_args.add_argument('-f', dest = 'function', help = "Function to call")
optional_args.add_argument('-it', dest = 'intree', help = "Input tree")
optional_args.add_argument('-ot', dest = 'outtree', help = "Output tree")
optional_args.add_argument('-br', dest = 'branch', help = "ID number for the branch of interest in the species tree")
optional_args.add_argument('-xml', dest = 'xml', help = "Path to XML file by Generax")
optional_args.add_argument('-map', dest = 'mapping', help = "Name of mapping file")
optional_args.add_argument('-imsa', dest = 'inmsa', help = "Input MSA file")
optional_args.add_argument('-omsa', dest = 'outmsa', help = "Output MSA file")
optional_args.add_argument('-csv', dest = 'Rcsv', help = "CSV file")

# Group all arguments in a list
args = parser.parse_args()

def treelabel(intree, outtree):
    # Label every internal branch of the tree with its number
    with open(intree, 'r') as in_tree:
        tree = in_tree.read()
    out_tree = Tree(tree, format=1)
    branch_num = 1
    for branch in out_tree.traverse():
        if not branch.is_leaf():
            branch.name = "spbranch"+str(branch_num)
        branch_num += 1
    out_tree.write(format=1, outfile=outtree)

def mapping(intree, mapping):
    # Create mapping file for generax using only the gene tree file, getting the species name from the leaves name
    with open(intree, 'r') as intree:
        intree = intree.read()
    tree = Tree(intree, format=1)
    for leaf in tree.iter_leaves():
        GAGAname = leaf.name.split('_')[2]
        with open(mapping, 'a') as mapfile:
            mapfile.write(GAGAname+':'+leaf.name+'\n')

def msalabel(inmsa, outmsa, mapping):
    # Label the sequences name in the alignment based on mapping file
    with open(mapping, "r") as mapfile:
        for line in mapfile:
            newname = line.split(":")[1][:-1]
            oldname = newname.split('_')
            oldname = '_'.join(oldname[:4])
            with open(inmsa, 'r') as in_msa:
                for record in SeqIO.parse(in_msa, "fasta"):
                    if record.name == oldname:
                        with open(outmsa, 'a') as out_msa:
                            out_msa.write(">" + newname + "\n")
                            out_msa.write(str(record.seq) + "\n")
                

def search(branch, Rcsv):
    clade = os.getcwd()
    clade = os.path.basename(clade)[8:]
    # Get gene tree branches id matching a specific branch of the species tree
    search_branch = rf'speciesLocation="spbranch{branch}"/>'
    in_rectree = False
    g_branches = []
    branch_list = []
    for file in os.listdir():
        xml = os.fsdecode(file)
        if xml.endswith(".xml"):
            sf = str(xml).split('_')[1]
            with open(xml, 'r') as xmlfile:
                for line in xmlfile:
                    if "<recGeneTree>" in line:
                        in_rectree = True
                    if in_rectree:
                        if "<name>"+sf in line:
                            g_branch = line.split('<name>')[1]
                            g_branch = g_branch.split('<')[0]
                            sp_branch = next(xmlfile)
                            sp_branch = next(xmlfile)
                            if search_branch in sp_branch: # define list of branches to focus on for gsea
                                branch_list.append(g_branch[:-1])
                        
    # build gmt and rnk file for gsea by going through csv and check if branch is in list
    # seperate in 2 lines: branch of interest (from list) and rest of the tree
    first_list = ["foreground_branch","foreground"]
    sec_list = ["background_branches","background"]
    gmt = clade + ".branch" + branch + ".gmt"
    rnk = clade + ".branch" + branch + ".rnk"
    open(rnk, 'w').close()
    open(gmt, 'w').close()
    with open(Rcsv, 'r') as csv_file:
        csv_data = csv.reader(csv_file)
        header = next(csv_data)
        for row in csv_data:
            sf = row[0]
            branch = row[1]
            pval = float(row[5])
            invert_pval = 1 - pval
            
            with open(rnk, 'a') as rnkf:
                rnkf.write(str(sf)+"_"+str(branch)+"\t"+str(pval)+"\n")
                
            search_pattern = sf + "_gbranch" + branch
            if search_pattern in branch_list:
                first_list.append(str(sf)+"_"+str(branch))
            else:
                sec_list.append(str(sf)+"_"+str(branch))
            
    data_list = [first_list,sec_list]
    with open(gmt, 'w', newline='') as gmtf:
        writer = csv.writer(gmtf, delimiter='\t')
        writer.writerows(data_list)
       

if __name__ == "__main__":
    if args.function == "mapping":
        mapping(args.intree, args.mapping)
    elif args.function == "search":
        search(args.branch, args.Rcsv)
    elif args.function == "treelabel":
        treelabel(args.intree, args.outtree)
    elif args.function == "msalabel":
        msalabel(args.inmsa, args.outmsa, args.mapping)
    else:
        print("ERROR: function '" + str(args.function) + "' doesn't exist")
        