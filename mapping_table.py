import os
import sys
from ete3 import Tree

with open("/mnt/c/Users/user/OneDrive/Bureau/Thesis/workspace/GAGA/GAGAsp.tree", 'r') as in_tree:
    tree = in_tree.read()
out_tree = Tree(tree, format=1)
dic_branch = {}
for branch_tree in out_tree.traverse():
    if "spbranch" in branch_tree.name:
        speciesbr = branch_tree.name.split('spbranch')[1]
    else:
        speciesbr = branch_tree.name
    dic_branch[speciesbr] = str(branch_tree.dist)

clade = os.path.basename(os.getcwd())
generax_dir = os.getcwd() + "/generax_" + clade + "/"
mapping = generax_dir + clade + ".csv"
godon_stats = clade + ".godonBS.all.fdr.csv"
mapping_stats = generax_dir + clade + ".mapping_table.branch_length.csv"

g_branches = []
branch_list = []

with open(mapping, 'w') as mapf:
    mapf.write("Clade,Subfamily,Gene branch,Gene branch length,Species branch,Species branch length\n")
for file_xml in os.listdir(generax_dir):
    xml = os.fsdecode(file_xml)
    if xml.endswith(".xml"):
        in_rectree = False
        sf_xml = xml.split('_')[1]
        print(sf_xml)
        xml = generax_dir + xml
        with open(xml, 'r') as xmlfile:
            for line_xml in xmlfile:
                if "<recGeneTree>" in line_xml:
                    in_rectree = True
                if in_rectree:
                    if '<name>' in line_xml and '_' in line_xml:
                        sp_branch = next(xmlfile)
                        sp_branch = next(xmlfile)
                        if 'species_0' in sp_branch:
                            continue
                        if 'speciesLocation=' in sp_branch:
                            print(line_xml)
                            g_branch = line_xml.split('<name>')[1]
                            g_branch = g_branch.split('<')[0]
                            
                            generaxtree = generax_dir+rf"generax.{sf_xml}.tree"
                            with open(generaxtree, 'r') as generax_tree:
                                generaxtree = generax_tree.read()
                            gene_tree = Tree(generaxtree, format=1)
                            for gene_branch in gene_tree.traverse():
                                if g_branch == gene_branch.name:
                                    g_branch_dist = str(gene_branch.dist)
                                    
                            g_branch = g_branch.split('_branch')[1]
                            sp_branch = sp_branch.split('"')[1]
                            if 'spbranch' in sp_branch:
                                sp_branch = sp_branch.split('spbranch')[1]
                            
                            with open(mapping, 'a') as mapf:
                                towrite = clade+","+sf_xml+","+g_branch+","+g_branch_dist+","+sp_branch+","+dic_branch[sp_branch]+"\n"
                                mapf.write(towrite)

with open(mapping_stats, 'w') as mapstatsf:
    mapstatsf.write("Clade,Subfamily,Gene branch,Gene branch length,Species branch,Species branch length,lnull,lalt,deltal2,pval,qval\n")

with open(mapping, 'r') as mapf:
    for line_map in mapf:
        clade = line_map.split(',')[0]
        sf_map = line_map.split(',')[1]
        gbranch_map = line_map.split(',')[2]
        gbranchlength_map = line_map.split(',')[3]
        spbranch_map = line_map.split(',')[4]
        spbranchlength_map = line_map.split(',')[5][:-1]
        
        with open(godon_stats, 'r') as godonstats:
            for line_godonstats in godonstats:
                sf_stats = line_godonstats.split(',')[0]
                gbranch_stats = line_godonstats.split(',')[1]
                lnull_stats = line_godonstats.split(',')[2]
                lalt_stats = line_godonstats.split(',')[3]
                deltal2_stats = line_godonstats.split(',')[4]
                pval_stats = line_godonstats.split(',')[5]
                qval_stats = line_godonstats.split(',')[6][:-1]
                
                if sf_map == sf_stats and gbranch_map == gbranch_stats:
                    with open(mapping_stats, 'a') as mapstatsf:
                        mapstatsf.write(clade+","+sf_map+","+gbranch_map+","+gbranchlength_map+","+spbranch_map+","+spbranchlength_map+","+lnull_stats+","+lalt_stats+","+deltal2_stats+","+pval_stats+","+qval_stats+"\n")
                        