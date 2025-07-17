#!/bin/bash

# Run as "bash PATH/TO/generax.sh" from within the clade directory

clade="$(basename "$PWD")"

rm -rf generax_$clade
mkdir generax_$clade
cp PATH/TO/GAGAsp.tree generax_$clade
echo "[FAMILIES]" > generax_$clade/families_$clade.txt

for tree in $(ls $PWD/*fdr.tree); do
	sf="$(cut -d'.' -f3 <<<"$(basename "$tree")")"
	
	cp $tree generax_$clade/generax.$sf.tree
	cp sf_$sf/guidance2/MSA.PRANK.$sf.cds.NNN0.*.trim.aln generax_$clade/
	
	python PATH/TO/generax.py -f mapping -it generax_$clade/generax.$sf.tree -map generax_$clade/mapping.$sf.link
	python PATH/TO/generax.py -f msalabel -imsa generax_$clade/MSA.PRANK.$sf.cds.NNN0.*.trim.aln -omsa generax_$clade/generax.$sf.fasta -map generax_$clade/mapping.$sf.link
	
	# Create family file
	echo "- family_$sf" >> generax_$clade/families_$clade.txt
	echo "starting_gene_tree = generax.$sf.tree" >> generax_$clade/families_$clade.txt
	echo "alignment = generax.$sf.fasta" >> generax_$clade/families_$clade.txt
	echo "mapping = mapping.$sf.link" >> generax_$clade/families_$clade.txt
	echo "subst_model = LG+I" >> generax_$clade/families_$clade.txt
done


cd generax_$clade/
rm -rf Generax MSA.PRANK* *xml slurm*
sbatch -n20 -N1 -J $clade.Generax -p privmanq --mem 10G --wrap \
"generax -f families_$clade.txt -s GAGAsp.tree --prune-species-tree -r UndatedDL --unrooted-gene-tree --reconcile \
&& cp GeneRax/reconciliations/*.xml ./
