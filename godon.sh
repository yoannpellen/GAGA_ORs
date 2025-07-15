#!/bin/bash

# Runs Godon for all orthogroups for which RAxML finished successfully
# See https://github.com/idavydov/godon-tutorial for more details on parameters
# Uses 40 CPUs
# Uses M0 model to estimate branch length
# 4 categories for codon rate variation
# Test all branches

clade=$(basename "$PWD")
for folder in $(ls -d $PWD/sf_*); do
	sf=$(basename "$folder" | cut -c 4-)
	jobname="${clade}_${sf}.Godon"
	# Skip running orthogroups
	if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "$jobname"; then 
		continue
	fi
	if [ -f "$folder/raxml/RAxML_bestTree.$sf" ] && [ ! -f "$folder/Godon_OK" ]; then
		cd $folder
		rm -rf godon && mkdir godon	
		sbatch -n40 -N1 -p preempt7d -J $jobname --wrap \
		"~/softwares/godon test BS \
		--procs 40 \
		--out $folder/godon/godonBS.$sf.out \
		--m0-tree \
		--ncat-codon-rate 4 \
		--all-branches \
		$folder/guidance2/MSA.PRANK.$sf.cds.NNN0.*.trim.aln \
		$folder/raxml/RAxML_bestTree.$sf \
		&& touch Godon_OK"
		cd ..
	fi
done
