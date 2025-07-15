#!/bin/bash

# Bash script to extract the values for the null and alternative models from godon's results
# Creates a csv file for all subfamilies, with as header orthogp / branch / lnull / lalt
# Runs the R script '~/scripts/branchSite.lrt.fdr.r' to do statistics and FDR correction
# FDR correction is set at 0.05, within the R script
# Outputs a tree file with the branches under positive selection marked with a '+' (to view in FigTree -> 'Branch Labels > Display > label')

clade=$(basename "$PWD")
for folder in $(ls -d $PWD/sf_*); do
	sf=$(basename "$folder" | cut -c 4-)
	if [ ! -f "$folder/Godon_OK" ]; then
		continue
	fi
	cd $folder/godon/
	# Preparing csv file from godon result file
	echo "sf,branch,lnull,lalt" > $clade.godonBS.$sf.out.csv
	egrep "Testing branch|lnL0=|lnL1=" godonBS.$sf.out >> $clade.godonBS.$sf.out.csv
	# Simple python command to handle the reformatting of the output from Godon (maybe move it into <godon_correction.py>)
	python <<PYTHON_COMMAND
import re
import os
sf = ""
godonfile = ""
for file in os.listdir():
	if file.endswith(".out.csv"):
		sf = file.split(".")[2]
		godonfile = file
with open(godonfile, "r") as godon:
	godon_content = godon.read()
godon_content = godon_content.replace("\n","")
godon_content = godon_content.replace("Testing branch ",f"\n{sf},")
godon_content = re.sub("/\d+lnL0=",",", godon_content)
godon_content = godon_content.replace(" lnL1=","")
with open(godonfile, "w") as godon_out:
	godon_out.write(godon_content + "\n")
PYTHON_COMMAND
	if [ $? -ne 0 ]; then 
		echo "ERROR: Could not create file $folder/godon/$clade.godonBS.$sf.out.csv" >&2 
		exit 1
	fi
	cd $folder/../
done

# Merge all subfamilies into one file for fdr correction
echo "sf,branch,lnull,lalt" > $PWD/$clade.godonBS.all.csv
for file in $PWD/sf_*/godon/$clade.godonBS.*.out.csv; do
	tail -n +2 $file >> $PWD/$clade.godonBS.all.csv # Copy $file content from the second line (without header)
	if [ $? -ne 0 ]; then 
		echo "ERROR: Could not add $file to $clade.godonBS.all.csv" >&2 
		exit 1
	fi
done

# Running R script to get the statistics
echo "yoann" | sudo -S -u yoann R --slave --args $PWD/$clade.godonBS.all.csv $PWD/$clade.godonBS.all.fdr.csv < /mnt/c/Users/user/OneDrive/Bureau/Thesis/workspace/GAGA/scripts/branchSite.lrt.fdr.r
if [ $? -ne 0 ]; then 
	echo "ERROR: Problem running R" >&2
	exit 1
fi

# Export result for each subfamily back to seperate files
clade=$(basename "$PWD")
for folder in $(ls -d $PWD/sf_*); do 
	sf=$(basename "$folder" | cut -c 4-)
	grep_pattern="$sf,"
	echo "sf,branch,lnull,lalt,deltal2,pval,qval" > $folder/godon/$clade.godonBS.$sf.fdr.csv
	grep "^$grep_pattern" $PWD/$clade.godonBS.all.fdr.csv >> $folder/godon/$clade.godonBS.$sf.fdr.csv
	# Creating output tree
	python /mnt/c/Users/user/OneDrive/Bureau/Thesis/workspace/GAGA/scripts/godon_correction.py $folder/godon/$clade.godonBS.$sf.fdr.csv $folder/godon/godonBS.$sf.out $folder/raxml/RAxML_bestTree.$sf $folder/$clade.godonBS.$sf.fdr.tree
	if [ $? -ne 0 ]; then 
		echo "ERROR getting tree for $sf" >&2 
		exit 1
	fi
done

cp sf_*/*fdr.tree .
