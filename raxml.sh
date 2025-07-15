#!/bin/bash

# Run RAxML for all subfamilies
# PROTCATLG model
# Random starting seed '123456'
# 100 bootstraps
# 40 CPUs

clade=$(basename $PWD)
rm -f RAxML_fail.log
for folder in $(ls -d $PWD/sf_*); do
	sf=$(basename "$folder" | cut -c 4-)
	jobname="${clade}_${sf}.RAxML"
	# Skip running orthogroups
	if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "$jobname"; then 
		continue
	fi
	if [ -f "$folder/guidance2/G_TRIM" ] && [ ! -f "$folder/guidance2/BAD_ALN" ] && [ ! -f "$folder/raxml/RAxML_bestTree.$sf" ]; then
		cd $folder
		sleep 1s
		sbatch -n40 -N1 -p preempt7d -J $jobname --wrap \
		"echo '$sf RAxML'
		module load raxml/1.2.2
		sleep 1s
		rm -rf raxml && mkdir raxml
		raxmlHPC-PTHREADS \
		-s $folder/guidance2/MSA.PRANK.$sf.aa.NNN0.*.trim.aln \
		-n $sf \
		-m PROTCATLG \
		-p 123456 \
		-N 100 \
		-f d \
		-T 40 \
		-w $folder/raxml
		sleep 1s
		if [ -f '$folder/raxml/RAxML_bestTree.$sf' ]; then
			rm -f raxml/*RUN*
		else
			echo -e '$sf\tFAILED' >> ../RAxML_fail.log
		fi"
		cd ..
	fi
done
