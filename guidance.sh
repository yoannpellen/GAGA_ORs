#!/bin/bash

# GAGA version of the Guidance2 script for the positive selection pipeline
# Each function successfully executed will create a file saying so, the next function can only be run if this file exists
# From the main working directory: bash PATH/TO/guidance.sh <function>
#	1st > align
#	2nd > mask
#	3rd > trim
#	4th > clean
# In case some subfamilies don't have the required 400 alignments, wait for all others to be done and run the 5th function
#	5th > fail

aligning() {
	# Run Guidance2 for each subfamily separately
	# codon alignment using PRANK
	# Set to use 20 CPUs
	clade=$(basename "$PWD")
	rm -f guidance_fail.txt
	for folder in $(ls -d $PWD/sf_*); do
		sf=$(basename "$folder" | cut -c 4-)
		jobname="${clade}_${sf}.G"
		# Skip running subfamilies
		if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "$jobname"; then 
			continue
		fi
		if [ ! -f "$folder/guidance2/ENDS_OK" ]; then
			rm -rf $folder/guidance2_err/ $folder/guidance2/ $folder/slurm*
			cd $folder
			sleep 1s
			sbatch -n20 -N1 -p preempt1d -J $jobname --wrap \
			"echo '$sf Guidance2 with PRANK'
			module load prank
			sleep 1s
			perl ~/softwares/guidance.v2.02/www/Guidance/guidance.pl \
			--seqfile $folder/$sf.cds.trim.fna \
			--seqType codon \
			--outDir $folder/guidance2 \
			--program GUIDANCE2 \
			--msaProgram PRANK \
			--proc_num 20
			sleep 1s
			if [ ! -f '$folder/guidance2/ENDS_OK' ]; then
				echo -e '$sf\tFAILED' >> ../guidance_fail.log
			fi"
			cd ..
		fi
	done
}

masking() {
	# Starting from masking at 0.8, loops with 0.1 steps to mask less than 20% of the alignment
	sf=$(basename "$PWD" | cut -c 4-)
	rm -f $PWD/guidance2/MSA.PRANK.$sf*NNN* $PWD/guidance2/G_CORREC*
	cp $PWD/guidance2/MSA.PRANK.aln.With_Names $PWD/guidance2/MSA.PRANK.$sf.cds.aln
	for ((i=8; i>=0; i--)); do
		if [[ "$i" -eq 0 ]]; then
			echo "The alignment of $sf has a bad score, no threshold masks less than 20%" > $PWD/guidance2/BAD_ALN
			break
		fi
		mask="0.$i"
		perl ~/softwares/guidance.v2.02/www/Guidance/maskLowScoreResidues.pl \
		$PWD/guidance2/MSA.PRANK.$sf.cds.aln \
		$PWD/guidance2/MSA.PRANK.Guidance2_res_pair_res.scr \
		$PWD/guidance2/MSA.PRANK.$sf.cds.NNN$mask.aln \
		$mask \
		nuc
		python ~/ORs/GAGA/PosSel/scripts/guidance_correction.py -f translate -nucl $PWD/guidance2/MSA.PRANK.$sf.cds.NNN$mask.aln -aa $PWD/guidance2/MSA.PRANK.$sf.aa.NNN$mask.aln
		masked=''
		masked=$(python ~/ORs/GAGA/PosSel/scripts/guidance_correction.py -f masking -sf $sf -aa $PWD/guidance2/MSA.PRANK.$sf.aa.NNN$mask.aln)
		if [[ "$masked" -lt 20 ]]; then
			echo "Final masking at $mask for subfamily $sf"
			break
		fi
		rm -f $PWD/guidance2/MSA.PRANK.$sf.*.NNN*.aln
	done
	
	if [ ! -f "$PWD/guidance2/BAD_ALN" ]; then
		# Checking for identical or fully masked sequences after translated to amino acids
		check_seq=''
		check_seq=$(python ~/ORs/GAGA/PosSel/scripts/guidance_correction.py -f seq_check -sf $sf -aa $PWD/guidance2/MSA.PRANK.$sf.aa.NNN$mask.aln -log $PWD/G_checkup.$sf.txt)
		echo "$check_seq" >> $PWD/G_checkup.$sf.txt
		touch $PWD/guidance2/G_CORREC
	fi
}

trimming() {
	# Delete sequences than can't be used after masking with Guidance2
	# Either 2 become identical, or 1 is completely masked
	clade=$(basename "$PWD")
	for folder in $(ls -d $PWD/sf_*); do
		sf=$(basename "$folder" | cut -c 4-)
		jobname="${clade}_${sf}.Gt"
		# Skip running subfamilies
		if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "$jobname"; then 
			continue
		fi
		if [ -f "$folder/guidance2/G_CORREC" ] && [ ! -f "$folder/guidance2/G_TRIM" ]; then
			cd $folder
			sleep 1s
			sbatch -n1 -N1 -p preempt1d -J $jobname --get-user-env --wrap \
			"python ~/ORs/GAGA/PosSel/scripts/guidance_correction.py -f trim -log G_checkup.$sf.txt \
			&& touch $folder/guidance2/G_TRIM"
		fi
		cd ..
	done
}

cleaning() {
	# Delete all useless files after Guidance2 succeeds
	# Can be runned even if all jobs aren't finished
	clade=$(basename "$PWD")
	for folder in $(ls -d $PWD/sf_*); do
		sf=$(basename "$folder" | cut -c 4-)
		jobname="${clade}_${sf}.Gc"
		if [ -f "$folder/guidance2/ENDS_OK" ]; then
			cd $folder/guidance2
			sbatch -n1 -N1 -p preempt1d -J $jobname --wrap \
			"rm -f COS.std *JalView* Seqs* *.tar.gz *html Sampled* *Without* *semphy.tree SampledOPVals.log ../slurm*.out"
			cd ../..
		fi
	done
}

failing() {
	clade=$(basename "$PWD")
	rm -f guidance_fail.txt
	for folder in $(ls -d $PWD/sf_*); do
		sf=$(basename "$folder" | cut -c 4-)
		jobname="${clade}_${sf}.Gf"
		# Skip running subfamilies
		if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "$jobname"; then 
			continue
		fi
		if [ ! -f "$folder/guidance2/ENDS_OK" ]; then
			cd $folder
			sbatch ~/ORs/GAGA/PosSel/scripts/guidance_fail.sh "$jobname"
			cd ..
		fi
	done
}

case $1 in

	align)
		aligning;;

	mask)
		# Easier to separate the loop and the call of the function through sbatch
		# The quote system in sbatch is giving me headaches
		clade=$(basename "$PWD")
		for folder in $(ls -d $PWD/sf_*); do
			sf=$(basename "$folder" | cut -c 4-)
			jobname="${clade}_${sf}.Gm"
		# Skip running orthogroups, whatever the function
			if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "$jobname"; then 
				continue
			fi
			if [ -f "$folder/guidance2/ENDS_OK" ] && [ ! -f "$folder/guidance2/G_CORREC" ]; then
				cd $folder
				sleep 1s
				sbatch -n1 -N1 -p preempt1d -J $jobname --wrap \
				"bash ~/ORs/GAGA/PosSel/scripts/guidance.sh masking"
				cd ..
			fi
		done;;

	masking)
		masking;;

	trim)
		trimming;;

	clean)
		cleaning;;

	fail)
		failing;;
	
esac
