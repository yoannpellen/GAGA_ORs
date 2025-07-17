# GAGA_ORs
Scripts for the pipeline of the GAGA ORs paper

All scripts have been specifically designed to work with the GAGA dataset and its name/file system. I have renamed all sequences in the original dataset as such: “AameOR190_A_GAGA-0001_Acromyrmex-ameliae”.
There’s still a need to enter some scripts and modify some values (list of species when preparing the data, thresholds for statistics, path to scripts/programs if needed, partitions to run on, etc…), but these corrections are minimal.
Start by creating the main directory where you want to work, its name will be used in some scripts as ‘clade’ to name some files and sub-directories. Then all scripts should be run from within that directory unless told otherwise.

# Preparing the data
prepare.py, 2 variables can be set by modifying the script:
GAGAdir => directory in which we have the fasta files of all ORs subfamilies. Emphasize on the fact that the files used as input are the subfamilies fasta (e.g. ‘Subfamily_*.cds.fasta’), and the script only takes from each the sequences belonging to the species we define in species. Also, because the 9-Exon has been divided (in 2 and then 6), make sure you’re only taking from one version of the division by renaming the files you don’t want to use (just add something after ‘.fasta’ to make the script skip the file). By default now, it takes the 9-Exon into 6 groups. Technically no need to change the path unless working with different subfamilies, default set as:
/data/home/privman/ypellen/ORs/GAGA_Final_OR_annots/GAGA_final_subfamilies/
species => list of species to study, identified with their GAGA id

There’s a few checks done to know which sequences are taken or not:
-	Must not have any internal stop codons
-	Its length must be a multiple of 3, as we want to build codon alignments it would be a problem for Guidance2 if it’s not divisible by 3
-	Removes final stop codon if it’s still present

Some subfamilies might have less than 4 sequences, which is the minimum for Guidance2 to work. In such cases, you can merge the small subfamily with the closest one based on this phylogeny GAGA_final_subfamilies.MAFFTeinsi.faa.tree.png.
You will end up with one directory for each subfamily, with a fasta file inside containing all sequences for the species you gave that are part of this subfamily.
To check the number of sequences: egrep -c ">" sf_*/*.cds.trim.fna

# Guidance2
guidance.sh and guidance_correction.py
 The script is called as: bash PATH/TO/guidance_multi.sh <function>
Masking and trimming rely on the python script, the path has to be updated in the bash script if you move it.
Note that you can run any function without risk of deleting anything because it always checks before if something is still running or if it has the necessary flag file.
The different functions are:
1.	align : actual run of Guidance2, building a codon alignment using PRANK. Guidance2 creates the file ENDS_OK in case of a success, used as a flag to run further functions.
2.	mask : masking the codon alignment based on the residues score from <subfamily>/guidance2/MSA.PRANK.Guidance2_res_pair_res.scr, with a loop on the threshold value, starting at 0.8, so that we don’t mask more than 20% of the alignment. If an alignment can’t meet the requirement of 20% masked residues maximum even with a threshold at 0.1, the orthogroup will be discarded from the analysis. With the python script, checks if some sequences become identical after masking, creates a file called G_checkup.$sf.txt with the % of masked residues and gaps in the alignment. Creates G_CORREC if it works, used as a flag to run further functions.
3.	trim : removing the bad sequences. Add the names of the sequences that have been removed to G_checkup.$sf.txt. Creates G_TRIM if it works, used as a flag to run further functions.
4.	clean : removing all intermediate files that are no longer useful to make some space, basically keeping the alignments and scoring files from Guidance2.

# RAxML
raxml.sh, nothing to change here.
RAxML is runned with the masked amino acid alignment created with Guidance2, using the PROTCATLG model and 123456 as seed for the initial tree, and 100 bootstraps.
It runs only if the flag file G_TRIM has been created.
If it runs successfully, the tree file RAxML_bestTree.$sf is created and used as a flag to run the following analysis.
It also removes all intermediate files to make some space.

# Godon
godon.sh, nothing to change here. Run ‘tail -n 1 sf_*/godon/godonBS.*.out’ to check if it ended OK, should finish with a “runtime” line.

When Godon is done, run godon_stats.Hive.sh to do the FDR correction and get a phylogenetic tree with branches under positive selection marked. It relies on branchSite.lrt.fdr.r for the FDR correction and statistical test (CSV file clade*.godonBS.all.fdr.csv in your main working directory).
ALL GODON STATS ONLY RUNS LOCALLY
Run godon_stats.sh, which uses godon_correction.py to build the gene trees with internal branch names numbered according to the Godon output. 

# Generax
If not already done, start by numbering the internal branches of the species tree, in_tree is your input species tree and out_tree is the name of the output tree:
python PATH/TO/generax.py -f treelabel -it <in_tree> -ot <out_tree>
The branches will be named ‘spbranchXX’ and can be viewed using FigTree  'Branch Labels > Display > label'.
The original GAGA species tree is GAGA_phylogeny_dated-1.png, then you have this version to see the numbers for the internal branches GAGA_sp_tree.branch_numbers.jpg.

Run generax.sh, which relies on generax.py, to prepare the files used as input by Generax.

When the Generax files are created, move everything back to Hive and run generax_Hive.sh. This just runs Generax, do it from within the main working directory.

When Generax is done, run mapping_table.py to build a table showing all statistics for each gene branch, and the species branch it’s mapped on.
Hypergeometric test
You can run it for as many branches as you want, then all the tests in one go with hypergeometric_test.py. Within it, you can change the threshold for the pvalue from Godon.
