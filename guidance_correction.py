# Utility python script for the <guidance_multi.sh> bash script
# Performs various operations helpful for the positive selection pipeline

import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import time

# Description for script whith '-h'
parser = argparse.ArgumentParser(description = "Different functions needed to deal with Guidance2")
optional_args = parser.add_argument_group(title = "Arguments (depending on the function called)")

# Define arguments with help message
optional_args.add_argument('-f', dest = 'function', help = "Function to call (rename, masking, seq_check, delete, translate, reverse)")
optional_args.add_argument('-sf', dest = 'sf', help = "Orthogroup")
optional_args.add_argument('-nucl', dest = 'nucl_aln', help = "Nucleic acid alignment")
optional_args.add_argument('-aa', dest = 'aa_aln', help = "Amino acid alignment")
optional_args.add_argument('-log', dest = 'logf', help = "Sequence check up file")

# Group all arguments in a list
args = parser.parse_args()

def rename(sf):
    out_name = "guidance2/MSA.PRANK." + sf + ".fullname.aln"
    with open("guidance2/MSA.PRANK.aln", "r") as msa:
        msa_data = msa.read()
    with open("guidance2/Seqs.Codes", "r") as seq_code:
        for line in seq_code:
            OR = line.split("\t")[0]
            OR = ">"+OR
            code = line.split("\t")[1]
            code = ">"+code
            msa_data = msa_data.replace(code, OR+"\n")
    with open(out_name, "w") as msa_out:
        msa_out.write(msa_data)

def delete():
    with open("sequences_checkup.txt", 'r') as log:
        for line in log:
            if line .startswith("Subfamily"):
                sf = line.split(":")[0][10:]
                mask = line.split(" ")[4][:-1]
                sf_file = "guidance2/MSA.PRANK." + sf + ".aa.NNN.aln"
                sf_file2 = "guidance2/MSA.PRANK." + sf + ".aa.NNN" + mask + ".trim.aln"
            elif line.startswith("   ->"):
                OR1 = line.split(" ")[4] + "_"
                OR2 = line.split(" ")[6] + "_"
                with open(sf_file, 'r') as sf_fasta:
                    for record in SeqIO.parse(sf_fasta, 'fasta'):
                        if record.id.startswith(OR1) or record.id.startswith(OR2):
                            os.system(rf"sed -i '/{record.id}/d' {sf_file2}")
                            os.system(rf"sed -i '/{record.seq}/d' {sf_file2}")

def reverse(aa_aln):
    # First, remove stop codon if present in nucleotide sequences
    os.remove("guidance2/Seqs.Orig_DNA.fas.FIXED")
    with open("guidance2/Seqs.Orig_DNA.fas", 'r') as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            if record.seq.endswith(("TAG", "TGA", "TAA")):
                record.seq = record.seq[:-3]
            with open("guidance2/Seqs.Orig_DNA.fas.FIXED", "a") as outfile:
                SeqIO.write(record, outfile, 'fasta')
    # Building alignment
    nucl_aln = aa_aln.replace(".aa.", ".cds.")
    open(nucl_aln, "w").close()
    with open(aa_aln, 'r') as aa_file:
        for aa_record in SeqIO.parse(aa_file, 'fasta'):
            with open("guidance2/Seqs.Orig_DNA.fas.FIXED", 'r') as nucl_source:
                for nucl_record in SeqIO.parse(nucl_source, 'fasta'):
                    if aa_record.id == nucl_record.id:
                        seq = str(nucl_record.seq)
                        nuclcodon = [seq[i:i+3] for i in range(0,len(seq),3)]
                        newcodon = []
                        aaflag = 0
                        for aa in aa_record.seq:
                            if aa == "-":
                                newcodon.append('---')
                            elif aa == "X":
                                newcodon.append('NNN')
                                aaflag += 1
                            else:
                                newcodon.append(nuclcodon[aaflag])
                                aaflag += 1
                        nucl_record.seq = ''.join(newcodon)
                        with open(nucl_aln, 'a') as nucl:
                            nucl.write(">" + nucl_record.id + "\n")
                            nucl.write(str(nucl_record.seq) + "\n")

def masking(aln, sf):
    # Checks the % of masked residues in the alignment
    percent_mask = 100
    masked_aa = 0
    total_aa = 0
    with open(aln, "r") as msa:
        for sequence in SeqIO.parse(msa, "fasta"):
            masked_aa += sequence.seq.count("X")
            total_aa += len(sequence.seq)
    percent_mask = masked_aa/total_aa
    print(int(percent_mask*100))

def translate(nucl, aa):
    # Translate nucleotide sequences into amino acids
    with open(nucl, 'r') as nucl_fasta:
        for sequence in SeqIO.parse(nucl_fasta, 'fasta'):
            seq_aa = Seq(str(sequence.seq))
            seq_aa = str(seq_aa.translate())
            with open(aa, 'a') as aa_fasta:
                aa_fasta.write(">" + sequence.id + "\n")
                aa_fasta.write(seq_aa + "\n")

def seq_check(aln, sf):
    # Checks sequences in the alignment right after masking
    # It opens the same alignment twice (msa and ref_msa), and compares the sequences one by one
    # To keep track of identical sequences, and make sure I don't have useless duplicates in the output, identical sequences names are stored in a list
    # for example, seq_1 and seq_5 are identical, I enter "seq_1seq_5" and "seq_5seq_1" in the list, this way I can skip sequences further in the alignment if they already had a match with a previous one
    # Outputs the % of masked residues and gaps for the calling bash script 
    id_list=[]
    masked_aa = 0
    gap_aa = 0
    total_aa = 0
    mask = aln.split("NNN")[1][:3]
    with open(aln, 'r') as msa:
        for msa_sequence in SeqIO.parse(msa, "fasta"):
            if all(aa in "-X" for aa in str(msa_sequence.seq)):
                print(rf"   -> all masked: {msa_sequence.id} ")
                continue
            with open(aln, "r") as ref_msa:
                for ref_sequence in SeqIO.parse(ref_msa, "fasta"):
                    if msa_sequence.seq == ref_sequence.seq:
                        merge_name = rf"{msa_sequence.id}{ref_sequence.id}"
                        merge_name_reverse = rf"{ref_sequence.id}{msa_sequence.id}"
                        if msa_sequence.id != ref_sequence.id and merge_name not in id_list and merge_name_reverse not in id_list:
                            # This 'print' is what's added to the Guidance checkup log file when there's duplicated sequences,
                            # as it's a simple print statement, it's returned in <$check_seq> in the bash script
                            print(rf"   -> same seq: {msa_sequence.id} and {ref_sequence.id} ")
                            id_list.extend([merge_name,merge_name_reverse])
            masked_aa += msa_sequence.seq.count("X")
            gap_aa += msa_sequence.seq.count("-")
            total_aa += len(msa_sequence.seq)
    percent_mask = masked_aa/total_aa
    percent_gap = gap_aa/total_aa
    outlog = rf"G_checkup.{sf}.txt"
    open(outlot, 'w').close()
    with open(outlog, 'a') as log:
        log.write(rf"{sf}: masking at {mask}")
        log.write("\n")
        log.write(rf"   {percent_mask:.0%} of the alignment is masked")
        log.write("\n")
        log.write(rf"   {percent_gap:.0%} of gaps in the alignment")
        log.write("\n")
    if len(id_list) == 0:
        print("")

def trim(logf):
    # Creates a trimmed msa without the identical sequences listed in the log file created by the 'seq_check' function
    with open(logf) as log:
        first_line = log.readline()
    sf = first_line.split(":")[0]
    mask = first_line.split(" ")[3][:-1]
    with open(logf, 'r') as log:
        logdata = log.read()
    orthoin_aa = "guidance2/MSA.PRANK." + sf + ".aa.NNN" + mask + ".aln"
    orthoin_cds = "guidance2/MSA.PRANK." + sf + ".cds.NNN" + mask + ".aln"
    orthoout_aa = "guidance2/MSA.PRANK." + sf + ".aa.NNN" + mask + ".trim.aln"
    open(orthoout_aa, 'w').close()
    orthoout_cds = "guidance2/MSA.PRANK." + sf + ".cds.NNN" + mask + ".trim.aln"
    open(orthoout_cds, 'w').close()
    with open(orthoin_aa, 'r') as in_aa:
        for sequence in SeqIO.parse(in_aa, 'fasta'):
            if sequence.id not in logdata:
                with open(orthoout_aa, 'a') as out_aa:
                    SeqIO.write(sequence, out_aa, 'fasta')
    with open(orthoin_cds, 'r') as in_cds:
        for sequence in SeqIO.parse(in_cds, 'fasta'):
            if sequence.id not in logdata:
                with open(orthoout_cds, 'a') as out_cds:
                    SeqIO.write(sequence, out_cds, 'fasta')


if __name__ == "__main__":
    
    if args.function == "rename":
        rename(args.sf)
    elif args.function == "delete":
        delete()
    elif args.function == "reverse":
        reverse(args.aa_aln)
    elif args.function == "masking":
        masking(args.aa_aln, args.sf)
    elif args.function == "seq_check":
        seq_check(args.aa_aln, args.sf)
    elif args.function == "trim":
        trim(args.logf)
    elif args.function == "translate":
        translate(args.nucl_aln, args.aa_aln)
    else:
        print("ERROR: function '" + args.function + "' doesn't exist")
