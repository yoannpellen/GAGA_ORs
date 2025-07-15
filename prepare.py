# Script to create a working directory for each subfamily for all selected species
# Run as 'python PATH/TO/prepare.py' from the clade directory

from Bio import SeqIO
import os
import shutil
import subprocess
import shlex

GAGAdir = "/lustre1/home/privman/ypellen/ORs/GAGA/GAGA_final_subfamilies/"
root = os.getcwd()
clade = os.path.basename(root)


# List of species to select based on their GAGA ID
#species = ["NCBI-0003","GAGA-0002","GAGA-0014","GAGA-0003","NCBI-0017","NCBI-0016","NCBI-0006","GAGA-0539","NCBI-0015"] # Clade 1 colony size
#species = ["GAGA-0495","GAGA-0359","GAGA-0502","GAGA-0494"] # Clade2 colony size
#species = ["GAGA-0366","GAGA-0524","GAGA-0024","GAGA-0177","GAGA-0332","GAGA-0341","GAGA-0364"] # Clade3 colony size New
#species = ["GAGA-0109","GAGA-0503","GAGA-0356","GAGA-0004"] # Clade4 colony size
#species = ["NCBI-0003","GAGA-0014","GAGA-0358","GAGA-0245","NCBI-0016","GAGA-0229","GAGA-0540","GAGA-0543","NCBI-0015"] # Clade 1 diet
#species = ["GAGA-0004","GAGA-0109","GAGA-0505","GAGA-0084","GAGA-0503","GAGA-0356","GAGA-0413"] # Clade 2 diet
#species = ["GAGA-0554","NCBI-0014","GAGA-0358","GAGA-0384","GAGA-0530"] # Clade 5 colony size
#species = ["NCBI-0009","GAGA-0379","GAGA-0528","GAGA-0363","OUT-0001","GAGA-0177","GAGA-0028","GAGA-0014","GAGA-0392","GAGA-0552"] # Formicoid branch
#species = ["GAGA-0024","GAGA-0026","GAGA-0275","GAGA-0198","GAGA-0332","GAGA-0341"] # Clade 3 diet
#species = ["GAGA-0087","GAGA-0028","GAGA-0401","GAGA-0114","NCBI-0012"] # Clade 4 diet
#species = ["GAGA-0510","GAGA-0099","GAGA-0395","GAGA-0378","GAGA-0533","GAGA-0382","GAGA-0520","GAGA-0405","GAGA-0515"] # clade 1 worker reproductive function
#species = ["GAGA-0002","NCBI-0003","NCBI-0004","GAGA-0538","GAGA-0229","GAGA-0530","GAGA-0384","GAGA-0358","NCBI-0014","GAGA-0554","NCBI-0002","GAGA-0454","GAGA-0245","GAGA-0346"] #clade 2 worker reproductive function
species = ["GAGA-0404","GAGA-0391","GAGA-0536","GAGA-0552","GAGA-0537","GAGA-0389","GAGA-0580","GAGA-0535"] # clade 3 worker reproductive function
#species = ["GAGA-0352","GAGA-0387","GAGA-0025","GAGA-0408","GAGA-0386","GAGA-0074","GAGA-0353","GAGA-0301","GAGA-0246"] # clade 4 worker reproductive function
#species = ["OUT-0001","NCBI-0008","GAGA-0359","GAGA-0502","GAGA-0494","GAGA-0354","GAGA-0304","GAGA-0055","NCBI-0005","GAGA-0200","GAGA-0374","GAGA-0360","GAGA-0336","GAGA-0187"] # clade 1 workers polymorphism
#species = ["GAGA-0177","GAGA-0366","GAGA-0020","GAGA-0026","GAGA-0275","GAGA-0198","GAGA-0332","GAGA-0341"] # clade 2 workers polymorphism
#species = ["GAGA-0530","GAGA-0384","GAGA-0358","NCBI-0014","GAGA-0554","NCBI-0002","GAGA-0454","GAGA-0245","GAGA-0346"] # clade 3 workers polymorphism
#species = ["GAGA-0003","GAGA-0002","GAGA-0014","NCBI-0004","NCBI-0016","GAGA-0538","GAGA-0541","GAGA-0229"] # clade 4 workers polymorphism
#species = ["GAGA-0378","GAGA-0382","GAGA-0579","GAGA-0533","GAGA-0520","GAGA-0085","GAGA-0333","GAGA-0405","GAGA-0515"] # clade 5 workers polymorphism


preplog = "prep.log"
open(preplog,'w').close()

for file in os.listdir(GAGAdir):
    GAGA_sf = os.fsdecode(file)
    if GAGA_sf.startswith("Subfamily_") and GAGA_sf.endswith(".cds.fasta"):
        sf = GAGA_sf[10:-10]
        with open(preplog, 'a') as log:
            log.write("Taking sequences from subfamily " + sf + "\n")
        sf_dir = os.getcwd() + "/sf_" + sf + "/"
        os.makedirs(os.path.dirname(sf_dir), exist_ok=True)
        sf_fasta = sf_dir + sf + ".cds.fna"
        GAGA_sf = GAGAdir + GAGA_sf
        with open(GAGA_sf, 'r') as GAGAsf:
            for sequence in SeqIO.parse(GAGAsf, "fasta"):
                GAGAid = sequence.id.split("_")[2]
                stop = False
                for sp in species:
                    if sp == GAGAid:
                    # Removing sequence if no divisible by 3
                        if len(sequence.seq)%3 != 0:
                            with open(preplog, 'a') as log:
                                log.write(str(sequence.id) + " discarded because not divisible by 3\n")
                            stop = True
                            continue
                        # Removing final stop codon
                        nucl_seq = str(sequence.seq)
                        nucl_seq = nucl_seq.upper()
                        if nucl_seq.endswith(("TAA", "TAG", "TGA")):
                            sequence.seq = sequence.seq[:-3]
                        nucl_seq = str(sequence.seq)
                        nucl_seq = nucl_seq.upper()
                        # Not taking sequences with internal stop codon
                        for i in range(0,len(nucl_seq),3):
                            codon = nucl_seq[i:i+3]
                            if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                                with open(preplog, 'a') as log:
                                    log.write(str(sequence.id) + " discarded because of internal stop codon\n")
                                stop = True
                                break
                        if stop == False:
                            with open(sf_fasta, 'a') as sffasta:
                                SeqIO.write(sequence, sffasta, 'fasta')
print("Done creating all subfamily files")

print("Checking for duplicated sequences and removing them")
listall = []
for path, dirs, files in os.walk(root):
    for file in files:
        if file.endswith(".cds.fna"):
            sf = file.split('.')[0]
            sf_fna = "sf_" + sf + "/" + file
            with open(sf_fna, 'r') as fna:
                id_list=[]
                for fna_seq in SeqIO.parse(fna, "fasta"):
                    with open(sf_fna, "r") as ref_fna:
                        for ref_seq in SeqIO.parse(ref_fna, "fasta"):
                            if fna_seq.seq == ref_seq.seq:
                                if fna_seq.id != ref_seq.id and fna_seq.id not in id_list and fna_seq.id not in listall:
                                    id_list.append(fna_seq.id)
                if id_list:
                    with open(preplog, 'a') as log:
                        log.write(sf + " -> identical sequences: " + str(id_list)[1:-1] + "\n")

with open(preplog, 'r') as log:
    logdata = log.read()
for path, dirs, files in os.walk(root):
    for file in files:
        if file.endswith(".cds.fna"):
            sf = file.split('.')[0]
            ref_fna = "sf_" + sf + "/" + sf + ".cds.fna"
            trim_fna = "sf_" + sf + "/" + sf + ".cds.trim.fna"
            open(trim_fna, 'w').close()
            with open(ref_fna, 'r') as ref:
                for sequence in SeqIO.parse(ref, 'fasta'):
                    if sequence.id not in logdata:
                        with open(trim_fna, 'a') as trim:
                            SeqIO.write(sequence, trim, 'fasta')

print("Done\nRun [egrep -c '>' sf_*/*.cds.trim.fna] for subfamilies count\n")

