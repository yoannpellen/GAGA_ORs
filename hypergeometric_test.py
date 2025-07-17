from scipy.stats import hypergeom
import argparse
import os
import sys


# Description for script whith '-h'
parser = argparse.ArgumentParser(description = "Script to run the hypergeometric test after dN/dS pipeline")
optional_args = parser.add_argument_group(title = "Arguments:")

# Define arguments with help message
optional_args.add_argument('-file', dest = 'inf', help = "Mapping file with all statistics.")
optional_args.add_argument('-clade', dest = 'clade', default = "all", help = "Clade to run (see first column). Default 'all'. To run only the 'colony size' for example, put '_cs'")
optional_args.add_argument('-sf', dest = 'sf', default = "all", help = "Default 'all'. To run only the '9-Exons', put '9E'")
optional_args.add_argument('-branch', dest = 'branch', help = "Species branch to focus on.To do several branches at once, give list separated by a comma (eg: 123,312,64)")
optional_args.add_argument('-pval', dest = 'pval', default = 0.05, help = "pvalue from Godon to consider positive selection. Default '0.05'.")

# Group all arguments in a list
args = parser.parse_args()

file = args.inf
clade = args.clade
subfamily = args.sf
branch = args.branch
pvalue = float(args.pval)

print(rf"Hypergeometric test using Godon pval < {pvalue} and accounting for gene tree branch length.")
print("Clade,Species branch,Sum of the length of all gene tree branches,Total number of positive gene tree branches,Sum of the length of gene tree branches mapped on the species branch of interest,Number of positive gene tree branches mapped on the species branch of interest,pval")


for br in branch.split(','):
    M=0 # population size (total number of gene branches)
    n=0 # number of successes in the population (branches under positive selection in all gene branches)
    N=0 # sample size (gene branches mapped on the species branch of interest)
    X=0 # number of drawn successes (branches under positive selection on the species branch of interest)
    gbranch_length = {}
    gbranch_count = {}
    cs = 0
    cg = 0
    with open(file, 'r') as statsf:
        next(statsf) # skip header
        for line in statsf:
            # pval => float(line.split(',')[9])
            # qval => float(line.split(',')[10][:-1])
            if (clade in line.split(',')[0] or clade == "all") and (subfamily == line.split(',')[1] or subfamily == "all"): # Clade and subfamily of interest
                M += float(line.split(',')[3])
                cs += 1
                if pvalue >= float(line.split(',')[9]): # Any positive gene branch
                    n += 1
                if br == line.split(',')[4]: # Species branch of interest
                    N += float(line.split(',')[3])
                    cg += 1
                    if pvalue >= float(line.split(',')[9]): # Positive gene branch on species branch of interest
                        X += 1
    
    N = round(N/(M/cs))
    M = round(M/(M/cs))
    n = round(n)
    X = round(X)
    hg_pval = hypergeom.pmf(X, M, n, N)
    # Print ready for csv format
    print(clade + "," + br + "," + str(M) + "," + str(n) + "," + str(N) + "," + str(X) + "," + str(hg_pval))
