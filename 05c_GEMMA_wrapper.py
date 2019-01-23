#runs a GWAS loop around given phenotypes; all input except for genotypes must be already in GEMMA format
import argparse
import numpy as np
import pandas as pd
import os
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
parser.add_argument("--relatedness", type = str, action = "store", dest = "relatedness", required = True, help = "Path to file containing relatedness matrix w/o IIDs for only individuals in analysis.")
parser.add_argument("--BIMBAM_prefix", type = str, action = "store", dest = "BIMBAM", default = "chr", required = False, help = "Prefix of BIMBAM files, not including 'BIMBAM/'. Default = 'chr'.")
parser.add_argument("--pheno", type = str, action = "store", dest = "pheno", required = True, help = "Path to file containing phenotypic information w/o IIDs for only individuals in analysis.")
parser.add_argument("--pheno_names", type = str, action = "store", dest = "pheno_names", required = True, help = "Path to file containing names of phenotypes.")
parser.add_argument("--covariates", type = str, action = "store", dest = "covariates", required = False, help = "Path to file containing covariates w/o IIDs for only individuals in analysis.")
parser.add_argument("--anno", type = str, action = "store", dest = "anno", required = False, help = "(Optional) Path to file containing the annotations, including 'anno'.")
parser.add_argument("--output", type = str, action = "store", dest = "output", required = False, default = "", help = "Name of output file")
args = parser.parse_args()

print("Reading input files.")
#following are to be used in GEMMA input
if args.anno is None:
    anno = " "
else:
    anno = " -a " + args.anno
if args.covariates is None:
    covariates_file = " "
else:
    covariates_file = " -c " + args.covariates + " "
pheno_file = args.pheno
pheno_names = args.pheno_names
relatedness = args.relatedness
output = args.output
BIMBAM_prefix = args.BIMBAM

'''
anno = " -a anno/anno"
covariates_file = " -c covar_woIID.txt "
pheno_file = "pheno_woIID.txt"
pheno_names = "pheno_names.txt"
relatedness = "relatedness_woIID.txt"
BIMBAM_prefix = "chr"
output = "AMR_"
'''

pheno_names = list(np.loadtxt(pheno_names, dtype = "str"))
pheno_nums = range(1, len(pheno_names) + 1) #cause GEMMA does things by indexes

#phenotype loop    
for pheno_num, pheno_name in zip(pheno_nums, pheno_names):
    print("Starting analyses on " + pheno_name + ".")
    for chr in range(1, 23):
        GEMMA_command = "gemma -g BIMBAM/" + BIMBAM_prefix + str(chr) + ".txt.gz -p " + pheno_file + " -n " + str(pheno_num) + anno + str(chr) + ".txt -k " + relatedness + covariates_file + " -lmm 4 -notsnp -o " + output + "_" + pheno_name
        os.system(GEMMA_command + " >> GEMMA_log.txt")
        print("Completed with chr. " + str(chr) + ".")
    print("Ending analyses on " + pheno_name + ".")
print("Analyses in all phenotypes is complete. Have a nice day :)!")
