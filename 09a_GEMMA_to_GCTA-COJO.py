#make GWAS output into GCTA-COJO format
#by Angela Andaleon (aandaleon@luc.edu)
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser() #I could do this in bash but flag passing in bash is terrible
parser.add_argument("--desc", type = str, action = "store", dest = "desc", required = False, help = "make GWAS output into GCTA-COJO format")
parser.add_argument("--fam", type = str, action = "store", dest = "fam", required = True, help = ".fam file path")
#parser.add_argument("--GWAS_prefix", type = str, action = "store", dest = "GWAS_prefix", required = True, help = "Prefix of GWAS results files (not including output/).")
#parser.add_argument("--output_prefix", type = str, action = "store", dest = "output_prefix", required = True, help = "Prefix of output file for GCTA-COJO input.")
parser.add_argument("--pheno_names", type = str, action = "store", dest = "pheno_names", required = False, default = "pheno_names.txt", help = "File containing pheno names. Default = 'pheno_names.txt'.")
args = parser.parse_args()

fam = args.fam
#GWAS_prefix = args.GWAS_prefix
#output_prefix = args.output_prefix
GWAS_prefix = ""
pheno_names = list(np.loadtxt(args.pheno_names, dtype = str))

'''
fam = "AMR.fam"
GWAS_prefix = "AMR_"
output_prefix = "AMR"
'''

num_ind = sum(1 for ind in open(fam, "r")) #https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
for pheno_name in pheno_names:
    run_awk = "cat output/" + GWAS_prefix + "" + pheno_name + "_chr*assoc.txt | awk ' {print $2\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8\"\\t\"$9\"\\t\"$12\"\\t\"" + str(num_ind) + "} ' > " + pheno_name + ".ma" #tabs are fun
    os.system(run_awk) #awk is nice and fast and so I don't have to load things; hooray bash for reorganizing data
    print("Output is in " + pheno_name + ".ma")
