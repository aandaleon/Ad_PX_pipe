#uses local ancestry as a dosage to run an ancestry-by-ancestry level admixture mapping analysis; all input except for 'genotypes' must be already in GEMMA format, and if some individuals and/or SNPs got removed, use 02_PrediXcan_dosages_to_GEMMA.py and 21_make_GEMMA_input.R to make input proper again
#by Angela Andaleon (aandaleon@luc.edu)
import argparse
import numpy as np
import pandas as pd
import os
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
#novel-ish part of using GEMMA
parser.add_argument("--snptable", type = str, action = "store", dest = "snptable", required = True, help = "Path to file containing the .csv output of 25_RFMix_loc_anc.py")
parser.add_argument("--refpop", type = str, action = "store", dest = "refpop", required = True, help = "Reference population (AFA or HIS)")

#established part of using GEMMA
parser.add_argument("--relatedness", type = str, action = "store", dest = "relatedness", required = False, default = "relatedness_woIID.txt", help = "Path to file containing relatedness matrix w/o IIDs for only individuals in analysis.")
parser.add_argument("--chr", type = str, action = "store", dest = "BIMBAM", required = False, help = "Path to file with BIMBAM-formatted genotypes.")
parser.add_argument("--pheno", type = str, action = "store", dest = "pheno", required = False, default = "pheno_woIID.txt", help = "Path to file containing phenotypic information w/o IIDs for only individuals in analysis.")
parser.add_argument("--covariates", type = str, action = "store", dest = "covariates", required = False, default = "GEMMA_covars.txt", help = "Path to file containing covariates w/o IIDs for only individuals in analysis.")
parser.add_argument("--anno", type = str, action = "store", dest = "anno", required = False, default = "anno/anno", help = "Path to file containing the annotations.")
#parser.add_argument("--output", type = str, action = "store", dest = "output", required = False, default = "", help = "Name of output file")
parser.add_argument("--pheno_names", type = str, action = "store", dest = "pheno_names", required = False, default = "pheno_names.txt", help = "File containing pheno names. Default = 'pheno_names.txt'.")
args = parser.parse_args()

print("Reading input files.")
loc_anc_cov = pd.read_csv(args.snptable, delimiter=',', encoding="utf-8-sig")
os.system("zcat BIMBAM/chr" + args.BIMBAM + ".txt.gz | awk '{ print $1, $2, $3 }' > loc_anc_output/SNPs_" + args.BIMBAM + ".txt") #there's not really a point to loading the entire BIMBAM if I'm just using the first three cols
SNPs = pd.read_csv("loc_anc_output/SNPs_" + args.BIMBAM + ".txt", delimiter=' ', encoding="utf-8-sig", header = None)

#following are just to be used in GEMMA input
if args.anno is None:
    anno = " "
else:
    anno = " -a " + args.anno + args.BIMBAM " "
if args.covariates is None:
    covariates_file = " "
else:
    covariates_file = " -c " + args.covariates + " "
pheno_file = args.pheno
pheno_names = list(np.loadtxt(args.pheno_names, dtype = str))
relatedness = args.relatedness
output = args.BIMBAM
refpop = args.refpop

'''
loc_anc_cov = pd.read_csv("loc_anc_input/chr22.csv", delimiter=',', encoding="utf-8-sig")
SNPs = pd.read_csv("loc_anc_input/SNPs_22.txt", delimiter=' ', encoding="utf-8-sig", header = None)
anno = " -a anno/anno22.txt "
covariates_file = " -c covariates_woIID.txt "
pheno_file = "pheno_woIID.txt"
relatedness = "relatedness.txt"
pheno_names = list(np.loadtxt("pheno_names.txt", dtype = str))
output = "chr22"
refpop = "HIS"
'''

print("Formatting input for processing.")
loc_anc_cov['IID'] = loc_anc_cov['IID'].str.replace(r'.*:', '')
inds = loc_anc_cov['IID'].tolist()
loc_anc_cov = loc_anc_cov.set_index('IID').transpose()
SNPs.columns = ['rs', 'A1', 'A0']
SNP_list = SNPs['rs'].tolist()
pheno = range(1, len(pheno_names) + 1) #format of phenotype file

print("Creating ancestry-specific dosage files.")
progress_landmarks_ind = np.linspace(0, len(inds), 21, dtype = int).tolist()
num_ind = 0
    
'''
Okay so for whoever bothers to read this
I tried making a large data frame and lists of lists with the dosages
However, this made the addition time for each individual raise incrementally
And that's not good for 12,000 people
So, I write each dosage as a line into a file individually so the addition time is constant
And then read that all back in and add the BIMBAM information
So I swear I tried to do it a fancier way but it's so slow with so many people
'''

if refpop == "HIS":
    IBS_file = open("BIMBAM/IBS" + output + ".txt", "a+")
    NAT_file = open("BIMBAM/NAT" + output + ".txt", "a+")
    YRI_file = open("BIMBAM/YRI" + output + ".txt", "a+")

elif refpop == "AFA":
    CEU_file = open("BIMBAM/CEU" + output + ".txt", "a+")
    YRI_file = open("BIMBAM/YRI" + output + ".txt", "a+")

for ind in inds: #iterate through cols
    ind = str(ind)
    ind_df = loc_anc_cov[[ind]] 
    if refpop == "HIS":
        ind_df['IBS'], ind_df['NAT'], ind_df['YRI'] = ind_df[ind].str.split('\t', 2).str #split each individual's column into 3
        ind_df = ind_df.drop(ind, axis = 1).transpose().applymap(str) #pull from one ancestry each
        
        IBS_list = ind_df.loc['IBS'].tolist() #assemble a BIMBAM file except it's local ancestries
        NAT_list = ind_df.loc['NAT'].tolist()
        YRI_list = ind_df.loc['YRI'].tolist()
       
        IBS_file.write("\t".join(IBS_list) + "\n")
        NAT_file.write("\t".join(NAT_list) + "\n")
        YRI_file.write("\t".join(YRI_list) + "\n")
    
    elif refpop == "AFA":
        ind_df['CEU'], ind_df['YRI'] = ind_df[ind].str.split('\t', 1).str
        ind_df = ind_df.drop(ind, axis = 1).transpose().applymap(str)
        
        CEU_list = ind_df.loc['CEU'].tolist() 
        YRI_list = ind_df.loc['YRI'].tolist()
       
        CEU_file.write("\t".join(CEU_list) + "\n")
        YRI_file.write("\t".join(YRI_list) + "\n")
    num_ind = num_ind + 1
    if num_ind in set(progress_landmarks_ind): #print progress by 5% increments
        progress = progress_landmarks_ind.index(num_ind)
        print("Individual ancestry dosage conversion is " + str(progress * 5) + "% complete.")
    
#so now you have dosages for each ancestry
if refpop == "HIS":
    IBS_file.close()
    NAT_file.close()
    YRI_file.close()
    
    #make into BIMBAM
    IBS_BIMBAM = pd.read_table("BIMBAM/IBS" + output + ".txt", sep = '\t', header = None).transpose()
    NAT_BIMBAM = pd.read_table("BIMBAM/NAT" + output + ".txt", sep = '\t', header = None).transpose()
    YRI_BIMBAM = pd.read_table("BIMBAM/YRI" + output + ".txt", sep = '\t', header = None).transpose()

    #add SNP info
    IBS_BIMBAM = pd.concat([SNPs, IBS_BIMBAM], axis = 1)
    NAT_BIMBAM = pd.concat([SNPs, NAT_BIMBAM], axis = 1)
    YRI_BIMBAM = pd.concat([SNPs, YRI_BIMBAM], axis = 1)
    
    #write to file
    IBS_BIMBAM.to_csv("BIMBAM/IBS" + output + ".txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")
    NAT_BIMBAM.to_csv("BIMBAM/NAT" + output + ".txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")
    YRI_BIMBAM.to_csv("BIMBAM/YRI" + output + ".txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")

elif refpop == "AFA":
    CEU_file.close()
    YRI_file.close()
    
    CEU_BIMBAM = pd.read_table("BIMBAM/CEU" + output + ".txt", sep = '\t', header = None).transpose()
    YRI_BIMBAM = pd.read_table("BIMBAM/YRI" + output + ".txt", sep = '\t', header = None).transpose()

    CEU_BIMBAM = pd.concat([SNPs, CEU_BIMBAM], axis = 1)
    YRI_BIMBAM = pd.concat([SNPs, YRI_BIMBAM], axis = 1)

    CEU_BIMBAM.to_csv("BIMBAM/CEU" + output + ".txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")
    YRI_BIMBAM.to_csv("BIMBAM/YRI" + output + ".txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")

#phenotype loop   
for pheno_num, pheno_name in zip(pheno, pheno_names):
    print("Starting analyses on " + pheno_name + ".")
    if refpop == "HIS":
        for pop in ['NAT', 'IBS', 'YRI']:
            GEMMA_command = "gemma -g BIMBAM/" + pop + output + ".txt.gz -p " + pheno_file + " -n " + str(pheno_num) + anno + " -k " + relatedness + covariates_file + " -lmm 4 -notsnp -o " + output + "_" + pheno_name + "_" + pop
            os.system(GEMMA_command + " >> GEMMA_log.txt")
    elif refpop == "AFA":
        for pop in ['CEU', 'YRI']:
            GEMMA_command = "gemma -g BIMBAM/" + pop + output + ".txt.gz -p " + pheno_file + " -n " + str(pheno_num) + anno + " -k " + relatedness + covariates_file + " -lmm 4 -notsnp -o " + output + "_" + pheno_name + "_" + pop
            os.system(GEMMA_command + " >> GEMMA_log.txt")
    print("Ending analyses on " + pheno_name + ".")

print("Removing extra files.")
os.system("rm -f loc_anc_output/SNPs_" + args.BIMBAM + ".txt") #cause I'm petty and GEMMA is annoying for not directing output
os.system("mv output/* loc_anc_output/")
os.system("rm -r output/")
if refpop == "HIS":
    os.system("rm -f BIMBAM/IBS" + output + ".txt.gz")
    os.system("rm -f BIMBAM/NAT" + output + ".txt.gz")
    os.system("rm -f BIMBAM/YRI" + output + ".txt.gz")
    os.system("rm -f BIMBAM/IBS" + output + ".txt")
    os.system("rm -f BIMBAM/NAT" + output + ".txt")
    os.system("rm -f BIMBAM/YRI" + output + ".txt")
elif refpop == "AFA":
    os.system("rm -f BIMBAM/CEU" + output + ".txt.gz")
    os.system("rm -f BIMBAM/YRI" + output + ".txt.gz")
    os.system("rm -f BIMBAM/CEU" + output + ".txt")
    os.system("rm -f BIMBAM/YRI" + output + ".txt")
print("Analyses in all phenotypes is complete. Have a nice day :)!")
