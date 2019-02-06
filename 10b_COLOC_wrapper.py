#Wrapper around COLOC translated from Bash because Python flags are superior
#@author: Angela Andaleon (aandaleon@luc.edu)

#make GWAS output into GCTA-COJO format
import argparse
import numpy as np
import os
#MESA and GTEx V6; this could be a dict but I already wrote it like this in bash soooooo
tissues = ["AFA", "AFHI", "ALL", "CAU", "HIS", "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood"]
tissues_sample_size= ["233", "578", "352", "585", "1163", "298", "298", "126", "197", "118", "285", "72", "72", "89", "103", "96", "92", "81", "81", "81", "82", "183", "114", "272", "124", "124", "127", "241", "218", "159", "159", "97", "278", "361", "256", "85", "149", "87", "87", "196", "302", "77", "89", "170", "157", "278", "278", "278", "338"]

parser = argparse.ArgumentParser() #I could do this in bash but flag passing in bash is terrible
parser.add_argument("--desc", type = str, action = "store", dest = "desc", required = False, help = "make GWAS output into GCTA-COJO format")
#parser.add_argument("--prefix", type = str, action = "store", dest = "prefix", required = True, help = "Prefix of input files (not including 'COLOC_input/').")
parser.add_argument("--GWAS_sample_size", type = str, action = "store", dest = "GWAS_sample_size", required = True, help = "Sample size of GWAS.")
parser.add_argument("--pheno_names", type = str, action = "store", dest = "pheno_names", required = False, default = "pheno_names.txt", help = "File containing pheno names. Default = 'pheno_names.txt'.")
args = parser.parse_args()

os.system("if cd summary-gwas-imputation; then git pull; else git clone https://github.com/hakyimlab/summary-gwas-imputation.git; fi") #contains COLOC wrapper, https://stackoverflow.com/questions/15602059/git-shortcut-to-pull-with-clone-if-no-local-there-yet
os.system("mkdir -p COLOC_results/")

#prefix = args.prefix
GWAS_sample_size = args.GWAS_sample_size
pheno_names = list(np.loadtxt(args.pheno_names, dtype = str))

'''
prefix = "AMR"
GWAS_sample_size = "347"
pheno_names = list(np.loadtxt("pheno_names.txt", dtype = str))
'''

print("Leave this script running overnight as well.")
for pheno_name in pheno_names:
    for tissue, tissue_sample_size in zip(tissues, tissues_sample_size):
        run_COLOC = "python3 summary-gwas-imputation/src/run_coloc.py -gwas COLOC_input/" + pheno_name + "_GWAS_" + tissue + ".txt.gz -gwas_sample_size " + GWAS_sample_size + " -eqtl COLOC_input/" + pheno_name + "_eQTL_" + tissue + ".txt.gz -eqtl_sample_size " + tissue_sample_size + " -output COLOC_results/" + pheno_name + "_" + tissue + ".txt.gz" #COLOC prints out too much so throw it into a file
        os.system(run_COLOC + " > /dev/null 2>&1") #this script is a noisy boi
        print("Completed running COLOC in " + tissue + " for " + pheno_name + ".")
print("Completed running COLOC. Have a good day :).")
