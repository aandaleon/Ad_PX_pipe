#convert predicted expressions into GEMMA-format pseudo-genotypes
#by Angela Andaleon (aandaleon@luc.edu)
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("--desc", type = str, action = "store", dest = "desc", required = False, help = "Runs a loop through all GTEx models and all MESA models through PrediXcan to make predicted expression files.")
#parser.add_argument("--pred_exp_prefix", type = str, action = "store", dest = "pred_exp_prefix", required = False, default = "pred_exp/", help = "Path to prefix of PrediXcan predicted expression input. Default = pred_exp/")
#parser.add_argument("--output_prefix", type = str, action = "store", dest = "output_prefix", required = False, default = "pred_exp_GEMMA/", help = "Path to GEMMA-stype predicted expression output. Default = pred_exp_GEMMA/")
args = parser.parse_args()

tissues = ["AFA", "AFHI", "ALL", "CAU", "HIS", "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood"]

#pred_exp_prefix = args.pred_exp_prefix
#output_prefix = args.output_prefix
pred_exp_prefix = "pred_exp/"
output_prefix = "pred_exp_GEMMA/"

'''
pred_exp_prefix = "pred_exp/"
output_prefix = "pred_exp_GEMMA/"
'''

os.system("mkdir -p " + output_prefix)
for tissue in tissues:
    PX_output = pd.read_csv(pred_exp_prefix + tissue + "_predicted_expression.txt", sep = "\t", engine = "python").transpose()
    GEMMA_input = PX_output.iloc[2:] #remove first two rows (FIDs and IIDs)
    GEMMA_input.insert(0, "A1", "NA")
    GEMMA_input.insert(1, "A0", "NA")
    GEMMA_input.insert(0, "gene", GEMMA_input.index)
    GEMMA_input.to_csv(output_prefix + tissue + ".txt", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3)
    print("Completed converting " + tissue + ".")




