#runs a loop through all GTEx models and all MESA models to run the --predicted expression flag in PrediXcan
#if you are using a new verison of GTEx, change the tissues list
#by Angela Andaleon (aandaleon@luc.edu)
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--desc", type = str, action = "store", dest = "desc", required = False, help = "Runs a loop through all GTEx models and all MESA models through PrediXcan to make predicted expression files.")
parser.add_argument("--PrediXcan_path", type = str, action = "store", dest = "PrediXcan_path", required = False, default = "/usr/local/bin/PrediXcan.py", help = "Path to PrediXcan executable. Default = /usr/local/bin/PrediXcan.py")
parser.add_argument("--dosages_path", type = str, action = "store", dest = "dosages_path", required = False, default = "dosages/", help = "Path to dosages. Default = dosages/")
parser.add_argument("--MESA_prefix", type = str, action = "store", dest = "MESA_prefix", required = False, default = "/home/lauren/files_for_revisions_plosgen/en_v7/dbs/", help = "Path to prefix of all MESA models. Default = /home/lauren/files_for_revisions_plosgen/en_v7/dbs/")
parser.add_argument("--MESA_suffix", type = str, action = "store", dest = "MESA_suffix", required = False, default = "_imputed_10_peer_3_pcs_2.db", help = "Path to suffix of all MESA models. Default = _imputed_10_peer_3.db")
parser.add_argument("--GTEx_prefix", type = str, action = "store", dest = "GTEx_prefix", required = False, default = "/home/wheelerlab3/Data/PrediXcan_db/GTEx-V6p-HapMap-2016-09-08/TW_", help = "Path to prefix of all GTEx models. Default = /home/wheelerlab3/Data/PrediXcan_db/GTEx-V6p-HapMap-2016-09-08/TW_")
parser.add_argument("--GTEx_suffix", type = str, action = "store", dest = "GTEx_suffix", required = False, default = "_0.5.db", help = "Path to suffix of all GTEx models. Default = _0.5.db")
parser.add_argument("--output_prefix", type = str, action = "store", dest = "output_prefix", required = False, default = "pred_exp/", help = "Path to prefix PrediXcan predicted expression output. Default = pred_exp/")
args = parser.parse_args()

pops = ["AFA", "AFHI", "ALL", "CAU", "HIS"] #all of Lauren's MESA models
tissues = ["Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood"] #this is set to GTEx V6 because I am familiar with these. Change as necessary.

PrediXcan_path = args.PrediXcan_path
dosages_path = args.dosages_path
MESA_prefix = args.MESA_prefix
MESA_suffix = args.MESA_suffix
GTEx_prefix = args.GTEx_prefix
GTEx_suffix = args.GTEx_suffix
output_prefix = args.output_prefix

print("Starting creation of predicted expressions.")
os.system("mkdir -p " + output_prefix) #make output folder if not already made
for pop in pops: #iterate through GTEx
    PX_command = PrediXcan_path + " --predict --dosages " + dosages_path + " --samples samples.txt --weights " + MESA_prefix + pop + MESA_suffix + " --output_prefix " + output_prefix + pop
    os.system(PX_command)
    print("Completed with " + pop + ".")
for tissue in tissues: #iterate through MESA
    PX_command = PrediXcan_path + " --predict --dosages " + dosages_path + " --samples samples.txt --weights " + GTEx_prefix + tissue + GTEx_suffix + " --output_prefix " + output_prefix + tissue
    os.system(PX_command)
    print("Completed with " + tissue + ".")
print("Completed with making predicted expressions.")
