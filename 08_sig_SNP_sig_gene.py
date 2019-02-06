#find significant SNPs and genes
#by Angela Andaleon (aandaleon@luc.edu)
import argparse
import numpy as np
tissues = ["AFA", "AFHI", "ALL", "CAU", "HIS", "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood"] #GTEx V6

parser = argparse.ArgumentParser()
parser.add_argument("--desc", type = str, action = "store", dest = "desc", required = False, help = "Finds significant SNPs and genes from GEMMA output.")
parser.add_argument("--SNP_sig", type = str, action = "store", dest = "SNP_sig", required = False, default = "5e-8", help = "Significance threshold for SNPs. Default = 5e-8")
parser.add_argument("--gene_sig", type = str, action = "store", dest = "gene_sig", required = False, default = "9.654e-6", help = "Significance threshold for genes. Default = 9.654-6")
#parser.add_argument("--input_prefix", type = str, action = "store", dest = "input_prefix", required = True, help = "Prefix for input, not including 'output/'.")
parser.add_argument("--pheno_names", type = str, action = "store", dest = "pheno_names", required = False, default = "pheno_names.txt", help = "File containing pheno names. Default = 'pheno_names.txt'.")
args = parser.parse_args()

SNP_sig = float(args.SNP_sig)
gene_sig = float(args.gene_sig)
#input_prefix = args.input_prefix
input_prefix = ""
pheno_names = list(np.loadtxt(args.pheno_names, dtype = str))

'''
SNP_sig = 5e-8
gene_sig = 9.654e-6
input_prefix = "AMR"
'''
 
for pheno_name in pheno_names:
    sig_snp = open("output/" + pheno_name + "_sig_snps.txt", "w")
    sig_snp.write("chr\trs\tps\tn_miss\tallele1\tallele0\taf\tbeta\tse\tl_remle\tl_mle\tp_wald\tp_lrt\tp_score\tpheno\n")
    sig_gene = open("output/" + pheno_name + "_sig_genes.txt", "w")
    sig_gene.write("chr\trs\tps\tn_miss\tallele1\tallele0\taf\tbeta\tse\tl_remle\tl_mle\tp_wald\tp_lrt\tp_score\ttissue\tpheno\n")
    for chrom in range(1, 23): #find sig snp
        for line in open("output/" + input_prefix + "" + pheno_name + "_chr" + str(chrom) + ".assoc.txt", "r"):
            arr = line.strip().split()
            (chr, rs, ps, n_miss, allele1, allele0, af, beta, se, l_remle, l_mle, p_wald, p_lrt, p_score) = arr[0:14]
            if line.startswith("chr"): #the header line messes things up so just skip it
                continue
            elif float(p_wald) < SNP_sig: #only sig. SNPs
                is_sig_snp = (chr + "\t" + rs + "\t" + ps + "\t" + n_miss + "\t" + allele1 + "\t" + allele0 + "\t" + af + "\t" + beta + "\t" + se + "\t" + l_remle + "\t" + l_mle + "\t" + p_wald + "\t" + p_lrt + "\t" + p_score + "\t" + pheno_name + "\n")
                sig_snp.write(is_sig_snp)

    for tissue in tissues: #find sig gene
        for line in open("output/" + input_prefix + "" + pheno_name + "_" + tissue + ".assoc.txt", "r"):
            arr = line.strip().split()
            (chr, rs, ps, n_miss, allele1, allele0, af, beta, se, l_remle, l_mle, p_wald, p_lrt, p_score) = arr[0:14]
            if line.startswith("chr"):
                continue
            elif float(p_wald) < gene_sig: #only sig. genes
                is_sig_gene = (chr + "\t" + rs + "\t" + ps + "\t" + n_miss + "\t" + allele1 + "\t" + allele0 + "\t" + af + "\t" + beta + "\t" + se + "\t" + l_remle + "\t" + l_mle + "\t" + p_wald + "\t" + p_lrt + "\t" + p_score + "\t" + tissue + "\t" + pheno_name + "\n")
                sig_gene.write(is_sig_gene)
                
sig_snp.close()
sig_gene.close()
print("Completed writing sig_SNPs.txt and sig_gene.txt for all phenotypes.")
