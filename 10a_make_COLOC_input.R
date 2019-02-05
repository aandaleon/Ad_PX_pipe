#converts GEMMA-style GWAS results and standard eQTL format to that suitable for COLOC
#the structure of this is weird b/c GTEx eQTL data is not separated by chr while MESA is
#by Angela Andaleon (aandaleon@luc.edu)
library(argparse)
library(data.table)
library(dplyr)
library(R.utils) #I forgot why I had this
"%&%" = function(a, b) paste(a, b, sep = "")

parser <- ArgumentParser()
parser$add_argument("--desc", help = "Perform a colocalization between GWAS results and eQTL data. This uses GTEx V6, so if you want V7, download those data and change paths accordingly.")
parser$add_argument("--ma_prefix", help = "Prefix of .ma file (the input of GCTA-COJO)")
parser$add_argument("--GWAS_prefix", help = "Prefix of GWAS results")
parser$add_argument("--sample_size", help = "Sample size of GWAS cohort")
parser$add_argument("--pheno_names", default = "pheno_names.txt", help = "File of pheno names. Defauly = 'pheno_names.txt'")
args <- parser$parse_args()

ma_prefix <- args$ma_prefix
GWAS_prefix <- args$GWAS_prefix
sample_size <- args$sample_size
pheno_names <- args$pheno_names

'
ma_prefix <- "AMR"
GWAS_prefix <- "AMR"
pheno_names <- "pheno_names.txt"
sample_size <- "347"
'

print("You're gonna want to leave this running overnight.")
pheno_names <- fread(pheno_names, header = F)$V1
system("mkdir -p COLOC_input/")

#GTEx frq
frq <- fread("/home/angela/Ad_PX_pipe_data/GTEx_WHLBLD.frq") #we're just gonna use WHLBLD dosages for everyone
bim <- fread("/home/angela/Ad_PX_pipe_data/GTEx_WHLBLD.bim") #for adding cpos
bim$cpos <- paste(bim$V1, bim$V4, sep = "_")
bim <- bim %>% dplyr::select(V2, cpos)
colnames(bim) <- c("SNP", "cpos") #GTEx just has to be special and use cpos instead of rsids
frq <- left_join(frq, bim, by = "SNP")
GTEx_frq <- frq %>% dplyr::select(SNP, MAF, cpos)

#general setup data
chrs <- c(1:22)
tissues <- c("AFA", "AFHI", "ALL", "CAU", "HIS", "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")             
tissues_sample_size <- c(233, 578, 352, 585, 1163, 298, 298, 126, 197, 118, 285, 72, 72, 89, 103, 96, 92, 81, 81, 81, 82, 183, 114, 272, 124, 124, 127, 241, 218, 159, 159, 97, 278, 361, 256, 85, 149, 87, 87, 196, 302, 77, 89, 170, 157, 278, 278, 278, 338)#R doesn't have dicts so we're doing it a slightly more ratchet way

for(tissue in 1:length(tissues)){ #read in tissue's .frq file for MAF
  if(tissue <= 5){ #if MESA
    frq <- fread("/home/angela/Ad_PX_pipe_data/MESA_frq/" %&% tissues[tissue] %&% ".frq", nThread = 30)
    frq <- frq %>% dplyr::select(SNP, MAF)
  }else{
    frq <- GTEx_frq
  }
  
  for(pheno_name in pheno_names){ #how do I structure this so I don't have to load GTEx multiple times  
    GWAS_SNPs <- fread(ma_prefix %&% "_" %&% pheno_name %&% ".ma", header = T)$rs
    eQTL_write <- data.frame(gene_id = character(), variant_id = character(), maf = numeric(), pval_nominal = numeric(), slope = numeric(), slope_se = numeric(), stringsAsFactors = F) 
    GWAS_write <- data.frame(panel_variant_id = character(), effect_size = numeric(), standard_error = numeric(), frequency = numeric(), sample_size = numeric(), stringsAsFactors = F) 
    
    if(tissue > 5){ #GTEx is not separated by chr
      meqtl <- fread("/home/angela/Ad_PX_pipe_data/GTEx_eQTL/" %&% tissues[tissue] %&% "_Analysis.v6p.all_snpgene_pairs.txt.gz", nThread = 30) #read in matrix eQTL results

      #why can't you just be in normal cpos format (force GTEx ids to be in normal cpos)
      meqtl$cpos <- gsub("^([^_]*_[^_]*)_.*$", "\\1", meqtl$variant_id) #https://stackoverflow.com/questions/7449564/regex-return-all-before-the-second-occurrence
      meQTL_for_COLOC <- left_join(meqtl, frq, by = "cpos") #add freq to COLOC input
      meQTL_for_COLOC <- meQTL_for_COLOC %>% dplyr::select(gene_id, SNP, MAF, pval_nominal, slope, slope_se) #subset to COLOC input
      meQTL_for_COLOC <- subset(meQTL_for_COLOC, SNP %in% GWAS_SNPs) #reduce to just the necessary SNPs
      meQTL_for_COLOC$n_samples <- tissues_sample_size[tissue]
      colnames(meQTL_for_COLOC) <- c("gene_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
      meQTL_for_COLOC <- meQTL_for_COLOC[complete.cases(meQTL_for_COLOC),]
      eQTL_write <- rbind(eQTL_write, meQTL_for_COLOC)
    }

    for(chr in chrs){ #yes triple loops are ratchet
      GEMMA_result <- fread("output/" %&% GWAS_prefix %&% "_" %&% pheno_name %&% "_chr" %&% chr %&% ".assoc.txt", header = F, nThread = 30)
      colnames(GEMMA_result) <- c("chr", "rs", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score")
      GEMMA_result$chr_pos <- paste(gsub("chr", "", GEMMA_result$chr), GEMMA_result$ps, sep = ":")
      GEMMA_for_COLOC <- GEMMA_result %>% dplyr::select(rs, beta, se, af) #subset to COLOC input
      GEMMA_for_COLOC$sample_size <- sample_size
      colnames(GEMMA_for_COLOC) <- c("panel_variant_id", "effect_size", "standard_error", "frequency", "sample_size")
      GEMMA_for_COLOC <- GEMMA_for_COLOC[complete.cases(GEMMA_for_COLOC),] #COLOC does not like missing values
      GEMMA_for_COLOC <- GEMMA_for_COLOC[2:nrow(GEMMA_for_COLOC),]
    
      if(tissue <= 5){ #MESA meQTL is by chr, then pop
        system("zcat -f /home/lauren/files_for_revisions_plosgen/meqtl_results/MESA/" %&% tissues[tissue] %&% "_Nk_10_PFs_chr" %&% chr %&% "pcs_3.meqtl.cis.* > COLOC_input/meQTL_input.txt") #fread doesn't seem to like wildcards so we're gonna do this the ugly way
        meqtl <- fread("COLOC_input/meQTL_input.txt", nThread = 40) #read in matrix eQTL results
        meqtl <- subset(meqtl, snps %in% GWAS_SNPs)
        meqtl$se <- meqtl$beta / meqtl$statistic #make your own standard error since it's not in the meQTL output
        meqtl$n_samples <- tissues_sample_size[tissue]
        meQTL_for_COLOC <- left_join(meqtl, frq, by = c("snps" = "SNP")) #add freq to COLOC input
        meQTL_for_COLOC <- meQTL_for_COLOC %>% dplyr::select(gene, snps, MAF, pvalue, beta, se) #subset to COLOC input
        colnames(meQTL_for_COLOC) <- c("gene_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
        meQTL_for_COLOC <- meQTL_for_COLOC[complete.cases(meQTL_for_COLOC),]
        eQTL_write <- rbind(eQTL_write, meQTL_for_COLOC)
      }
      GWAS_write <- rbind(GWAS_write, GEMMA_for_COLOC)
    }
  
    snps_in_both <- intersect(GWAS_write$panel_variant_id, eQTL_write$variant_id) #is there a better way to do this? Probably. Do I feel like figuring it out? Nah.
    snps_in_all <- intersect(snps_in_both, GWAS_SNPs) #only keep overlapping SNPs
    GWAS_write <- subset(GWAS_write, panel_variant_id %in% snps_in_all)
    eQTL_write <- subset(eQTL_write, variant_id %in% snps_in_all)
    eQTL_write <- eQTL_write[order(eQTL_write$gene_id),] #results are weird when not ordered
  
    fwrite(eQTL_write, "COLOC_input/" %&% GWAS_prefix %&% "_" %&% pheno_name %&% "_eQTL_" %&% tissues[tissue] %&% ".txt", quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)
    gzip("COLOC_input/" %&% GWAS_prefix %&% "_" %&% pheno_name %&% "_eQTL_" %&% tissues[tissue] %&% ".txt", destname = "COLOC_input/" %&% GWAS_prefix %&% "_" %&% pheno_name %&% "_eQTL_" %&% tissues[tissue] %&% ".txt.gz", overwrite = T) #script may only take .gz values so can't hurt to be too careful
    fwrite(GWAS_write, "COLOC_input/" %&% GWAS_prefix %&% "_" %&% pheno_name %&% "_GWAS_" %&% tissues[tissue] %&% ".txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
    gzip("COLOC_input/" %&% GWAS_prefix %&% "_" %&% pheno_name %&% "_GWAS_" %&% tissues[tissue] %&% ".txt", "COLOC_input/" %&% GWAS_prefix %&% "_" %&% pheno_name %&% "_GWAS_" %&% tissues[tissue] %&% ".txt.gz", overwrite = T)
    print("Completed with " %&% tissues[tissue] %&% ", for " %&% pheno_name %&%".")
  }
}
print("Completed making input for COLOC.")
