#backward elimination of significant genes to determine which ones are independent
#if anyone tries to debug this... I'm sorry...
library(argparse)
library(data.table)
library(MASS)
library(dplyr)
library(tidyverse)
"%&%" = function(a,b) paste(a,b,sep="")

parser <- ArgumentParser()
parser$add_argument("--desc", help = "Perform a colocalization between GWAS results and eQTL data. This uses GTEx V6, so if you want V7, download those data and change paths accordingly.")
parser$add_argument("--sig_gene", help = "Path to significant gene file")
parser$add_argument("--pheno", help = "Path to phenotype file WITH IDS")
parser$add_argument("--pred_exp_prefix", help = "Prefix of PrediXcan-GEMMA results")
parser$add_argument("--pheno_name", help = "Name of pheno to test")
#parser$add_argument("--sample_size", help = "Sample size of GWAS cohort")
args <- parser$parse_args()

sig_gene <- args$sig_gene
pheno <- args$pheno
pred_exp_prefix <- args$pred_exp_prefix
pheno_name <- args$pheno_name

'
sig_gene <- "output/AMR_sig_genes.txt"
pheno <- "pheno_wIID.txt"
pred_exp_prefix <- "AMR_"
pheno_name <- "pheno"
"%&%" = function(a,b) paste(a,b,sep="")
'

sig_gene <- fread(sig_gene, header = T)
pheno <- fread(pheno, header = T)
keep_pheno <- c("FID", "IID", pheno_name)
pheno <- pheno %>% dplyr::select(keep_pheno)
colnames(pheno) <- c("FID", "IID", "pheno")

#convert tissues to tissue abbreviations
tissues <- fread("/home/angela/tiss_abbreviations_no_TW.txt")
tissues_abbv_df <- tissues
tiss_abbv <- tissues$tiss
tissues <- tissues$tissue #it would be nice if R had dicts

#prep data
BP_Chrome <- read.table("/home/angela/px_yri_chol/PrediXcan/BP_Chrome.txt", header = T, stringsAsFactors = F)
BP_Chrome$gene <- gsub("\\..*", "", BP_Chrome$gene)
back_elim <- pheno[, 2] #extract only IID
sig_gene$rs <- gsub("\\..*", "", sig_gene$rs)
sig_gene <- left_join(sig_gene, BP_Chrome, by = c("rs" = "gene"))
sig_gene$gene_tiss_combos <- paste(sig_gene$tissue, sig_gene$gene_name, sep = "_")
gene_tiss_combos <- sig_gene$gene_tiss_combos
gene_tiss_combos <- c("IID", gene_tiss_combos)
sig_gene$gene_tiss_combos <- NULL

#get table of all predicted expressions
for(tiss_index in 1:length(tissues)){
  print("Running analyses on " %&% tiss_abbv[tiss_index] %&% ".")
  tiss <- fread('pred_exp/' %&% tissues[tiss_index] %&% '_predicted_expression.txt', header = T, nThread = 30) #load in predicted expression file
  tiss$FID <- NULL #remove FID
  
  #convert current gene IDs to gene names
  genes_in_pred_exp <- as.data.frame(colnames(tiss)[2:length(colnames(tiss))])
  colnames(genes_in_pred_exp) <- "gene"
  genes_in_pred_exp$gene <- gsub("\\..*", "", genes_in_pred_exp$gene)
  genes_in_pred_exp <- left_join(genes_in_pred_exp, BP_Chrome, by = "gene")
  genes_in_pred_exp <- genes_in_pred_exp %>% dplyr::select(gene, gene_name)
  if(tissues[tiss_index] == "Cells_EBV-transformed_lymphocytes"){
    genes_in_pred_exp$gene_name <- paste("Cells_EBV_transformed_lymphocytes", genes_in_pred_exp$gene_name, sep = "_") #the dash causes downstream problems
  }else{
    genes_in_pred_exp$gene_name <- paste(tissues[tiss_index], genes_in_pred_exp$gene_name, sep = "_")
  }
  genes_in_pred_exp <- genes_in_pred_exp$gene_name
  genes_in_pred_exp <- c("IID", genes_in_pred_exp)
  colnames(tiss) <- genes_in_pred_exp
  tiss <- subset(tiss, select = which(!duplicated(names(tiss))))  #remove NA's
  tiss_gene_in_tiss <- intersect(genes_in_pred_exp, gene_tiss_combos)
  for_back_elim <- tiss %>% dplyr::select(tiss_gene_in_tiss) #select IIDs and intersection of tiss-gene combos
  back_elim <- left_join(back_elim, for_back_elim, by = "IID") #add to list to back-elim
}

back_elim <- left_join(back_elim, pheno, by = "IID")

#make clusters of genes
gene_clusters <- sig_gene %>% dplyr::select(CHR, rs, tissue)
colnames(gene_clusters) <- c("chr", "gene", "tissue")
gene_clusters$gene <- gsub("\\..*", "", gene_clusters$gene) #remove all after .
gene_clusters <- gene_clusters[order(gene_clusters$chr),] 
gene_clusters$first_occur <- paste(gene_clusters$chr, gene_clusters$gene, sep = "_") #instead of checking across three columns for unique, just check one
gene_clusters <- gene_clusters[match(unique(gene_clusters$gene), gene_clusters$gene),]
gene_clusters$first_occur <- NULL
gene_clusters <- aggregate(data = gene_clusters, gene ~ chr, FUN=paste) #https://stackoverflow.com/questions/38125125/combine-rows-which-have-same-value-in-two-columns-r

all_tiss_gene_results <- data.frame(Estimate = numeric(), SE = numeric(), T = numeric(), P = numeric(), tiss_gene = character(), stringsAsFactors = F) 
back_elim_results <- data.frame(Estimate = numeric(), SE = numeric(), T = numeric(), P = numeric(), tiss_genes = character(), stringsAsFactors = F) 

for(chr in 1:nrow(gene_clusters)){
  #pheno_chr <- (gene_clusters$pheno_chr)[1]
  chr_num <- as.integer(gene_clusters[chr, 1])
  genes_to_test <- unlist(gene_clusters[chr, 2], recursive = F) #https://stackoverflow.com/questions/11048208/converting-a-list-of-lists-in-a-single-list
  tissues_to_test <- subset(sig_gene, rs %in% genes_to_test & CHR == as.numeric(chr_num))
  tissue_gene <- paste(tissues_to_test$tissue, tissues_to_test$gene, sep = "_")
  tissue_gene <- c("IID", tissue_gene)
  #there was probably a MUCH less convoluted way to do this but here we are
  
  print("Started making models for chr. " %&% chr %&% ".")
  tiss_gene_to_keep <- c(intersect(colnames(back_elim), tissue_gene))
  back_elim_cluster <- back_elim %>% dplyr::select(tiss_gene_to_keep)
  back_elim_cluster <- back_elim_cluster[complete.cases(back_elim_cluster),]
  back_elim_cluster <- cbind(pheno, back_elim_cluster)
  back_elim_cluster$FID <- NULL
  back_elim_cluster$IID <- NULL
  back_elim_cluster$IID <- NULL #why are there two of these
  predictor_genes <- colnames(back_elim_cluster)[2:length(colnames(back_elim_cluster))] #https://stackoverflow.com/questions/5251507/how-to-succinctly-write-a-formula-with-many-variables-from-a-data-frame
  fmla <- paste("pheno ~ ", paste(predictor_genes, collapse= "+"))
  fmla <- as.formula(fmla)
  
  #run full model
  all_tiss_gene <- lm(fmla, data = back_elim_cluster) 
  all_tiss_gene_info <- data.frame(summary(all_tiss_gene)[4])
  colnames(all_tiss_gene_info) <- c("Estimate", "SE", "T", "P")
  all_tiss_gene_info$chr <- chr_num
  all_tiss_gene_info$tiss_gene <- rownames(all_tiss_gene_info)
  all_tiss_gene_results <- rbind(all_tiss_gene_results, all_tiss_gene_info)
  print("Finished making full model for " %&% chr_num %&% ".")
  
  #using stepAIC, is there a better option somewhere?
  back_elim_complete <- stepAIC(all_tiss_gene, direction = "backward", trace = FALSE) #perform backward analysis on full model
  back_elim_complete_info <- data.frame(summary(back_elim_complete)[4])
  colnames(back_elim_complete_info) <- c("Estimate", "SE", "T", "P")
  back_elim_complete_info$chr <- chr_num
  back_elim_complete_info$tiss_gene <- rownames(back_elim_complete_info)
  back_elim_results <- rbind(back_elim_results, back_elim_complete_info)
  print("Finished making backward-eliminated model for " %&% chr_num %&% ".")
}

gene_genename <- sig_gene %>% dplyr::select(rs, gene_name)
back_elim_results[back_elim_results == "(Intercept)"] <- NA
back_elim_results <- back_elim_results[complete.cases(back_elim_results),]
back_elim_results <- back_elim_results %>% separate(tiss_gene, into = c("tiss", "gene"), sep="_(?=[^_]+$)") #https://stackoverflow.com/questions/50518137/separate-a-column-into-2-columns-at-the-last-underscore-in-r
back_elim_results <- back_elim_results %>% dplyr::select(chr, tiss, gene, P)
colnames(back_elim_results)[3] <- "gene_name"
back_elim_results <- left_join(back_elim_results, BP_Chrome, by = "gene_name")
back_elim_results <- back_elim_results[order(back_elim_results$CHR, back_elim_results$BP, back_elim_results$gene_name),]
back_elim_results <- back_elim_results %>% dplyr::select(chr, BP, gene_name, tiss, P)
fwrite(back_elim_results, "back_elim_results.csv", row.names = F, col.names = T, sep = ",", quote = F, na = NA)
print("Results are in back_elim_results.csv")