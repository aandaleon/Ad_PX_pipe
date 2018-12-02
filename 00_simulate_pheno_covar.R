library(argparse)
library(data.table)
library(dplyr)
"%&%" = function(a, b) paste(a, b, sep = "")

parser <- ArgumentParser()
parser$add_argument("--desc", help = "This script simulates phenotypes and covariates to use with the sample data set, AMR.")
parser$add_argument("--bfile", help = "PLINK-format genotypes")
args <- parser$parse_args()

bfile_path <- args$bfile
#bfile_path <- "AMR"

fam <- fread(bfile_path %&% ".fam")
fam <- fam %>% select(V1, V2)
colnames(fam) <- c("FID", "IID")

#generate random pheno and covar
set.seed(1000)
fam$pheno <- rnorm(mean = 0, sd = 1, n = nrow(fam))
fam$covar <- rbinom(nrow(fam), 1, 0.25)

#write pheno file w/ and w/o IDs
pheno <- fam %>% select(FID, IID, pheno)
fwrite(pheno, "pheno_wIID.txt", row.names = F, col.names = T, na = "NA", quote = F, sep = "\t")
pheno$FID <- NULL
pheno$IID <- NULL
fwrite(pheno, "pheno_woIID.txt", row.names = F, col.names = F, na = "NA", quote = F, sep = "\t")

#write covar file w/ and w/o IDs
covar <- fam %>% select(FID, IID, covar)
fwrite(covar, "covar_wIID.txt", row.names = F, col.names = T, na = "NA", quote = F, sep = "\t")
covar$FID <- 1 #GEMMA requires this as an intercept term
covar$IID <- NULL
fwrite(covar, "covar_woIID.txt", row.names = F, col.names = F, na = "NA", quote = F, sep = "\t")
