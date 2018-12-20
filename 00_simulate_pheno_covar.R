library(argparse)
library(data.table)
library(dplyr)
"%&%" = function(a, b) paste(a, b, sep = "")

parser <- ArgumentParser()
parser$add_argument("--desc", help = "This script simulates phenotypes and covariates to use with the sample data set, AMR. Note: the first column of the covariates file must be a column of 1's.")
parser$add_argument("--bfile", help = "PLINK-format genotypes")
args <- parser$parse_args()

bfile_path <- args$bfile
#bfile_path <- "AMR"

fam <- fread(bfile_path %&% ".fam")
fam <- fam %>% select(V1, V2) #keep only FID and IID columns
colnames(fam) <- c("FID", "IID")

#generate random pheno and covar
set.seed(1000)
fam$pheno <- rnorm(mean = 0, sd = 1, n = nrow(fam)) #produce a random continuous variable (ex. total cholesterol levels)
fam$covar <- rbinom(nrow(fam), 1, 0.25) #produce a random discrete binary variable for covariates (ex. not on lipid medication/on lipid medication)

#write pheno file w/ and w/o IDs
pheno <- fam %>% select(FID, IID, pheno) #give each individual a random phenotype
fwrite(pheno, "pheno_wIID.txt", row.names = F, col.names = T, na = "NA", quote = F, sep = "\t") #for readability purposes
pheno$FID <- NULL #GEMMA requires only columns of phenotypes and doesn't take IDs
pheno$IID <- NULL
fwrite(pheno, "pheno_woIID.txt", row.names = F, col.names = F, na = "NA", quote = F, sep = "\t") #for GEMMA purposes, which doesn't like column names

#write covar file w/ and w/o IDs
covar <- fam %>% select(FID, IID, covar)
fwrite(covar, "covar_wIID.txt", row.names = F, col.names = T, na = "NA", quote = F, sep = "\t") #for readability purposes
covar$FID <- 1 #GEMMA requires this as an intercept term
covar$IID <- NULL #GEMMA requires only columns of covariates and doesn't take IDs
fwrite(covar, "covar_woIID.txt", row.names = F, col.names = F, na = "NA", quote = F, sep = "\t") #for GEMMA purposes, which doesn't like column names
