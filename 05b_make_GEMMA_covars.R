#Makes a new covariance file for GEMMA based on known covariates and KING PCs.
#by Angela Andaleon (aandaleon@luc.edu)
library(argparse)
library(data.table)
library(dplyr)
"%&%" = function(a, b) paste(a, b, sep = "")

parser <- ArgumentParser()
parser$add_argument("--desc", help = "This script makes covariance file from known covariates and KING PCs.")
parser$add_argument("--covar", help = "Covariance file (w/o IDs) produced in step 0")
#parser$add_argument("--pcs_file", help = "Principal components file produced in step 2. Default = kingpc.ped", default = "kingpc.ped")
parser$add_argument("--pcs_num", help = "Number of principal components to include in covariates file. Default = 5.", default = "5")
parser$add_argument("--output", help = "Name of output file. Default = GEMMA_covars.txt", default = "GEMMA_covars.txt")
args <- parser$parse_args()

covar_file_name <- args$covar
#pcs_file_name <- args$pcs_file
pcs_file_name <- "kingpc.ped"
pcs_num <- as.integer(args$pcs_num)
output <- args$output

'
covar_file_name <- "covar_woIID.txt"
pcs_file_name <- "kingpc.ped"
pcs_num <- 5
output <- "GEMMA_covars.txt"
'

KING_pcs <- fread(pcs_file_name)
old_covars <- fread(covar_file_name)
if(pcs_num == 0){
  new_covars <- old_covars #no change in covars
}else{
  pcs_to_keep <- KING_pcs[, 7:(6 + pcs_num)] #keep specified number of PCs from kingpc.ped
  new_covars <- cbind(old_covars, pcs_to_keep) #make old covariates and pcs into one file
}
fwrite(new_covars, output, sep = "\t", row.names = F, col.names = F, na = "NA", quote = F) #write to file for direct GEMMA input

print("Completed making new covariate file.")
