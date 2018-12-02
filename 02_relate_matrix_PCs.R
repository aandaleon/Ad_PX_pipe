library(argparse)
library(data.table)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(OmicKriging)
"%&%" = function(a, b) paste(a, b, sep = "")

parser <- ArgumentParser()
parser$add_argument("--desc", help = "This script takes in the path to your PLINK-format genotype files and calculates  principal components and a relationship matrix for use in GEMMA. You can calculate these with a subset of SNPs (ex. chr22) for speed.")
parser$add_argument("--king", help = "Path to KING executable. Default = /home/angela/px_his_chol/KING_LAPACK/king-offline", default = "/home/angela/px_his_chol/KING_LAPACK/king-offline")
parser$add_argument("--bfile", help = "PLINK-format genotypes")
args <- parser$parse_args()

king_path <- args$king
bfile_path <- args$bfile

'
king_path <- "/home/angela/px_his_chol/KING_LAPACK/king-offline"
bfile_path <- "AMR_chr22"
'

print("Calculating relationship matrix in GEMMA and GCTA format.")
#run KING --kinship
KING_kinship <- king_path %&% " -b " %&% bfile_path %&% ".bed --kinship --prefix relatedness"
system(KING_kinship)

snpgdsBED2GDS(bed.fn = bfile_path %&% ".bed", bim.fn = bfile_path %&% ".bim", fam.fn = bfile_path %&% ".fam", out.gdsfn = bfile_path %&% ".gds", family = T)
geno <- GdsGenotypeReader(filename = bfile_path %&% ".gds")
genoData <- GenotypeData(geno)
iids <- getScanID(genoData)

#if .kin does or doesn't have entries
kin <- fread("relatedness.kin", header = T)
if(nrow(kin) > 0){
  KINGmat <- king2mat(file.kin0 = "relatedness.kin0", file.kin = "relatedness.kin", iids = iids)
}else{
  KINGmat <- king2mat(file.kin0 = "relatedness.kin0", iids = iids)
}

#convert to matrix format
KINGmat2 <- as.data.frame(KINGmat)
colnames(KINGmat2) <- iids
rownames(KINGmat2) <- iids
fwrite(KINGmat2, "relatedness_wIID.txt", sep = "\t", row.names = T, col.names = T, nThread = 40)
fwrite(KINGmat2, "relatedness_woIID.txt", sep = "\t", row.names = F, col.names = F, nThread = 40)

#for GCTA b/c they don't like neg. values
KINGmat3 <- KINGmat2
KINGmat3[KINGmat3 < 0] <- 0
fwrite(KINGmat3, "relatedness_wIID_noNeg.txt", sep = "\t", row.names = T, col.names = T, nThread = 40)
fwrite(KINGmat3, "relatedness_woIID_noNeg.txt", sep = "\t", row.names = F, col.names = F, nThread = 40)

#then GCTA style
KINGmat3 <- as.matrix(KINGmat3)
rownames(KINGmat3) <- iids
colnames(KINGmat3) <- iids
write_GRMBin(KINGmat3, prefix = "relatedness")
print("Completed creating relatedness matrix. Now calculating PCs.")

#calculate PCs
KING_pca <- king_path %&% " --pca -b " %&% bfile_path %&% ".bed"
system(KING_pca)
print("Completed calculating PCs, saved into the files kingpc.dat and kingpc.ped.")
