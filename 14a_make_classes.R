#makes .classes file for RFMix from HAPI-UR
#by Angela Andaleon (aandaleon@luc.edu)
"%&%" = function(a,b) paste(a,b,sep="")
library(dplyr)
library(data.table)
args <- commandArgs(trailingOnly = T)
phind_file_name <- args[1] #name of .phind file that HAPI-UR output
pop_file_name <- args[2] #name of pop file with three columns: FID, IID, and pop
#test_pop <- args[3] #code of population to be tested

#phind_file_name <- "haplotypes/chr22.phind"
#pop_file_name <- "HIS"

#parse .phind file
phind <- fread(phind_file_name, header = F)
phind$V2 <- NULL
phind$V3 <- NULL
phind$FID <- gsub("\\:.*", "", phind$V1) #separate FID and IID columns
phind$IID <- gsub(".*:\\s*|_.*", "", phind$V1)
phind$V1 <- NULL

#add pops
pop_file <- fread("/home/angela/Ad_PX_pipe_data/RFMix/RefPop/" %&% pop_file_name %&% "_pop.txt", header = F)
colnames(pop_file) <- c("FID", "IID", "pop")
ordered_phind <- left_join(phind, pop_file, by = c("FID", "IID"))
ordered_phind$pop <- as.numeric(factor(ordered_phind$pop)) #convert pops to ints
ordered_phind[is.na(ordered_phind)] <- 0
#ordered_phind$pop[ordered_phind$pop == max(ordered_phind$pop)] <- "0" #assign 0 to test pop (which should be the max if ordered properly)
write(paste(as.character(ordered_phind$pop), collapse = " "), "RFMix.classes")

