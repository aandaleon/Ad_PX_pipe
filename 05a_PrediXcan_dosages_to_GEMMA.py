#Takes information from PrediXcan-style dosages to make SNP annotation and BIMBAM files for GEMMA
#by Angela Andaleon (aandaleon@luc.edu)
import argparse
import gzip
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("--desc", type = str, action = "store", dest = "desc", required = False, help = "Takes information from PrediXcan-style dosages to make SNP annotation and BIMBAM files for GEMMA.")
parser.add_argument("--dosage_path", type = str, action = "store", dest = "dosage_path", required = True, help = "Path to folder containing dosages and samples.txt")
parser.add_argument("--chr", type = str, action = "store", dest = "chr", required = False, help = "Path to chromosome to analyze. If no input, analyzes all 22 pairs of chromosomes.")
parser.add_argument("--dosage_suffix", type = str, action = "store", dest = "dosage_suffix", default = ".maf0.01.r20.8.dosage.txt.gz", help = "Suffix of dosages. Default = .maf0.01.r20.8.dosage.txt.gz")
args = parser.parse_args()

print("Reading input files.")
dosage_path = args.dosage_path
dosage_suffix = args.dosage_suffix
BIMBAM_path = "BIMBAM/" #this was going to be a user option but I didn't see a point
anno_path = "anno/"

#make BIMBAM and anno folders if they don't already exist
os.system("mkdir -p BIMBAM")
os.system("mkdir -p anno")

if args.chr is None: #if no argument given to run a certain chr, run all chr (I don't think this flag works... might remove)
    chrs_to_test = range(1, 23)
if args.chr is not None:
    chrs_to_test = range(int(args.chr), (int(args.chr) + 1))

'''
dosage_path = "dosages/"
'''

dosage_samples = pd.read_csv(dosage_path + "samples.txt", delim_whitespace = True, header = None) #read in samples.txt file from PrediXcan-style dosages
dosage_samples = dosage_samples[[1]] #take only IIDs

print("Starting conversion from PLINK dosage to GEMMA input BIMBAM and anno.")
for i in chrs_to_test: #for each chr
    anno = open(anno_path + "anno" + str(i) + ".txt", "w") #start an annotation file
    BIMBAM = gzip.open(BIMBAM_path + "chr" + str(i) + ".txt.gz", "wb") #start a genotype file
    for line in gzip.open(dosage_path + "chr" + str(i) + dosage_suffix, "rb"): #open PrediXcan-style dosage file
        arr = line.strip().split() #split each line of PrediXcan-style dosage file into an array
        (chr, rs, bp, A1, A2, MAF) = arr[0:6] #extract non-genotype information into tuple
        if len(A1) < 2 and len(A2) < 2: #keep only alleles, no indels
           dosages = arr[6:] #pull out dosages into array
           dosages = [str(dosage) for dosage in dosages] #convert dosages all into str
           dosages_str = '\t'.join(dosages) #separate dosages by delimiter and convert to str
           BIMBAM_format = (rs + "\t" + A1 + "\t" + A2 + "\t" + dosages_str + "\n") #write genotype to BIMBAM file in proper format
           BIMBAM.write(BIMBAM_format)
           anno.write(rs + "\t" + bp + "\t" + str(i) + "\n") #write annotations to anno file in proper format
    BIMBAM.close()
    anno.close()
    print("Completed with chr " + str(i) + ".")
print("Conversion from PLINK dosage to GEMMA input has been completed. Have a nice day!")
