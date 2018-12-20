#converts from BIMBAM format (comma delimited) to PrediXcan dosage format
from __future__ import division #because python 2 doesn't like real numbers
import argparse #allows for input of command line arguments
parser = argparse.ArgumentParser() #open the parser
parser.add_argument("--anno", type = str, action = "store", dest = "anno", required = True, help = "Input annotation file")
parser.add_argument("--BIMBAM", type = str, action = "store", dest = "BIMBAM", required = True, help = "Input BIMBAM file")
parser.add_argument("--output", type = str, action = "store", dest = "output", required = False, default = "dosages.txt", help = "Output dosage file")
args = parser.parse_args() #parse given arguments

anno = open(args.anno, "r") #open anno file
BIMBAM = open(args.BIMBAM, "r") #open BIMBAM file
output = open(args.output, "w") #open output file

'''
BIMBAM = open("YRI.TGP_and_imputed.rmBAD.20130526.geno.TEST_10000", "r")
anno = open("YRI.TGP_and_imputed.rmBAD.20130526.snp.info.TEST_10000", "r")
output = open("dosages.txt", "w")
'''

for anno_line, BIMBAM_line in zip(anno, BIMBAM): #go through each line in anno and BIMBAM (should be the same)
    a_line = anno_line.strip().split(",") #split anno line into list
    (rs_anno, bp, chr) = a_line[0:3] #assign each part of anno line to a variable
    B_line = BIMBAM_line.strip().split(",") #split BIMBAM line into list
    (rs_BIMBAM, A1, A2) = B_line[0:3] #assign first part of BIMBAM into variable
    dosages = B_line[3:] #assign last part of BIMBAM to a list
    float_dosages = [float(dosage) for dosage in dosages] #convert each dosage into floats
    MAF = sum(float_dosages)/(len(dosages) * 2) #calculate minor allele frequency (I hope the math is correct)
    if MAF > 0.5: #if minor allele frequency is the other allele
        MAF = 1 - MAF #take the complement (is that the right use of the word complement?)
    MAF = str(MAF) #convert to string b/c Python doesn't like when you concatenate strings and numbers later
    dosages_str = "\t".join(dosages) #create string of dosages separate by tabs
    if rs_anno == rs_BIMBAM: #if the rs in anno and the rs in BIMBAM match, which they always should
        dosage_format = (chr + "\t" + rs_anno + "\t" + bp + "\t" + A1 + "\t" + A2 + "\t" + MAF + "\t" + dosages_str + "\n") #make all variables into PrediXcan dosage format
        output.write(dosage_format) #write to output file

anno.close() #close anno file
BIMBAM.close() #close BIMBAM file
output.close() #close output file
