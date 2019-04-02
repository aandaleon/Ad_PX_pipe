#this file can be run to take user input and pass it to the querying file
import argparse
import glob
import numpy as np
import os
import pandas as pd

#if the python script is run without any flags, just output genenames, cv_R2_avg, rsid, and weights (most common things we use)
parser = argparse.ArgumentParser()

#string inputs
#if none, query all available
parser.add_argument("--db", type = str, action = "store", dest = "db", required = False, help = ".db you want to query.") #"db/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db"
parser.add_argument("--genes", type = str, action = "store", dest = "genes", required = False, help = "File containing genes (Ensembl IDs).") #"practice_python_queries/genenames.txt"
parser.add_argument("--genenames", type = str, action = "store", dest = "genenames", required = False, help = "File containing gene names.") #"practice_python_queries/genenames.txt"

#boolean inputs
#EXTRA
parser.add_argument("--n.snps.in.model", action = "store_true", dest = "n_snps_in_model", default = False, help = "Output the number of SNPs within the cis window that have non-zero weights, as found by elastic-net.")
parser.add_argument("--test_R2_avg", action = "store_true", dest = "test_R2_avg", default = False, help = "Output the average coefficient of determination when predicting values of the hold out fold during nested cross validation.")
parser.add_argument("--cv_R2_avg", action = "store_true", dest = "cv_R2_avg", default = False, help = "Output the average coefficient of determination for each of the hold out folds when cross-validation was performed on the entire data set.")
parser.add_argument("--rho_avg", action = "store_true", dest = "rho_avg", default = False, help = "Output the average correlation between predicted and observed on the hold out folds when doing nested cross-validation.")
parser.add_argument("--rho_zscore", action = "store_true", dest = "rho_zscore", default = False, help = "Output the transformation of rho_avg into a z-score using Stouffer's Method.")
parser.add_argument("--pred.perf.R2", action = "store_true", dest = "pred_perf_R2", default = False, help = "Output the rho_avg squared.")
parser.add_argument("--pred.perf.pval", action = "store_true", dest = "pred_perf_pval", default = False, help = "Output the p-value for rho_zscore.")

#WEIGHTS (this table depends on the EXTRA table)
parser.add_argument("--rsid", action = "store_true", dest = "rsid", default = False, help = "Output the rsids in the models of queried genes.")
parser.add_argument("--varID", action = "store_true", dest = "varID", default = False, help = "Output the variant IDs in the models of queried genes. These are string labels of the format chromosome_position_allele1_allele2_build. All varIDs are from build 37 of the HUman Reference Genome.")
parser.add_argument("--ref_allele", action = "store_true", dest = "ref_allele", default = False, help = "Output the reference alleles of the SNPs in the models of the queried genes.")
parser.add_argument("--eff_allele", action = "store_true", dest = "eff_allele", default = False, help = "Output the effect alleles of the SNPs in the models of the queried genes.")
parser.add_argument("--weight", action = "store_true", dest = "weight", default = False, help = "Output the weights for the SNPs that are used to calculate predicted expression for the gene. In predicting the expression for the gene, the weight is multiplied by the count, or estimated count, of the effect allele in individual. This value is added to all other weighted SNPs in the model.")

#SAMPLE INFO
parser.add_argument("--n_samples", action = "store_true", dest = "n_samples", default = False, help = "Output the number of samples used the make the .db file.")
parser.add_argument("--population", action = "store_true", dest = "population", default = False, help = "Output the population studied.")
parser.add_argument("--tissue", action = "store_true", dest = "tissue", default = False, help = "Output the tissue or MESA population from which RNA was sequenced.")

#THRESHOLDS (may mess around with)
parser.add_argument("--test_R2_avg_thres", type = float, dest = "test_R2_avg_thres", default = 0, help = "Restrict the test_R2_avg to values above this threshold.")
parser.add_argument("--cv_R2_avg_thres", type = float, dest = "cv_R2_avg_thres", default = 0, help = "Restrict the cv_R2_avg to values above this threshold.")
parser.add_argument("--rho_avg_thres", type = float, dest = "rho_avg_thres", default = 0, help = "Restrict the rho_avg to values above this threshold.")
parser.add_argument("--pred.perf.R2_thres", type = float, dest = "pred_perf_R2_thres", default = 0, help = "Restrict the test_R2_avg to values above this threshold.")
parser.add_argument("--pred.perf.pval_thres", type = float, dest = "pred_perf_pval_thres", default = 0, help = "Restrict the pred_perf_pval to values below this threshold.")

args = parser.parse_args() #then pass these arguments to further things


###EXTRA
extra_flags = [] #store the flags the user passes
if args.n_snps_in_model:
    extra_flags.append("n.snps.in.model")
if args.test_R2_avg:
    extra_flags.append("test_R2_avg")
if args.cv_R2_avg:
    extra_flags.append("cv_R2_avg")
if args.rho_avg:
    extra_flags.append("rho_avg")
if args.rho_zscore:
    extra_flags.append("rho_zscore")
if args.pred_perf_R2:
    extra_flags.append("pred.perf.R2")
if args.pred_perf_pval:
    extra_flags.append("pred.perf.pval")
print(extra_flags)
  
###WEIGHTS
weights_query = False #default to not query weights
weights_flags = []
if args.rsid or args.varID or args.ref_allele or args.eff_allele or args.weight:
    weights_query = True #you will be querying the weights matrix
if args.rsid:
    weights_flags.append("rsid")
if args.varID:
    weights_flags.append("varID")
if args.ref_allele:
    weights_flags.append("ref_allele")
if args.eff_allele:
    weights_flags.append("eff_allele")
if args.weight:
    weights_flags.append("weight")
print(weights_flags)

###SAMPLE INFO
sample_info_flags = []
if args.n_samples:
    sample_info_flags.append("n_samples")
if args.population:
    sample_info_flags.append("population")
if args.tissue:
    sample_info_flags.append("tissue")
print(sample_info_flags)

if len(extra_flags) + len(weights_flags) + len(sample_info_flags) == 0:
    print("The user has passed no flags. Do a thing.")

args_db = args.db #.endswith doesn't like arguments
if args_db.endswith(".db"): #its a single .db file
    dbs = list(args.db)
else: #its (I assume) a folder
    if args_db.endswith("/"):
        folder_name = args_db
    else:
        folder_name = args_db + "/"
    dbs = []
    for file in glob.glob(".db"): #find files in a folder - https://stackoverflow.com/questions/3964681/find-all-files-in-a-directory-with-extension-txt-in-python
        dbs.append(folder_name + file)
    print(dbs)
    

'''
So the flags the user wants are stored in:
    extra
    weights
    sample info
'''





########
#CARLEE#
########

###QUERYING 
#NOTE: we'll just query EVERYTHING (aka use *) in SQL b/c I'm more comfortable parsing in pandas than in SQL
#for db/, have the user direct to a FOLDER of db files
  #OR have two db flags - one that iterates through a folder of .dbs and one that is for a specific .db file
    #if ends w/ .db, it's a db (make a list of one)
    #else find all the .dbs in that folder
    #Folder or single file, that's what this if/else hierarhcy does 
#have giant list to make into list of lists
#iterate through each db file
  #db file will be a column name in final output
  #everything is encased in one big loop around a list of db files to parse

  #if user did pass genes or genenames
    #translates genenames to genes
    #store list in genes variable
  #else
    #select gene from extra
    #make all genes into a list

  #iterate through each gene (see python practice example query stuff)
    #if the user chose any general flags (anything from EXTRA, WEIGHTS, or SAMPLE INFO)

      #if chose anything from extra
        #select n.snps.in.model,test_R2_avg,etc. from extra where gene is iterated gene
          #it'll output a tuple
          #store all these in variables to be used later on (at least the gene)
      
      #repeat these lines for the sample info flags

      #if user chooses any weights flags
        #select rsid,var_ID,etc. from weights where gene is iterated gene
          #it'll output a tuple
          #store all these in variables to be used later on
          #then append all information (all info in weights, all info in extra, all info in sample info flags) to list of lists

  #else (user chooses no flags)
    #get list of all genes in model
    #iterate through each gene
      #pull out cv_R2_avg
      #search in weights for all of them
        #pull rsid and weights

########        
#SHREYA#
########

###PARSING QUERY OUTPUT
#alright so we're done querying shiz and we got a big boi list of lists
#convert list of lists into dataframe

#if user has flags 
  #get everything "true" the user wants
    #always include genename
    #ex. if they want cv_r2_avg, rho_avg, and rsid, make a list of ["genename", "cv_r2_avg", "rho_avg", "rsid"]
    
#if user has no flags
  #genename, cv_r2_avg, rs, and weights
#give column names
  #(this will all be in the same order so we just need to figure out what this order is)
#only pull columns of what the user wants

#restrict to thresholds the user wants (see threshold flags)
  #ex. if only want cv_r2_avg > 0.1, only keep those
   #(this will involve using a bunch of .loc)
  
#delete duplicate rows
#print to csv
  #remove indexes cause they're annoying
  
#additional ideas:
  #pull what the genes have been implicated in the GWAS catalog
  
  
