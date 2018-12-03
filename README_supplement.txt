This is a supplement to Ad_PX_pipe that better explains what each step and each script does.
  *I'm going to assume whoever's reading this at least knows what a GWAS is and has a vague understanding of PrediXcan

*The test data we're going to use is 1000G AMR, a cohort of 347 individuals from the Americas.
*Individuals include Mexican Ancestry from Los Angeles, USA (MXL); Puerto Ricans from Puerto Rico (PUR); Colombians from Medellin, Colombia (CLM); and Peruvians from Lima, Peru (PEL)
*These correspond nicely with my study of the Hispanic Community Health Study/Study of Latinos, as these individuals also have multi-continental ancestry
  *Also, HCHS is too large to play well as test data
*I subset these data to only include 100k SNPs so it's faster for the purpose of serving as test data
*See: https://www.researchgate.net/profile/Taras_Oleksyk/publication/239524406/figure/fig1/AS:298625939853314@1448209380432/a-Individual-ancestry-proportions-in-the-1000-Genomes-CLM-MXL-and-PUR-populations.png

0. Produce phenotypes and covariates (ex. medicines) in R (for test data only)
INPUT: PLINK binary format genotypes (--bfile)
OUTPUT: two covariate files (one w/ IDs, one w/o) and two phenotype files (one w/ IDs, one w/o)
  *NOTE: in a real data set, fill all missing values with "NA" or you're going to have a bad time
  *1000 Genomes data doesn't have real phenotype data available
  *We will simulate a random phenotype and a random covariate
	*Covariates are used to prevent inflation of results by non-genotype confounders, such as cholesterol-lowering medications in a study on lipid traits
    *Phenotypes are randomly selected from a normal distribution (rnorm) and covariates are randomly selected from a binomial distribution (binom)
  *We have two sets of each output file because GEMMA only takes files without headers and IDs. Though we only have one covariate and one phenotype, in real data, it may be hard to keep track of all the phenotypes and covariates floating around, so we save a w/ ID column file for human readability

1. Perform quality control in PLINK using gwasqc_pipeline
  *Documentation for this step is available at: https://github.com/WheelerLab/gwasqc_pipeline
  
2. Calculate principal components and a relationship matrix in KING
INPUT: PLINK binary format genotypes (--bfile)
OUTPUT: four relatedness matrices (w/ ID w/ negs., w/o ID w/ negs., w/ ID w/o negs, and w/o ID w/o negs.), GCTA-style relatedness file, and a file of top 20 principal componenents
  *The software KING is optimized to calculate principal components in admixed populations, so we are using that instead of PLINK or SMARTPCA
  *We first calculate relatedness between all the individuals in our population and convert it to a matrix format for GEMMA, as well as GCTA-COJO format for later downstream
    *GCTA requires a separate format, and also does not accept negative values of relatedness that KING outputs
  *We then calculate PCs, which is output as kingpc.ped

3. Impute data with the Michigan imputation server
  *See PowerPoint within repository for better documentation

4. Convert imputed data to PrediXcan dosage 
OUTPUT: 22 dosage files (genotypes) and samples.txt (IDs) file 
  *This python script takes one UMich imputed vcf files as input, removes SNPs with R2<0.8 and MAF>0.01 (options to change), finds the rsID for each SNP, and makes output files:
  *We use a slightly different approach for our data because I don't want to make someone go through imputation just for test data
  
5. Perform a genome-wide association study in GEMMA

5a. Convert from PrediXcan dosage format to BIMBAM format (GEMMA genotype imput)
INPUT: path to dosage folder (--dosage-path), suffix of dosages (--dosage_suffix)
OUTPUT: BIMBAM-format genotypes in a BIMBAM/ folder, and annotation files in an anno/ folder
  *GEMMA requires yet another genotype format called BIMBAM, which is similar to PrediXcan dosage format
  *This script reorganizes the PrediXcan dosages into two that are part of GEMMA: BIMBAM (genotype format) and anno (SNP location information)

5b. Make covariance files from known covariates and KING PCs
INPUT: covariates file (without IDs), principal components file, number of principal components to keep, and output name
OUTPUT: combined covariates file in GEMMA format
  *This simple script combines the principal components and covariances files make earlier into a format appropriate for GEMMA

5c. Run all chrs in a loop
INPUT: prefix of genotype files (-g), phenotype file (w/o IDs) (-p), prefix of annotation files (-a), relatedness matrix (w/o IDs) (-k), covariates file (-c), and output prefix (-o)
OUTPUT: association files in output/ folder
  *This simple script is just a wrapper around GEMMA for 22 autosomal chromosomes

6. Calculate predicted gene expressions in PrediXcan using GTEx and MESA models
INPUT: path to dosages (--dosages_path) and prefix of output, which should be a folder name (--output_prefix)
OUTPUT: 49 predicted expression files 
  *This script is yet another wrapper around PrediXcan to use all current available models (44 GTEx, 5 MESA)
  *There is also an option to change file paths of where the models are coming from. Currently, I set GTEx to V6.  
  
7. Perform an imputed transcriptome-based association study in GEMMA

7a. Convert predicted expression to GEMMA-style pseudo-genotypes
INPUT: prefix of PrediXcan predicted expression files (--pred_exp_prefix) and prefix of output, which should be a folder (--output_prefix)
OUTPUT: GEMMA-style pseudo-genotypes, similar to the BIMBAM files earlier
  *This script takes in the just-produced predicted gene expressions and converts them to GEMMA-style pseudo-genotypes in similar format to SNPs 

7b. Run all pops. and tissues in a loop
INPUT: see 5c, but also run the flag "-h", standing for "hacked" GEMMA
OUTPUT: see 5c
  *This is the same as 5c, but using genes instead of SNPs
  
8. Find significant SNPs from GWAS and significant genes from PrediXcan
INPUT: prefix of input (except for "output/") and optional parameters for significance of SNPs and GEMMA
OUTPUT: significant SNP file and significant gene file
  *Goes through all GEMMA results and outputs the SNPs or genes that are below the input thresholds. Defaults are 5e-8 (SNP) and 9.654e-6 (gene)
  
  
  
  
  
  
  
  
  
  
  
  
  