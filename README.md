# Ad_PX_pipe
This repository reorganizes and restructures scripts from ["Analysis of the genetic architecture and predicted gene expression of lipid traits in Hispanic cohorts"](https://github.com/WheelerLab/px_his_chol) to be more user-friendly. We will be using genotypes from 1000 Genomes American superpopulation and simulating phenotypes and covariances in R. For exact details on the inner workings of each script, use the `--help` flag. All paths to softwares are defaulted to those on wheelerlab3, with genotypes available at `/home/angela/Ad_PX_pipe/AMR`.

00. Produce phenotypes and covariates (ex. medicines) in R (for test data only)

    * `Rscript 00_simulate_pheno_covar.R --bfile AMR`
    * Note: may have to cheat and change this later to force some to be significant

01. Perform quality control in [PLINK](https://www.cog-genomics.org/plink/1.9/filter) using [gwasqc_pipeline](https://github.com/WheelerLab/gwasqc_pipeline)
    
    * Test data has already been filtered

02. Calculate principal components and a relationship matrix in [KING](http://people.virginia.edu/~wc9c/KING/manual.html)

    * `Rscript 02_relate_matrix_PCs.R --bfile AMR_chr22`
 
03. Impute data with the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)

    * Test data is already imputed
    * Follow instructions in 03_Michigan_Imputation_Server.pptx

04. Convert imputed data to PrediXcan dosages using [this script](https://github.com/WheelerLab/Imputation/blob/master/UMich_vcf2pxfixCAAPA.py).

    * Since test data is already imputed in PLINK format, used [this script](https://github.com/hakyimlab/PrediXcan/blob/master/Software/convert_plink_to_dosage.py) instead

05. Perform a genome wide association study in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)
    * Convert from PrediXcan dosage format to BIMBAM format (GEMMA genotype input)
      * 05a. `python 05a_PrediXcan_dosages_to_GEMMA.py --dosage_path dosages/ --dosage_suffix .txt.gz`
    * Make covariance file from known covariates and KING PCs
      * 05b. `Rscript 05b_make_GEMMA_covars.R --covar covar_woIID.txt`

06. Calculate predicted gene expressions in [PrediXcan](https://github.com/hakyimlab/PrediXcan) using [GTEx](http://predictdb.org/) and [MESA](https://github.com/aandaleon/DivPop) models

07. Perform an imputed transcriptome-based association study in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)

08. Find significant SNPs from GWAS and significant genes from PrediXcan in R

09. Calculate independent significant SNPs in a joint analysis in [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO)

10. Perform colocalization between GWAS results and eQTL data using [COLOC](https://cran.r-project.org/web/packages/coloc/coloc.pdf) in a [COLOC wrapper](https://github.com/hakyimlab/summary-gwas-imputation)

11. Perform backward elimination modeling of all significant genes in R

12. Make reference populations for use in local ancestry inference in [PLINK](https://www.cog-genomics.org/plink/1.9/data)

13. Infer haplotypes with [HAPI-UR](https://code.google.com/archive/p/hapi-ur/)

14. Calculate local ancestry inference in [RFMix](https://sites.google.com/site/rfmixlocalancestryinference/)

15. Perform admixture mapping in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)


