# Ad_PX_pipe
This repository reorganizes and restructures scripts from ["Analysis of the genetic architecture and predicted gene expression of lipid traits in Hispanic cohorts"](https://github.com/WheelerLab/px_his_chol) to be more user-friendly. We will be using genotypes from 1000 Genomes American superpopulation and simulating phenotypes and covariances in R. For exact details on the inner workings of each script, use the `--help` flag. All paths to softwares are defaulted to those on wheelerlab3, with genotypes available at `/home/angela/Ad_PX_pipe/AMR`.

00. Produce phenotypes and covariates (ex. medicines) in R (for test data only)

    * `Rscript 00_simulate_pheno_covar.R --bfile AMR`
    * Note: may have to cheat and change this later to force some to be significant

01. Perform quality control in [PLINK](https://www.cog-genomics.org/plink/1.9/filter) using [gwasqc_pipeline](https://github.com/WheelerLab/gwasqc_pipeline)

02. Calculate principal components and a relationship matrix in [KING](http://people.virginia.edu/~wc9c/KING/manual.html)

    * `Rscript 02_relate_matrix_PCs.R --bfile AMR_chr22`
 
03. Impute data with the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)

04. Convert imputed data to PrediXcan dosages using [this script](https://github.com/WheelerLab/Imputation/blob/master/UMich_vcf2pxfixCAAPA.py).

04. Perform a genome wide association study in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)

05. Calculate predicted gene expressions in [PrediXcan](https://github.com/hakyimlab/PrediXcan) using [GTEx](http://predictdb.org/) and [MESA](https://github.com/aandaleon/DivPop) models

06. Perform an imputed transcriptome-based association study in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)

07. Find significant SNPs from GWAS and significant genes from PrediXcan in R

08. Calculate independent significant SNPs in a joint analysis in [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO)

09. Perform colocalization between GWAS results and eQTL data using [COLOC](https://cran.r-project.org/web/packages/coloc/coloc.pdf) in a [COLOC wrapper](https://github.com/hakyimlab/summary-gwas-imputation)

10. Perform backward elimination modeling of all significant genes in R

11. Make reference populations for use in local ancestry inference in [PLINK](https://www.cog-genomics.org/plink/1.9/data)

12. Infer haplotypes with [HAPI-UR](https://code.google.com/archive/p/hapi-ur/)

13. Calculate local ancestry inference in [RFMix](https://sites.google.com/site/rfmixlocalancestryinference/)

14. Perform admixture mapping in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)

