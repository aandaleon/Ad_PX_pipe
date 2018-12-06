# Ad_PX_pipe
This repository reorganizes and restructures scripts from ["Analysis of the genetic architecture and predicted gene expression of lipid traits in Hispanic cohorts"](https://github.com/WheelerLab/px_his_chol) to be more user-friendly. All paths to softwares are defaulted to those on wheelerlab3, with genotypes available at `/home/angela/Ad_PX_pipe/AMR`, and it is expected that all scripts are run from the same directory. For exact details on the inner workings of each script, use the `--help` flag or see the manual. We will be using genotypes from 1000 Genomes American superpopulation and simulating phenotypes and covariances in R. The input population does not need to be admixed, but local ancestry steps (12-15) are for use mainly in admixed (African-American, Hispanic) populations. You should be able to copy the genotypes and copy and paste all the commands from this README and everything should run... but will you understand it?

For much more detail on the process of everything in here, please see the manual, I spent a lot of time writing it :).

00. Produce phenotypes and covariates (ex. medicines) in R (for test data only)
    * `Rscript 00_simulate_pheno_covar.R --bfile AMR`

#For use in all populations:

01. Perform quality control in [PLINK](https://www.cog-genomics.org/plink/1.9/filter) using [gwasqc_pipeline](https://github.com/WheelerLab/gwasqc_pipeline)
    * Test data has already been filtered

02. Calculate principal components and a relationship matrix in [KING](http://people.virginia.edu/~wc9c/KING/manual.html)
    * `plink --bfile AMR --chr 22 --make-bed --out AMR_chr22; Rscript 02_relate_matrix_PCs.R --bfile AMR_chr22`
 
03. Impute data with the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)
    * Test data is already imputed
    * Follow instructions in 03_Michigan_Imputation_Server.pptx

04. Convert imputed data to PrediXcan dosages using [this script](https://github.com/WheelerLab/Imputation/blob/master/UMich_vcf2pxfixCAAPA.py).

    * Since test data is already imputed in PLINK format, used [this script](https://github.com/hakyimlab/PrediXcan/blob/master/Software/convert_plink_to_dosage.py) instead
      * `mkdir dosages/; cd dosages/; wget https://raw.githubusercontent.com/hakyimlab/PrediXcan/master/Software/convert_plink_to_dosage.py; awk '{ print $1"\t"$2 }' ../AMR.fam > samples.txt; python convert_plink_to_dosage.py -b ../AMR -p /usr/local/bin/plink; cd ..`

05. Perform a genome wide association study in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)
    * a. Convert from PrediXcan dosage format to BIMBAM format (GEMMA genotype input)
      * `python 05a_PrediXcan_dosages_to_GEMMA.py --dosage_path dosages/ --dosage_suffix .txt.gz`
    * b. Make covariance file from known covariates and KING PCs
      * `Rscript 05b_make_GEMMA_covars.R --covar covar_woIID.txt --pcs_file kingpc.ped --pcs_num 5 --output GEMMA_covars.txt`
    * c. Run all chrs. in a loop
      * `bash 05c_GEMMA_loop.sh -g BIMBAM/chr -p pheno_woIID.txt -a anno/anno -k relatedness_woIID.txt -c GEMMA_covars.txt -o AMR_`
      * Note to self: make easier to iterate through more than one pheno

06. Calculate predicted gene expressions in [PrediXcan](https://github.com/hakyimlab/PrediXcan) using [GTEx](http://predictdb.org/) and [MESA](https://github.com/aandaleon/DivPop) models
    * `python 06_make_pred_exp.py --dosages_path dosages/ --output_prefix pred_exp/`

07. Perform an imputed transcriptome-based association study in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)
    * a. Convert predicted expression to GEMMA-style pseudo-genotypes
      * `python 07a_convert_PrediXcan_to_GEMMA.py --pred_exp_prefix pred_exp/ --output_prefix pred_exp_GEMMA/`
    * b. Run all pops. and tissues in a loop
      * `bash 05c_GEMMA_loop.sh -g pred_exp_GEMMA/ -p pheno_woIID.txt -a anno/anno -k relatedness_woIID.txt -c GEMMA_covars.txt -o AMR_ -h`    

08. Find significant SNPs from GWAS and significant genes from PrediXcan in R
    * `python 08_sig_SNP_sig_gene.py --SNP_sig 5e-4 --gene_sig 0.05 --input_prefix AMR_`
      * NOTE: I set these thresholds arbitrarilty low so we have significant genes to with with later. Usually, use the defaults, 5e-8 and 9.654e-6, respectively.

09. Calculate independent significant SNPs in a joint analysis in [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO)
    * a. Make GWAS output into GCTA-COJO format 
      * `python 09a_GEMMA_to_GCTA-COJO.py --fam AMR.fam --GWAS_prefix AMR_ --output_prefix AMR`
    * b. Run GCTA (this is not a script, this is actual gcta)
      * `gcta64 --cojo-file AMR.ma --cojo-slct --cojo-p 5e-4 --bfile AMR --cojo-actual-geno --out AMR`
        * Again, set P arbitrarily low for example; set to 5e-8 for real analyses

10. Perform colocalization between GWAS results and eQTL data using [COLOC](https://cran.r-project.org/web/packages/coloc/coloc.pdf) in a [COLOC wrapper](https://github.com/hakyimlab/summary-gwas-imputation)
    * a. Put GWAS and eQTL data into COLOC input format (takes a while depending on number of SNPs in GWAS)
      * `Rscript 10a_make_COLOC_input.R --ma AMR.ma --GWAS_prefix AMR_ --sample_size 347`
    * b. Run COLOC wrapper (first argument sample size, second argument GWAS prefix)
      * `bash 10b_run_COLOC.sh 347 AMR_`

11. Perform backward elimination modeling of all significant genes to find independent signals in R
    * `Rscript 11_back_elim.R --sig_gene output/AMR_sig_genes.txt --pheno pheno_wIID.txt --pred_exp_prefix AMR_ --pheno_name pheno`

#For use in admixed populations:

12. Make reference populations for use in local ancestry inference in [PLINK](https://www.cog-genomics.org/plink/1.9/data)
    * Note: may just make these ahead of time and have users choose

13. Infer haplotypes with [HAPI-UR](https://code.google.com/archive/p/hapi-ur/)

14. Calculate local ancestry inference in [RFMix](https://sites.google.com/site/rfmixlocalancestryinference/)

15. Perform admixture mapping in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)


