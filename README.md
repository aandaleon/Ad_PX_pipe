# Ad_PX_pipe
This repository reorganizes and restructures scripts from ["Genetically regulated gene expression underlies lipid traits in Hispanic cohorts"](https://github.com/WheelerLab/px_his_chol) to be more user-friendly. All paths to softwares are defaulted to those on wheelerlab3, and it is expected that all scripts are run from the same directory. For exact details on the inner workings of each script, use the `--help` flag or see the manual. We will be using genotypes from 1000 Genomes American superpopulation (`AMR.bed/.bim/.fam`) and simulating phenotypes and covariances in R. The input population does not need to be admixed, but local ancestry steps (12-15) are for use mainly in admixed (African-American, Hispanic) populations. You should be able to copy and paste all the commands from this README and everything should run... but will you understand it?

For much more detail on the process of everything in here, please see the manual, I spent a lot of time writing it :).

00. Produce phenotypes and covariates (ex. medicines) in R (for test data only)
    * `Rscript 00_simulate_pheno_covar.R --bfile AMR`

### For use in all populations:

01. Perform quality control in [PLINK](https://www.cog-genomics.org/plink/1.9/filter) using [gwasqc_pipeline](https://github.com/WheelerLab/gwasqc_pipeline)
    * Test data has already been filtered

02. Calculate principal components and a relationship matrix in [KING](http://people.virginia.edu/~wc9c/KING/manual.html)
    * `plink --bfile AMR --chr 22 --make-bed --out AMR_chr22; Rscript 02_relate_matrix_PCs.R --bfile AMR_chr22`
 
03. Impute data with the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)
    * Test data is already imputed
    * Follow instructions in 03_Michigan_Imputation_Server.pptx

04. Convert imputed data to PrediXcan dosages using [this script](https://github.com/WheelerLab/Imputation/blob/master/UMich_vcf2pxfixCAAPA.py).

    * Since test data is already imputed in PLINK format, run this (for example data only)
      * `mkdir dosages/; cd dosages/; wget https://raw.githubusercontent.com/hakyimlab/PrediXcan/master/Software/convert_plink_to_dosage.py; awk '{ print $1"\t"$2 }' ../AMR.fam > samples.txt; python convert_plink_to_dosage.py -b ../AMR -p /usr/local/bin/plink; cd ..; echo "Completed converting dosages."`

05. Perform a genome wide association study in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)
    * a. Convert from PrediXcan dosage format to BIMBAM format (GEMMA genotype input)
      * `python 05a_PrediXcan_dosages_to_GEMMA.py --dosage_path dosages/ --dosage_suffix .txt.gz`
    * b. Make covariance file from known covariates and KING PCs
      * `Rscript 05b_make_GEMMA_covars.R --covar covar_woIID.txt --pcs_file kingpc.ped --pcs_num 5 --output GEMMA_covars.txt`
    * c. Run all chrs. in a loop
      * `python 05c_GEMMA_wrapper.py --relatedness relatedness_woIID.txt --geno_prefix BIMBAM/chr --pheno pheno_woIID.txt --pheno_names pheno_names.txt --covariates GEMMA_covars.txt --anno anno/anno --output AMR --GWAS`

06. Calculate predicted gene expressions in [PrediXcan](https://github.com/hakyimlab/PrediXcan) using [GTEx](http://predictdb.org/) and [MESA](https://github.com/aandaleon/DivPop) models
    * `python 06_make_pred_exp.py --dosages_path dosages/ --output_prefix pred_exp/`

07. Perform an imputed transcriptome-based association study in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)
    * a. Convert predicted expression to GEMMA-style pseudo-genotypes
      * `python 07a_convert_PrediXcan_to_GEMMA.py --pred_exp_prefix pred_exp/ --output_prefix pred_exp_GEMMA/`
    * b. Run all pops. and tissues in a loop
      * `python 05c_GEMMA_wrapper.py --relatedness relatedness_woIID.txt --geno_prefix pred_exp_GEMMA/ --pheno pheno_woIID.txt --pheno_names pheno_names.txt --covariates GEMMA_covars.txt --anno anno/anno --output AMR --pred_exp`    

08. Find significant SNPs from GWAS and significant genes from PrediXcan in R
    * `python 08_sig_SNP_sig_gene.py --SNP_sig 5e-4 --gene_sig 0.05 --input_prefix AMR --pheno_names pheno_names.txt`
      * NOTE: I set these thresholds arbitrarilty low so we have significant genes to with with later. Usually, use the defaults, 5e-8 and 9.654e-6, respectively.

09. Calculate independent significant SNPs in a joint analysis in [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO)
    * a. Make GWAS output into GCTA-COJO format 
      * `python 09a_GEMMA_to_GCTA-COJO.py --fam AMR.fam --GWAS_prefix AMR --output_prefix AMR --pheno_names pheno_names.txt`
    * b. Run GCTA (this is not a script, this is actual GCTA)
      * `gcta64 --cojo-file AMR_pheno1.ma --cojo-slct --cojo-p 5e-4 --bfile AMR --cojo-actual-geno --out AMR_pheno1; gcta64 --cojo-file AMR_pheno2.ma --cojo-slct --cojo-p 5e-4 --bfile AMR --cojo-actual-geno --out AMR_pheno2`
        * Again, set P arbitrarily low for example; set to 5e-8 for real analyses

10. Perform colocalization between GWAS results and eQTL data using [COLOC](https://cran.r-project.org/web/packages/coloc/coloc.pdf) in a [COLOC wrapper](https://github.com/hakyimlab/summary-gwas-imputation)
    * a. Put GWAS and eQTL data into COLOC input format (takes a while due to size of eQTL files)
      * `Rscript 10a_make_COLOC_input.R --ma_prefix AMR --GWAS_prefix AMR --sample_size 347 --pheno_names pheno_names.txt`
    * b. Run COLOC wrapper (takes a while due to size of eQTL files)
      * `python 10b_COLOC_wrapper.py --prefix AMR --GWAS_sample_size 347 --pheno_names pheno_names.txt`

11. Perform backward elimination modeling of all significant genes to find independent signals in R
    * `Rscript 11_back_elim.R --sig_gene output/AMR_sig_genes.txt --pheno pheno_wIID.txt --pheno_name pheno1`

### For use in admixed populations:

12. Prepare genotype data with [PLINK](https://www.cog-genomics.org/plink/1.9/data)
    * a. Merge genotypes with a reference populations (options are African-American (AFA) and Hispanic (HIS)), restrict to only SNPs in the study population
      * `plink --bfile AMR --bmerge /home/angela/Ad_PX_pipe_data/RFMix/RefPop/HIS --make-bed --out AMR_w_ref --extract AMR.bim; plink --bfile AMR --exclude AMR_w_ref-merge.missnp --make-bed --out AMR_for_merge; plink --bfile /home/angela/Ad_PX_pipe_data/RFMix/RefPop/HIS --exclude AMR_w_ref-merge.missnp --make-bed --out HIS_for_merge --extract AMR.bim; plink --bfile AMR_for_merge --bmerge HIS_for_merge --make-bed --out AMR_w_ref --allow-no-sex`
    * b. Order genotype file by population and add cM positions
      * `cat HIS_for_merge.fam AMR_for_merge.fam > merged.fam; plink --bfile AMR_w_ref --indiv-sort file merged.fam  --cm-map /home/angela/Ad_PX_pipe_data/cm_map/genetic_map_chr@_combined_b37.txt --make-bed --out AMR_w_ref_ordered`
    * c. Split input by chr
      * `mkdir -p merged_w_ref/; for i in {1..22}; do plink --bfile AMR_w_ref_ordered --chr ${i} --make-bed --out merged_w_ref/chr${i}; done`

13. Phase haplotypes with [HAPI-UR](https://code.google.com/archive/p/hapi-ur/)
   * `mkdir -p haplotypes/; for i in {1..22}; do /home/angela/Ad_PX_pipe_data/HAPI-UR/hapi-ur -p merged_w_ref/chr${i} -w 64 -o haplotypes/chr${i}; done`

14. Calculate local ancestry inference in [RFMix](https://sites.google.com/site/rfmixlocalancestryinference/)
   * a. Make additional files for RFMix input
   
    
15. Perform admixture mapping in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)

