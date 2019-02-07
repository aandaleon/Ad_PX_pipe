# Ad_PX_pipe
This repository reorganizes and restructures scripts from [Genetically regulated gene expression underlies lipid traits in Hispanic cohorts](https://github.com/WheelerLab/px_his_chol) to be more user-friendly. All data input (PLINK .bed/.bim/.fam, covariates w/ and w/o IIDs, phenotypes w/ and w/o IIDs, and phenotype names) must be **exactly** in the same format as the example data, and all missing data must be replaced with "NA". All paths to softwares are defaulted to those on wheelerlab3, and it is expected that all scripts are run from the same directory. For exact details on the inner workings of each script, use the `--help` flag or see the [wiki](https://github.com/aandaleon/Ad_PX_pipe/wiki). We will be using genotypes from 1000 Genomes American superpopulation (`AMR.bed/.bim/.fam`) and simulating phenotypes and covariances in R. The input population does not need to be admixed, but local ancestry steps (12-15) are for use mainly in admixed (African-American, Hispanic) populations. You should be able to copy and paste all the commands from this README and everything should run... but will you understand it?

For much more detail on the process of everything in here, please see the [wiki](https://github.com/aandaleon/Ad_PX_pipe/wiki), I spent a lot of time writing it :).

00. Produce phenotypes and covariates (ex. medicines) in R (for test data only)
    * `Rscript 00_simulate_pheno_covar.R --bfile AMR`

### For use in all populations:

01. Perform quality control in [PLINK](https://www.cog-genomics.org/plink/1.9/filter) using [gwasqc_pipeline](https://github.com/WheelerLab/gwasqc_pipeline)
    * Test data has already been filtered

02. Calculate principal components and a relationship matrix in [KING](http://people.virginia.edu/~wc9c/KING/manual.html)
    * `Rscript 02_relate_matrix_PCs.R --bfile AMR`
 
03. Impute data with the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!)
    * Test data is already imputed
    * Follow instructions in 03_Michigan_Imputation_Server.pptx

04. Convert imputed data to PrediXcan dosages using [this script](https://github.com/WheelerLab/Imputation/blob/master/UMich_vcf2pxfixCAAPA.py).

    * Since test data is already imputed in PLINK format, run this (for example data only)
      * `mkdir dosages/; cd dosages/; wget https://raw.githubusercontent.com/hakyimlab/PrediXcan/master/Software/convert_plink_to_dosage.py; awk '{ print $1"\t"$2 }' ../AMR.fam > samples.txt; python convert_plink_to_dosage.py -b ../AMR -p /usr/local/bin/plink; cd ..; echo "Completed converting dosages."`

05. Perform a genome wide association study in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)
    * a. Convert from PrediXcan dosage format to BIMBAM format (GEMMA genotype input)
      * `python 05a_PrediXcan_dosages_to_GEMMA.py --dosage_suffix .txt.gz`
    * b. Make covariance file from known covariates and KING PCs
      * `Rscript 05b_make_GEMMA_covars.R --pcs_num 5`
    * c. Run all chrs. in a loop
      * `python 05c_GEMMA_wrapper.py --GWAS`

06. Calculate predicted gene expressions in [PrediXcan](https://github.com/hakyimlab/PrediXcan) using [GTEx](http://predictdb.org/) and [MESA](https://github.com/aandaleon/DivPop) models
    * `python 06_make_pred_exp.py`

07. Perform an imputed transcriptome-based association study in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)
    * a. Convert predicted expression to GEMMA-style pseudo-genotypes
      * `python 07a_convert_PrediXcan_to_GEMMA.py`
    * b. Run all pops. and tissues in a loop
      * `python 05c_GEMMA_wrapper.py --pred_exp`    

08. Find significant SNPs from GWAS and significant genes from PrediXcan in R
    * `python 08_sig_SNP_sig_gene.py --SNP_sig 5e-6 --gene_sig 5e-3`
      * NOTE: I set these thresholds arbitrarilty low so we have significant genes to with with later. Usually, use the defaults, 5e-8 and 9.654e-6, respectively.

09. Calculate independent significant SNPs in a joint analysis in [GCTA-COJO](https://cnsgenomics.com/software/gcta/#COJO)
    * a. Make GWAS output into GCTA-COJO format 
      * `python 09a_GEMMA_to_GCTA-COJO.py --fam AMR.fam`
    * b. Run GCTA (see [manual](https://cnsgenomics.com/software/gcta/#COJO))
      * `gcta64 --cojo-file pheno1.ma --cojo-slct --cojo-p 5e-6 --bfile AMR --cojo-actual-geno --out pheno1; gcta64 --cojo-file pheno2.ma --cojo-slct --cojo-p 5e-6 --bfile AMR --cojo-actual-geno --out pheno2`
        * Again, set P arbitrarily low for example; set to 5e-8 for real analyses

10. Perform colocalization between GWAS results and eQTL data using [COLOC](https://cran.r-project.org/web/packages/coloc/coloc.pdf) in a [COLOC wrapper](https://github.com/hakyimlab/summary-gwas-imputation)
    * a. Put GWAS and eQTL data into COLOC input format (takes a while due to size of eQTL files)
      * `Rscript 10a_make_COLOC_input.R --sample_size 320`
    * b. Run COLOC wrapper (takes a while due to size of eQTL files)
      * `python 10b_COLOC_wrapper.py --sample_size 320`

11. Perform backward elimination modeling of all significant genes to find independent signals in R
    * `Rscript 11_back_elim.R --sig_gene output/pheno1_sig_genes.txt --pheno pheno_wIID.txt --pheno_name pheno1`

### For use in admixed populations:

12. Prepare genotype data with [PLINK](https://www.cog-genomics.org/plink/1.9/data)
    * a. Merge genotypes with a reference populations (options are African-American (AFA) and Hispanic (HIS)), restrict to only SNPs in the study population
      * `plink --bfile AMR --bmerge /home/angela/Ad_PX_pipe_data/RFMix/RefPop/HIS --make-bed --out AMR_w_ref --extract AMR.bim; plink --bfile AMR --exclude AMR_w_ref-merge.missnp --make-bed --out AMR_for_merge; plink --bfile /home/angela/Ad_PX_pipe_data/RFMix/RefPop/HIS --exclude AMR_w_ref-merge.missnp --make-bed --out HIS_for_merge --extract AMR.bim; plink --bfile AMR_for_merge --bmerge HIS_for_merge --make-bed --out AMR_w_ref --allow-no-sex`
    * b. Order genotype file by population and add cM positions
      * `cat HIS_for_merge.fam AMR_for_merge.fam > merged.fam; plink --bfile AMR_w_ref --indiv-sort file merged.fam  --cm-map /home/angela/Ad_PX_pipe_data/cm_map/genetic_map_chr@_combined_b37.txt --make-bed --out AMR_w_ref_ordered`
    * c. Split input by chr
      * `mkdir -p merged_w_ref/; for i in {1..22}; do plink --bfile AMR_w_ref_ordered --chr ${i} --make-bed --out merged_w_ref/chr${i}; done`

13. Phase haplotypes with [HAPI-UR](https://code.google.com/archive/p/hapi-ur/). Change value of `-w` depending on parameters described in the [manual](https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/hapi-ur/hapi-ur-manual-09_27_2012.pdf)
    * `mkdir -p haplotypes/; for i in {1..22}; do /home/angela/Ad_PX_pipe_data/HAPI-UR/hapi-ur -p merged_w_ref/chr${i} -w 64 -o haplotypes/chr${i}; done`

14. Calculate local ancestry inference in [RFMix](https://sites.google.com/site/rfmixlocalancestryinference/)
    * a. Make additional files for RFMix input. When making classes, choose a reference population (options are African-American (AFA) and Hispanic (HIS)).
      * `for i in {1..22}; do awk '{print $3}' haplotypes/chr${i}.phsnp > haplotypes/chr${i}.snp_locations; done; Rscript 14a_make_classes.R haplotypes/chr22.phind HIS`
      * Note: the sample data is rather small, so some chromosomes may not output
    * b. Run RFMix to estimate local ancestry for all individuals. This process takes a very long time, so we will only run chr. 22 for this example. Push chrs. to multiple cores when running real data.
      * `mkdir -p RFMix/; cd /home/angela/Ad_PX_pipe_data/RFMix/; python RunRFMix.py -e 2 -w 0.2 --num-threads 10 --use-reference-panels-in-EM --forward-backward  PopPhased /home/angela/Ad_PX_pipe/haplotypes/chr22.phgeno /home/angela/Ad_PX_pipe/RFMix.classes  /home/angela/Ad_PX_pipe/haplotypes/chr22.snp_locations -o /home/angela/Ad_PX_pipe/RFMix/chr22.rfmix`
    
15. Perform admixture mapping in [GEMMA](http://www.xzlab.org/software/GEMMAmanual.pdf)
    * a. Convert RFMix output into an intermediate for GEMMA input downstream
      * `python 15a_RFMix_to_BIMBAM.py --Viterbi RFMix/chr22.rfmix.2.Viterbi.txt --phind haplotypes/chr22.phind --fam AMR.fam --phsnp haplotypes/chr22.phsnp --output_prefix chr22`
    * b. Run phenotype associations with each local ancestry in GEMMA
      * `python 15b_loc_anc_wrapper.py --snptable loc_anc_input/chr22.csv --refpop HIS --chr 22`
