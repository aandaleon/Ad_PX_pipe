This optional subfolder creates plots from the results of scripts from Ad_PX_pipe. For plots related to quality control, please refer to [Ryan's GWAS QC pipeline](https://github.com/RyanSchu/gwasqc_pipeline/wiki).

In this lab, we also use many external plotting tools, such as:
* [LocusZoom](http://locuszoom.org/)
* [Geography of Genetic Variants Browser](https://popgen.uchicago.edu/ggv/?data=%221000genomes%22&chr=11&pos=116663707)

The plots in this folder include:

* P1. Parallel coordinates plot, with principal components on the x-axis and eigenvalues on the y-axis for the first 5 PCs calculated by KING
  * `Rscript P1_rotation_plot.R`
  
![](https://i.imgur.com/ThDW8b4.png)

* qqplot
* manhattan plot - SNP
* manhattan plot - gene
* predicted expression vs. phenotype
* admixture mapping (AFA and HIS separate?)


