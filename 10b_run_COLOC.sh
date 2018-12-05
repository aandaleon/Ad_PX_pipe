#wrapper around COLOC
#$1 is sample size (required)
#$2 is prefix (not required)
#pass an argument of sample size
if cd summary-gwas-imputation; then git pull; else git clone https://github.com/hakyimlab/summary-gwas-imputation.git; fi #contains COLOC wrapper, https://stackoverflow.com/questions/15602059/git-shortcut-to-pull-with-clone-if-no-local-there-yet
mkdir -p COLOC_results

tissues=(AFA AFHI ALL CAU HIS Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Breast_Mammary_Tissue Cells_EBV-transformed_lymphocytes Cells_Transformed_fibroblasts Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Liver Lung Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis Thyroid Uterus Vagina Whole_Blood)
tissues_sample_size=(233 578 352 585 1163 298 298 126 197 118 285 72 72 89 103 96 92 81 81 81 82 183 114 272 124 124 127 241 218 159 159 97 278 361 256 85 149 87 87 196 302 77 89 170 157 278 278 278 338)

for tiss in ${!tissues[*]};
do
  /usr/bin/python3 summary-gwas-imputation/src/run_coloc.py -gwas COLOC_input/"$2"GWAS_${tissues[$tiss]}.txt.gz -gwas_sample_size "$1" -eqtl COLOC_input/"$2"eQTL_${tissues[$tiss]}.txt.gz -eqtl_sample_size ${tissues_sample_size[$tiss]} -output COLOC_results/"$2"${tissues[$tiss]}.txt.gz &> COLOC.log #COLOC prints out too much so throw it into a file
  echo Completed with ${tissues[$tiss]}
done
