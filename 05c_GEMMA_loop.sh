#!/bin/bash
#I basically copied the structure of the GWAS QC pipeline, thanks Ryan - https://github.com/WheelerLab/gwasqc_pipeline/blob/master/shellscripts/01MissingnessFiltering

GEMMAPath=/usr/local/bin/gemma
hackedDefault=False #set to normal GWAS GEMMA
hacked=False
while :
do
    case "$1" in
      -g) #genotype file prefix
	        GenoFile="$2"
	        shift 2
	        ;;
      -p) #phenotype file (w/o IDs)
	        PhenoFile="$2"  
	        shift 2
	        ;;
      -a) #annotation file prefix
        	AnnoFile="$2"
	        shift 2
	        ;;
      -k) #relatedness file (w/o IDs)
      	  	RelFile="$2"
	        shift 2
	        ;;
      -c) #covariates file (w/o IDs)
      	 	CovFile="$2"
	        shift 2
	        ;;
      -o) #output prefix
      	 	Prefix="$2"
	        shift 2
	        ;;
      -h) #switch to "hacked" mode
		hacked=True
		shift 1
		;;
      -gp) #path to GEMMA
		GEMMAPath="$2"
		shift 2
		;;
      -*) #unknown 
      		echo "Error: Unknown option: $1" >&2
	        exit 1
	        ;;
      *)  # No more options
         	shift
	        break
	        ;;
     esac #why is bash syntax so weird
done

case "${hacked:=$hackedDefault}" in
	True) #for "hacked" tissue GEMMA
		declare -a tissues=("AFA" "AFHI" "ALL" "CAU" "HIS" "Adipose_Subcutaneous" "Adipose_Visceral_Omentum" "Adrenal_Gland" "Artery_Aorta" "Artery_Coronary" "Artery_Tibial" "Brain_Anterior_cingulate_cortex_BA24" "Brain_Caudate_basal_ganglia" "Brain_Cerebellar_Hemisphere" "Brain_Cortex" "Brain_Frontal_Cortex_BA9" "Brain_Hippocampus" "Brain_Hypothalamus" "Brain_Nucleus_accumbens_basal_ganglia" "Brain_Putamen_basal_ganglia" "Breast_Mammary_Tissue" "Cells_EBV-transformed_lymphocytes" "Cells_Transformed_fibroblasts" "Colon_Sigmoid" "Colon_Transverse" "Esophagus_Gastroesophageal_Junction" "Esophagus_Mucosa" "Esophagus_Muscularis" "Heart_Atrial_Appendage" "Heart_Left_Ventricle" "Liver" "Lung" "Muscle_Skeletal" "Nerve_Tibial" "Ovary" "Pancreas" "Pituitary" "Prostate" "Skin_Not_Sun_Exposed_Suprapubic" "Skin_Sun_Exposed_Lower_leg" "Small_Intestine_Terminal_Ileum" "Spleen" "Stomach" "Testis" "Thyroid" "Uterus" "Vagina" "Whole_Blood") 
		for tissue in "${tissues[@]}"; #iterate through list of tissues
		do
			echo "${GEMMAPath}" -g "${GenoFile}""${tissue}".txt -p "${PhenoFile}" -k "${RelFile}" -c "${CovFile}" -lmm 4 -notsnp -o "${Prefix}""${tissue}"
			"${GEMMAPath}" -g "${GenoFile}""${tissue}".txt -p "${PhenoFile}" -k "${RelFile}" -c "${CovFile}" -lmm 4 -notsnp -o "${Prefix}""${tissue}"
		done
		;;
	False) #for regular GEMMA
		for chr in {1..22}; #iterate through each chromosome
		do
			echo "${GEMMAPath}" -g "${GenoFile}""${chr}".txt.gz -p "${PhenoFile}" -a "${AnnoFile}" -k "${RelFile}" -c "${CovFile}" -lmm 4 -o "${Prefix}"chr"${chr}"
			"${GEMMAPath}" -g "${GenoFile}""${chr}".txt.gz -p "${PhenoFile}" -a "${AnnoFile}""${chr}".txt -k "${RelFile}" -c "${CovFile}" -lmm 4 -o "${Prefix}"chr"${chr}"
		done
		;;
esac
