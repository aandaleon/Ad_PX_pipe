#!/bin/bash
GEMMAPath=/usr/local/bin
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
      -*) #unknown 
      		echo "Error: Unknown option: $1" >&2
	        exit 1
	        ;;
      *)  # No more options
         	shift
	        break
	        ;;
     esac
done

#Run gemma loop
for chr in {1..22};
do
	echo "${GEMMAPath}"/gemma -g "${GenoFile}""${chr}".txt.gz -p "${PhenoFile}" -a "${AnnoFile}" -k "${RelFile}" -c "${CovFile}" -lmm 4 -o "${Prefix}"chr"${chr}"
	"${GEMMAPath}"/gemma -g "${GenoFile}""${chr}".txt.gz -p "${PhenoFile}" -a "${AnnoFile}""${chr}".txt -k "${RelFile}" -c "${CovFile}" -lmm 4 -o "${Prefix}"chr"${chr}"
done