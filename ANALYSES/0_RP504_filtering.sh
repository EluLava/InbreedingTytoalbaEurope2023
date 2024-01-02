#!/bin/bash

#save path shortcut
myPATH=/work/FAC/FBM/DEE/jgoudet/barn_owl/elavanc1/ChapterIII/data/

mkdir -p ${myPATH}/1.1_VCFs

#First filter is to remove the two indivduals excluded from REF panel
bcftools view -S ${myPATH}/inputFILES/RP502_Libnames.list -O z -o ${myPATH}/1.1_VCFs/RP502_Libnames_TF1_Mask_indDP.vcf.gz ${myPATH}/1.1_VCFs/RP504_Libnames_TF1_Mask_indDP.vcf.gz

#Second filter is for bi-allelic SNPs (min and max numbers of alleles = 2) only and MAC 3
bcftools view -v snps -m 2 -M 2 -c 3:minor --threads 10 -O z -o ${myPATH}/1.1_VCFs/RP502_Libnames_TF1_Mask_indDP_biall_mac3.vcf.gz ${myPATH}/1.1_VCFs/RP502_Libnames_TF1_Mask_indDP.vcf.gz

#Third filter is on missing data (only sites will less than 10% of missing data are kept
bcftools view -i 'F_PASS(GT!="mis") > 0.9' --threads 10 -O z -o ${myPATH}/1.1_VCFs/RP502_Libnames_TF1_Mask_indDP_biall_mac3_missing90.vcf.gz ${myPATH}/1.1_VCFs/RP502_Libnames_TF1_Mask_indDP_biall_mac3.vcf.gz

#Then we want to rename INDIVIDUALs so that they have the FINAL name (!= Libname)
bcftools reheader -s ${myPATH}/inputFILES/RP502_LibNames_NewNames.list -o ${myPATH}/1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90.vcf.gz ${myPATH}/1.1_VCFs/RP502_Libnames_TF1_Mask_indDP_biall_mac3_missing90.vcf.gz

#index the final filtered file
bcftools index ${myPATH}/1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90.vcf.gz