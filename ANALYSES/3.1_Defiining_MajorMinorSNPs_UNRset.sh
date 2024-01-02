#! /bin/bash

cd ../data

mkdir -p 3_SNP_EFF

#### We need to filter for the unrelated ONLY to estimate alleleic frequencies ####

#Get a VCF with only the unrelated individuals
bcftools view -S <(awk 'NR > 1 {print $2}' ./inputFILES/RPUNRELATED_LibNames_NewNames.list) -Oz -o ./1.1_VCFs/RP187UNR_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10.vcf.gz \
./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10.vcf.gz

bcftools index ./1.1_VCFs/RP187UNR_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10.vcf.gz

#stopped here

#100 bootstraps to select 6 swiss individuals
for i  in {1..100}; do

    #sample 70 individuals from the list of UNRELATED swiss that we will NOT subsample
    shuf -n 70 <(awk 'NR > 1 {print $0}' ./inputFILES/RP76_CH_UNRELATED_Libnames_Newnames.txt) > ./inputFILES/tmp_70_UNRCH.list

    #subsample the VCF
    bcftools view -S ^<(awk '{print $2}' ./inputFILES/tmp_70_UNRCH.list) -Oz -o ./1.1_VCFs/TMP_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10.vcf.gz \
    ./1.1_VCFs/RP187UNR_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10.vcf.gz

    bcftools index ./1.1_VCFs/TMP_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10.vcf.gz

    #Extract AF
    bcftools query -f '%REF\t%ALT\t%AF\n' ./1.1_VCFs/TMP_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10.vcf.gz > ./3_SNP_EFF/TMP_REF_ALT_AF.txt

    #Retrieve MAJOR allele
    awk '{if($3 > 0.5){print $1}else{print $2}}' ./3_SNP_EFF/TMP_REF_ALT_AF.txt > ./3_SNP_EFF/TMP_Bootstrap_${i}_MAJOR_ALL.txt

done

#rm tmp list of indv file
rm ./inputFILES/tmp_70_UNRCH.list
#rm tmp vcf
rm ./1.1_VCFs/TMP_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10.vcf.gz
#rm AF file
rm ./3_SNP_EFF/TMP_REF_ALT_AF.txt

#### Get major allele in most case

R --vanilla << EOF

    library(foreach)
    library(doParallel)
    cl = 40
    registerDoParallel(cl)

    #read first bootstraps
    dta = read.table("./3_SNP_EFF/TMP_Bootstrap_1_MAJOR_ALL.txt", h = F)

    #loop through the rest of the 100 boostraps
    dta = foreach(i = 1:100, .combine=cbind) %dopar% {

        #read the file
        dta.sub = read.table(paste0("./3_SNP_EFF/TMP_Bootstrap_", i, "_MAJOR_ALL.txt"), h = F)

        #merge with previous df
        return(dta.sub)

    }

    #Create output dataframe
    MAJall = as.data.frame(x = rep(NA, nrow(dta)))
    colnames(MAJall) = "MAJORall"

    #Extract allele most selected as major per line
    MAJall[,1] = apply(dta, 1, function(x){names(sort(table(x), decreasing = T)[1])})

    #Loop through the rows
    output = foreach(i=1:nrow(dta)) %dopar% {

        #Extract allele most selected as major per line
        return(names(sort(table(unlist(dta[i,])), decreasing = T)[1]))

    }

    MAJall[,1] = unlist(output)

    #save the true major allele
    write.table(MAJall, "./3_SNP_EFF/TRUE_MAJOR_ALL.txt", quote = F, col.names = F, row.names = F)

EOF

#### RECODE 0 and 1 allele in VCF ####

#Of course we first need to update SNPID info, we'll juste use CHROM_POS
bcftools query -f '%CHROM\t%POS\n' ./1.1_VCFs/RP187UNR_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10.vcf.gz | awk '{print $1"\t"$2"\t"$1"_"$2}' > ./3_SNP_EFF/CHROM_POS_SNPsIDs.txt

#The merge with MAJ ALL
R --vanilla << EOF

    dta1 = read.table("./3_SNP_EFF/TRUE_MAJOR_ALL.txt", h = F)
    dta2 = read.table("./3_SNP_EFF/CHROM_POS_SNPsIDs.txt", h = F)

    #merge
    dta = cbind(dta2, dta1)

    #save
    write.table(dta[,3:4], "./3_SNP_EFF/MAJOR_ALL_RECODING.txt", quote = F, col.names = F, row.names = F)

EOF

#rm tmp files
rm ./3_SNP_EFF/TRUE_MAJOR_ALL.txt

#zip the annotation file
bgzip ./3_SNP_EFF/CHROM_POS_SNPsIDs.txt
#tabix the annotated file
tabix -p vcf ./3_SNP_EFF/CHROM_POS_SNPsIDs.txt.gz

#add SNPID into the VCF
bcftools annotate -c CHROM,POS,ID -a ./3_SNP_EFF/CHROM_POS_SNPsIDs.txt.gz -Oz -o ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_withSNPsIDs.vcf.gz \
./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10.vcf.gz

#recoding allelles MAJ MIN
plink --vcf ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_withSNPsIDs.vcf.gz --recode vcf-iid --reference-allele ./3_SNP_EFF/MAJOR_ALL_RECODING.txt \
--allow-extra-chr --out ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot

#zip vcf
bcftools view -Oz -o ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot.vcf.gz ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot.vcf
bcftools index ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot.vcf.gz

#rm unzipped
rm ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot.vcf
rm ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot.log
rm ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot.nosex


