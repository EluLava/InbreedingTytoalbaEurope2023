#! /bin/bash

cd ../data/3_SNP_EFF
mkdir -p SNPeff_Output/

##### 0 : Configure our Genome ######

# snpEff configure new genome
# https://pcingola.github.io/SnpEff/se_build_db/

# Download all mandatory files
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/691/265/GCF_018691265.1_T.alba_DEE_v4.0/GCF_018691265.1_T.alba_DEE_v4.0_genomic.gtf.gz
mv GCF_018691265.1_T.alba_DEE_v4.0_genomic.gtf.gz genes.gtf.gz
gzip -kd genes.gtf.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/691/265/GCF_018691265.1_T.alba_DEE_v4.0/GCF_018691265.1_T.alba_DEE_v4.0_genomic.fna.gz
mv GCF_018691265.1_T.alba_DEE_v4.0_genomic.fna.gz T.alba_DEE_v4.0.fa.gz
gzip -kd T.alba_DEE_v4.0.fa.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/691/265/GCF_018691265.1_T.alba_DEE_v4.0/GCF_018691265.1_T.alba_DEE_v4.0_rna_from_genomic.fna.gz
mv GCF_018691265.1_T.alba_DEE_v4.0_rna_from_genomic.fna.gz cds.fa.gz
gzip -kd cds.fa.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/691/265/GCF_018691265.1_T.alba_DEE_v4.0/GCF_018691265.1_T.alba_DEE_v4.0_protein.faa.gz
mv GCF_018691265.1_T.alba_DEE_v4.0_protein.faa.gz protein.fa.gz
gunzip -kd protein.fa.gzgunzip -kd protein.fa.gz

# Build the DB for T.alba
module load gcc/10.4.0 openjdk/17.0.3_7 snpeff/5.1d
snpeff build -gtf22 -v T.alba_DEE_v4.0 -c snpEff.config -noCheckProtein

##### 1 : run snpEff #####

# transform VCF rename CHR like annotation
awk '!/#/ {print $1,$7}' ./T.alba_DEE_v4.0/GCF_018691265.1_T.alba_DEE_v4.0_assembly_report.txt > $$
bcftools annotate --rename-chrs $$ ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot.vcf.gz \
> ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed.vcf

cd ./3_SNP_EFF

# Run snpEff
snpEff eff -v T.alba_DEE_v4.0 -c ./snpEff.config -canon ../1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed.vcf \
> ../1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated.vcf

cat ../1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated.vcf | vcfEffOnePerLine.pl > \
../1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated_OneEffPerLine.vcf

#rename CHRs
awk '!/#/ {print $7,$1}' ./T.alba_DEE_v4.0/GCF_018691265.1_T.alba_DEE_v4.0_assembly_report.txt > $$
bcftools annotate --rename-chrs $$ ../1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated_OneEffPerLine.vcf > \
../1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated_OneEffPerLine_renamed.vcf

# get the summary table
snpSift extractFields ../1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated_OneEffPerLine_renamed.vcf \
CHROM POS "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].TRID" \
> ./SNPeff_Output/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated_OneEffPerLine_renamed_BialelicOnly.txt

#In this summary table, we have some SNPs with more than one line (if they have different impacts on different genes) --> We want to keep only the biggest impact

#We have a problem with the output of snp eff because it DOES NOT have constant nb of columns --> just keep the first 5 for the UNIQUE R section
awk '{OFS=" "; for(i=1;i<=5;i++){printf $i"\t"};printf "\n"}' ./SNPeff_Output/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated_OneEffPerLine_renamed_BialelicOnly.txt > \
./SNPeff_Output/TMP_RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated_OneEffPerLine_renamed_BialelicOnly_FIRST5columns.txt

#Now since we added +10 in position for positions with 0 recombination detected we miss the "ANN[*].EFFECT" for these regions, we can just take the previous one
awk '{OFS=" ";if(NF == 5){print $0; previousEff=$4}else{print $1"\t"$2"\t"$3"\t"previousEff"\t"$4}}' ./SNPeff_Output/TMP_RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated_OneEffPerLine_renamed_BialelicOnly_FIRST9columns.txt \
> ./SNPeff_Output/TMP_RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_COMPLETE_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated_OneEffPerLine_renamed_BialelicOnly_FIRST5columnsforRinput.txt

cd ..

#yes the names of VCFs became ridiculously long ...

R --vanilla << EOF

    library(foreach)
    library(doParallel)
    cl = 40
    registerDoParallel(cl)

    #read the table
    dta = read.table("./3_SNP_EFF/SNPeff_Output/TMP_RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_COMPLETE_GenPOSplus10_MajorMinorUNRBoot_CHRrenamed_annotated_OneEffPerLine_renamed_BialelicOnly_FIRST5columnsforRinput.txt", header = T)

    #Add one column where we'll paste SS and POS so that we have one unique identifier per SNP
    dta = as.data.frame(cbind(SNPid=vector(mode="character", length = nrow(dta)), dta))
    #paste CHR and POS
    dta[,1] = paste(dta[,2],dta[,3],sep="_")

    #Get nb of lines per SNP ID
    nbSNPID = table(dta[,1])

    #select SNP IDs which only appear once
    uniqSNP = names(nbSNPID[nbSNPID == 1])
    #Select SNPs which appear more than once
    severalSNP = names(nbSNPID[nbSNPID != 1])

    #Add all uniq identifier to dtaUNIQUE
    dtaUNIQUE = dta[dta[,1] %in% uniqSNP,]

    #Loop through the multiple SNP IDs
    foreachOUTPUT = foreach(snp=1:length(severalSNP), .combine=rbind) %dopar% {

        #print SNP / length(snp) to get an idea where we are
        print(paste0("Processing SNP ", snp, " out of ", length(severalSNP)))

        #extract lines with this SNP effect
        dtasub = dta[dta[,1] == severalSNP[snp],]


        #select SNP with higher effect

        #If only one effect present, randomly select a line
        if(length(unique(dtasub[,6])) == 1){

            #randomly select one
            effectline = dtasub[sample(1:nrow(dtasub),1),]

        #If not just one effect, if "HIGH" effect is present, take that one
        } else if("HIGH" %in% dtasub[,6]){

            #Count how many high we have, if one take this one, else randomly select one
            if(length(dtasub[dtasub[,6] == "HIGH",6]) == 1){

                #select the HIGH row
                effectline = dtasub[dtasub[,6] == "HIGH",]

            } else {

                #randomly select one line with HIGH
                effectline = dtasub[sample(which(dtasub[6] == "HIGH"),1),]

            }

        #If not just one effect, and no HIGH effect is present, take MODERATE effect
        } else if("MODERATE" %in% dtasub[,6]){

            #Count how many moderate we have, if one take this one, else randomly select one
            if(length(dtasub[dtasub[,6] == "MODERATE",6]) == 1){

                #select the MODERATE row
                effectline = dtasub[dtasub[,6] == "MODERATE",]

            } else {

                #randomly select one line with MODERATE
                effectline = dtasub[sample(which(dtasub[6] == "MODERATE"),1),]

            }

        #If not just one effect, and no HIGH nor MODERATE effect is present, take LOW effect
        } else if("LOW" %in% dtasub[,6]){

            #Count how many low we have, if one take this one, else randomly select one
            if(length(dtasub[dtasub[,6] == "LOW",6]) == 1){

                #select the LOW row
                effectline = dtasub[dtasub[,6] == "LOW",]

            } else {

                #randomly select one line with LOW
                effectline = dtasub[sample(which(dtasub[6] == "LOW"),1),]

            }

	}

        #Add the effect line
        return(effectline)

    }

    #merge single SNPs and foreach output
    dtaUNIQUE.2 = rbind(dtaUNIQUE, foreachOUTPUT)

    #Read the list of SNPs we used for ROHs analyses
    dtaCHRPOS = read.table("./3_ROHs/RZooRoH_SNPs_CHR_POSGEN.txt", h = T)

    #Subsample from dtaUNIQUE only the SNPs in the dtaCHRPOS
    dtaUNIQUE.3 = dtaUNIQUE.2[dtaUNIQUE.2[,1] %in% dtaCHRPOS[,3],]

    #Sort the output
    dtaUNIQUE.4 = dtaUNIQUE.3[match(dtaCHRPOS[,3],dtaUNIQUE.3[,1]),]

    #Save the UNIQUE dta
    write.table(unique(dtaUNIQUE.4),"./3_SNP_EFF/SNPeff_Output/RP502_Libnames_TF1_Mask_indDP_missing90_MAJORasREFUNRBoot_BialelicOnly_CURATED.txt", quote = F, col.names = T, row.names = F)

    #Save R env for debuf in case of problems
    save.image("./3_SNP_EFF/SNPeff_Output/REnv_CURATING_MAJORasREFUNRBoot_in_case.RData")

    #read the annotation file from ENSEMBL (for subsampling CDS regions only)
    dtaANN = read.table("./3_SNP_EFF/data/T.alba_DEE_v4.0/TytoAlbaAnnotationGENEsSSnames.txt", h = T)
    #subsample only protein coding genes
    dtaANNcoding = dtaANN[dtaANN[,4] == "protein-coding",]

    #subsample effects (to only get EXOME effects)

    #create output file
    GenesEffectFile = as.data.frame(matrix(nco = ncol(dtaUNIQUE.4), nrow = 0))
    colnames(GenesEffectFile) = colnames(dtaUNIQUE.4)

    #Loop through genes to extract variants ony in protein coding sequence
    for(rw in 1:nrow(dtaANNcoding)){

        #Extract variant within the gene
        variants = dtaUNIQUE.4[((dtaUNIQUE.4[,2] == dtaANNcoding[rw,1]) & (dtaUNIQUE.4[,3] > dtaANNcoding[rw,2]) & (dtaUNIQUE.4[,3] < dtaANNcoding[rw,3])),]

        #Add this to exome file
        GenesEffectFile = rbind(GenesEffectFile, variants)

    }

    #Exclude intron
    ExomeEffectFile = GenesEffectFile[GenesEffectFile[,5] != "intron_variant",]
    #Only keep intron (for R')
    IntronEffectFile = GenesEffectFile[GenesEffectFile[,5] == "intron_variant",]

    #write both
    write.table(ExomeEffectFile, "./3_SNP_EFF/SNPeff_Output/EXOME_RP502_Libnames_TF1_Mask_indDP_missing90_MAJORasREFUNRBoot_BialelicOnly_CURATED.txt", quote = F, col.names = T, row.names = F)
    write.table(IntronEffectFile, "./3_SNP_EFF/SNPeff_Output/INTRONs_RP502_Libnames_TF1_Mask_indDP_missing90_MAJORasREFUNRBoot_BialelicOnly_CURATED.txt", quote = F, col.names = T, row.names = F)

EOF

