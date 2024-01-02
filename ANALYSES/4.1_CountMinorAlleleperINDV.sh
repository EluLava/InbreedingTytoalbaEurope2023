#! /bin/bash

cd ../data

#from vcf to bed for dosage conversion
plink --vcf ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot.vcf.gz --make-bed --allow-extra-chr --keep-allele-order --out \
./1.2_BED/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot

#Open R
R --vanilla << EOF

    library(gaston)
    library(hierfstat)
    library(SNPRelate)

    #from bed to gds
    snpgdsBED2GDS("./1.2_BED/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot",
		  out.gdsfn = "1.2_BED/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot.gds", cvt.chr = "char", cvt.snpid = "auto")

    #dosage data
    genos = snpgdsGetGeno('./1.2_BED/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_GenPOSplus10_MajorMinorUNRBoot.gds')

    #read indvs list
    INDVs = as.vector(read.table("6_MinorAllelesCount/INDVs.list", h = F)[,1])
    #read SNPS list
    SNPs = as.vector(read.table("6_MinorAllelesCount/SNPs.list", h = F)[,1])

    #col and row names in genos matrix
    rownames(genos) = INDVs
    colnames(genos) = SNPs

    #read the effect file
    EffectFile = read.table("./3_SNP_EFF/SNPeff_Output/EXOME_RP502_Libnames_TF1_Mask_indDP_missing90_MAJORasREFUNRBoot_BialelicOnly_CURATED.txt", h = T)

    #change SNPid with GenPOS
    #EffectFile[,1] = paste0(EffectFile[,2], "_", EffectFile[,3])

    #subsample dosage columns so that it matches exome SNPs
    genos.2 = genos[,colnames(genos) %in% EffectFile[,1]]

    #sanity check
    unique(EffectFile[,1] == colnames(genos.2))

    #for each neutral minor variant, count the number of minor all per INDV
    NeutralSum = rowSums(genos.2[,colnames(genos.2) %in% EffectFile[EffectFile[,6] == "MODIFIER",1]], na.rm = T)
    #for each highly deleterious minor variant, count the number of minor all per INDV
    HighDeletSum = rowSums(genos.2[,colnames(genos.2) %in% EffectFile[EffectFile[,6] == "HIGH",1]], na.rm = T)
    #for each moderately deleterious minor variant, count the number of minor all per INDV
    ModerDeletSum = rowSums(genos.2[,colnames(genos.2) %in% EffectFile[EffectFile[,6] == "MODERATE",1]], na.rm = T)
    #for each lowly deleterious minor variant, count the number of minor all per INDV
    LowDeletSum = rowSums(genos.2[,colnames(genos.2) %in% EffectFile[EffectFile[,6] == "LOW",1]], na.rm = T)

    #merge
    data = merge(NeutralSum, HighDeletSum, by = 0)
    colnames(data) = c("INDVs", "Neutral", "HighDelet")
    data.2 = merge(ModerDeletSum, LowDeletSum, by = 0)
    colnames(data.2) = c("INDVs", "ModerDelet", "LowDelet")
    data.3 = merge(data, data.2, by = "INDVs")
    data.3[,1] = as.character(data.3[,1])

    #save the counts
    write.table(data.3, "./Analyses/CountMinorAllelesperCategoryperINDV_EXOME_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

    #Now count how many are in homozygous form !

    #for each neutral minor variant, count the number of minor all per INDV
    NeutralSum = rowSums(genos.2[,colnames(genos.2) %in% EffectFile[EffectFile[,6] == "MODIFIER",1]] == 2, na.rm = T)
    #for each highly deleterious minor variant, count the number of minor all per INDV
    HighDeletSum = rowSums(genos.2[,colnames(genos.2) %in% EffectFile[EffectFile[,6] == "HIGH",1]] == 2, na.rm = T)
    #for each moderately deleterious minor variant, count the number of minor all per INDV
    ModerDeletSum = rowSums(genos.2[,colnames(genos.2) %in% EffectFile[EffectFile[,6] == "MODERATE",1]] == 2, na.rm = T)
    #for each lowly deleterious minor variant, count the number of minor all per INDV
    LowDeletSum = rowSums(genos.2[,colnames(genos.2) %in% EffectFile[EffectFile[,6] == "LOW",1]] == 2, na.rm = T)

    #merge
    data = merge(NeutralSum, HighDeletSum, by = 0)
    colnames(data) = c("INDVs", "Neutral", "HighDelet")
    data.2 = merge(ModerDeletSum, LowDeletSum, by = 0)
    colnames(data.2) = c("INDVs", "ModerDelet", "LowDelet")
    data.homozyg = merge(data, data.2, by = "INDVs")
    data.homozyg[,1] = as.character(data.homozyg[,1])

    #save the counts
    write.table(data.homozyg, "./Analyses/CountHomozygMinorAllelesperCategoryperINDV_EXOME_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

EOF
