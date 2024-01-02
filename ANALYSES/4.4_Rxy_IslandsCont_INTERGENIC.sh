#! /bin/bash

cd ../data

#Open R
R --vanilla << EOF

    library(gaston)
    library(hierfstat)
    library(SNPRelate)

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
    EffectFile_introns = read.table("./2_SNP_EFF/SNPeff_Output/INTRONs_RP502_Libnames_TF1_Mask_indDP_missing90_MAJORasREFUNRBoot_BialelicOnly_CURATED.txt", h = T)
    EffectFile_exons = read.table("./2_SNP_EFF/SNPeff_Output/EXOME_RP502_Libnames_TF1_Mask_indDP_missing90_MAJORasREFUNRBoot_BialelicOnly_CURATED.txt", h = T)

    #merge both
    EffectFile = rbind(EffectFile_introns, EffectFile_exons)

    #subset variants NOT in exome from dosage matrix
    genos.2 = genos[,!(colnames(genos) %in% EffectFile[,1])]

    #INDvs POPs
    dtaF = read.table("./Analyses/GeneticPositionsFHBDs_FHOM_FROH_FASs_FPED_FAsUNR_RP502.txt", h = T)[,c(1,7)]
    colnames(dtaF) = c("INDVs", "Population")

    #add Continental info
    dtaF = as.data.frame(cbind(dtaF, Continent=vector(length = nrow(dtaF))))
    dtaF[,3][dtaF[,2] %in% c("CH","DK","FR","GE","GR","IS","IT","MA","PT","SB")] = 1
    dtaF[,3][dtaF[,2] %in% c("AE","CO","CT","CY","EC","GB","IO","IR","WC")] = 2
    #Pass Continent into factor
    dtaF[,3] = factor(dtaF[,3], levels = c(1,2))


    #define Rxy function
    Rxyfun = function(dosage, PopIDs){

    #takes as arguments:
    #dosage = the dosage matrix, individuals in rows and SNPs in columns; rownames and colnames MUST be defined
    #PopIDs = dataframe describing the individuals populations: two columns: INDV name and either one or two: INDV name MUST be the same as rownames from the dosage matrix

    #Then separate it for both populations
    dos.sub.popx = dosage[rownames(dosage) %in% PopIDs[PopIDs[,2] == 1, 1],]
    dos.sub.popy = dosage[rownames(dosage) %in% PopIDs[PopIDs[,2] == 2, 1],]

    #get the fraction of minor allele in pop1 per locus
    fracdeletpopx = colSums(dos.sub.popx, na.rm = T)/(2*nrow(dos.sub.popx))
    #get fraction major allele in pop2 per locus
    fracnondeletpopy = 1 - (colSums(dos.sub.popy, na.rm = T)/(2*nrow(dos.sub.popy)))
    #get Lxy: sum of multiplication of both
    Lxy = sum((fracdeletpopx * fracnondeletpopy))

    #same for pop2:
    #get the fraction of minor allele in pop2 per locus
    fracdeletpopy = colSums(dos.sub.popy, na.rm = T)/(2*nrow(dos.sub.popy))
    #get fraction major allele in pop1 per locus
    fracnondeletpopx = 1 - colSums(dos.sub.popx, na.rm = T)/(2*nrow(dos.sub.popx))
    #get Lyx: sum of multiplication of both
    Lyx = sum((fracdeletpopy * fracnondeletpopx))

    #Ratio
    Rxy = Lxy/Lyx

    #return ratio
    return(Rxy)

    }


    #Get the full Rxy
    RxyMeanNeutral = Rxyfun(dosage = genos.2, PopIDs = dtaF[,c(1,3)])

    #merge the data
    output = as.data.frame(cbind(RxyMeanSYN = RxyMeanNeutral))
    #save
    write.table(output, "./6_MinorAllelesCount/stats/EXOME_Meanvalues_Rxy_INTERGENICSNPs_islandvsCont_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

    #Get the 100 jackknife

    #create the 100 intervals
    intervals = cut((1:ncol(genos.2)),100)

    #create empty df for output
    output = as.data.frame(matrix(nrow = 0, ncol = 2))
    colnames(output) = c("interval", "RxyMeanSYN")

    #Loop through intervals
    for(int in 1:100){

	#subset the dosage matrix WITHOUT the interval
	dos.sub = genos.2[,which(intervals != levels(intervals)[int])]
	#same for the effect for each variant
	EffectFile.sub = EffectFile[which(intervals != levels(intervals)[int]),]

	#rerun the four Rxy
	RxyMeanNeutral.sub = Rxyfun(dosage = dos.sub, PopIDs = dtaF[,c(1,3)])

	#merge the fourth stats and intervals
	line = as.data.frame(cbind(interval = int, RxyMeanSYN = RxyMeanNeutral.sub))

	#merge
	output = rbind(output, line)

    }

    #save the jackife results
    write.table(output, "./6_MinorAllelesCount/stats/EXOME_Jackknife_Rxy_INTERGENICSNPs_islandvsCont_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

    ## Now about LxySQUARED

    #define Rxy function
    Rxy2fun = function(dosage, PopIDs){

    #takes as arguments:
    #dosage = the dosage matrix, individuals in rows and SNPs in columns; rownames and colnames MUST be defined
    #PopIDs = dataframe describing the individuals populations: two columns: INDV name and either one or two: INDV name MUST be the same as rownames from the dosage matrix

    #Then separate it for both populations
    dos.sub.popx = dosage[rownames(dosage) %in% PopIDs[PopIDs[,2] == 1, 1],]
    dos.sub.popy = dosage[rownames(dosage) %in% PopIDs[PopIDs[,2] == 2, 1],]

    dx = colSums(dos.sub.popx, na.rm = T)
    nx = 2*nrow(dos.sub.popx)
    dy = colSums(dos.sub.popy, na.rm = T)
    ny = 2*nrow(dos.sub.popy)

    #get the fraction of minor allele in pop1 per locus
    fracdeletpopx = (2*dx*(nx - dx))/(nx*(nx - 1))
    #get fraction major allele in pop2 per locus
    fracnondeletpopy = 1 - ((2*dy*(ny - dy))/(ny*(ny - 1)))
    #get Lxy: sum of multiplication of both
    L2xy = sum((fracdeletpopx * fracnondeletpopy))

    #same for pop2

    #get the fraction of minor allele in pop1 per locus
    fracdeletpopy = (2*dy*(ny - dy))/(ny*(ny - 1))
    #get fraction major allele in pop2 per locus
    fracnondeletpopx = 1 - ((2*dx*(nx - dx))/(nx*(nx - 1)))
    #get Lxy: sum of multiplication of both
    L2yx = sum((fracdeletpopy * fracnondeletpopx))

    #Ratio
    R2xy = L2xy/L2yx

    #return ratio
    return(R2xy)

    }

    #Get the full Rxy
    R2xyMeanNeutral = Rxy2fun(dosage = genos.2, PopIDs = dtaF[,c(1,3)])

    #merge the data
    output = as.data.frame(cbind(R2xyMeanSYN = R2xyMeanNeutral))
    #save
    write.table(output, "./6_MinorAllelesCount/stats/EXOME_Meanvalues_R2xy_INTERGENICSNPs_islandvsCont_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

    #Get the 100 jackknife

    #create empty df for output
    output = as.data.frame(matrix(nrow = 0, ncol = 2))
    colnames(output) = c("interval", "R2xyMeanSYN")

    #Loop through intervals
    for(int in 1:100){

        #subset the dosage matrix WITHOUT the interval
        dos.sub = genos.2[,which(intervals != levels(intervals)[int])]
        #same for the effect for each variant
        EffectFile.sub = EffectFile[which(intervals != levels(intervals)[int]),]

        #rerun the four Rxy
        RxyMeanNeutral.sub = Rxy2fun(dosage = dos.sub, effectDF = EffectFile.sub[,c(1,5)], varianType = "intron_variant", PopIDs = dtaF[,c(1,3)])

        #merge the fourth stats and intervals
        line = as.data.frame(cbind(interval = int, R2xyMeanSYN = RxyMeanNeutral.sub))

        #merge
	output = rbind(output, line)

    }

    #save the jackife results
    write.table(output, "./6_MinorAllelesCount/stats/EXOME_Jackknife_R2xy_INTERGENICSNPs_islandvsCont_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

EOF
