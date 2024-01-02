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
    EffectFile = read.table("./3_SNP_EFF/SNPeff_Output/EXOME_RP502_Libnames_TF1_Mask_indDP_missing90_MAJORasREFUNRBoot_BialelicOnly_CURATED.txt", h = T)

    #change SNPid with GenPOS
    #EffectFile[,1] = paste0(EffectFile[,2], "_", EffectFile[,3])

    #sanity check
    #EffectFile[,1] == colnames(genos)

    #subset variants in exome from dosage matrix
    genos.2 = genos[,colnames(genos) %in% EffectFile[,1]]

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
    Rxyfun = function(dosage, effectDF, varianType, PopIDs){

    #takes as arguments:
    #dosage = the dosage matrix, individuals in rows and SNPs in columns; rownames and colnames MUST be defined
    #effectDF = dataframe describing the effect of each variant: two columns: SNP identifier & effect; SNP identifier MUST be the same as column names of the dosage matrix
    #varianType = type of mutation, character string;
    #PopIDs = dataframe describing the individuals populations: two columns: INDV name and either one or two: INDV name MUST be the same as rownames from the dosage matrix

    #First subset the data for the type of variant
    dos.sub = dosage[,colnames(dosage) %in% effectDF[effectDF[,2] == varianType, 1]]

    #Then separate it for both populations
    dos.sub.popx = dos.sub[rownames(dos.sub) %in% PopIDs[PopIDs[,2] == 1, 1],]
    dos.sub.popy = dos.sub[rownames(dos.sub) %in% PopIDs[PopIDs[,2] == 2, 1],]

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
    RxyMeanNeutral = Rxyfun(dosage = genos.2, effectDF = EffectFile[,c(1,6)], varianType = "MODIFIER", PopIDs = dtaF[,c(1,3)])
    RxyMeanLow = Rxyfun(dosage = genos.2, effectDF = EffectFile[,c(1,6)], varianType = "LOW", PopIDs = dtaF[,c(1,3)])
    RxyMeanMod = Rxyfun(dosage = genos.2, effectDF = EffectFile[,c(1,6)], varianType = "MODERATE", PopIDs = dtaF[,c(1,3)])
    RxyMeanHigh = Rxyfun(dosage = genos.2, effectDF = EffectFile[,c(1,6)], varianType = "HIGH", PopIDs = dtaF[,c(1,3)])

    #merge the data
    output = as.data.frame(cbind(RxyMeanNeutral = RxyMeanNeutral, RxyMeanLow = RxyMeanLow, RxyMeanMod = RxyMeanMod, RxyMeanHigh = RxyMeanHigh))
    #save
    write.table(output, "./6_MinorAllelesCount/stats/EXOME_Meanvalues_Rxy_islandvsCont_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

    #Get the 100 jackknife

    #create the 100 intervals
    intervals = cut((1:ncol(genos.2)),100)

    #create empty df for output
    output = as.data.frame(matrix(nrow = 0, ncol = 5))
    colnames(output) = c("interval", "RxyMeanNeutral", "RxyMeanLow", "RxyMeanMod", "RxyMeanHigh")

    #Loop through intervals
    for(int in 1:100){

	#subset the dosage matrix WITHOUT the interval
	dos.sub = genos.2[,which(intervals != levels(intervals)[int])]
	#same for the effect for each variant
	EffectFile.sub = EffectFile[which(intervals != levels(intervals)[int]),]

	#rerun the four Rxy
	RxyMeanNeutral.sub = Rxyfun(dosage = dos.sub, effectDF = EffectFile.sub[,c(1,6)], varianType = "MODIFIER", PopIDs = dtaF[,c(1,3)])
	RxyMeanLow.sub = Rxyfun(dosage = dos.sub, effectDF = EffectFile.sub[,c(1,6)], varianType = "LOW", PopIDs = dtaF[,c(1,3)])
	RxyMeanMod.sub = Rxyfun(dosage = dos.sub, effectDF = EffectFile.sub[,c(1,6)], varianType = "MODERATE", PopIDs = dtaF[,c(1,3)])
	RxyMeanHigh.sub = Rxyfun(dosage = dos.sub, effectDF = EffectFile.sub[,c(1,6)], varianType = "HIGH", PopIDs = dtaF[,c(1,3)])

	#merge the fourth stats and intervals
	line = as.data.frame(cbind(interval = int, RxyMeanNeutral = RxyMeanNeutral.sub, RxyMeanLow = RxyMeanLow.sub, RxyMeanMod = RxyMeanMod.sub, RxyMeanHigh = RxyMeanHigh.sub))

	#merge
	output = rbind(output, line)

    }

    #save the jackife results
    write.table(output, "./6_MinorAllelesCount/stats/EXOME_Jackknife_Rxy_islandvsCont_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

    ## Now about LxySQUARED

    #define Rxy function
    Rxy2fun = function(dosage, effectDF, varianType, PopIDs){

    #takes as arguments:
    #dosage = the dosage matrix, individuals in rows and SNPs in columns; rownames and colnames MUST be defined
    #effectDF = dataframe describing the effect of each variant: two columns: SNP identifier & effect; SNP identifier MUST be the same as column names of the dosage matrix
    #varianType = type of mutation, character string;
    #PopIDs = dataframe describing the individuals populations: two columns: INDV name and either one or two: INDV name MUST be the same as rownames from the dosage matrix

    #First subset the data for the type of variant
    dos.sub = dosage[,colnames(dosage) %in% effectDF[effectDF[,2] == varianType, 1]]

    #Then separate it for both populations
    dos.sub.popx = dos.sub[rownames(dos.sub) %in% PopIDs[PopIDs[,2] == 1, 1],]
    dos.sub.popy = dos.sub[rownames(dos.sub) %in% PopIDs[PopIDs[,2] == 2, 1],]

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
    R2xyMeanNeutral = Rxy2fun(dosage = genos.2, effectDF = EffectFile[,c(1,6)], varianType = "MODIFIER", PopIDs = dtaF[,c(1,3)])
    R2xyMeanLow = Rxy2fun(dosage = genos.2, effectDF = EffectFile[,c(1,6)], varianType = "LOW", PopIDs = dtaF[,c(1,3)])
    R2xyMeanMod = Rxy2fun(dosage = genos.2, effectDF = EffectFile[,c(1,6)], varianType = "MODERATE", PopIDs = dtaF[,c(1,3)])
    R2xyMeanHigh = Rxy2fun(dosage = genos.2, effectDF = EffectFile[,c(1,6)], varianType = "HIGH", PopIDs = dtaF[,c(1,3)])

    #merge the data
    output = as.data.frame(cbind(R2xyMeanNeutral = R2xyMeanNeutral, R2xyMeanLow = R2xyMeanLow, R2xyMeanMod = R2xyMeanMod, R2xyMeanHigh = R2xyMeanHigh))
    #save
    write.table(output, "./6_MinorAllelesCount/stats/EXOME_Meanvalues_R2xy_islandvsCont_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

    #Get the 100 jackknife

    #create empty df for output
    output = as.data.frame(matrix(nrow = 0, ncol = 5))
    colnames(output) = c("interval", "R2xyMeanNeutral", "R2xyMeanLow", "R2xyMeanMod", "R2xyMeanHigh")

    #Loop through intervals
    for(int in 1:100){

        #subset the dosage matrix WITHOUT the interval
        dos.sub = genos.2[,which(intervals != levels(intervals)[int])]
        #same for the effect for each variant
        EffectFile.sub = EffectFile[which(intervals != levels(intervals)[int]),]

        #rerun the four Rxy
        RxyMeanNeutral.sub = Rxy2fun(dosage = dos.sub, effectDF = EffectFile.sub[,c(1,6)], varianType = "MODIFIER", PopIDs = dtaF[,c(1,3)])
        RxyMeanLow.sub = Rxy2fun(dosage = dos.sub, effectDF = EffectFile.sub[,c(1,6)], varianType = "LOW", PopIDs = dtaF[,c(1,3)])
        RxyMeanMod.sub = Rxy2fun(dosage = dos.sub, effectDF = EffectFile.sub[,c(1,6)], varianType = "MODERATE", PopIDs = dtaF[,c(1,3)])
        RxyMeanHigh.sub = Rxy2fun(dosage = dos.sub, effectDF = EffectFile.sub[,c(1,6)], varianType = "HIGH", PopIDs = dtaF[,c(1,3)])

        #merge the fourth stats and intervals
        line = as.data.frame(cbind(interval = int, R2xyMeanNeutral = RxyMeanNeutral.sub, R2xyMeanLow = RxyMeanLow.sub, R2xyMeanMod = RxyMeanMod.sub, R2xyMeanHigh = RxyMeanHigh.sub))

        #merge
	output = rbind(output, line)

    }

    #save the jackife results
    write.table(output, "./6_MinorAllelesCount/stats/EXOME_Jackknife_R2xy_islandvsCont_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

EOF
