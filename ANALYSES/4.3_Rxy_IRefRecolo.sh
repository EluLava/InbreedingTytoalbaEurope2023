#! /bin/bash

#SBATCH --mail-user eleonore.lavanchy@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --time 11:59:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 64G
#SBATCH --job-name RXYstats
#SBATCH -o STD/%x_%j.stdout
#SBTACH -e STD/%x_%j.stderr
#SBATCH --export=NONE
#SBATCH --account jgoudet_barn_owl

set -e
set -x

cd /users/elavanc1/BARNOWL/ChapterIII/data

module load gcc r

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
    EffectFile = read.table("./2_SNP_EFF/SNPeff_Output/EXOME_RP502_Libnames_TF1_Mask_indDP_missing90_MAJORasREFUNRBoot_BialelicOnly_CURATED.txt", h = T)

    #change SNPid with GenPOS
    #EffectFile[,1] = paste0(EffectFile[,2], "_", EffectFile[,3])

    #sanity check
    #EffectFile[,1] == colnames(genos)

    genos.2 = genos[,colnames(genos) %in% EffectFile[,1]]

    #INDvs POPs
    dtaF = read.table("./Analyses/GeneticPositionsFHBDs_FHOM_FROH_FASs_FPED_FAsUNR_RP502.txt", h = T)[,c(1,7)]
    colnames(dtaF) = c("INDVs", "Population")

    #add Continental info
    dtaF = as.data.frame(cbind(dtaF, Refugium=vector(length = nrow(dtaF))))
    dtaF[,3][dtaF[,2] %in% c("PT","IT","IS","MA","GR")] = 1
    dtaF[,3][dtaF[,2] %in% c("CH","DK","FR","GE","SB")] = 2
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
    dos.sub.pop1 = dos.sub[rownames(dos.sub) %in% PopIDs[PopIDs[,2] == 1, 1],]
    dos.sub.pop2 = dos.sub[rownames(dos.sub) %in% PopIDs[PopIDs[,2] == 2, 1],]

    #get the fraction of minor allele in pop1 per locus
    fracdeletpop1 = colSums(dos.sub.pop1, na.rm = T)/(2*nrow(dos.sub.pop1))
    #get fraction major allele in pop2 per locus
    fracnondeletpop2 = 1 - colSums(dos.sub.pop2, na.rm = T)/(2*nrow(dos.sub.pop2))
    #get Lxy: sum of multiplication of both
    Lxy = sum((fracdeletpop1 * fracnondeletpop2))

    #same for pop2:
    #get the fraction of minor allele in pop2 per locus
    fracdeletpop2 = colSums(dos.sub.pop2, na.rm = T)/(2*nrow(dos.sub.pop2))
    #get fraction major allele in pop1 per locus
    fracnondeletpop1 = 1 - colSums(dos.sub.pop1, na.rm = T)/(2*nrow(dos.sub.pop1))
    #get Lyx: sum of multiplication of both
    Lyx = sum((fracdeletpop2 * fracnondeletpop1))

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
    write.table(output, "./6_MinorAllelesCount/stats/EXOME2_Meanvalues_Rxy_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

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
    write.table(output, "./6_MinorAllelesCount/stats/EXOME2_Jackknife_Rxy_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

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
    dos.sub.pop1 = dos.sub[rownames(dos.sub) %in% PopIDs[PopIDs[,2] == 1, 1],]
    dos.sub.pop2 = dos.sub[rownames(dos.sub) %in% PopIDs[PopIDs[,2] == 2, 1],]

    dx = colSums(dos.sub.pop1, na.rm = T)
    nx = 2*nrow(dos.sub.pop1)
    dy = colSums(dos.sub.pop2, na.rm = T)
    ny = 2*nrow(dos.sub.pop2)

    #get the fraction of minor allele in pop1 per locus
    fracdeletpop1 = ((2*dx*(nx - dx))/(nx*(nx - 1)))
    #get fraction major allele in pop2 per locus
    fracnondeletpop2 = 1 - ((2*dy*(ny - dy))/(ny*(ny - 1)))
    #get Lxy: sum of multiplication of both
    L2xy = sum((fracdeletpop1 * fracnondeletpop2))

    #same for pop2

    #get the fraction of minor allele in pop1 per locus
    fracdeletpop1 = ((2*dy*(ny - dy))/(ny*(ny - 1)))
    #get fraction major allele in pop2 per locus
    fracnondeletpop2 = 1 - ((2*dx*(nx - dx))/(nx*(nx - 1)))
    #get Lxy: sum of multiplication of both
    L2yx = sum((fracdeletpop1 * fracnondeletpop2))

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
    write.table(output, "./6_MinorAllelesCount/stats/EXOME2_Meanvalues_R2xy_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

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
    write.table(output, "./6_MinorAllelesCount/stats/EXOME2_Jackknife_R2xy_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", quote = F, col.names = T, row.names = F)

EOF
