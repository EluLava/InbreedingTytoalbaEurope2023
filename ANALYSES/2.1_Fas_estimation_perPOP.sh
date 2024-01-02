#! /bin/bash

cd ../data

R --vanilla << EOF

    library(gaston)
    library(hierfstat)
    library(JGTeach)

    #source for fas funnction
    source("../scripts/functions/functions.R")

    #read the bed matrix
    bed = readVCF("./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90.vcf.gz")

    #Update the famID column so that it contains populations, we'll use FHBD dataframe created before
    dta = read.table("./Analyses/GeneticPositionsFHBDs_RP502.txt", header = T)

    #Update famid with the corresponding Pooulation value in dta
    bed@ped[,1] = dta[,5][match(bed@ped[,2], dta[,1])]

    #Get GENERAL FAS
    FasALL = get.fasT(bed, nb.cores = 1)

    #Save that into a dataframe with new column for FASpop that we'll fill per pop
    dtaFas = as.data.frame(cbind(INDVs = bed@ped[,2], Population = bed@ped[,1], FASmetapop = FasALL, FASpop=vector(mode = "numeric", length = length(FasALL))))

    #Get Fas per population

    #Loop through populations and getFas (also for the swiss so we can compare)
    for(pop in unique(bed@ped[,1])){

	   #subset bedmatrix
	   bedpop = select.inds(bed, famid == pop)
	   #rm monomorphic SNPs
	   bedpop.snps = select.snps(bedpop, maf > 0)
    
	   #Get Fas
	   Faspop = get.fasT(bedpop.snps, nb.cores = 1)
    
	   #merge with sample IDs
	   dtaFaspop = as.data.frame(cbind(INDVs = bedpop.snps@ped[,2], FASpop = Faspop))
    
	   #Fill the df we created before
	   dtaFas[match(dtaFaspop[,1],dtaFas[,1]),4] = dtaFaspop[,2]
    }

    #FOR SWISS, we need to use another function to get MiiT adn MbT for both the entire swiss pop and the unrelated
    pop = "CH"

    #subset bedmatrix
    bedpop = select.inds(bed, famid == pop)
    #rm monomorphic SNPs
    bedpop.snps = select.snps(bedpop, maf > 0)

    #get MiiT and MbT for the entire population
    ALLCHmatchingvalues = get.Matchingvalues(bed = bedpop.snps, nb.cores=1)

    #read the list of unrelated individuals
    unrelated = read.table("./inputFILES/RPUNRELATED_LibNames_NewNames.list", header = T)

    #subset the bed matrix so that it only contains unrelated swiss individuals
    bedpopUNR = select.inds(bedpop.snps, id %in% unrelated[,2])
    #rm monomorphic SNPs
    bedpopUNR.snps = select.snps(bedpopUNR, maf > 0)

    #get MiiT and MbT for the UNRELATED swiss
    UNRCHmatchingvalues = get.Matchingvalues(bed = bedpopUNR.snps, nb.cores=1)

    #Actually calculate the ENTIRE swiss Fas by using MiiT from the entire population (ALLCHmatchingvalues) and MbT from the set of UNRELATED (UNRCHmatchingvalues)
    FaspopCH = (((ALLCHmatchingvalues[[1]])*2-1)-UNRCHmatchingvalues[[2]])/(1-UNRCHmatchingvalues[[2]])

    #Create a new column to the dataframe
    dta.2 = as.data.frame(cbind(dtaFas, FASpopCH = vector(mode = "numeric", length = nrow(dta))))

    #Fill the other pops --> same as FASpop
    dta.2[,5][dta.2[,2] != "CH"] = dta.2[,4][dta.2[,2] != "CH"]

    #Fill the swiss individuals
    dta.2[,5][match(names(FaspopCH), dta.2[,1])] = FaspopCH

    #Merge the df with other Fs df
    dta.3 = merge(dta, dta.2, by = "INDVs")

    #Write the df
    write.table(dta.3, "./Analyses/GeneticPositionsFHBDs_FASs_RP502.txt", quote = F, row.names = F, col.names = T)

    #trim bed matrix for UNRELATED individuals (all pops for figure S4)
    bedUNR = select.inds(bed, id %in% UNR[,2])

    #filter monomorphic SNPs out
    bedUNR.snps = select.snps(bedUNR, maf > 0)
    
    #Get Fas per population
    #Save that into a dataframe with new column for FASpop that we'll fill per pop
    dtaFasUNR = as.data.frame(cbind(INDVs=bed@ped[,2], Population=bed@ped[,1], FASpopUNR=vector(mode = "numeric", length = length(nrow(dta)))))
    dtaFasUNR[,3] = as.numeric(NA)
    
    #Loop through populations and getFas (also for the swiss so we can compare)
    for(pop in unique(bed@ped[,1])){
        
        #subset bedmatrix
        bedpop = select.inds(bedUNR.snps, famid == pop)
        
        #rm monomorphic SNPs
        bedpop.snps = select.snps(bedpop, maf > 0)
        
        #Get Fas
        Faspop = get.fasT(bedpop.snps, nb.cores = 1)
        
        #merge with sample IDs
        dtaFaspop = as.data.frame(cbind(INDVs=bedpop.snps@ped[,2], FASpopUNR=Faspop))
        dtaFaspop[,2] = as.numeric(dtaFaspop[,2])
        
        #Fill the df we created before
        dtaFasUNR[match(dtaFaspop[,1],dtaFasUNR[,1]),3] = dtaFaspop[,2]
    }

    #Merge the df with other Fs df
    dta.4 = merge(dta.3, dtaFasUNR[,c(1,3)], by = "INDVs")
    
    #Write the df
    write.table(dta.4, "./Analyses/GeneticPositionsFHBDs_FASs_FAsUNR_RP502.txt", quote = F, row.names = F, col.names = T)

EOF
