#! /bin/bash

R --vanilla <<EOF

    library("RZooRoH")

    #read sample IDs
    sampleIDs = as.vector(read.table("./inputFILES/RP502_NewNames.list")[,1])
    #Read NB SNPs per interval
    SNPsNB = as.vector(read.table("./inputFILES/SNPsnb_per_SupScaff.list", h = F)[,2])

    #Create empty dataframe that we'll fill with superscaffolds, ncol = 16: nb classes (13 HBD and 1 non-HBD) + 2 (CHR + POS)
    HBDprob = as.data.frame(matrix(ncol = 16, nrow = 0))
    colnames(HBDprob) = c("CHROM","POS",paste0("class",seq(1,14)))

    #Read the SuperScaff R ENV files
    for(ss in SuperScaffolds){

        print(paste0("Starting SuperScaff ", ss))

        #Load the R file
        load(list.files(path="./1.2_HBDsegments", pattern=paste0("EntireRsession_GeneticPosplus10_Model13HBDclasses_ss_",ss,".RData"), full.names = T))

        #Create empty dataframe that we'll fill with each individuals
		HBDprobint = as.data.frame(matrix(ncol = 16, nrow = SNPsNB[ss]))
		colnames(HBDprobint) = colnames(HBDprob)
		#fill CHR column
		HBDprobclass[,1] = ss
		#fill POS column
		HBDprobclass[,2] = data_Rohs@bp

        #Loop through HBD classes
        for(class in 1:14){

	    #create empty dataframe that we'll fill with all INDIVIDUALS so we can average among individuals later on
	    HBDprobclass = as.data.frame(matrix(ncol = 502, nrow = SNPsNB[ss]))

	    #Loop throgh individuals to extract their prob for this specific class
	    for(indv in 1:502){

		#Extract values from model
		HBDprobclass[,indv] = loc_mod@hbdp[[indv]][class,]

	    }

	    #Get the mean per SNP and fill the superscaff df
	    HBDprobint[,class+2] = apply(HBDprobclass, 1, mean)

        }

        #merge the interval HBD prob with the overall HBD prob
        HBDprob = rbind(HBDprob, HBDprobint)

    }

    #Write the final df
    write.table(HBDprob, "./Analyses/HBD_prob_per_SNP_perHBDclass_RP502_RecMapIncluded.txt", sep = "\t", row.names = F, col.names = T, quote = F)

EOF

