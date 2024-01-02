#! /bin/bash

mkdir -p Analyses/

#Get nb of SNPs per Super-Scaffolds (for F weighted mean among CHRs)
for ss in $(cat ./inputFILES/SuperScaffolds.list); do echo '${ss}' $(bcftools view -H ./1.1_VCFs/RP502_Libnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}.vcf.gz | wc -l) >> ./1.1_VCFs/SNPsnb_per_SupScaff.list; done

#Open R
R --vanilla <<EOF

    library("RZooRoH")

    #read sample names
    sampleIDs = as.vector(read.table("./inputFILES/RP502_Newnames.list")[,1])

    #read super scaffolds names
    SuperScaffolds = as.vector(read.table("./inputFILES/SuperScaffolds.list")[,1])

    #### F ####

    #Create DF with just sample IDs that we'll complete with each F super-scaffold
    F_HBD_512genint = as.data.frame(sampleIDs)

    #create the list of dataframes that we'll fill later (for CUMULATIVE estimation of FHBD)
    Ftofill = list(as.data.frame(matrix(ncol = 14, nrow=502)),as.data.frame(matrix(ncol = 14, nrow=502)),as.data.frame(matrix(ncol = 14, nrow=502)),as.data.frame(matrix(ncol = 14, nrow=502)),as.data.frame(matrix(ncol = 14, nrow=502)), 
               as.data.frame(matrix(ncol = 14, nrow=502)),as.data.frame(matrix(ncol = 14, nrow=502)),as.data.frame(matrix(ncol = 14, nrow=502)),as.data.frame(matrix(ncol = 14, nrow=502)),as.data.frame(matrix(ncol = 14, nrow=502)),
               as.data.frame(matrix(ncol = 14, nrow=502)),as.data.frame(matrix(ncol = 14, nrow=502)),as.data.frame(matrix(ncol = 14, nrow=502)))

    #The different HBD classes we used
    classes = c(2,4,8,16,32,64,128,256,512,1024,2048,4096,8192)
    names(Ftofill) = classes

    #First column of FtoFill is INDVs
    for(i in 1:13){Ftofill[[i]][,1] = sampleIDs}

    #Read the SuperScaffolds R ENV files and extract F
    for(ss in SuperScaffolds){

        print(paste0("Starting SuperScaff ", ss))

        #Load the R file
        load(list.files(path="./1.2_HBDsegments", pattern=paste0("EntireRsession_GeneticPosplus10_Model13HBDclasses_ss_",ss,".RData"), full.names = T))

        #Estimate Fgen 512
	    Fgen512int = cumhbd(zres=loc_mod, T = 1024)
	    #merge Fgen512int with the full dataframe
	    F_HBD_512genint = as.data.frame(cbind(F_HBD_512genint, Fgen512int))
	    #Set new column name
	    colnames(F_HBD_512genint)[ncol(F_HBD_512genint)] = paste0("FHBD_int_",int)

        #Estimate CUMULATIVE FHBD per class
        #loop through classes
        for(class in 1:13){
            #merge Fclass with the full dataframe
            Ftofill[[class]][,(int + 1)] = as.vector(cumhbd(zres=loc_mod, T = classes[class]*2))
        }


    }

    #Read the file with nb of SNPs per SuperScaffold (for F weigthed mean)
    SNPsNB = as.vector(read.table("./1.1_VCFs/SNPsnb_per_SupScaff.list", h = F)[,2])

    #Create new vector with weighted means (per intervals)
    F_HBD512gen = as.data.frame(cbind(INDVs=F_HBD_512genint[,1], FHBD512gen=as.numeric(apply(F_HBD_512genint[,2:11],1 ,function(x){mean(x, weigths = SNPsNB)}))))

    colnames(F_HBD512gen) = c("INDVs","FHBD512genREC")

    #create empty new dataframe that we'll fill with mean among int
    F_HBD_CUM = as.data.frame(matrix(nrow = 502, ncol = 14))

    #first column is INDV names
    F_HBD_CUM[,1] = sampleIDs

    #get the mean per class
    for(i in 1:13){F_HBD_CUM[,(i + 1)] = as.numeric(apply(Ftofill[[i]][,2:11], 1, function(x){mean(x, weigths = SNPsNB)}))}

    #colnames new df
    colnames(F_HBD_CUM) = c("INDVs",paste0("FHBD_CUM_", classes/2, "gen"))

    #merge both F est
    F_HBD = merge(F_HBD512gen, F_HBD_CUM, by = "INDVs")

    #read sample descriptions INTERESTING columns (Pops, sex, etc.)
    sampleDESC = read.csv("./inputFILES/Metadata_3kOwls-metadata_highcov_cramStats_v1_17_04_23.csv")[,c(2:4,6)]
    colnames(sampleDESC) = c("INDVs","Libnames","Population","GeneticSex")

    #Merge F_HBD with sample DESC
    F_HBD_RP502 = merge(F_HBD, sampleDESC, by = "INDVs")

    #save the FHBD "final" table
    write.table(F_HBD_RP502, "./Analyses/GeneticPositionsFHBDs_RP502.txt", sep = "\t", quote = F, col.names = T, row.names = F)

EOF
