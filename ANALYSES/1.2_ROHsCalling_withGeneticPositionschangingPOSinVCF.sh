#! /bin/bash

mkdir -p 2.2_HBDsegments

#First step is to create a file with correct format for changing POS field in VCF
awk 'NR > 1 {print $1"\t"$2"\t"$3}' ../data/2.1_RecMaps/PhysicalandGeneticPositionsPLUSTEN_ss_${ss}.txt| bgzip > ../data/1.1_VCFs/TMPfromPHYStoGENpos_RP502_Libnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}.tab.gz

#index
tabix -p vcf ../data/1.1_VCFs/TMPfromPHYStoGENpos_RP502_Libnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}.tab.gz

#replace POS in VCF
bcftools annotate -a ../data/1.1_VCFs/TMPfromPHYStoGENpos_RP502_Libnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}.tab.gz -c CHROM,POS,~POS -O z -o ../data/1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}_GenPOSplus10.vcf.gz \
../data/1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}.vcf.gz

bcftools index ../data/1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}_GenPOSplus10.vcf.gz

#rm the TMP file
rm ../data/1.1_VCFs/TMPfromPHYStoGENpos_RP502_Libnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}.tab.gz*

#Convert from VCF to OXford format
plink --vcf ../data/1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}_GenPOSplus10.vcf.gz \
--recode oxford --allow-extra-chr --out ../data/1.1_VCFs/TMP_OX_RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}_GenPOSplus10

#Open R
	R --vanilla <<EOF

		library(RZooRoH)
		library(doParallel)
		library(foreach)
		cl <- 10
		registerDoParallel(cl)

		# read the OX file
		data_Rohs = zoodata("../data/1.1_VCFs/TMP_OX_RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90_ss_${ss}_GenPOSplus10.gen", samplefile = "../data/inputFILES/RP502_NewNames.list", zformat = "gp")

		#### RUN THE MODEL

		# Create Mdodel: 13 HBD and 1 non-HBD class.es
		Mod <- zoomodel(K=14, krates=c(2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,8192), err = 3.76e-5)

		# Run the model
		loc_mod <- zoorun(Mod, data_Rohs, localhbd = TRUE, nT = cl)

		# Extract table with all detected segments
		loc_table <- loc_mod@hbdseg

	        #We need to update the information about INDVs NAMES
        	loc_table[,1] = data_Rohs@sample_ids[loc_table[,1]]

	        #Loop through CHR to fill the columns
        	for(chr in 1:length(data_Rohs@chrnames)){
            	    #Fill the CHROM column
            	    loc_table[,2][loc_table[,2] == chr] = data_Rohs@chrnames[chr]
        	}

		#Save this file
		write.table(loc_table, "../data/2.2_HBDsegments/RP502_allHBDsegments_GeneticPosplus10_ss_${ss}.hom", quote = F, col.names = T, row.names = F)

		#Save the table with F per HBD (and non-HBD) class (actually we won't use this)
		#write.table(loc_mod@realized, "../data/2.2_HBDsegments/RP502_FHBD_GeneticPosplus10_ss_${ss}.hom.indv", quote = F, col.names = T, row.names = F)

		#Just for this time we'll save the entire SESSION
		save.image(file = "../data/2.2_HBDsegments/EntireRsession_GeneticPosplus10_Model13HBDclasses_ss_${ss}.RData")

EOF

#rm TMP OX file
rm ../data/1.1_VCFs/TMP_OX_RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90ss_${ss}.*
