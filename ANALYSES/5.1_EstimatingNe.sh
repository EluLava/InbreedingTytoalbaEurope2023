#! /bin/bash

R --vanilla << EOF

    library(gaston)
    library(hierfstat)
    library(SNPRelate)

    #Read the bed matrix
    bed = read.bed.matrix("./1.2_BED/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90")

    #Update the famID column so that it contains populations, we'll use FHBD daatframe created before
    dta = read.table("./Analyses/GeneticPositionsFHBDs_RP502.txt", header = T)

    #read SuperScaffolds list
    SuperScaffolds = as.vector(read.table("./inputFILES/SuperScaffolds.list")[,1])

    #Update famid with the corresponding Pooulation value in dta
    bed@ped[,1] = dta[,6][match(bed@ped[,2], dta[,1])]

    #Add subsequent SNP IDs
    bed@snps[,2] = paste(bed@snps[,1], bed@snps[,4], sep = ".")

    #Create empty dataframe with ss start and stop positions of the geomic chunk we'll use for boostrapping
    chunksBOOT = as.data.frame(matrix(ncol = 6, nrow = 0))
    colnames(chunksBOOT) = c("SuperScaffold", "start", "stop", "SNP1", "SNP2", "nbSNPs")

    #Divide genome in 1MB regions
    #Loop through chromosomes
    for(ss in SuperScaffolds){

	print(paste0("Starting ss "), ss, " !!!!")

        #extract length of the ss
        len = max(bed@snps['pos'][bed@snps['chr'] == ss]) - min(bed@snps['pos'][bed@snps['chr'] == ss])

    	if(len >= 1e6) {

            #divide the ss into 100KB chunks
            startchunk = seq(1,(len-1e6),1e6)
            stopchunks = startchunk + 1e6

	    SNP1 = vector(mode = "character", length = 0)
            SNP2 = vector(mode = "character", length = 0)
	    nbSNPs = vector(mode = "character", length = 0)

            #loop through start
            for(win in 1:length(startchunk)){

		#get nb SNPs
		nbSNPs.win = length(bed@snps['id'][(bed@snps['chr'] == ss) & (bed@snps['pos'] > startchunk[win]) & (bed@snps['pos'] < stopchunks[win])])

		#check that we have at least one 2 SNPs
		if(nbSNPs.win > 1){

		    SNP1.win = bed@snps['id'][(bed@snps['chr'] == ss) & (bed@snps['pos'] > startchunk[win]) & (bed@snps['pos'] < stopchunks[win])][1]
		    SNP2.win = bed@snps['id'][(bed@snps['chr'] == ss) & (bed@snps['pos'] > startchunk[win]) & (bed@snps['pos'] < stopchunks[win])][length(bed@snps['id'][(bed@snps['chr'] == ss) & (bed@snps['pos'] > startchunk[win]) & (bed@snps['pos'] < stopchunks[win])])]

		} else {SNP1.win = NA; SNP2.win = NA}

	 	    SNP1 = append(SNP1, SNP1.win)
		    SNP2 = append(SNP2, SNP2.win)
		    nbSNPs = append(nbSNPs, nbSNPs.win)
	    }

            #create output df
            linetoapp = as.data.frame(matrix(nrow = length(startchunk), ncol = 6))
            colnames(linetoapp) = colnames(chunksBOOT)
            linetoapp[,1] = rep(ss, length(startchunk))
            linetoapp[,2] = startchunk
            linetoapp[,3] = stopchunks
            linetoapp[,4] = SNP1
            linetoapp[,5] = SNP2
            linetoapp[,6] = nbSNPs

            #add output
	    	chunksBOOT = rbind(chunksBOOT, linetoapp)

        }

        #add output
        #chunksBOOT = rbind(chunksBOOT, linetoapp)
    }

    #as numeric start ad stop positions
    chunksBOOT[,2] = as.numeric(chunksBOOT[,2])
    chunksBOOT[,3] = as.numeric(chunksBOOT[,3])

    #save tmp chunkboot
    #write.table(chunksBOOT, "./8_NeEst/ChunkBOOTswindows.txt", quote = F, col.names = T, row.names = F)
    #read tmp chunkboots
    #chunksBOOT = read.table("./8_NeEst/ChunkBOOTswindows.txt", h = T)

    #rm windows with NA as SNPs
    chunksBOOT = chunksBOOT[(!is.na(chunksBOOT['SNP1'])) & (!is.na(chunksBOOT['SNP2'])),]

    #Create output table with pop, pi and nb snps per window !
    PiTable = as.data.frame(matrix(nrow = (nrow(chunksBOOT)*length(unique(bed@ped[,1]))), ncol = (4)))
    colnames(PiTable) = c("window", "Pi", "nbsnps", "pop")

    #Fill window column
    PiTable[,1] = rep(c(1:nrow(chunksBOOT)), length(unique(bed@ped[,1])))
    #Fill population column
    PiTable[,4] = rep(unique(bed@ped[,1]), each = nrow(chunksBOOT))

    #Loop through populations
    for(pop in unique(bed@ped[,1])){

        print(paste0("Starting pop ", pop))
        #Subsample pop
        bedpop = select.inds(bed, famid == pop)

		#write a tmp bed matrix to extract dosage for pi and nb snps estimation
		#save the bed
		write.bed.matrix(x = bedpop, basename = paste0("/scratch/elavanc1/ROHsRP502/bedMatperNe_",pop))
		#from bed to gds
		snpgdsBED2GDS(bed.fn=paste0("/scratch/elavanc1/ROHsRP502/bedMatperNe_", pop, ".bed"), fam.fn = paste0("/scratch/elavanc1/ROHsRP502/bedMatperNe_", pop, ".fam"), bim.fn = paste0("/scratch/elavanc1/ROHsRP502/bedMatperNe_", pop, ".bim"),
				out.gdsfn = paste0("/scratch/elavanc1/ROHsRP502/gdsMatperNe_", pop, ".gds"), cvt.chr = "char")
	
		#bed2gds
		genos = snpgdsGetGeno(paste0('/scratch/elavanc1/ROHsRP502/gdsMatperNe_',pop,'.gds'))
	
		#set col and row names
		indvid = unlist(bedpop@ped['id'])
		names(indvid) = NULL
		rownames(genos) = indvid
		snpid = unlist(bedpop@snps['id'])
		names(snpid) = NULL
		colnames(genos) = snpid
	
		#create the TMP dataframe in case of timeout
		TMPpiTable = PiTable[PiTable['pop'] == pop,]
	
		#Loop through chunks
		for(chunk in 1:nrow(chunksBOOT)){
	
		    #get SNP1
		    SNP1 = chunksBOOT[chunk, 'SNP1']
		    #get SNP2
    	        SNP2 = chunksBOOT[chunk, 'SNP2']
	
		    #Get genos columns which corresponds to SNP1 and SNP2
		    col1 = which(colnames(genos) == SNP1)
		    col2 = which(colnames(genos) == SNP2)
	
		    #sanity check with nb of SNPs in chunksBOOT (plus one because oth SNP are included !)
		    if(((col2 - col1) + 1) != chunksBOOT[chunk, 'nbSNPs']){stop("not the same number of SNPs in chunksBOT and sub dosage matrix")}
	
		    #subsample the dosage matrix
		    genos.chunk = genos[,col1:col2]
	
		    #get nb of variants
		    nbvar = ncol(genos.chunk)
		    #get pi
		    pi = pi.dosage(genos.chunk, L = 1e6)
	
		    #Add this to the outputdataframe full
		    PiTable[,'Pi'][(PiTable['window'] == chunk) & (PiTable['pop'] == pop)] = pi
    	        PiTable[,'nbsnps'][(PiTable['window'] == chunk) & (PiTable['pop'] == pop)] = nbvar
	
		    #Add this to the tmp pop datframe in case of timeout
		    TMPpiTable[,'Pi'][(TMPpiTable['window'] == chunk) & (TMPpiTable['pop'] == pop)] = pi
		    TMPpiTable[,'nbsnps'][(TMPpiTable['window'] == chunk) & (TMPpiTable['pop'] == pop)] = nbvar
	
		}
	
		#in case of time out, save the popdta
		write.table(TMPpiTable, paste0("/scratch/elavanc1/ROHsRP502/PiTable_",pop,".txt"), quote = F, col.names = T, row.names = F)

    }

    #write output
    write.table(PiTable, "./8_NeEst/PiTableperwindow.txt", quote = F, col.names = T, row.names = F)

    #PiTable = read.table("./8_NeEst/PiTableperwindow.txt", h = T)

    #set mut rate
    mu = 4.6e-9

    #list of unique windows
    windows = unique(PiTable[,'window'])

    #create output
    NeEstTable = as.data.frame(matrix(nrow = 0, ncol = 3))
    colnames(NeEstTable) = c("Boot", "Pop", "NeEst")

    #loop through bootstraps
    for(boot in 1:1000){

	print(paste0("Starting bootstrap ", boot, "/1000 !!!"))

	#randomly sample windows with replacment
	samp = sample(windows, length(windows), replace = T)

	#Loop through populations
	for(pop in unique(PiTable[,'pop'])){

	    #sample pis and windows
	    pipop = PiTable[(PiTable['pop'] == pop),]
	    pipopboot = pipop[samp,]

	    #get sum pi nb markers
	    piboot = sum(pipopboot[,'nbsnps'] * pipopboot[,'Pi'])/sum(pipopboot[,'nbsnps'])

	    #get Ne
	    Neest = piboot/mu

	    #create output
	    lineout = as.data.frame(cbind(boot, pop, Neest))
	    colnames(lineout) = colnames(NeEstTable)

	    #Add to final output
	    NeEstTable = rbind(lineout, NeEstTable)
	}

    }

    NeEstTable[,'NeEst'] = as.numeric(NeEstTable[,'NeEst'])

    #save output
    write.table(NeEstTable, "./8_NeEst/NeEstperPop.txt", quote = F, col.names = T, row.names = F)

EOF
