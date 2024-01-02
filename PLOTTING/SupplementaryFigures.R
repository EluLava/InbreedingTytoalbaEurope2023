##################################################################################################################################
#################################################### SUPPEMENTARY PLOTSsssss #####################################################
##################################################################################################################################

setwd("~/PhD/Barn_owl/ChapterIV/data/Analyses/")

library(vioplot)

################################################
##### CUMULATIVE FHBD ISLANDS vs CONTINENT #####
################################################

dtaFCUM = read.table("./GeneticPositions_CUMULATIVE_FHBDs.txt", header = T)
#Create new column for ISLAND VS CONTINENT
dtaFCUM = as.data.frame(cbind(dtaFCUM, Continent=vector(length = nrow(dtaFCUM))))
dtaFCUM$Continent[dtaFCUM$Population %in% c("CH","DK","FR","GE","GR","IS","IT","MA","PT","SB")] = "continent"
dtaFCUM$Continent[dtaFCUM$Population %in% c("AE","CO","CT","CY","EC","GB","IO","IR","WC")] = "island"
#Pass Continent into factor
dtaFCUM$Continent = factor(dtaFCUM$Continent, levels = c("continent","island"))

#change data frame form
dtaFCUM.2 = as.data.frame(matrix(nrow=0, ncol=6))
colnames(dtaFCUM.2) = c("INDVs", "Generation", "CumFHBD", "Population", "GeneticSex", "Continent")

#loop through dtaFCUM to fill new data
for(gen in c(1,2,4,8,16,32,64,128,256,512,1024,2048,4096)){
  
  #extract corresponding F
  newlines = dtaFCUM[names(dtaFCUM) %in% c("INDVs",paste0("FHBD_CUM_", gen, "gen"),"Population","GeneticSex","Continent")]
  
  #new dataframe newlines with gen and correct order
  newlines.2 = as.data.frame(cbind(newlines$INDVs, Generation = rep(gen, nrow(newlines)), CumFHBD = newlines[,1], newlines[,3:5]))
  
  #merge
  dtaFCUM.2 = rbind(dtaFCUM.2, newlines.2)
  
}

labelsspaceCUM = c("1g\n[50]","2g\n[25]","4g\n[12.5]","8g\n[6.25]","16g\n[3.13]","32g\n[1.56]","64g\n[0.78]",
                   "128g\n[0.39]","256g\n[0.20]","512g\n[0.10]","1024g\n[0.05]","2048g\n[0.02]","4096g\n[0.01]")

pdf(file = "./PLOTs/SMPLOTS/FIGS1_CUMULATIVEFHBD500.pdf", width = 15, height = 10)

par(mar=c(10,10,0,0))

vioplot(dtaFCUM.2$CumFHBD ~ dtaFCUM.2$Continent + dtaFCUM.2$Generation, xlab = NA, ylab = NA, yaxt = 'n', ylim = c(0,0.45),
        col = rep(c("#741558", "#48989e"),13), frame.plot=FALSE, xaxt = 'n')
#Add x axis
axis(1, at = seq(1.5,25.5, length = 13), labels = labelsspaceCUM, cex.axis = 1.5, padj = 1)
#Add x axis name
mtext(1,text = "# generations back to the coalescence event\n[expected mean HBD segments length in Mb]", at = 13, line = 9, cex = 2)
#Add y axis
axis(2, at = seq(0,0.4,0.1), labels = seq(0,0.4,0.1), cex.axis = 2, las = 1, hadj = 1.25)
#Add y axis name
mtext(2, text = expression(F[HBD]), at = 0.225, line = 8, cex = 2)

dev.off()



###############################
##### FHBD per POPULATION #####
###############################

pdf(file = "./PLOTs/SMPLOTS/FIGS2_FHBD512genperPOP.pdf", width = 15, height = 5)

par(mar=c(7.5,7.5,0,0))

vioplot(dtaF$FHBD512genREC ~ dtaF$Population, xlab = NA, ylab = NA, yaxt = 'n', ylim = c(0,0.4),
        col = c(rep("#741558",10),rep("#48989e",9)), frame.plot=FALSE, xaxt = 'n')
#Add x axis
axis(1, at = 1:19, labels = levels(dtaF$Population), cex.axis = 1.5, padj = 1.25)
#Add x axis name
mtext(1,text = "Population", at = 9.5, line = 5, cex = 2)
#Add y axis
axis(2, at = seq(0,0.4,0.1), labels = seq(0,0.4,0.1), cex.axis = 2, las = 1, hadj = 1.25)
#Add y axis name
mtext(2, text = expression(F[HBD]), at = 0.2, line = 5, cex = 2)

dev.off()


###############################################
##### CONTINENTAL REFUGIUM vs RECOLONIZED #####
###############################################

dtaF = read.table("./GeneticPositionsFHBDs_FHOM_FROH_FASs_FPED_FAsUNR_RP502.txt", header = T)
#rm duplictaed columns
dtaF = dtaF[,c(1:8,10:12,15:22)]
colnames(dtaF)[c(2,7,8)] = c("Libnames","Population","GeneticSex")
#Pass Population into factor
dtaF$Population = factor(dtaF$Population, levels = c("CH","DK","FR","IS","PT","IT","MA","GE","SB","GR","AE","IO","CO","CT","CY","EC","WC","GB","IR"))
#Create new column for ISLAND VS CONTINENT
dtaF = as.data.frame(cbind(dtaF, Continent=vector(length = nrow(dtaF))))
dtaF$Continent[dtaF$Population %in% c("CH","DK","FR","GE","GR","IS","IT","MA","PT","SB")] = "continent"
dtaF$Continent[dtaF$Population %in% c("AE","CO","CT","CY","EC","GB","IO","IR","WC")] = "island"
#Pass Continent into factor
dtaF$Continent = factor(dtaF$Continent, levels = c("continent","island"))
#Create new column for REFUGIUM VS RECOLONIZED
dtaF = as.data.frame(cbind(dtaF, Refugium=vector(length = nrow(dtaF))))
dtaF$Refugium[dtaF$Population %in% c("PT","IT","IS","MA","GR")] = "refugium"
dtaF$Refugium[dtaF$Population %in% c("CH","DK","FR","GE","SB")] = "recolonized"
#Pass Refugium into factor
dtaF$Refugium = factor(dtaF$Refugium, levels = c("refugium","recolonized"))

#### SROH vs NROH ####

#HBD seg distrbution
dtaDIST = read.table("./RP502_HBDsegments_Viterbi_RecMapIncluded.txt", header = T)
samples = as.vector(read.table("../RP502_Libnames.list")[,1])

#Add Pop and Continent column from dtaF
dtaDIST.2 = merge(dtaDIST, dtaF[,c(2,7,8,20,21)], by = "Libnames")

#Create empty dataframe
dtaSROHNROH = as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(dtaSROHNROH) = c("INDVs","SROH","NROH")
#Loop through individuals
for(indv in unique(dtaF$INDVs)){
  
  linea=as.data.frame(matrix(nrow = 1, ncol = 3))
  colnames(linea) = colnames(dtaSROHNROH)
  
  linea[1,1] = indv
  linea[1,2] = sum(dtaDIST.2$length[dtaDIST.2$INDVs == indv])
  linea[1,3] = nrow(dtaDIST.2[dtaDIST.2$INDVs == indv,])    
  
  dtaSROHNROH = rbind(dtaSROHNROH, linea)    
}

dtaSROHNROH$SROH = as.integer(dtaSROHNROH$SROH)
dtaSROHNROH$NROH = as.integer(dtaSROHNROH$NROH)
#Add pop and rest of info
dtaSROHNROH.2 = merge(dtaSROHNROH, dtaF[,c(1,2,7,8,20,21)], by = "INDVs")

#### without the less than 10 HBD classes

dtaDIST.4 = dtaDIST.2[!(dtaDIST.2$HBDclass %in% c(11,12,13)),]

#Create empty dataframe
dtaSROHNROH.5 = as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(dtaSROHNROH.5) = c("INDVs","SROH","NROH")
#Loop through individuals
for(indv in unique(dtaF$INDVs)){
  
  linea=as.data.frame(matrix(nrow = 1, ncol = 3))
  colnames(linea) = colnames(dtaSROHNROH.5)
  
  linea[1,1] = indv
  linea[1,2] = sum(dtaDIST.4$length[dtaDIST.4$INDVs == indv])
  linea[1,3] = nrow(dtaDIST.4[dtaDIST.4$INDVs == indv,])    
  
  dtaSROHNROH.5 = rbind(dtaSROHNROH.5, linea)    
}

dtaSROHNROH.5$SROH = as.integer(dtaSROHNROH.5$SROH)
dtaSROHNROH.5$NROH = as.integer(dtaSROHNROH.5$NROH)

#Add pop and rest of info
dtaSROHNROH.6 = merge(dtaSROHNROH.5, dtaF[,c(1,2,7,8,20,21)], by = "INDVs")


#### HBD DIST ####

#Create empty df to fill
meanSUMofLENgthperClass = as.data.frame(matrix(nrow = 13, ncol = 25))
colnames(meanSUMofLENgthperClass) = c("HBDclass","ALL","Continent","Islands", "Refugium","Recolonized",as.character(unique(dtaDIST.2$Population)))
#Fill the HDB class column
meanSUMofLENgthperClass$HBDclass = 1:13

## OVER ALL

# Loop through classes
for(class in 1:13){
  
  #Subset the dataset
  dtasub = dtaDIST.2[dtaDIST.2$HBDclass == class,]
  
  #Create empt vector thatb we'll fill with INDV sum
  lengthssum = vector(mode = "numeric", length = 0)
  
  #Loop through individuals to get sum of lengths
  for(indv in unique(dtaDIST.2$INDVs)){
    
    #Estimate the sum
    sumindv = sum(dtasub$length[dtasub$INDVs == indv])
    
    #Add this to the vector
    lengthssum = append(lengthssum, sumindv)
    
  }
  
  #Fill the output dataframe
  meanSUMofLENgthperClass$ALL[meanSUMofLENgthperClass$HBDclass == class] = mean(lengthssum)
  
}

## CONTINENT

# Loop through classes
for(class in 1:13){
  
  #Subset the dataset
  dtasub = dtaDIST.2[dtaDIST.2$HBDclass == class & dtaDIST.2$Continent == "continent",]
  
  #Create empt vector thatb we'll fill with INDV sum
  lengthssum = vector(mode = "numeric", length = 0)
  
  #Loop through individuals to get sum of lengths
  for(indv in unique(dtaDIST.2$INDVs[dtaDIST.2$Continent == "continent"])){
    
    #Estimate the sum
    sumindv = sum(dtasub$length[dtasub$INDVs == indv])
    
    #Add this to the vector
    lengthssum = append(lengthssum, sumindv)
    
  }
  
  #Fill the output dataframe
  meanSUMofLENgthperClass$Continent[meanSUMofLENgthperClass$HBDclass == class] = mean(lengthssum)
  
}

## ISLANDs

# Loop through classes
for(class in 1:13){
  
  #Subset the dataset
  dtasub = dtaDIST.2[dtaDIST.2$HBDclass == class & dtaDIST.2$Continent == "island",]
  
  #Create empt vector thatb we'll fill with INDV sum
  lengthssum = vector(mode = "numeric", length = 0)
  
  #Loop through individuals to get sum of lengths
  for(indv in unique(dtaDIST.2$INDVs[dtaDIST.2$Continent == "island"])){
    
    #Estimate the sum
    sumindv = sum(dtasub$length[dtasub$INDVs == indv])
    
    #Add this to the vector
    lengthssum = append(lengthssum, sumindv)
    
  }
  
  #Fill the output dataframe
  meanSUMofLENgthperClass$Islands[meanSUMofLENgthperClass$HBDclass == class] = mean(lengthssum)
  
}

## REFUGIUM

# Loop through classes
for(class in 1:13){
  
  #Subset the dataset
  dtasub = dtaDIST.2[(dtaDIST.2$HBDclass == class) & (dtaDIST.2$Refugium == "refugium") & (!is.na(dtaDIST.2$Refugium)),]
  
  #Create empt vector thatb we'll fill with INDV sum
  lengthssum = vector(mode = "numeric", length = 0)
  
  #Loop through individuals to get sum of lengths
  for(indv in unique(dtaDIST.2$INDVs[(dtaDIST.2$Refugium == "refugium") & (!is.na(dtaDIST.2$Refugium))])){
    
    #Estimate the sum
    sumindv = sum(dtasub$length[dtasub$INDVs == indv])
    
    #Add this to the vector
    lengthssum = append(lengthssum, sumindv)
    
  }
  
  #Fill the output dataframe
  meanSUMofLENgthperClass$Refugium[meanSUMofLENgthperClass$HBDclass == class] = mean(lengthssum)
  
}

## RECOLONIZED

# Loop through classes
for(class in 1:13){
  
  #Subset the dataset
  dtasub = dtaDIST.2[(dtaDIST.2$HBDclass == class) & (dtaDIST.2$Refugium == "recolonized") & (!is.na(dtaDIST.2$Refugium)),]
  
  #Create empt vector thatb we'll fill with INDV sum
  lengthssum = vector(mode = "numeric", length = 0)
  
  #Loop through individuals to get sum of lengths
  for(indv in unique(dtaDIST.2$INDVs[(dtaDIST.2$Refugium == "recolonized") & (!is.na(dtaDIST.2$Refugium))])){
    
    #Estimate the sum
    sumindv = sum(dtasub$length[dtasub$INDVs == indv])
    
    #Add this to the vector
    lengthssum = append(lengthssum, sumindv)
    
  }
  
  #Fill the output dataframe
  meanSUMofLENgthperClass$Recolonized[meanSUMofLENgthperClass$HBDclass == class] = mean(lengthssum)
  
}

## perPOP

#Loop through populations
for(pop in as.character(unique(dtaDIST.2$Population))){
  
  # Loop through classes
  for(class in 1:13){
    
    #Subset the dataset
    dtasub = dtaDIST.2[dtaDIST.2$HBDclass == class & dtaDIST.2$Population == pop,]
    
    #Create empt vector thatb we'll fill with INDV sum
    lengthssum = vector(mode = "numeric", length = 0)
    
    #Loop through individuals to get sum of lengths
    for(indv in unique(dtaDIST.2$INDVs[dtaDIST.2$Population == pop])){
      
      #Estimate the sum
      sumindv = sum(dtasub$length[dtasub$INDVs == indv])
      
      #Add this to the vector
      lengthssum = append(lengthssum, sumindv)
      
    }
    #Fill the output dataframe
    meanSUMofLENgthperClass[,colnames(meanSUMofLENgthperClass) == pop][meanSUMofLENgthperClass$HBDclass == class] = mean(lengthssum)
    
  }
}

#### per CONTINENT

#Create new dataframe
barplotSistCONT = as.data.frame(matrix(nrow = 26, ncol = 3))
colnames(barplotSistCONT) = c("HBDclass","Continent","meansumlength")
#Fill class column
barplotSistCONT$HBDclass = rep(1:13,2)
#Fill Continent column
barplotSistCONT$Continent = rep(c("continent","island"), each = 13)
#Fill meansumlength columns
barplotSistCONT$meansumlength = c(meanSUMofLENgthperClass$Continent, meanSUMofLENgthperClass$Islands)
#Only recent
barplotSistCONTRECENT = barplotSistCONT[barplotSistCONT$HBDclass %in% c(1:10),]

#### per REFUGIUM

#Create new dataframe
barplotSistREF = as.data.frame(matrix(nrow = 26, ncol = 3))
colnames(barplotSistREF) = c("HBDclass","Refugium","meansumlength")
#Fill class column
barplotSistREF$HBDclass = rep(1:13,2)
#Fill Continent column
barplotSistREF$Refugium = rep(c("refugium","recolonized"), each = 13)
#Fill meansumlength columns
barplotSistREF$meansumlength = c(meanSUMofLENgthperClass$Refugium, meanSUMofLENgthperClass$Recolonized)
#Only recent
barplotSistREFRECENT = barplotSistREF[barplotSistREF$HBDclass %in% c(1:10),]

#### per POPULATION

#Create new dataframe
barplotSistPOP = as.data.frame(matrix(nrow = 247, ncol = 3))
colnames(barplotSistPOP) = c("HBDclass","Population","meansumlength")
#Fill class column
barplotSistPOP$HBDclass = rep(1:13,19)
#Fill Continent column
barplotSistPOP$Population = rep(colnames(meanSUMofLENgthperClass)[5:23], each = 13)
#Fill meansumlength columns
barplotSistPOP$meansumlength = unlist(c(meanSUMofLENgthperClass[,5:23]))
#Only recent
barplotSistPOPRECENT = barplotSistPOP[barplotSistPOP$HBDclass %in% c(1:10),]


#### PLOT ####

pdf(file = "./PLOTs/FIG1RecMapgoddD_HBD500gen_Fhbd_FasFhbd_SrohNroh_HBDdist_REFvsRECOL.pdf", width = 21, height = 15)

layout(mat = matrix(c(rep(0,7),0,1,0,2,2,2,0,rep(0,7),0,3,0,4,0,5,0,rep(0,7)), 5, 7, byrow = T), widths = c(.2,1,.15,.58,.02,.87,0), heights = c(.075,1,.15,1,.25))

#### PANEL A: FHBD 500 GEN CONTINENT vs ISLANDS ####

vioplot(dtaF$FHBD512genREC ~ dtaF$Refugium, xlab = NA, ylab = NA, yaxt = 'n', ylim = c(0,0.4),
        col = c("#ffbfee", "#91faff"), frame.plot = FALSE, xaxt = 'n')
#Add x axis
axis(1, at = 1:2, labels = c("refugium\nn = 44", "recolonized\nn = 372"), cex.axis = 2, padj = 1, cex.axis = 2.5)
#Add x axis name
mtext(1, text = "Population type", at = 1.5, line = 9, cex = 2.5)
#Add y axis
axis(2, at = seq(0,0.4,0.1), labels = seq(0,0.4,0.1), cex.axis = 3, las = 1, hadj = 1.25)
#Add y axis name
mtext(2, text = expression(F[HBD]), at = 0.2, line = 9, cex = 2.5)
#Add panel letter
mtext(side = 3, text = "A", at = -0.135, line = 3, cex = 4)


#### PANEL B: FAS vs FHBD

plot(dtaF$FHBD512genREC[!is.na(dtaF$Refugium)] ~ dtaF$FASpopCH[!is.na(dtaF$Refugium)], ylab = NA, xlab = NA, pch = c(3,4,8,15:18,25,11,13)[dtaF$Population[dtaF$Population %in% levels(dtaF$Population)[1:10]]],
     col = c("#ffbfee", "#91faff")[dtaF$Refugium[!is.na(dtaF$Refugium)]], ylim = c(0,0.35), xlim = c(-0.2,0.5), xaxt = 'n', yaxt = 'n', frame.plot = F, cex = 1.5)
abline(0,1)
#x axis
axis(side = 1, at = seq(-20,40,10)/100, labels = seq(-20,40,10)/100, cex.axis = 2.5 , padj = 1.5)
mtext(side = 1, text = expression(F[AS]), cex = 2.5, line = 9, at = 0.1)
#y axis
axis(side = 2, at = seq(0,35,5)/100, labels = seq(0,35,5)/100, cex.axis = 2.5, hadj = 1.5, las = 1)
mtext(side = 2, text = expression(F[HBD]), cex = 2.5, line = 9)
#Add legend
legend(bty = "n", x = 0.45, y = 0.25, legend = levels(dtaF$Population)[1:10], pch = c(3,4,8,15:18,25,11,13),
       col = c(rep("#91faff",3),rep("#ffbfee",4),rep("#91faff",2),"#ffbfee"), cex = 1.75)
#Add panel letter
mtext(side = 3, text = "B", at = -0.35, line = 3, cex = 4)

#### PANEL C: SROH vs NROH ####

plot(dtaSROHNROH.6$NROH[!is.na(dtaSROHNROH.6$Refugium)] ~ dtaSROHNROH.6$SROH[!is.na(dtaSROHNROH.6$Refugium)], ylab = NA, xlab = NA, pch = c(3,4,8,15:18,25,11,13)[dtaSROHNROH.6$Population[dtaSROHNROH.6$Population %in% levels(dtaSROHNROH.6$Population)[1:10]]], col = c("#ffbfee", "#91faff")[dtaSROHNROH.6$Refugium[!is.na(dtaSROHNROH.6$Refugium)]],
     xlim = c(0, 5e8), xaxt = 'n', yaxt = 'n', frame.plot = F, cex = 2)
#x axis
axis(side = 1, at = seq(0,4e8,1e8), labels = c(expression("0"),expression("1 x 10"^"8"),expression("2 x 10"^"8"),expression("3 x 10"^"8"),expression("4 x 10"^"8")), cex.axis = 2.5, padj = 1.5)
mtext(side = 1, text = expression(S["HBD segments [bp]"]), cex = 2.5, line = 10.5, at = 2e8)
#y axis
axis(side = 2, at = seq(0,600,200), labels = seq(0,600,200), cex.axis = 2.5, hadj = 1.5, las = 1)
mtext(side = 2, text = expression(N["HBD segments"]), cex = 2.5, line = 9, at = 300)
#Add legend
legend(bty = "n", x = 4.5e8, y = 700, legend = levels(dtaF$Population)[1:10], pch = c(3,4,8,15:18,25,11,13),
       col = c(rep("#91faff",3),rep("#ffbfee",4),rep("#91faff",2),"#ffbfee"), cex = 1.75)
#Add panel letter
mtext(side = 3, text = "C", at = -1.5e8, line = 3, cex = 4)

#### PANEL D: HBD segments distributions ####

#in here is actually split into two plots
#split dataframe into twogroups of classes
barplotSistREFRECENTsmall = barplotSistREFRECENT[barplotSistREFRECENT$HBDclass <= 4,]
barplotSistREFRECENTlarge = barplotSistREFRECENT[barplotSistREFRECENT$HBDclass > 4,]
barplotSistREFRECENTsmall$Refugium = factor(barplotSistREFRECENTsmall$Refugium, levels = levels(dtaF$Refugium))
barplotSistREFRECENTlarge$Refugium = factor(barplotSistREFRECENTlarge$Refugium, levels = levels(dtaF$Refugium))

#Two sets of labels, once with both KB and MB for precise dating but not space, one with only Mb for more space
labelsprecise = c("1g\n[50Mb]","2g\n[25Mb]","4g\n[12.5Mb]","8g\n[6.25Mb]","16g\n[3.13Mb]","32g\n[1.56Mb]","64g\n[781.25Kb]",
                  "128g\n[390.63Kb]","256g\n[195.31Kb]","512g\n[97.66Kb]")
labelsspace = c("1g\n[50]","2g\n[25]","4g\n[12.5]","8g\n[6.25]","16g\n[3.13]","32g\n[1.56]","64g\n[0.78]",
                "128g\n[0.39]","256g\n[0.20]","512g\n[0.10]")

#First four classes
barplot(formula = (meansumlength)/1000000 ~ Refugium + HBDclass, data = barplotSistREFRECENTsmall, beside = T, axes = F, col = c("#ffbfee", "#91faff"),
        xlab = NA, names.arg = NA, ylab = NA)
#add x axis
axis(1, at = seq(2,11,3), labels = labelsspace[1:4], line = 1, padj = 1, cex.axis = 2.45)
#Add y axis
axis(2, cex.axis = 2.5, las = 1, hadj = 1.25)
#Add axis title
mtext(side = 2, text = "Mean sum of length [Mb]", line = 9, cex = 2)
#Add panel letter
mtext(side = 3, text = "D", at = -4.5, line = 3, cex = 4)

#Last 6 classes
barplot(formula = (meansumlength)/1000000 ~ Refugium + HBDclass, data = barplotSistREFRECENTlarge, beside = T, axes = F, col = c("#ffbfee", "#91faff"),
        xlab = NA, names.arg = NA, ylab = NA)
#add x axis
axis(1, at = seq(2,17,3), labels = labelsspace[5:10], line = 1, padj = 1, cex.axis = 2.45)
#add axis title
mtext(side = 1, text = "# generations back to the coalescence event\n[expected mean HBD segments length in Mb]",
      line = 15, cex = 2, at = 2)
#Add y axis
axis(2, cex.axis = 2.5, las = 1, hadj = 1.5)

dev.off()

#### WILCOXON TEST FHBD + mean REFUGIUM vs RECOLONIZED

wilcox.test(dtaF$FHBD512genREC[!is.na(dtaF$Refugium)] ~ dtaF$Refugium[!is.na(dtaF$Refugium)])
mean(dtaF$FHBD512genREC[dtaF$Refugium == "refugium" & !is.na(dtaF$Refugium)])
mean(dtaF$FHBD512genREC[dtaF$Refugium == "recolonized" & !is.na(dtaF$Refugium)])



###############################################
##### FHBD vs Fas trimmed for relatedness #####
###############################################

#We get list of UNRELATED

UNRELATED = read.table("../RPUNRELATED_LibNames_NewNames.list", header = T)

#then we subsample dtaF with only UNR
dtaFUNR = dtaF[dtaF$INDVs %in% UNRELATED$Newnames,]

pdf(file = "./PLOTs/SMPLOTS/FIGS2_REC_FHBD500genvsFas187UNRELATED.pdf", width = 20, height = 15)

par(mar=c(10,13,0,0))

plot(dtaFUNR$FHBD512genREC ~ dtaFUNR$FASpopUNR, ylab = NA, xlab = NA, pch = c(3,4,8,15:18,25,11,13,3,4,8,15:18,25,13)[dtaFUNR$Population], col = c("#741558","#48989e")[dtaFUNR$Continent],
     ylim = c(0,0.35), xlim = c(-0.2,0.5), xaxt = 'n', yaxt = 'n', frame.plot = F, cex = 1.5)
abline(0,1)
#x axis
axis(side = 1, at = seq(-20,40,10)/100, labels = seq(-20,40,10)/100, cex.axis = 2.5 , padj = 1.5)
mtext(side = 1, text = expression(F[AS]), cex = 3, line = 8, at = 0.1)
#y axis
axis(side = 2, at = seq(0,35,5)/100, labels = seq(0,35,5)/100, cex.axis = 2.5, hadj = 1.5, las = 1)
mtext(side = 2, text = expression(F[HBD]), cex = 3, line = 9)
#Add legend
legend(bty = "n", x = 0.45, y = 0.3, legend = levels(dtaFUNR$Population), pch = c(3,4,8,15:18,25,11,13,3,4,8,15:18,25,13),
       col = c(rep("#741558",10),rep("#48989e",9)), cex = 2.5)

dev.off()



#####################################################
########### FHBD VS FPED Swiss population ###########
#####################################################

pdf(file = "./PLOTs/SMPLOTS/FIGS3RecMap_FHBD500genvsFped.pdf", width = 15, height = 15)

par(mar=c(10,10,0,0))

plot(dtaF$FHBD500genREC[dtaF$Population == "CH"] ~ dtaF$FPED[dtaF$Population == "CH"], xlab = NA, ylab = NA, yaxt = 'n', ylim = c(0,0.4),
     xlim = c(0,0.4), col = "#741558", frame.plot=FALSE, xaxt = 'n', cex = 2, pch = 20)
#one to one line
abline(0,1,lwd = 2)
#Add x axis
axis(1, at = seq(0,0.4,0.1), labels = seq(0,0.4,0.1), cex.axis = 2, padj = 1.5)
#Add x axis name
mtext(1,text = expression(F[PED]), at = 0.2, line = 7, cex = 3)
#Add y axis
axis(2, at = seq(0,0.4,0.1), labels = seq(0,0.4,0.1), cex.axis = 2, las = 1, hadj = 1.5)
#Add y axis name
mtext(2, text = expression(F[HBD]), at = 0.2, line = 7, cex = 3)

dev.off()



##########################################################
########### HBD segments distributions per POP ###########
##########################################################


barplotSistPOPRECENT$Population = factor(barplotSistPOPRECENT$Population, levels = levels(dtaF$Population))

popcolors = c("#e10d47","#f6a98b","#f75802","#b6a4a3","#7a1816","#ad1baf","#bf7c82","#ced8de","#132428","#efefa1","#5397dc","#1e6682","#48465b","#d0f2fc","#b548f7","#a2f59a","#fbd20b","#3cd02c","#c9c3bd")

pdf(file = "./PLOTs/SMPLOTS/FIGS4_REC_HBDsegDistperPop.pdf", width = 35, height = 17.5)

par(mar=c(15,10,2,2))

#Two sets of labels, once with both KB and MB for precise dating but not space, one with only Mb for more space
labelsprecise = c("1g\n[50Mb]","2g\n[25Mb]","4g\n[12.5Mb]","8g\n[6.25Mb]","16g\n[3.13Mb]","32g\n[1.56Mb]","64g\n[781.25Kb]",
                  "128g\n[390.63Kb]","256g\n[195.31Kb]","512g\n[97.66Kb]")
labelsspace = c("1g\n[50]","2g\n[25]","4g\n[12.5]","8g\n[6.25]","16g\n[3.13]","32g\n[1.56]","64g\n[0.78]",
                "128g\n[0.39]","256g\n[0.20]","512g\n[0.01]")

#barplot(formula = (meansumlength)/1000000 ~ Population + HBDclass, data = barplotSistPOPRECENT, beside = T, axes = F, col = c(rep("#741558",10),rep("#48989e",9)),
#        xlab = NA, names.arg = NA, ylab = NA, ylim = c(0,100))
barplot(formula = (meansumlength)/1000000 ~ Population + HBDclass, data = barplotSistPOPRECENT, beside = T, axes = F, col = popcolors,
        xlab = NA, names.arg = NA, ylab = NA, ylim = c(0,120), xlim = c(0,215))

#add x axis
axis(1, at = (seq(10,200,20)), labels = labelsspace, line = 1, padj = 1, cex.axis = 2.45)
#add axis title
mtext(side = 1, text = "# generations back to the coalescence event\n[expected mean HBD segments length in Mb]",
      line = 13, cex = 3)
#Add y axis
axis(2, cex.axis = 2.5, las = 1, hadj = 1.5)
#Add axis title
mtext(side = 2, text = "Mean sum of length [Mb]", line = 8, cex = 3)
legend(x = 205, y = 110, legend = levels(barplotSistPOPRECENT$Population), fill = popcolors, bty = 'n', cex = 3)

dev.off()


#######################################################################################
##### COUNT MINOR & HOMOZYGOUS MINOR ALLELES / NB SEG SITES ISLANDS vs CONTINENTS #####
#######################################################################################

library(vioplot)

setwd("~/PhD/Barn_owl/ChapterIV/data/Analyses/")

#### PLOT FIG 3 ####

dtaNB = read.table("./CountMinorAllelesperCategoryperINDV_EXOME_MajorMinorUNR.txt", header = T)
dtaF = read.table("./GeneticPositionsFHBDs_FHOM_FROH_FASs_FPED_FAsUNR_RP502.txt", header = T)
#rm duplictaed columns
dtaF = dtaF[,c(1:8,10:12,15:22)]
colnames(dtaF)[c(2,7,8)] = c("Libnames","Population","GeneticSex")
#Pass Population into factor
dtaF$Population = factor(dtaF$Population, levels = c("CH","DK","FR","IS","PT","IT","MA","GE","SB","GR","AE","IO","CO","CT","CY","EC","WC","GB","IR"))
#Create new column for ISLAND VS CONTINENT
dtaF = as.data.frame(cbind(dtaF, Continent=vector(length = nrow(dtaF))))
dtaF$Continent[dtaF$Population %in% c("CH","DK","FR","GE","GR","IS","IT","MA","PT","SB")] = "continent"
dtaF$Continent[dtaF$Population %in% c("AE","CO","CT","CY","EC","GB","IO","IR","WC")] = "island"
#Pass Continent into factor
dtaF$Continent = factor(dtaF$Continent, levels = c("continent","island"))

#merge that with the actual F
dtaNB.2 = merge(dtaNB, dtaF[,c(1,2,5,7,8,20)], by = "INDVs")
dtaNB.2$Continent = factor(dtaNB.2$Continent, levels = c("continent","island"))

#read the nb SNPs per INDV
dtaNBsnps = read.table("./NBsnpsPerINDV.txt", h = T)
colnames(dtaNBsnps) = c("INDVs", "SNPsNB")
#merge with dtaNB
dtaNB.2 = merge(dtaNB.2, dtaNBsnps, by = "INDVs")

#HOMOZYG
dtaNBHomozyg = read.table("./CountHomozygMinorAllelesperCategoryperINDV_EXOME_MajorMinorUNR.txt", header = T)
#merge that with the actual F
dtaNBHomozyg.2 = merge(dtaNBHomozyg, dtaF[,c(1,2,5,7,8,20)], by = "INDVs")
dtaNBHomozyg.2$Continent = factor(dtaNBHomozyg.2$Continent, levels = c("continent","island"))

#merge with dtaNBsnps
dtaNBHomozyg.2 = merge(dtaNBHomozyg.2, dtaNBsnps, by = "INDVs")

#plot
pdf(file = "./PLOTs/SMPLOTS/FIGS7_EXOME_MinCountANDHomozygUNRSegSitesCORRECTED.pdf", width = 20, height = 30)

layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5),0,5,0,6,0,rep(0,5),0,7,0,8,0,rep(0,5)), 9, 5, byrow = T), widths = c(.35,1,.2,1,.3), heights = c(.3,1,.3,1,1,1,.3,1,.2))

#normal margins
par(mar = c(5,5,5,5))

#MINOR ALLELE COUNT
#neutral
vioplot(dtaNB.2$Neutral/dtaNB.2$SNPsNB ~ dtaNB.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = NA, padj = 1.75, cex.axis = 5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of minor alleles"),
      cex = 2.5, line = 15)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Neutral", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "A", cex = 6, line = 6, at = -.5)

#lowly delet
vioplot(dtaNB.2$LowDelet/dtaNB.2$SNPsNB ~ dtaNB.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = NA, padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Lowly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "B", cex = 6, line = 6, at = -.5)

#mildly delet
vioplot(dtaNB.2$ModerDelet/dtaNB.2$SNPsNB ~ dtaNB.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNB.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of minor alleles"),
      cex = 2.5, line = 15)
#add panel name
mtext(side = 3, text = "Mildly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "C", cex = 6, line = 6, at = -.5)

#highly delet
vioplot(dtaNB.2$HighDelet/dtaNB.2$SNPsNB ~ dtaNB.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNB.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Highly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "D", cex = 6, line = 6, at = -.5)

#HOMOZYG MINOR ALLELE COUNT
#neutral
vioplot(dtaNBHomozyg.2$Neutral/dtaNBHomozyg.2$SNPsNB ~ dtaNBHomozyg.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of homozygous\nminor alleles"),
      cex = 2.5, line = 15)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Neutral", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "E", cex = 6, line = 6, at = -.5)

#lowly delet
vioplot(dtaNBHomozyg.2$LowDelet/dtaNBHomozyg.2$SNPsNB ~ dtaNBHomozyg.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Lowly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "F", cex = 6, line = 6, at = -.5)

#mildly delet
vioplot(dtaNBHomozyg.2$ModerDelet/dtaNBHomozyg.2$SNPsNB ~ dtaNBHomozyg.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of homozygous\nminor alleles"),
      cex = 2.5, line = 15)
#add panel name
mtext(side = 3, text = "Mildly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "G", cex = 6, line = 6, at = -.5)

#highly delet
vioplot(dtaNBHomozyg.2$HighDelet/dtaNBHomozyg.2$SNPsNB ~ dtaNBHomozyg.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Highly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "H", cex = 6, line = 6, at = -.5)

dev.off()


#Wilcoxon tests
wilcox.test(dtaNB.2$Neutral/dtaNB.2$SNPsNB ~ dtaNB.2$Continent)
wilcox.test(dtaNB.2$LowDelet/dtaNB.2$SNPsNB ~ dtaNB.2$Continent)
wilcox.test(dtaNB.2$ModerDelet/dtaNB.2$SNPsNB ~ dtaNB.2$Continent)
wilcox.test(dtaNB.2$HighDelet/dtaNB.2$SNPsNB ~ dtaNB.2$Continent)

wilcox.test(dtaNBHomozyg.2$Neutral/dtaNBHomozyg.2$SNPsNB ~ dtaNBHomozyg.2$Continent)
wilcox.test(dtaNBHomozyg.2$LowDelet/dtaNBHomozyg.2$SNPsNB ~ dtaNBHomozyg.2$Continent)
wilcox.test(dtaNBHomozyg.2$ModerDelet/dtaNBHomozyg.2$SNPsNB ~ dtaNBHomozyg.2$Continent)
wilcox.test(dtaNBHomozyg.2$HighDelet/dtaNBHomozyg.2$SNPsNB ~ dtaNBHomozyg.2$Continent)




########################################################################################
##### COUNT MINOR & HOMOZYGOUS MINOR ALLELES / NB SEG SITES REFUGIA vs RECOLONIZED #####
########################################################################################

library(vioplot)

setwd("~/PhD/Barn_owl/ChapterIV/data/Analyses/")

#### PLOT FIG 3 ####

dtaNB = read.table("./CountMinorAllelesperCategoryperINDV_EXOME_MajorMinorUNR.txt", header = T)
dtaF = read.table("./GeneticPositionsFHBDs_FHOM_FROH_FASs_FPED_FAsUNR_RP502.txt", header = T)
#rm duplictaed columns
dtaF = dtaF[,c(1:8,10:12,15:22)]
colnames(dtaF)[c(2,7,8)] = c("Libnames","Population","GeneticSex")
#Pass Population into factor
dtaF$Population = factor(dtaF$Population, levels = c("CH","DK","FR","IS","PT","IT","MA","GE","SB","GR","AE","IO","CO","CT","CY","EC","WC","GB","IR"))
#Create new column for ISLAND VS CONTINENT
dtaF = as.data.frame(cbind(dtaF, Continent=vector(length = nrow(dtaF))))
dtaF$Continent[dtaF$Population %in% c("CH","DK","FR","GE","GR","IS","IT","MA","PT","SB")] = "continent"
dtaF$Continent[dtaF$Population %in% c("AE","CO","CT","CY","EC","GB","IO","IR","WC")] = "island"
#Pass Continent into factor
dtaF$Continent = factor(dtaF$Continent, levels = c("continent","island"))

#merge that with the actual F
dtaNB.2 = merge(dtaNB, dtaF[,c(1,2,5,7,8,20)], by = "INDVs")
dtaNB.2$Continent = factor(dtaNB.2$Continent, levels = c("continent","island"))

#Create new column, refigium vs recolonized
dtaNB.3 = as.data.frame(cbind(dtaNB.2, Refugium=vector(length = nrow(dtaNB.2))))
dtaNB.3$Refugium[dtaNB.3$Population %in% c("PT","IT","IS","MA","GR")] = "refugium"
dtaNB.3$Refugium[dtaNB.3$Population %in% c("CH","DK","FR","GE","SB")] = "recolonized"
#Pass Refugium into factor
dtaNB.3$Refugium = factor(dtaNB.3$Refugium, levels = c("refugium","recolonized"))

#read the nb SNPs per INDV
dtaNBsnps = read.table("./NBsnpsPerINDV.txt", h = T)
colnames(dtaNBsnps) = c("INDVs", "SNPsNB")
#merge with dtaNB
dtaNB.3 = merge(dtaNB.3, dtaNBsnps, by = "INDVs")

#HOMOZYG
dtaNBHomozyg = read.table("./CountHomozygMinorAllelesperCategoryperINDV_EXOME_MajorMinorUNR.txt", header = T)
#merge that with the actual F
dtaNBHomozyg.2 = merge(dtaNBHomozyg, dtaF[,c(1,2,5,7,8,20)], by = "INDVs")
#Create new column, refigium vs recolonized
dtaNBHomozyg.3 = as.data.frame(cbind(dtaNBHomozyg.2, Refugium=vector(length = nrow(dtaNBHomozyg.2))))
dtaNBHomozyg.3$Refugium[dtaNBHomozyg.3$Population %in% c("PT","IT","IS","MA","GR")] = "refugium"
dtaNBHomozyg.3$Refugium[dtaNBHomozyg.3$Population %in% c("CH","DK","FR","GE","SB")] = "recolonized"
#Pass Refugium into factor
dtaNBHomozyg.3$Refugium = factor(dtaNBHomozyg.3$Refugium, levels = c("refugium","recolonized"))
#merge with dtaNB
dtaNBHomozyg.3 = merge(dtaNBHomozyg.3, dtaNBsnps, by = "INDVs")


#plot
pdf(file = "./PLOTs/SMPLOTS/FIGS8_EXOME_MinCountANDHomozygUNRSegSitesCORRECTED_RefRecolo.pdf", width = 22, height = 30)

layout(mat = matrix(c(rep(0,5),0,1,0,2,0,rep(0,5),0,3,0,4,0,rep(0,5),0,5,0,6,0,rep(0,5),0,7,0,8,0,rep(0,5)), 9, 5, byrow = T), widths = c(.35,1,.2,1,.3), heights = c(.3,1,.3,1,1,1,.3,1,.2))

#normal margins
par(mar = c(5,5,5,5))

#MINOR ALLELE COUNT
#neutral
vioplot(dtaNB.3$Neutral/dtaNB.3$SNPsNB ~ dtaNB.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = NA, padj = 1.75, cex.axis = 5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of minor alleles"),
      cex = 2.5, line = 15)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Neutral", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "A", cex = 6, line = 6, at = 0)

#lowly delet
vioplot(dtaNB.3$LowDelet/dtaNB.3$SNPsNB ~ dtaNB.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = NA, padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Lowly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "B", cex = 6, line = 6, at = 0)

#mildly delet
vioplot(dtaNB.3$ModerDelet/dtaNB.3$SNPsNB ~ dtaNB.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNB.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of minor alleles"),
      cex = 2.5, line = 15)
#add panel name
mtext(side = 3, text = "Mildly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "C", cex = 6, line = 6, at = 0)

#highly delet
vioplot(dtaNB.3$HighDelet/dtaNB.3$SNPsNB ~ dtaNB.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNB.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Highly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "D", cex = 6, line = 6, at = 0)

#HOMOZYG MINOR ALLELE COUNT
#neutral
vioplot(dtaNBHomozyg.3$Neutral/dtaNBHomozyg.3$SNPsNB ~ dtaNBHomozyg.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of homozygous\nminor alleles"),
      cex = 2.5, line = 15)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Neutral", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "D", cex = 6, line = 6, at = 0)

#lowly delet
vioplot(dtaNBHomozyg.3$LowDelet/dtaNBHomozyg.3$SNPsNB ~ dtaNBHomozyg.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Lowly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "E", cex = 6, line = 6, at = 0)

#mildly delet
vioplot(dtaNBHomozyg.3$ModerDelet/dtaNBHomozyg.3$SNPsNB ~ dtaNBHomozyg.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of homozygous\nminor alleles"),
      cex = 2.5, line = 15)
#add panel name
mtext(side = 3, text = "Mildly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "E", cex = 6, line = 6, at = 0)

#highly delet
vioplot(dtaNBHomozyg.3$HighDelet/dtaNBHomozyg.3$SNPsNB ~ dtaNBHomozyg.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Highly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "F", cex = 6, line = 6, at = 0)

dev.off()


#Wilcoxon tests
wilcox.test(dtaNB.3$Neutral/dtaNB.3$SNPsNB ~ dtaNB.3$Refugium)
wilcox.test(dtaNB.3$LowDelet/dtaNB.3$SNPsNB ~ dtaNB.3$Refugium)
wilcox.test(dtaNB.3$ModerDelet/dtaNB.3$SNPsNB ~ dtaNB.3$Refugium)
wilcox.test(dtaNB.3$HighDelet/dtaNB.3$SNPsNB ~ dtaNB.3$Refugium)

wilcox.test(dtaNBHomozyg.3$Neutral/dtaNBHomozyg.3$SNPsNB ~ dtaNBHomozyg.3$Refugium)
wilcox.test(dtaNBHomozyg.3$LowDelet/dtaNBHomozyg.3$SNPsNB ~ dtaNBHomozyg.3$Refugium)
wilcox.test(dtaNBHomozyg.3$ModerDelet/dtaNBHomozyg.3$SNPsNB ~ dtaNBHomozyg.3$Refugium)
wilcox.test(dtaNBHomozyg.3$HighDelet/dtaNBHomozyg.3$SNPsNB ~ dtaNBHomozyg.3$Refugium)






