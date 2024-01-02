setwd("~/PhD/Barn_owl/ChapterIV/data/Analyses/")

library(vioplot)

##############################################################################################
########################################## FIGURE 1 ##########################################
##############################################################################################


#### PLOT FIG S1 MAP ####

library(mapdata)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(sf)
library(ggplot2)
library(ggmap)
library(tmap)
library(rworldxtra)
library(sp)
library(OpenStreetMap)

#############################

#Read data
metadata <- read.table('../metadata_with_lat_lon.tab',sep='\t',h=T)
#rm USA and SG indvs
metadata <- metadata[! metadata$Population %in% c('USA','SGP'),]

colnames(metadata)[2] = "INDVs"

#merge with cont vs pop info
metadata = merge(metadata, dtaF[,c(1,16)], by = "INDVs")

xy <- cbind.data.frame(metadata$lon, metadata$lat)
spdf <- SpatialPointsDataFrame(xy, metadata[,c(2,4,6,7,8,12,13,14,15)])
sf <- st_as_sf(spdf,'+proj=moll' )

#making map 
data('countriesHigh')
margin_box <- function(box, margin=0){
  box[[1]] <- box[[1]] - margin
  box[[2]] <- box[[2]] - margin
  box[[3]] <- box[[3]] + margin
  box[[4]] <- box[[4]] + margin
  return(box) 
}

box <- st_bbox(sf)
box2 <- margin_box(box, margin=4)

meta <- metadata[! is.na(metadata$lat) & ! is.na(metadata$lon),]
boxlc <- st_bbox(sf)
box2lc <- margin_box(boxlc, margin=4)

register_google(key='AIzaSyA0C1kC_AXDZsXTRvLQzim450yz4l8-994') # PUT KEY HERE 
tmp <-  get_googlemap(center = c(6.591467,46.551096 ),zoom=9,maptype='satellite')
ggmap(tmp)

names(box2) <- c('left','bottom','right','top')
tmp <-  get_map(box2,zoom=3,scale=4,maptype='satellite',source='google')
pdf('./PLOTS/FIG1_SamplesMap.pdf',height=5,width=5)
par(mar = c(0,0,0,0))
ggmap(tmp) + scale_x_continuous(limits=box2[c('left','right')])+
  scale_y_continuous(limits=box2[c('bottom','top')])+ 
  geom_point(data=metadata, aes(lon,lat),
             size=1.5, colour = c("white","black")[metadata$Continent], fill = c("#741558","#48989e")[metadata$Continent], shape = 21) + xlab('Longitude') + 
  ylab('Latitude')

dev.off()

##############################################################################################
########################################## FIGURE 2 ##########################################
##############################################################################################

#### FHBD ####

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



#### PLOT FIG 1 ISLANDS vs CONTINENT ####

pdf(file = "./PLOTs/FIG2_FHBD512gen_Fhbd_FasFhbd_SrohNroh_HBDdist.pdf", width = 21, height = 15)

layout(mat = matrix(c(rep(0,7),0,1,0,2,2,2,0,rep(0,7),0,3,0,4,0,5,0,rep(0,7)), 5, 7, byrow = T), widths = c(.2,1,.15,.58,.02,.87,0), heights = c(.075,1,.15,1,.25))

#### PANEL A: FHBD 500 GEN CONTINENT vs ISLANDS ####

vioplot(dtaF$FHBD512genREC ~ dtaF$Continent, xlab = NA, ylab = NA, yaxt = 'n', ylim = c(0,0.4),
        col = c("#741558","#48989e"), frame.plot = FALSE, xaxt = 'n')
#Add x axis
axis(1, at = 1:2, labels = c("continent\nn = 416", "island\nn = 86"), cex.axis = 2, padj = 1, cex.axis = 2.5)
#Add x axis name
mtext(1, text = "Population type", at = 1.5, line = 9, cex = 2.5)
#Add y axis
axis(2, at = seq(0,0.4,0.1), labels = seq(0,0.4,0.1), cex.axis = 3, las = 1, hadj = 1.25)
#Add y axis name
mtext(2, text = expression(F[HBD]), at = 0.2, line = 9, cex = 2.5)
#Add panel letter
mtext(side = 3, text = "A", at = -0.135, line = 3, cex = 4)


#### PANEL B: FAS vs FHBD

plot(dtaF$FHBD512genREC ~ dtaF$FASpopCH, ylab = NA, xlab = NA, pch = c(3,4,8,15:18,25,11,13,3,4,8,15:18,25,13)[dtaF$Population], col = c("#741558","#48989e")[dtaF$Continent],
     ylim = c(0,0.35), xlim = c(-0.2,0.5), xaxt = 'n', yaxt = 'n', frame.plot = F, cex = 1.5)
abline(0,1)
#x axis
axis(side = 1, at = seq(-20,40,10)/100, labels = seq(-20,40,10)/100, cex.axis = 2.5 , padj = 1.5)
mtext(side = 1, text = expression(F[AS]), cex = 2.5, line = 9, at = 0.1)
#y axis
axis(side = 2, at = seq(0,35,5)/100, labels = seq(0,35,5)/100, cex.axis = 2.5, hadj = 1.5, las = 1)
mtext(side = 2, text = expression(F[HBD]), cex = 2.5, line = 9)
#Add legend
legend(bty = "n", x = 0.45, y = 0.35, legend = levels(dtaF$Population), pch = c(3,4,8,15:18,25,11,13,3,4,8,15:18,25,13),
       col = c(rep("#741558",10),rep("#48989e",9)), cex = 1.75)
#Add panel letter
mtext(side = 3, text = "B", at = -0.35, line = 3, cex = 4)

#### PANEL C: SROH vs NROH ####

plot(dtaSROHNROH.6$NROH ~ dtaSROHNROH.6$SROH, ylab = NA, xlab = NA, pch = c(3,4,8,15:18,25,11,13,3,4,8,15:18,25,13)[dtaSROHNROH.6$Population], col = c("#741558","#48989e")[dtaSROHNROH.6$Continent],
     xlim = c(0, 5e8), xaxt = 'n', yaxt = 'n', frame.plot = F, cex = 2)
#x axis
axis(side = 1, at = seq(0,4e8,1e8), labels = c(expression("0"),expression("1 x 10"^"8"),expression("2 x 10"^"8"),expression("3 x 10"^"8"),expression("4 x 10"^"8")), cex.axis = 2.5, padj = 1.5)
mtext(side = 1, text = expression(S["HBD segments [bp]"]), cex = 2.5, line = 10.5, at = 2e8)
#y axis
axis(side = 2, at = seq(0,600,200), labels = seq(0,600,200), cex.axis = 2.5, hadj = 1.5, las = 1)
mtext(side = 2, text = expression(N["HBD segments"]), cex = 2.5, line = 9, at = 300)
#Add legend
legend(bty = "n", x = 4.5e8, y = 700, legend = levels(dtaF$Population), pch = c(3,4,8,15:18,25,11,13,3,4,8,15:18,25,13),
       col = c(rep("#741558",10),rep("#48989e",9)), cex = 1.75)
#Add panel letter
mtext(side = 3, text = "C", at = -1.5e8, line = 3, cex = 4)

#### PANEL D: HBD segments distributions ####

#in here is actually split into two plots
#split dataframe into twogroups of classes
barplotSistCONTRECENTsmall = barplotSistCONTRECENT[barplotSistCONTRECENT$HBDclass <= 4,]
barplotSistCONTRECENTlarge = barplotSistCONTRECENT[barplotSistCONTRECENT$HBDclass > 4,]

#Two sets of labels, once with both KB and MB for precise dating but not space, one with only Mb for more space
labelsprecise = c("1g\n[50Mb]","2g\n[25Mb]","4g\n[12.5Mb]","8g\n[6.25Mb]","16g\n[3.13Mb]","32g\n[1.56Mb]","64g\n[781.25Kb]",
                  "128g\n[390.63Kb]","256g\n[195.31Kb]","512g\n[97.66Kb]")
labelsspace = c("1g\n[50]","2g\n[25]","4g\n[12.5]","8g\n[6.25]","16g\n[3.13]","32g\n[1.56]","64g\n[0.78]",
                "128g\n[0.39]","256g\n[0.20]","512g\n[0.10]")

#First four classes
barplot(formula = (meansumlength)/1000000 ~ Continent + HBDclass, data = barplotSistCONTRECENTsmall, beside = T, axes = F, col = c("#741558","#48989e"),
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
barplot(formula = (meansumlength)/1000000 ~ Continent + HBDclass, data = barplotSistCONTRECENTlarge, beside = T, axes = F, col = c("#741558","#48989e"),
        xlab = NA, names.arg = NA, ylab = NA)
#add x axis
axis(1, at = seq(2,17,3), labels = labelsspace[5:10], line = 1, padj = 1, cex.axis = 2.45)
#add axis title
mtext(side = 1, text = "# generations back to the coalescence event\n[expected mean HBD segments length in Mb]",
      line = 15, cex = 2, at = 2)
#Add y axis
axis(2, cex.axis = 2.5, las = 1, hadj = 1.5)

dev.off()

#### WILCOXON TEST FHBD + mean CONTINENT vs ISLAND

wilcox.test(dtaF$FHBD512genREC ~ dtaF$Continent)
mean(dtaF$FHBD512genREC[dtaF$Continent == "island"])
mean(dtaF$FHBD512genREC[dtaF$Continent == "continent"])


##############################################################################################
########################################## FIGURE 3 ##########################################
##############################################################################################

setwd("~/PhD/Barn_owl/ChapterIV/data/Analyses/")

#Package for scanning windows
#library(devtools)
#install_github('tavareshugo/windowscanr')
library(windowscanr)
library(dplyr)
library(tidyverse)
library(viridis)
library(ggplotify)
library(ggplot2)

HBDprob = read.table("./HBDprob_512G_meanAmongINDVsRecMapIncluded.txt", header = T)

HBDprobperSS = aggregate(HBDprob$pIBD, by = list(HBDprob$CHROM), FUN = mean)
colnames(HBDprobperSS) = c("Super.Scaffold","HBDavProb")

## some scaffolds are really shitty (either too small or devoid of SNPs I guess so I'll just focus on LARGE SCaffolds (defined as > 10K SNPs))

#Empty vector with Nb SNPs per ss
SNPsperSS = vector(mode = "numeric", length = length(unique(HBDprob$CHROM)))
names(SNPsperSS) = unique(HBDprob$CHROM)

for (ss in unique(HBDprob$CHROM)) {
  SNPsperSS[names(SNPsperSS) == ss] = nrow(HBDprob[HBDprob$CHROM == ss,])
}

#How many ss with less than 10K SNPs
SSlarge = names(SNPsperSS[SNPsperSS > 10000])

running_roh_100KB_cM <- winScan(x = HBDprob,
                                groups = "CHROM",
                                position = "POScM50000",
                                values = "pIBD",
                                win_size = 100000, #100KB
                                win_step = 20000, #20KB
                                funs = c("mean"))

#Get rid of windows where we have no SNPs
running_roh_100KB_cM %>% filter(pIBD_n > 0) -> running_roh_100KB_cM_p

#RM SS with less than 10K SNPs
running_roh_100KB_cM_p_largeSS = running_roh_100KB_cM_p[running_roh_100KB_cM_p$CHROM %in% SSlarge, ]

#Pass SS to factor and order per SS size (for better plot)
#Read SS sizes
SSsizes = read.table("../SSlengths.txt")

SSordered = SSsizes$V1[order(SSsizes$V2, decreasing = T)]
SSorderedlarge = SSordered[SSordered %in% SSlarge]

#get SS per size
running_roh_100KB_cM_p_largeSS$CHROM = factor(running_roh_100KB_cM_p_largeSS$CHROM, levels = SSorderedlarge)

#plot

#Colors
fill_cols <- viridis(20, option = "A")
qn <- seq(0, 1, length.out=length(fill_cols))

p1 <- ggplot(running_roh_100KB_cM_p_largeSS, aes(x = win_start, y = 0.5, fill = pIBD_mean)) + 
  #geom_tile(color = "grey", size = 0) +
  geom_tile() +
  theme_minimal(base_size = 13, base_family = "Helvetica") + 
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), 
                     breaks = seq(0, 100000000, by = 10000000),
                     labels = as.character(seq(0, 100, 10))) +
  ylab("Chromosome") +
  scale_fill_gradientn("Probability of belonging to an HBD segment",
                       colors = rev(fill_cols), values = qn,
                       breaks = c(0.1,0.2,0.3,0.4),
                       labels = c(0.1,0.2,0.3,0.4)) +
  facet_grid(CHROM~., switch="both") +
  xlab("Position [Mb]") +
  theme(panel.spacing.y=unit(0.1, "lines"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(margin=margin(t=5), size = 20),
        axis.title.y = element_text(margin=margin(r=5), size = 20),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 15),
        axis.ticks.x = element_line(linewidth = 0.5),
        #axis.title.y = element_blank(),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.805,0.13),
        legend.direction = "horizontal",
        strip.text.y.left = element_text(size = 15, angle = 0),
        #  axis.line.x = element_line(size = 0.3)
  ) +
  guides(fill = guide_colourbar(title.position = "bottom" ,
                                barwidth = 10.95, barheight = 0.5))
p1

#### "LEGEND"

x <- running_roh_100KB_bp.2_p_largeSS$pIBD_mean
y <- density(x, n = 2^12)
p2 <- ggplot(data.frame(x = y$x, y = y$y), aes(x, y)) + 
  geom_line() + 
  theme_minimal(base_size = 13, base_family = "Helvetica") + 
  geom_segment(aes(xend = x, yend = 0, color = x)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0.1, 0.2, 0.3, 0.4), labels = c(0.1, 0.2, 0.3, 0.4)) +
  scale_color_gradientn("Probability of belonging to an HBD segment",
                        colors = rev(fill_cols), values = qn,
                        breaks = c(0.1,0.2, 0.3, 0.4)) +
  theme(legend.position = "none",
        plot.margin = margin(1, 1, 1, 1, unit = "cm"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.text.x = element_text(color = "black", size = 25, margin = margin(t=15)),
        axis.text.y = element_text(color = "black", size = 25, margin = margin(t=15)),
        axis.title.x = element_text(margin=margin(t=25), size = 30),
        axis.title.y = element_text(margin=margin(r=25), size = 30)) + 
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Density") +
  xlab("Probability of belonging\nto an HBD segment")

p2
ggsave("./PLOTS/FIG3_ROHsDesnity_HeatMapperSS_allPops_DENSITY.pdf", p2, width = 10, height = 10)

##### Probability to be HBD according to nb genes #####

#Read annotation table
dtaANN = read.csv("../Annotation/tyto_alba_annotationGENES.csv")

#we need to change SS assembly names to oiur SS/CHR names
repASS = read.table("../Annotation/GCF_018691265.1_T.alba_DEE_v4.0_assembly_report.txt")[,c(1,7)]
colnames(repASS) = c("Super.Scaffold","Annotation.Genomic.Range.Accession")

dtaANN.2 = merge(dtaANN, repASS, by = "Annotation.Genomic.Range.Accession")

#read ss lengths
SSlengths = read.table("../SSlengths.txt")
colnames(SSlengths) = c("Super.Scaffold","lengths.bp")

nbGenesperSS = as.data.frame(cbind(Super.Scaffold=names(table(dtaANN.2$Super.Scaffold)),NBgenese=table(dtaANN.2$Super.Scaffold)))

dtagenesperss = merge(nbGenesperSS, SSlengths, by = "Super.Scaffold")
dtagenesperss$NBgenese = as.integer(dtagenesperss$NBgenese)

#Merge HBD prob with nb genes per ss
#dtagenesperss.2 = merge(dtagenesperss, HBDprobperSS, by = "Super.Scaffold")
#dtagenesperss.3 = as.data.frame(cbind(dtagenesperss.2, NbSNPsper1Mb = ((dtagenesperss.2$NBgenese/dtagenesperss.2$lengths.bp)*1000000)))

colnames(running_roh_100KB_bp.2_p_largeSS)[1] = "Super.Scaffold"

#Same for nb genes per windows but we first need to get mean position of each gene (rather than start or stop)
dtaANN.3 = as.data.frame(cbind(dtaANN.2, MEANPOS = (dtaANN.2$Annotation.Genomic.Range.Start + ((dtaANN.2$Annotation.Genomic.Range.Stop - dtaANN.2$Annotation.Genomic.Range.Start)/2))))

#ordered (just in case)
dtaANN.3.ord = dtaANN.3[order(dtaANN.3$MEANPOS),]

dtaANN.4.ord = dtaANN.3.ord[!is.na(dtaANN.3.ord$Gene.Group.Identifier) & dtaANN.3.ord$Gene.Type == "PROTEIN_CODING",]

#sliding windows
running_nbgenes.2_100KB <- winScan(x = dtaANN.4.ord,
                                   groups = "Super.Scaffold",
                                   position = "MEANPOS",
                                   values = "NCBI.GeneID",
                                   win_size = 100000, #100KB
                                   win_step = 20000, #20KB
                                   funs = c("length"))

colnames(running_roh_100KB_bp.2_p_largeSS)[1] = "Super.Scaffold"

dta100kbhndnbg.2 = merge(running_roh_100KB_bp.2_p_largeSS, running_nbgenes.2_100KB, by = c("Super.Scaffold","win_mid","win_start","win_end"), all.x = T)
#dta100kbhndnbg.2$NCBIGeneID_length NAs simply means no genes were present for this window in the annotation --> goes to 0
dta100kbhndnbg.2$NCBI.GeneID_length[is.na(dta100kbhndnbg.2$NCBI.GeneID_length)] = 0

#The plotttttt

#cumulative dirtbution to get quantile value of each value in dta100kbhndnbg.2$HBDavProb_mean (for col)
my.ecdf <- ecdf(dta100kbhndnbg.2$pIBD_mean)(dta100kbhndnbg.2$pIBD_mean)

pdf("./PLOTS/FIG3_RecMap_ROHsDensity_NBgenes.pdf", width = 10, height = 10)

par(mar=c(10,10,.5,.2))
plot(x = dta100kbhndnbg.2$NCBI.GeneID_length, y = dta100kbhndnbg.2$pIBD_mean, pch = 20,
     xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NA, ylab = NA, xlim = c(0,16), ylim = c(0,.26), col = rev(fill_cols)[round(quantile(x = dta100kbhndnbg.2$pIBD_mean, probs = my.ecdf)*75)+1])
axis(1, at = seq(0, 15, 5), padj = 1.25, cex.axis = 2)
mtext(1, text = "# genes", at = 8, line = 7, cex = 3)
axis(2, at = seq(0, 0.25, 0.05), labels = seq(0, 0.25, 0.05), las = 1, hadj = 1.25, cex.axis = 2)
mtext(2,text = expression(p[HBD]), at = 0.125, line = 7, cex = 3)
mtext(3,text = "C", at = -4.25, line = -2.75, cex = 5)

dev.off()

##############################################################################################
########################################## FIGURE 4 ##########################################
##############################################################################################

#### CONTINENTS vs ISLANDS ####

library(vioplot)

setwd("~/PhD/Barn_owl/ChapterIV/data/Analyses/")

#### PLOT FIG 4 ####

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

#Read Rxy files OVERALL
Rxy = read.table("./EXOME2_Meanvalues_Rxy_islandvsCont_MajorMinorUNR.txt", h = T)
#Read Rxy file SYNONYMOUS ONLY
RxySYN = read.table("EXOME_Meanvalues_Rxy_INTERGENICSNPs_islandvsCont_MajorMinorUNR.txt", h = T)

#Estimate R'
Rxy = Rxy/as.numeric(RxySYN[1,1])

#Read Rxy files OVERALL for SE
Rxyjackknive = read.table("./EXOME2_Jackknife_Rxy_islandvsCont_MajorMinorUNR.txt", h = T)
#Read Rxy files SYNONYMOUS ONLY for SE
RxyjackkniveSYN = read.table("EXOME_Jackknife_Rxy_INTERGENICSNPs_islandvsCont_MajorMinorUNR.txt", h = T)

#estimate R'
Rxyjackknive[,2:5] = Rxyjackknive[,2:5]/RxyjackkniveSYN[,2]

#Estimate standard error with the Jackknive: Sqrt( (n-1)/n 2; (T; - Tava)**2)
SeRxyNeutral = sqrt(((nrow(Rxyjackknive) - 1)/nrow(Rxyjackknive)) * sum((Rxyjackknive$RxyMeanNeutral - Rxy$RxyMeanNeutral)^2))
SeRxyLow = sqrt(((nrow(Rxyjackknive) - 1)/nrow(Rxyjackknive)) * sum((Rxyjackknive$RxyMeanLow - Rxy$RxyMeanLow)^2))
SeRxyMod = sqrt(((nrow(Rxyjackknive) - 1)/nrow(Rxyjackknive)) * sum((Rxyjackknive$RxyMeanMod - Rxy$RxyMeanMod)^2))
SeRxyHigh = sqrt(((nrow(Rxyjackknive) - 1)/nrow(Rxyjackknive)) * sum((Rxyjackknive$RxyMeanHigh - Rxy$RxyMeanHigh)^2))

#Read R2xy files OVERALL
R2xy = read.table("./EXOME2_Meanvalues_R2xy_islandvsCont_MajorMinorUNR.txt", h = T)
#Read R2xy file SYNONYMOUS ONLY
R2xySYN = read.table("EXOME_Meanvalues_R2xy_INTERGENICSNPs_islandvsCont_MajorMinorUNR.txt", h = T)

#Estimate R2'
R2xy = R2xy/as.numeric(R2xySYN[1,1])

#Read R2xy files OVERALL for SE
R2xyjackknive = read.table("./EXOME2_Jackknife_R2xy_islandvsCont_MajorMinorUNR.txt", h = T)
#Read R2xy files SYNONYMOUS ONLY for SE
R2xyjackkniveSYN = read.table("EXOME_Jackknife_R2xy_INTERGENICSNPs_islandvsCont_MajorMinorUNR.txt", h = T)

#estimate R2'
R2xyjackknive[,2:5] = R2xyjackknive[,2:5]/R2xyjackkniveSYN[,2]

#Estimate standard error with the Jackknive: Sqrt( (n-1)/n 2; (T; - Tava)**2)
SeR2xyNeutral = sqrt(((nrow(R2xyjackknive) - 1)/nrow(R2xyjackknive)) * sum((R2xyjackknive$R2xyMeanNeutral - R2xy$R2xyMeanNeutral)^2))
SeR2xyLow = sqrt(((nrow(R2xyjackknive) - 1)/nrow(R2xyjackknive)) * sum((R2xyjackknive$R2xyMeanLow - R2xy$R2xyMeanLow)^2))
SeR2xyMod = sqrt(((nrow(R2xyjackknive) - 1)/nrow(R2xyjackknive)) * sum((R2xyjackknive$R2xyMeanMod - R2xy$R2xyMeanMod)^2))
SeR2xyHigh = sqrt(((nrow(R2xyjackknive) - 1)/nrow(R2xyjackknive)) * sum((R2xyjackknive$R2xyMeanHigh - R2xy$R2xyMeanHigh)^2))


#plot
pdf(file = "./PLOTs/FIG4_EXOME_MinCountANDHomozygUNR_RprimexyINTERGENIC.pdf", width = 38, height = 27)

layout(mat = matrix(c(rep(0,7),0,1,0,2,0,5,0,rep(0,7),0,3,0,4,0,5,0,rep(0,7),0,6,0,7,0,10,0,rep(0,7),0,8,0,9,0,10,0,rep(0,7)), 9, 7, byrow = T), widths = c(.35,1,.2,1,.3,1,.02), heights = c(.3,1,.3,1,1,1,.3,1,.2))

#normal margins
par(mar = c(5,5,5,5))

#MINOR ALLELE COUNT
#neutral
vioplot(dtaNB.2$Neutral ~ dtaNB.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = NA, padj = 1.75, cex.axis = 5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of minor alleles"),
      cex = 2.5, line = 20)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Neutral", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "A", cex = 6, line = 8, at = 0)

#lowly delet
vioplot(dtaNB.2$LowDelet ~ dtaNB.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = NA, padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Lowly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "B", cex = 6, line = 8, at = 0)

#mildly delet
vioplot(dtaNB.2$ModerDelet ~ dtaNB.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNB.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of minor alleles"),
      cex = 2.5, line = 20)
#add panel name
mtext(side = 3, text = "Mildly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "C", cex = 6, line = 8, at = 0)

#highly delet
vioplot(dtaNB.2$HighDelet ~ dtaNB.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNB.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Highly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "D", cex = 6, line = 8, at = 0)

#Add margins for smaller plots
par(mar = c(20,10,10,10))

#Rxy
plot(1 ~ Rxy$RxyMeanNeutral, xlim = c(0.8, 1.2), ylim = c(0,4), type = 'n', xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NA, ylab = NA)
#axes
axis(1, at = seq(0.8, 1.2, .1), padj = 1.75, cex.axis = 5, lwd = 2)
mtext(text = expression("Ratio continent/islands"), side = 1, line = 17, cex = 3)
#axis(2, at = 4:1, labels = c("Neutral", "Lowly\ndeleterious", "Moderately\ndeleterious", "Highly\ndeleterious"), las = 1, hadj = 1.25, cex.axis = 4, lwd = 0, tick = 3)
mtext(text = "Neutral", side = 2, at = 4, line = 1, cex = 3, las = 2)
mtext(text = "Lowly\ndeleterious", side = 2, at = 3, line = 1, cex = 3, las = 2)
mtext(text = "Moderately\ndeleterious", side = 2, at = 2, line = 1, cex = 3, las = 2)
mtext(text = "Highly\ndeleterious", side = 2, at = 1, line = 1, cex = 3, las = 2)
#Add the zero dashed line
abline(v = 1, col = "darkgrey", lty = 2, lwd = 8)
#add panel letter
mtext(side = 3, text = "E", cex = 6, line = 8, at = 0.6)
#Add the points +- se
#Neutral
points(x = Rxy$RxyMeanNeutral, y = 4, pch = 20, cex = 8)
segments(x0 = (Rxy$RxyMeanNeutral - SeRxyNeutral), x1 = (Rxy$RxyMeanNeutral + SeRxyNeutral), y0 = 4, y1 = 4, lwd = 8)
segments(x0 = (Rxy$RxyMeanNeutral - SeRxyNeutral), x1 = (Rxy$RxyMeanNeutral - SeRxyNeutral), y0 = 3.9, y1 = 4.1, lwd = 8)
segments(x0 = (Rxy$RxyMeanNeutral + SeRxyNeutral), x1 = (Rxy$RxyMeanNeutral + SeRxyNeutral), y0 = 3.9, y1 = 4.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((Rxy$RxyMeanNeutral - SeRxyNeutral) > 1){text("*", x = (Rxy$RxyMeanNeutral + SeRxyNeutral + 0.02), y = 4, cex = 4)}
if((Rxy$RxyMeanNeutral + SeRxyNeutral) < 1){text("*", x = (Rxy$RxyMeanNeutral - SeRxyNeutral - 0.02), y = 4, cex = 4)}
#Low
points(x = Rxy$RxyMeanLow, y = 3, pch = 20, cex = 8)
segments(x0 = (Rxy$RxyMeanLow - SeRxyLow), x1 = (Rxy$RxyMeanLow + SeRxyLow), y0 = 3, y1 = 3, lwd = 8)
segments(x0 = (Rxy$RxyMeanLow - SeRxyLow), x1 = (Rxy$RxyMeanLow - SeRxyLow), y0 = 2.9, y1 = 3.1, lwd = 8)
segments(x0 = (Rxy$RxyMeanLow + SeRxyLow), x1 = (Rxy$RxyMeanLow + SeRxyLow), y0 = 2.9, y1 = 3.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((Rxy$RxyMeanLow - SeRxyLow) > 1){text("*", x = (Rxy$RxyMeanLow + SeRxyLow + 0.02), y = 3, cex = 4)}
if((Rxy$RxyMeanLow + SeRxyLow) < 1){text("*", x = (Rxy$RxyMeanLow - SeRxyLow - 0.02), y = 3, cex = 4)}
#Moderate
points(x = Rxy$RxyMeanMod, y = 2, pch = 20, cex = 8)
segments(x0 = (Rxy$RxyMeanMod - SeRxyMod), x1 = (Rxy$RxyMeanMod + SeRxyMod), y0 = 2, y1 = 2, lwd = 8)
segments(x0 = (Rxy$RxyMeanMod - SeRxyMod), x1 = (Rxy$RxyMeanMod - SeRxyMod), y0 = 1.9, y1 = 2.1, lwd = 8)
segments(x0 = (Rxy$RxyMeanMod + SeRxyMod), x1 = (Rxy$RxyMeanMod + SeRxyMod), y0 = 1.9, y1 = 2.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((Rxy$RxyMeanMod - SeRxyMod) > 1){text("*", x = (Rxy$RxyMeanMod + SeRxyMod + 0.02), y = 2, cex = 4)}
if((Rxy$RxyMeanMod + SeRxyMod) < 1){text("*", x = (Rxy$RxyMeanMod - SeRxyMod - 0.02), y = 2, cex = 4)}
#High
points(x = Rxy$RxyMeanHigh, y = 1, pch = 20, cex = 8)
segments(x0 = (Rxy$RxyMeanHigh - SeRxyHigh), x1 = (Rxy$RxyMeanHigh + SeRxyHigh), y0 = 1, y1 = 1, lwd = 8)
segments(x0 = (Rxy$RxyMeanHigh - SeRxyHigh), x1 = (Rxy$RxyMeanHigh - SeRxyHigh), y0 = .9, y1 = 1.1, lwd = 8)
segments(x0 = (Rxy$RxyMeanHigh + SeRxyHigh), x1 = (Rxy$RxyMeanHigh + SeRxyHigh), y0 = .9, y1 = 1.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((Rxy$RxyMeanHigh - SeRxyHigh) > 1){text("*", x = (Rxy$RxyMeanHigh + SeRxyHigh + 0.02), y = 1, cex = 4)}
if((Rxy$RxyMeanHigh + SeRxyHigh) < 1){text("*", x = (Rxy$RxyMeanHigh - SeRxyHigh - 0.02), y = 1, cex = 4)}

#back to normal margins
par(mar = c(5,5,5,5))

#HOMOZYG MINOR ALLELE COUNT
#neutral
vioplot(dtaNBHomozyg.2$Neutral ~ dtaNBHomozyg.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of homozygous\nminor alleles"),
      cex = 2.5, line = 20)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Neutral", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "F", cex = 6, line = 8, at = 0)

#lowly delet
vioplot(dtaNBHomozyg.2$LowDelet ~ dtaNBHomozyg.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Lowly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "G", cex = 6, line = 8, at = 0)

#mildly delet
vioplot(dtaNBHomozyg.2$ModerDelet ~ dtaNBHomozyg.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of homozygous\nminor alleles"),
      cex = 2.5, line = 20)
#add panel name
mtext(side = 3, text = "Mildly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "H", cex = 6, line = 8, at = 0)

#highly delet
vioplot(dtaNBHomozyg.2$HighDelet ~ dtaNBHomozyg.2$Continent, col = c("#741558", "#48989e"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.2$Continent), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Highly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "I", cex = 6, line = 8, at = 0)

#Add margins for smaller plots
par(mar = c(20,10,10,10))

#R2xy
plot(1 ~ R2xy$R2xyMeanNeutral, xlim = c(0.8, 1.2), ylim = c(0,4), type = 'n', xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NA, ylab = NA)
#axes
axis(1, at = seq(0.8, 1.2, .1), padj = 1.75, cex.axis = 5, lwd = 2)
mtext(text = expression(Ratio^{2} ~ "continent/islands"), side = 1, line = 17, cex = 3)
#axis(2, at = 4:1, labels = c("Neutral", "Lowly\ndeleterious", "Moderately\ndeleterious", "Highly\ndeleterious"), las = 1, hadj = 1.25, cex.axis = 4, lwd = 0, tick = 3)
mtext(text = "Neutral", side = 2, at = 4, line = 1, cex = 3, las = 2)
mtext(text = "Lowly\ndeleterious", side = 2, at = 3, line = 1, cex = 3, las = 2)
mtext(text = "Moderately\ndeleterious", side = 2, at = 2, line = 1, cex = 3, las = 2)
mtext(text = "Highly\ndeleterious", side = 2, at = 1, line = 1, cex = 3, las = 2)
#Add the zero dashed line
abline(v = 1, col = "darkgrey", lty = 2, lwd = 8)
#add panel letter
mtext(side = 3, text = "J", cex = 6, line = 8, at = 0.6)
#Add the points +- se
#Neutral
points(x = R2xy$R2xyMeanNeutral, y = 4, pch = 20, cex = 8)
segments(x0 = (R2xy$R2xyMeanNeutral - SeR2xyNeutral), x1 = (R2xy$R2xyMeanNeutral + SeR2xyNeutral), y0 = 4, y1 = 4, lwd = 8)
segments(x0 = (R2xy$R2xyMeanNeutral - SeR2xyNeutral), x1 = (R2xy$R2xyMeanNeutral - SeR2xyNeutral), y0 = 3.9, y1 = 4.1, lwd = 8)
segments(x0 = (R2xy$R2xyMeanNeutral + SeR2xyNeutral), x1 = (R2xy$R2xyMeanNeutral + SeR2xyNeutral), y0 = 3.9, y1 = 4.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((R2xy$R2xyMeanNeutral - SeR2xyNeutral) > 1){text("*", x = (R2xy$R2xyMeanNeutral + SeR2xyNeutral + 0.02), y = 4, cex = 4)}
if((R2xy$R2xyMeanNeutral + SeR2xyNeutral) < 1){text("*", x = (R2xy$R2xyMeanNeutral - SeR2xyNeutral - 0.02), y = 4, cex = 4)}
#Low
points(x = R2xy$R2xyMeanLow, y = 3, pch = 20, cex = 8)
segments(x0 = (R2xy$R2xyMeanLow - SeR2xyLow), x1 = (R2xy$R2xyMeanLow + SeR2xyLow), y0 = 3, y1 = 3, lwd = 8)
segments(x0 = (R2xy$R2xyMeanLow - SeR2xyLow), x1 = (R2xy$R2xyMeanLow - SeR2xyLow), y0 = 2.9, y1 = 3.1, lwd = 8)
segments(x0 = (R2xy$R2xyMeanLow + SeR2xyLow), x1 = (R2xy$R2xyMeanLow + SeR2xyLow), y0 = 2.9, y1 = 3.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((R2xy$R2xyMeanLow - SeR2xyLow) > 1){text("*", x = (R2xy$R2xyMeanLow + SeR2xyLow + 0.02), y = 3, cex = 4)}
if((R2xy$R2xyMeanLow + SeR2xyLow) < 1){text("*", x = (R2xy$R2xyMeanLow - SeR2xyLow - 0.02), y = 3, cex = 4)}
#Moderate
points(x = R2xy$R2xyMeanMod, y = 2, pch = 20, cex = 8)
segments(x0 = (R2xy$R2xyMeanMod - SeR2xyMod), x1 = (R2xy$R2xyMeanMod + SeR2xyMod), y0 = 2, y1 = 2, lwd = 8)
segments(x0 = (R2xy$R2xyMeanMod - SeR2xyMod), x1 = (R2xy$R2xyMeanMod - SeR2xyMod), y0 = 1.9, y1 = 2.1, lwd = 8)
segments(x0 = (R2xy$R2xyMeanMod + SeR2xyMod), x1 = (R2xy$R2xyMeanMod + SeR2xyMod), y0 = 1.9, y1 = 2.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((R2xy$R2xyMeanMod - SeR2xyMod) > 1){text("*", x = (R2xy$R2xyMeanMod + SeR2xyMod + 0.02), y = 2, cex = 4)}
if((R2xy$R2xyMeanMod + SeR2xyMod) < 1){text("*", x = (R2xy$R2xyMeanMod - SeR2xyMod - 0.02), y = 2, cex = 4)}
#High
points(x = R2xy$R2xyMeanHigh, y = 1, pch = 20, cex = 8)
segments(x0 = (R2xy$R2xyMeanHigh - SeR2xyHigh), x1 = (R2xy$R2xyMeanHigh + SeR2xyHigh), y0 = 1, y1 = 1, lwd = 8)
segments(x0 = (R2xy$R2xyMeanHigh - SeR2xyHigh), x1 = (R2xy$R2xyMeanHigh - SeR2xyHigh), y0 = .9, y1 = 1.1, lwd = 8)
segments(x0 = (R2xy$R2xyMeanHigh + SeR2xyHigh), x1 = (R2xy$R2xyMeanHigh + SeR2xyHigh), y0 = .9, y1 = 1.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((R2xy$R2xyMeanHigh - SeR2xyHigh) > 1){text("*", x = (R2xy$R2xyMeanHigh + SeR2xyHigh + 0.02), y = 1, cex = 4)}
if((R2xy$R2xyMeanHigh + SeR2xyHigh) < 1){text("*", x = (R2xy$R2xyMeanHigh - SeR2xyHigh - 0.02), y = 1, cex = 4)}

dev.off()


#Wilcoxon tests SINGLE COPY
wilcox.test(dtaNB.2$Neutral ~ dtaNB.2$Continent)
wilcox.test(dtaNB.2$LowDelet ~ dtaNB.2$Continent)
wilcox.test(dtaNB.2$ModerDelet ~ dtaNB.2$Continent)
wilcox.test(dtaNB.2$HighDelet ~ dtaNB.2$Continent)

#Effect sizes SINGLE COPY
effectSizeNeut = rstatix::wilcox_effsize(data = dtaNB.2, formula = Neutral ~ Continent)
effectSizeLow = rstatix::wilcox_effsize(data = dtaNB.2, formula = LowDelet ~ Continent)
effectSizeMid = rstatix::wilcox_effsize(data = dtaNB.2, formula = ModerDelet ~ Continent)
effectSizeHig = rstatix::wilcox_effsize(data = dtaNB.2, formula = HighDelet ~ Continent)

#Wilcoxon tests HOMOZYGOUS
wilcox.test(dtaNBHomozyg.2$Neutral ~ dtaNBHomozyg.2$Continent)
wilcox.test(dtaNBHomozyg.2$LowDelet ~ dtaNBHomozyg.2$Continent)
wilcox.test(dtaNBHomozyg.2$ModerDelet ~ dtaNBHomozyg.2$Continent)
wilcox.test(dtaNBHomozyg.2$HighDelet ~ dtaNBHomozyg.2$Continent)

#Effect sizes HOMOZYGOUS
effectSizeNeut = rstatix::wilcox_effsize(data = dtaNBHomozyg.2, formula = Neutral ~ Continent)
effectSizeLow = rstatix::wilcox_effsize(data = dtaNBHomozyg.2, formula = LowDelet ~ Continent)
effectSizeMid = rstatix::wilcox_effsize(data = dtaNBHomozyg.2, formula = ModerDelet ~ Continent)
effectSizeHig = rstatix::wilcox_effsize(data = dtaNBHomozyg.2, formula = HighDelet ~ Continent)


##############################################################################################
########################################## FIGURE 5 ##########################################
##############################################################################################

#### REFUGIUM vs RECOLONIZED ####

library(vioplot)

setwd("~/PhD/Barn_owl/ChapterIV/data/Analyses/")

#### PLOT FIG 5 ####

library(vioplot)

setwd("~/PhD/Barn_owl/ChapterIV/data/Analyses/")

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

#Read Rxy file
Rxy = read.table("./EXOME2_Meanvalues_Rxy_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", h = T)
#Read Rxy file SYNONYMOUS ONLY
RxySYN = read.table("EXOME_Meanvalues_Rxy_INTERGENICSNPs_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", h = T)

#Estimate R'
Rxy = Rxy/as.numeric(RxySYN[1,1])

#Read Rxy files OVERALL for SE
Rxyjackknive = read.table("./EXOME2_Jackknife_Rxy_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", h = T)
#Read Rxy files SYNONYMOUS ONLY for SE
RxyjackkniveSYN = read.table("EXOME_Jackknife_Rxy_INTERGENICSNPs_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", h = T)

#estimate R'
Rxyjackknive[,2:5] = Rxyjackknive[,2:5]/RxyjackkniveSYN[,2]

#Estimate standard error with the Jackknive: Sqrt( (n-1)/n 2; (T; - Tava)**2)
SeRxyNeutral = sqrt(((nrow(Rxyjackknive) - 1)/nrow(Rxyjackknive)) * sum((Rxyjackknive$RxyMeanNeutral - Rxy$RxyMeanNeutral)^2))
SeRxyLow = sqrt(((nrow(Rxyjackknive) - 1)/nrow(Rxyjackknive)) * sum((Rxyjackknive$RxyMeanLow - Rxy$RxyMeanLow)^2))
SeRxyMod = sqrt(((nrow(Rxyjackknive) - 1)/nrow(Rxyjackknive)) * sum((Rxyjackknive$RxyMeanMod - Rxy$RxyMeanMod)^2))
SeRxyHigh = sqrt(((nrow(Rxyjackknive) - 1)/nrow(Rxyjackknive)) * sum((Rxyjackknive$RxyMeanHigh - Rxy$RxyMeanHigh)^2))

#Read R2xy files OVERALL
R2xy = read.table("./EXOME2_Meanvalues_R2xy_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", h = T)
#Read R2xy file SYNONYMOUS ONLY
R2xySYN = read.table("EXOME_Meanvalues_R2xy_INTERGENICSNPs_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", h = T)

#Estimate R2'
R2xy = R2xy/as.numeric(R2xySYN[1,1])

#Read R2xy files OVERALL for SE
R2xyjackknive = read.table("./EXOME2_Jackknife_R2xy_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", h = T)
#Read R2xy files SYNONYMOUS ONLY for SE
R2xyjackkniveSYN = read.table("EXOME_Jackknife_R2xy_INTERGENICSNPs_RefugiumvsnonRefugiumnoISLnoGBIR_MajorMinorUNR.txt", h = T)

#estimate R2'
R2xyjackknive[,2:5] = R2xyjackknive[,2:5]/R2xyjackkniveSYN[,2]

#Estimate standard error with the Jackknive: Sqrt( (n-1)/n 2; (T; - Tava)**2)
SeR2xyNeutral = sqrt(((nrow(R2xyjackknive) - 1)/nrow(R2xyjackknive)) * sum((R2xyjackknive$R2xyMeanNeutral - R2xy$R2xyMeanNeutral)^2))
SeR2xyLow = sqrt(((nrow(R2xyjackknive) - 1)/nrow(R2xyjackknive)) * sum((R2xyjackknive$R2xyMeanLow - R2xy$R2xyMeanLow)^2))
SeR2xyMod = sqrt(((nrow(R2xyjackknive) - 1)/nrow(R2xyjackknive)) * sum((R2xyjackknive$R2xyMeanMod - R2xy$R2xyMeanMod)^2))
SeR2xyHigh = sqrt(((nrow(R2xyjackknive) - 1)/nrow(R2xyjackknive)) * sum((R2xyjackknive$R2xyMeanHigh - R2xy$R2xyMeanHigh)^2))


#plot
pdf(file = "./PLOTs/FIG5_EXOME_MinCountANDHomozygUNR_REFnsnonREFnoISLnoIRGB_RprimexyINTERGENIC.pdf", width = 38, height = 27)

layout(mat = matrix(c(rep(0,7),0,1,0,2,0,5,0,rep(0,7),0,3,0,4,0,5,0,rep(0,7),0,6,0,7,0,10,0,rep(0,7),0,8,0,9,0,10,0,rep(0,7)), 9, 7, byrow = T), widths = c(.35,1,.2,1,.3,1,.02), heights = c(.3,1,.3,1,1,1,.3,1,.2))

#normal margins
par(mar = c(5,5,5,5))

#MINOR ALLELE COUNT
#neutral
vioplot(dtaNB.3$Neutral ~ dtaNB.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = NA, padj = 1.75, cex.axis = 5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of minor alleles"),
      cex = 2.5, line = 20)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Neutral", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "A", cex = 6, line = 8, at = 0)

#lowly delet
vioplot(dtaNB.3$LowDelet ~ dtaNB.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = NA, padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Lowly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "B", cex = 6, line = 8, at = 0)

#mildly delet
vioplot(dtaNB.3$ModerDelet ~ dtaNB.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNB.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of minor alleles"),
      cex = 2.5, line = 20)
#add panel name
mtext(side = 3, text = "Mildly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "C", cex = 6, line = 8, at = 0)

#highly delet
vioplot(dtaNB.3$HighDelet ~ dtaNB.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNB.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Highly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "D", cex = 6, line = 8, at = 0)

#Add margins for smaller plots
par(mar = c(20,10,10,10))

#Rxy
plot(1 ~ Rxy$RxyMeanNeutral, xlim = c(0.8, 1.2), ylim = c(0,4), type = 'n', xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NA, ylab = NA)
#axes
axis(1, at = seq(0.8, 1.2, .1), padj = 1.75, cex.axis = 5, lwd = 2)
mtext(text = expression("Ratio continent/islands"), side = 1, line = 17, cex = 3)
#axis(2, at = 4:1, labels = c("Neutral", "Lowly\ndeleterious", "Moderately\ndeleterious", "Highly\ndeleterious"), las = 1, hadj = 1.25, cex.axis = 4, lwd = 0, tick = 3)
mtext(text = "Neutral", side = 2, at = 4, line = 1, cex = 3, las = 2)
mtext(text = "Lowly\ndeleterious", side = 2, at = 3, line = 1, cex = 3, las = 2)
mtext(text = "Moderately\ndeleterious", side = 2, at = 2, line = 1, cex = 3, las = 2)
mtext(text = "Highly\ndeleterious", side = 2, at = 1, line = 1, cex = 3, las = 2)
#Add the zero dashed line
abline(v = 1, col = "darkgrey", lty = 2, lwd = 8)
#add panel letter
mtext(side = 3, text = "E", cex = 6, line = 8, at = 0.6)
#Add the points +- se
#Neutral
points(x = Rxy$RxyMeanNeutral, y = 4, pch = 20, cex = 8)
segments(x0 = (Rxy$RxyMeanNeutral - SeRxyNeutral), x1 = (Rxy$RxyMeanNeutral + SeRxyNeutral), y0 = 4, y1 = 4, lwd = 8)
segments(x0 = (Rxy$RxyMeanNeutral - SeRxyNeutral), x1 = (Rxy$RxyMeanNeutral - SeRxyNeutral), y0 = 3.9, y1 = 4.1, lwd = 8)
segments(x0 = (Rxy$RxyMeanNeutral + SeRxyNeutral), x1 = (Rxy$RxyMeanNeutral + SeRxyNeutral), y0 = 3.9, y1 = 4.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((Rxy$RxyMeanNeutral - SeRxyNeutral) > 1){text("*", x = (Rxy$RxyMeanNeutral + SeRxyNeutral + 0.02), y = 4, cex = 4)}
if((Rxy$RxyMeanNeutral + SeRxyNeutral) < 1){text("*", x = (Rxy$RxyMeanNeutral - SeRxyNeutral - 0.02), y = 4, cex = 4)}
#Low
points(x = Rxy$RxyMeanLow, y = 3, pch = 20, cex = 8)
segments(x0 = (Rxy$RxyMeanLow - SeRxyLow), x1 = (Rxy$RxyMeanLow + SeRxyLow), y0 = 3, y1 = 3, lwd = 8)
segments(x0 = (Rxy$RxyMeanLow - SeRxyLow), x1 = (Rxy$RxyMeanLow - SeRxyLow), y0 = 2.9, y1 = 3.1, lwd = 8)
segments(x0 = (Rxy$RxyMeanLow + SeRxyLow), x1 = (Rxy$RxyMeanLow + SeRxyLow), y0 = 2.9, y1 = 3.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((Rxy$RxyMeanLow - SeRxyLow) > 1){text("*", x = (Rxy$RxyMeanLow + SeRxyLow + 0.02), y = 3, cex = 4)}
if((Rxy$RxyMeanLow + SeRxyLow) < 1){text("*", x = (Rxy$RxyMeanLow - SeRxyLow - 0.02), y = 3, cex = 4)}
#Moderate
points(x = Rxy$RxyMeanMod, y = 2, pch = 20, cex = 8)
segments(x0 = (Rxy$RxyMeanMod - SeRxyMod), x1 = (Rxy$RxyMeanMod + SeRxyMod), y0 = 2, y1 = 2, lwd = 8)
segments(x0 = (Rxy$RxyMeanMod - SeRxyMod), x1 = (Rxy$RxyMeanMod - SeRxyMod), y0 = 1.9, y1 = 2.1, lwd = 8)
segments(x0 = (Rxy$RxyMeanMod + SeRxyMod), x1 = (Rxy$RxyMeanMod + SeRxyMod), y0 = 1.9, y1 = 2.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((Rxy$RxyMeanMod - SeRxyMod) > 1){text("*", x = (Rxy$RxyMeanMod + SeRxyMod + 0.02), y = 2, cex = 4)}
if((Rxy$RxyMeanMod + SeRxyMod) < 1){text("*", x = (Rxy$RxyMeanMod - SeRxyMod - 0.02), y = 2, cex = 4)}
#High
points(x = Rxy$RxyMeanHigh, y = 1, pch = 20, cex = 8)
segments(x0 = (Rxy$RxyMeanHigh - SeRxyHigh), x1 = (Rxy$RxyMeanHigh + SeRxyHigh), y0 = 1, y1 = 1, lwd = 8)
segments(x0 = (Rxy$RxyMeanHigh - SeRxyHigh), x1 = (Rxy$RxyMeanHigh - SeRxyHigh), y0 = .9, y1 = 1.1, lwd = 8)
segments(x0 = (Rxy$RxyMeanHigh + SeRxyHigh), x1 = (Rxy$RxyMeanHigh + SeRxyHigh), y0 = .9, y1 = 1.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((Rxy$RxyMeanHigh - SeRxyHigh) > 1){text("*", x = (Rxy$RxyMeanHigh + SeRxyHigh + 0.02), y = 1, cex = 4)}
if((Rxy$RxyMeanHigh + SeRxyHigh) < 1){text("*", x = (Rxy$RxyMeanHigh - SeRxyHigh - 0.02), y = 1, cex = 4)}

#back to normal margins
par(mar = c(5,5,5,5))

#HOMOZYG MINOR ALLELE COUNT
#neutral
vioplot(dtaNBHomozyg.3$Neutral ~ dtaNBHomozyg.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of homozygous\nminor alleles"),
      cex = 2.5, line = 20)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Neutral", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "F", cex = 6, line = 8, at = 0)

#lowly delet
vioplot(dtaNBHomozyg.3$LowDelet ~ dtaNBHomozyg.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
#axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Lowly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "G", cex = 6, line = 8, at = 0)

#mildly delet
vioplot(dtaNBHomozyg.3$ModerDelet ~ dtaNBHomozyg.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#Add y axis title
mtext(side = 2, text =
        expression("Count of homozygous\nminor alleles"),
      cex = 2.5, line = 20)
#add panel name
mtext(side = 3, text = "Mildly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "H", cex = 6, line = 8, at = 0)

#highly delet
vioplot(dtaNBHomozyg.3$HighDelet ~ dtaNBHomozyg.3$Refugium, col = c("#ffbfee", "#91faff"), xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NULL, ylab = NULL)
axis(1, at = seq(1,2), labels = levels(dtaNBHomozyg.3$Refugium), padj = 1.75, cex.axis = 5, lwd = 2)
axis(2, las = 1, hadj = 1.25, cex.axis = 3.5, lwd = 2)
#add panel name
mtext(side = 3, text = "Highly deleterious", cex = 3, line = 8)
#add panel letter
mtext(side = 3, text = "I", cex = 6, line = 8, at = 0)

#Add margins for smaller plots
par(mar = c(20,10,10,10))

#R2xy
plot(1 ~ R2xy$R2xyMeanNeutral, xlim = c(0.8, 1.2), ylim = c(0,4), type = 'n', xaxt = 'n', yaxt = 'n', frame.plot = F, xlab = NA, ylab = NA)
#axes
axis(1, at = seq(0.8, 1.2, .1), padj = 1.75, cex.axis = 5, lwd = 2)
mtext(text = expression(Ratio^{2} ~ "continent/islands"), side = 1, line = 17, cex = 3)
#axis(2, at = 4:1, labels = c("Neutral", "Lowly\ndeleterious", "Moderately\ndeleterious", "Highly\ndeleterious"), las = 1, hadj = 1.25, cex.axis = 4, lwd = 0, tick = 3)
mtext(text = "Neutral", side = 2, at = 4, line = 1, cex = 3, las = 2)
mtext(text = "Lowly\ndeleterious", side = 2, at = 3, line = 1, cex = 3, las = 2)
mtext(text = "Moderately\ndeleterious", side = 2, at = 2, line = 1, cex = 3, las = 2)
mtext(text = "Highly\ndeleterious", side = 2, at = 1, line = 1, cex = 3, las = 2)
#Add the zero dashed line
abline(v = 1, col = "darkgrey", lty = 2, lwd = 8)
#add panel letter
mtext(side = 3, text = "J", cex = 6, line = 8, at = 0.6)
#Add the points +- se
#Neutral
points(x = R2xy$R2xyMeanNeutral, y = 4, pch = 20, cex = 8)
segments(x0 = (R2xy$R2xyMeanNeutral - SeR2xyNeutral), x1 = (R2xy$R2xyMeanNeutral + SeR2xyNeutral), y0 = 4, y1 = 4, lwd = 8)
segments(x0 = (R2xy$R2xyMeanNeutral - SeR2xyNeutral), x1 = (R2xy$R2xyMeanNeutral - SeR2xyNeutral), y0 = 3.9, y1 = 4.1, lwd = 8)
segments(x0 = (R2xy$R2xyMeanNeutral + SeR2xyNeutral), x1 = (R2xy$R2xyMeanNeutral + SeR2xyNeutral), y0 = 3.9, y1 = 4.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((R2xy$R2xyMeanNeutral - SeR2xyNeutral) > 1){text("*", x = (R2xy$R2xyMeanNeutral + SeR2xyNeutral + 0.02), y = 4, cex = 4)}
if((R2xy$R2xyMeanNeutral + SeR2xyNeutral) < 1){text("*", x = (R2xy$R2xyMeanNeutral - SeR2xyNeutral - 0.02), y = 4, cex = 4)}
#Low
points(x = R2xy$R2xyMeanLow, y = 3, pch = 20, cex = 8)
segments(x0 = (R2xy$R2xyMeanLow - SeR2xyLow), x1 = (R2xy$R2xyMeanLow + SeR2xyLow), y0 = 3, y1 = 3, lwd = 8)
segments(x0 = (R2xy$R2xyMeanLow - SeR2xyLow), x1 = (R2xy$R2xyMeanLow - SeR2xyLow), y0 = 2.9, y1 = 3.1, lwd = 8)
segments(x0 = (R2xy$R2xyMeanLow + SeR2xyLow), x1 = (R2xy$R2xyMeanLow + SeR2xyLow), y0 = 2.9, y1 = 3.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((R2xy$R2xyMeanLow - SeR2xyLow) > 1){text("*", x = (R2xy$R2xyMeanLow + SeR2xyLow + 0.02), y = 3, cex = 4)}
if((R2xy$R2xyMeanLow + SeR2xyLow) < 1){text("*", x = (R2xy$R2xyMeanLow - SeR2xyLow - 0.02), y = 3, cex = 4)}
#Moderate
points(x = R2xy$R2xyMeanMod, y = 2, pch = 20, cex = 8)
segments(x0 = (R2xy$R2xyMeanMod - SeR2xyMod), x1 = (R2xy$R2xyMeanMod + SeR2xyMod), y0 = 2, y1 = 2, lwd = 8)
segments(x0 = (R2xy$R2xyMeanMod - SeR2xyMod), x1 = (R2xy$R2xyMeanMod - SeR2xyMod), y0 = 1.9, y1 = 2.1, lwd = 8)
segments(x0 = (R2xy$R2xyMeanMod + SeR2xyMod), x1 = (R2xy$R2xyMeanMod + SeR2xyMod), y0 = 1.9, y1 = 2.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((R2xy$R2xyMeanMod - SeR2xyMod) > 1){text("*", x = (R2xy$R2xyMeanMod + SeR2xyMod + 0.02), y = 2, cex = 4)}
if((R2xy$R2xyMeanMod + SeR2xyMod) < 1){text("*", x = (R2xy$R2xyMeanMod - SeR2xyMod - 0.02), y = 2, cex = 4)}
#High
points(x = R2xy$R2xyMeanHigh, y = 1, pch = 20, cex = 8)
segments(x0 = (R2xy$R2xyMeanHigh - SeR2xyHigh), x1 = (R2xy$R2xyMeanHigh + SeR2xyHigh), y0 = 1, y1 = 1, lwd = 8)
segments(x0 = (R2xy$R2xyMeanHigh - SeR2xyHigh), x1 = (R2xy$R2xyMeanHigh - SeR2xyHigh), y0 = .9, y1 = 1.1, lwd = 8)
segments(x0 = (R2xy$R2xyMeanHigh + SeR2xyHigh), x1 = (R2xy$R2xyMeanHigh + SeR2xyHigh), y0 = .9, y1 = 1.1, lwd = 8)
#if SE does not overlap 1, add a star !
if((R2xy$R2xyMeanHigh - SeR2xyHigh) > 1){text("*", x = (R2xy$R2xyMeanHigh + SeR2xyHigh + 0.02), y = 1, cex = 4)}
if((R2xy$R2xyMeanHigh + SeR2xyHigh) < 1){text("*", x = (R2xy$R2xyMeanHigh - SeR2xyHigh - 0.02), y = 1, cex = 4)}

dev.off()


#Wilcoxon tests SINGLE COPY
wilcox.test(dtaNB.3$Neutral ~ dtaNB.3$Refugium)
wilcox.test(dtaNB.3$LowDelet ~ dtaNB.3$Refugium)
wilcox.test(dtaNB.3$ModerDelet ~ dtaNB.3$Refugium)
wilcox.test(dtaNB.3$HighDelet ~ dtaNB.3$Refugium)

#Effect sizes SINGLE COPY
effectSizeNeut = rstatix::wilcox_effsize(data = dtaNB.3, formula = Neutral ~ Refugium)
effectSizeLow = rstatix::wilcox_effsize(data = dtaNB.3, formula = LowDelet ~ Refugium)
effectSizeMid = rstatix::wilcox_effsize(data = dtaNB.3, formula = ModerDelet ~ Refugium)
effectSizeHig = rstatix::wilcox_effsize(data = dtaNB.3, formula = HighDelet ~ Refugium)

#Wilcoxon tests HOMOZYGOUS
wilcox.test(dtaNBHomozyg.3$Neutral ~ dtaNBHomozyg.3$Refugium)
wilcox.test(dtaNBHomozyg.3$LowDelet ~ dtaNBHomozyg.3$Refugium)
wilcox.test(dtaNBHomozyg.3$ModerDelet ~ dtaNBHomozyg.3$Refugium)
wilcox.test(dtaNBHomozyg.3$HighDelet ~ dtaNBHomozyg.3$Refugium)

#Effect sizes HOMOZYGOUS
effectSizeNeut = rstatix::wilcox_effsize(data = dtaNBHomozyg.3, formula = Neutral ~ Refugium)
effectSizeLow = rstatix::wilcox_effsize(data = dtaNBHomozyg.3, formula = LowDelet ~ Refugium)
effectSizeMid = rstatix::wilcox_effsize(data = dtaNBHomozyg.3, formula = ModerDelet ~ Refugium)
effectSizeHig = rstatix::wilcox_effsize(data = dtaNBHomozyg.3, formula = HighDelet ~ Refugium)
