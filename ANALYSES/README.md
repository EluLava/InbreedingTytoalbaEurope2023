# ANALYSES

All the code is present in this folder:

- 0_RP504_filtering.sh: contains the filtering steps of the Reference panels (removing weird individuals, MAC 3 and missingness 10\%).

- 1.1_CreatingRecMap_RZoofriendlyplus10.sh: creates text files with RzooRoH-friendly positions including the regional recombination rates !

- 1.2_ROHsCalling_withGeneticPositionschangingPOSinVCF.sh: changes physical positions into genomic positions (estimated with the previous script 1.1) and calls HBD segments with RZooROH ! (RZooRoH uses tremendous amount of memory, we had to split the genome into super-scaffolds to run this script).

- 1.3_HBDsegments_FandDistEstimation.sh:  Estimates FHBD and cumulatie FHBD.

- 1.4_pHBD_perSite.sh: Estimate probability (among individuals) that a SNP is HBD for each of the different HBD classes.

- 2.1_Fas_estimation_perPOP.sh: estimates Fas general,  per population and with mean allele-sharing taken from an unrelated set of individuals for the Swiss, and also only with the unrelated set of individuals (all pops)

- 3.1_Defiining_MajorMinorSNPs_UNRset.sh: get the major and minor allele for each site according to a set of unrelated (swiss) and create a new VCF with 0 as MAJOR allele in this set and 1 as minor allele in this set!

- 3.2_RunningSNPEff.sh: all the code to run SNPEff.

- 4.1_CountMinorAlleleperINDV.sh: count number of minor alleles and homozygous minor alleles per individual.

- 4.2_Rxy_IslandsCont.sh: Estimates $R_{XY}$ and $R^{2}_{XY}$ ratios continental versus islands populations.

- 4.3_Rxy_IRefRecolo.sh:  Estimates $R_{XY}$ and $R^{2}_{XY}$ ratios continental refugium versus recolonized populations.

- 4.4_Rxy_IslandsCont_INTERGENIC.sh: Estimates $R_{XY}$ and $R^{2}_{XY}$ ratios continental versus islands populations on intergenic SNPs only.

- 4.5_Rxy_IRefRecolo_INTERGENIC.sh: Estimates $R_{XY}$ and $R^{2}_{XY}$ ratios continental refugium versus recolonized populations on intergenic SNPs only.

- 5.1_EstimatingNe.sh: Estimate pi (per population) per 1Mb windows and then Ne per population with 1,000 bootstraps.

- 6.1_CountingNBSegSiteperINDV.sh: count the number of segregating sites per individuals (figures S7 \& S8).


### functions folder

contains the code with the R functions used at different steps of the analyses !

### inputFILES folder

contains the input files needed for several scripts to run (such as list of individual names, list of super-scaffolds, list of unrelated individuals, etc.).
