cd ..

#subset the VCF with only one indv and rm monomorphic SNPs
parallel -j 20 "bcftools view -s {} -c 1:minor -H -Oz ./1.1_VCFs/RP502_FINALnames_TF1_Mask_indDP_biall_mac3_missing90.vcf.gz | wc -l > ./6_MinorAllelesCount/tmpVCFs/SNPsCount_{}.txt" ::: $(cat ./inputFILES/RP502_NewNames.list)

#create header
echo -e "INDV SNPsNB" > ./Analyses/NBsnpsPerINDV.txt

#merge all outputs
for file in $(find ./6_MinorAllelesCount/tmpVCFs/ -name "SNPsCount*.txt"); do

    #Extract INDV name
    INDV=$(basename -s ".txt" ${file} | cut -d'_' -f2)

    echo "${INDV} " $(cat ${file}) >> ./Analyses/NBsnpsPerINDV.txt

done

#rm all files
rm ./6_MinorAllelesCount/tmpVCFs/SNPsCount*.txt