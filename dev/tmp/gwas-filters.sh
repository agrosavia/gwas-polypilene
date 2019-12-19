#!/bin/bash

#-------------------------------------------------------------

GENO=$1

mkdir -p tmp
CMM="ln -f -s $PWD/$1.ped tmp/"
echo $CMM; eval $CMM
CMM="ln -f -s $PWD/$1.map tmp/"
echo $CMM; eval $CMM

echo ""
echo "-----------------------------------------------------"
echo " Convert text to binary format"
echo " Filter missingness per sample"
echo " Filter SNPs with a low minor allele frequency (MAF)"
echo " Filter SNPs which are not in Hardy-Weinberg equilibrium (HWE)."
echo " Filter individuals closely related (Cryptic relatedness)"
echo "-----------------------------------------------------"
CMM1="plink --file $GENO --mind 0.1 --maf 0.01 --hwe 1e-10 --genome --min 0.0001 --out $GENO-QC --make-bed"
CMM="plink --file $GENO --mind 0.1 --maf 0.01 --hwe 1e-10 --out $GENO-QC --make-bed"
echo $CMM; eval $CMM
plink --bfile $GENO-QC --recode tab --out $GENO-QC

echo "Copy links of filtered plink files to main dir"
CMM="cp $PWD/$GENO-QC.ped out/filtered-plink-genotype.ped"
echo $CMM; eval $CMM
CMM="cp $PWD/$GENO-QC.map out/filtered-plink-genotype.map"
echo $CMM; eval $CMM

echo ""
echo "-----------------------------------------------------"
echo " Get final markers and individuals"
echo "-----------------------------------------------------"
CMM="cut -f 2 $GENO-QC.ped > out/filtered-names-samples.tbl"
echo -e "\n>>> " $CMM "\n"
eval $CMM

CMM="cut -f 2 $GENO-QC.map > out/filtered-names-markers.tbl"
echo -e "\n>>> " $CMM "\n"
eval $CMM


