#!/bin/bash

GENO=$1
PHENO=$2

cmm="plink --file $1 --make-bed --out plink1"
echo $cmm
eval $cmm

cmm="plink2 --bfile plink1 --make-pgen --out plink2"
echo $cmm
eval $cmm


cmm="plink2 --pfile plink2 --make-king-table --out plink3-kingship.table"
echo $cmm
eval $cmm

# It is conventional to use a cutoff of ~0.354 (the geometric mean of 0.5 and 0.25) 
# to screen for monozygotic twins and duplicate samples, ~0.177 to add first-degree relations, etc.
cmm="plink2 --pfile plink2 --king-cutoff 0.177 --out plink3-kingship"
echo $cmm
eval $cmm
