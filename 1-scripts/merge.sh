#!/bin/bash

export PATH="/storage/work/fjz5078/Courses/STAT555/project/1-scripts/bedtools2/bin:$PATH"
export PATH="/storage/work/fjz5078/Courses/STAT555/project/1-scripts:$PATH"

cd /storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq

prefixes=("./CMP/ENCFF343PTQ" "./CMP/ENCFF832UUS" "./Erythroblast/ENCFF181AMY" "./Erythroblast/ENCFF616EWK" "./HSC/ENCFF662DYG" "./HSC/ENCFF255IVU")

# merge replicates' bed files and regions

cat "${prefixes[0]}.bed" "${prefixes[1]}.bed" | sort -k1,1 -k2,2n | bedtools merge > ./CMP/CMP_merged.bed
cat "${prefixes[2]}.bed" "${prefixes[3]}.bed" | sort -k1,1 -k2,2n | bedtools merge > ./Erythroblast/ERY_merged.bed
cat "${prefixes[4]}.bed" "${prefixes[5]}.bed" | sort -k1,1 -k2,2n | bedtools merge > ./HSC/HSC_merged.bed


mkdir -p CMP_vs_ERY
cd CMP_vs_ERY

cat ../CMP/CMP_merged.bed ../Erythroblast/ERY_merged.bed | sort -k1,1 -k2,2n | bedtools merge > CMP_vs_ERY.bed
cat CMP_vs_ERY.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' > CMP_vs_ERY.withpkname.bed

cd ..

mkdir -p HSC_vs_CMP
cd HSC_vs_CMP

cat ../HSC/HSC_merged.bed ../CMP/CMP_merged.bed | sort -k1,1 -k2,2n | bedtools merge > HSC_vs_CMP.bed
cat HSC_vs_CMP.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' > HSC_vs_CMP.withpkname.bed

