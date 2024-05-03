#!/bin/bash

export PATH="/storage/work/fjz5078/Courses/STAT555/project/1-scripts/bedtools2/bin:$PATH"
export PATH="/storage/work/fjz5078/Courses/STAT555/project/1-scripts:$PATH"

cd /storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq/Hierachy

prefixes=("../CMP/ENCFF343PTQ" "../CMP/ENCFF832UUS" "../Erythroblast/ENCFF181AMY" "../Erythroblast/ENCFF616EWK" "../HSC/ENCFF662DYG" "../HSC/ENCFF255IVU")

cat ../CMP/CMP_merged.bed ../Erythroblast/ERY_merged.bed ../HSC/HSC_merged.bed | sort -k1,1 -k2,2n | bedtools merge > CMP_ERY_HSC.bed
cat CMP_ERY_HSC.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' > CMP_ERY_HSC.withpkname.bed

bigWigAverageOverBed "${prefixes[0]}.bw" CMP_ERY_HSC.withpkname.bed "CMP1.tab"
bigWigAverageOverBed "${prefixes[1]}.bw" CMP_ERY_HSC.withpkname.bed "CMP2.tab"

bigWigAverageOverBed "${prefixes[2]}.bw" CMP_ERY_HSC.withpkname.bed "ERY1.tab"
bigWigAverageOverBed "${prefixes[3]}.bw" CMP_ERY_HSC.withpkname.bed "ERY2.tab"

bigWigAverageOverBed "${prefixes[4]}.bw" CMP_ERY_HSC.withpkname.bed "HSC1.tab"
bigWigAverageOverBed "${prefixes[5]}.bw" CMP_ERY_HSC.withpkname.bed "HSC2.tab"