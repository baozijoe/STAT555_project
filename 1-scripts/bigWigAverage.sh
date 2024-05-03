#!/bin/bash

export PATH="/storage/work/fjz5078/Courses/STAT555/project/1-scripts/bedtools2/bin:$PATH"
export PATH="/storage/work/fjz5078/Courses/STAT555/project/1-scripts:$PATH"

cd /storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq

prefixes=("../CMP/ENCFF343PTQ" "../CMP/ENCFF832UUS" "../Erythroblast/ENCFF181AMY" "../Erythroblast/ENCFF616EWK" "../HSC/ENCFF662DYG" "../HSC/ENCFF255IVU")

cd CMP_vs_ERY
    
bigWigAverageOverBed "${prefixes[0]}.bw" CMP_vs_ERY.withpkname.bed "CMP1.tab"
bigWigAverageOverBed "${prefixes[1]}.bw" CMP_vs_ERY.withpkname.bed "CMP2.tab"
bigWigAverageOverBed "${prefixes[2]}.bw" CMP_vs_ERY.withpkname.bed "ERY1.tab"
bigWigAverageOverBed "${prefixes[3]}.bw" CMP_vs_ERY.withpkname.bed "ERY2.tab"

cd ..

cd HSC_vs_CMP

bigWigAverageOverBed "${prefixes[0]}.bw" HSC_vs_CMP.withpkname.bed "CMP1.tab"
bigWigAverageOverBed "${prefixes[1]}.bw" HSC_vs_CMP.withpkname.bed "CMP2.tab"
bigWigAverageOverBed "${prefixes[4]}.bw" HSC_vs_CMP.withpkname.bed "HSC1.tab"
bigWigAverageOverBed "${prefixes[5]}.bw" HSC_vs_CMP.withpkname.bed "HSC2.tab"
