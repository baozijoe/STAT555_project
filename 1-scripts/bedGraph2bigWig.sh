#!/bin/bash

export PATH="/storage/work/fjz5078/Courses/STAT555/project/1-scripts/bedtools2/bin:$PATH"
export PATH="/storage/work/fjz5078/Courses/STAT555/project/1-scripts:$PATH"

cd /storage/work/fjz5078/Courses/STAT555/project/0-rawdata/ATAC-seq

prefixes=("./CMP/ENCFF343PTQ" "./CMP/ENCFF832UUS" "./Erythroblast/ENCFF181AMY" "./Erythroblast/ENCFF616EWK" "./HSC/ENCFF662DYG" "./HSC/ENCFF255IVU")

for prefix in "${prefixes[@]}"; do
    # file name
    file="${prefix}.bed"

    # check existence
    if [ -f "$file" ]; then
        
        cut -f 1,2,3,4 "$file" | sort -k 1,1 -k 2,2n > "${prefix}.bedGraph"
     
        cat "${prefix}.bedGraph" | awk -F '\t' -v OFS='\t' '{print $1,$2,$3-1,$4}' > "${prefix}.bedGraph.tmp1"

        bedtools sort -i "${prefix}.bedGraph.tmp1" | bedtools merge -i - -c 4 -o mean > "${prefix}.bedGraph.tmp2"

        cat "${prefix}.bedGraph.tmp2" | awk -F '\t' -v OFS='\t' '{print $1,$2,$3+1,$4}' > "${prefix}.merge.bedGraph"
        
        bedGraphToBigWig "${prefix}.merge.bedGraph" mm10.chrom.sizes "${prefix}.bw"

        rm "${prefix}.bedGraph.tmp1" "${prefix}.bedGraph.tmp2"
    else
        echo "File not found: $file"
    fi
done
