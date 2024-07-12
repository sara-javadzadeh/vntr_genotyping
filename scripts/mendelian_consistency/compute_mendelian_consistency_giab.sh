#!/bin/bash

vcfs=""

#conda activate bioinfo
for sample in HG002 HG003 HG004; do
    for file in $(ls ../../logs/logs_genotype_giab_hifi_integrated_pacbcloud/${sample}/combined/${sample}*.vcf | grep -v sorted | grep -v merged_trio); do
        #echo "working on $file"
        sorted_file=$(echo $file | sed -e 's/\.vcf/_sorted\.vcf/g')
        if [ ! -e $sorted_file.gz.tbi ]; then
            echo "sorting file $file"
            bcftools sort -o ${sorted_file} $file
            echo "bgzip file $sorted_file"
            bgzip -c $sorted_file > $sorted_file.gz
            echo "index file $sorted_file.gz"
            bcftools index -t $sorted_file.gz
        fi
        vcfs="$sorted_file.gz,$vcfs"
    done
done
vcfs=$(echo $vcfs | sed -e 's/,$//g')
#conda activate trtools
merged_file="../../logs/logs_genotype_giab_hifi_integrated_pacbcloud/merged_trio_10k"
if [ ! -e ${merged_file}.vcf ]; then
    echo "Merging VCF files with mergeSTR"
    mergeSTR --vcfs $vcfs --out $merged_file
fi
vcf_path=$merged_file.vcf
python compute_mendelian_consistency.py $vcf_path
