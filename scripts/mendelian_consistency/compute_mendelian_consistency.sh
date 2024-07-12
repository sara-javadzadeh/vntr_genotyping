#!/bin/bash

vcfs=""

# Print pedigree sample ids
if [ ! -e pedigree_hprc.txt ]; then
    sample_ids=$(ls logs_illumina_wgs_5k_short_txt/ | grep HG | grep -v error)
    for sample_id in $sample_ids; do
        parents=$(ls logs_illumina_wgs_5k_short_txt/$sample_id/log_*.final.cram.log | \
                  grep -v log_${sample_id}.final.cram.log | \
                  sed -e 's/.final.cram.log//g' | \
                  sed -e 's/log_//g')
        parents=$(basename --multiple $parents)
        echo $sample_id $sample_id $parents
    done
fi

#conda activate bioinfo
#for sample in HG03516  HG01106  HG01978  HG02886  HG00735  HG01928  HG02630  HG02622  HG00621\
#              HG02148  HG03453  HG00741  HG01952  HG02717  HG00673  HG01891  HG01361  HG00438; do
for sample in HG00438 HG00621 HG00673 HG00735 HG00741 HG01106 HG01123 HG01175 HG01258 HG01358\
              HG01361 HG01891 HG01928 HG01952 HG01978 HG02148 HG02257 HG02486 HG02559 HG02572\
              HG02622 HG02630 HG02717 HG02886 HG03453 HG03516 HG03540 HG03579; do


    for file in $(ls logs_illumina_wgs_5k_short/${sample}/*.vcf | grep -v sorted | grep -v merged_trio); do
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
merged_file="logs_illumina_wgs_5k_short/merged_trio"
if [ ! -e ${merged_file}.vcf ]; then
    echo "Merging VCF files with mergeSTR"
    mergeSTR --vcfs $vcfs --out $merged_file
fi
vcf_path=$merged_file.vcf
python compute_mendelian_consistency.py $vcf_path
