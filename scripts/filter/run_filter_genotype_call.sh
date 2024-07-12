#!/bin/bash

output="output_disease_all.txt"
rm -f $output; touch $output
samples_dir="/ribosome/projects/saraj/vntr/data/pangenome_project/pacbio_hifi_mapped_reads"

for sample_id in $(ls $samples_dir | grep HG) ; do
    python filter_genotype_call.py $sample_id &>> $output
done
#for sample_id in HG001; do #HG00621 HG00673 HG01175 HG01891 HG02148 HG02572 HG02886 HG03579; do
#    python filter_genotype_call.py $sample_id &>> $output
#done
