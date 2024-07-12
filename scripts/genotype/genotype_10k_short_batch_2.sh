#!/bin/bash

#PBS -q home
#PBS -N advntr_COH
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -o output_file
#PBS -e error_file
#PBS -V
#PBS -m ae

log_dir="/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_10k_short"
samples_dir="/ribosome/projects/saraj/vntr/data/pangenome_project/pacbio_hifi_mapped_reads"

# Doesn't include added disease associated VNTRs.
#all_vntr_db="/pedigree2/projects/jonghun/filtering_similar_vntrs/hg38_VNTRs_by_TRF.db"
# Most updated VNTR DB for now.
#all_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/combined_trf_hg38/hg38_VNTRs_by_TRF.db"
#long_vntrs="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.db"
short_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.db"

# Create a directory for the logs and outputs.
mkdir -p $log_dir

# Raise the number of allowed open files by the OS. Default is 1024.
# This is to avoid the error: IOError: [Errno 24] Too many open files
#ulimit -n 2048


for input_sample_dir in $(cat $samples_dir/sample_ids.txt | tail -n 8); do
    #HG03540 HG02622 HG01258 HG02559 HG02257 HG02148 HG02717 HG02630 HG02486 HG01175 HG01123 HG01891; do
    # Remove previous output file.
    output="$log_dir/output_per_sample_${input_sample_dir}.txt"
    rm -f $output; touch $output
    for input_bam_file in $(ls $samples_dir/$input_sample_dir/*.bam); do
        echo "----------------" >> $output
        echo "working with input file $input_bam_file" >> $output
        # To avoid the error about too many open files, increase ulimit from the default value of 1024.
        ulimit -n 4096
        /usr/bin/time -v advntr genotype --alignment_file $input_bam_file --working_directory $log_dir -m $short_vntr_db --pacbio --log_pacbio_reads &>> $output
    done
done
