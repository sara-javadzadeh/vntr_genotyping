#!/bin/bash

#PBS -q home
#PBS -N advntr_COH
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -o output_file
#PBS -e error_file
#PBS -V
#PBS -m ae

log_dir="/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_high_coverage_subset_w_advntr_1.5.0"
output="$log_dir/output.txt"
samples_dir="/ribosome/projects/saraj/vntr/data/pangenome_project/pacbio_hifi_mapped_reads"

# We work with two databases for now. Use one at each run.
short_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.db"
long_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.db"

# Target set of VNTR ids.
#target_genotype="target_vntrs/disease_associated_long_vntrs.txt"
#target_genotype="target_vntrs/disease_associated_short_vntrs.txt"
target_genotype="target_vntrs/short_vntrs_high_coverage.txt"

# Create a directory for the logs and outputs.
mkdir -p $log_dir
# Remove previous output file.
rm -f $output; touch $output

# Get the VNTR IDs from the corresponding file to pass on the the genotype command.
vids=$(cat $target_genotype | tr '\n' ',' | sed -E 's/,$//g')

for input_sample_dir in $(ls $samples_dir | grep HG); do
    for input_bam_file in $(ls $samples_dir/$input_sample_dir/*.bam); do
        echo "----------------" >> $output
        echo "working with input file $input_bam_file" >> $output
        /usr/bin/time -v advntr genotype --alignment_file $input_bam_file --working_directory $log_dir -m $short_vntr_db --pacbio --log_pacbio_reads -vid $vids &>> $output
    done
done
