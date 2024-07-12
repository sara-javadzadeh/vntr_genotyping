#!/bin/bash



genotype_log_dir="/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_ensemble_validation_set"
aln_log_dir="/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_pairwise_align_ensemble_validation_set"
command_output="$aln_log_dir/run_aln_command_output.txt"

# We work with two databases for now. Use one at each run.
all_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/combined_trf_hg38/hg38_VNTRs_by_TRF.db"
vntr_db=$all_vntr_db

# Path to run the pairwise align script from.
adnvtr_path="/nucleus/projects/saraj/vntr/sources/COH_analysis/scripts/advntr-1.5.0/adVNTR/advntr"

# Create a directory for the logs and outputs.
mkdir -p $aln_log_dir
# Remove previous output file.
#rm -f $command_output; touch $command_output
for log_file_from_advntr in $(ls -t ${genotype_log_dir}/*.log | grep "HG001"); do
    for vid in 290964 4237 674126 915594 123860 515774 525610 30693 564306; do
        echo "----------------" >> $command_output
        echo "working with input file $log_file_from_advntr and vntr id $vid" >> $command_output
        log_file_base_name=$(basename ${log_file_from_advntr})
        aln_output_file="${aln_log_dir}/${log_file_base_name}_vid_${vid}.txt"
        /usr/bin/time -v python ${adnvtr_path}/pairwise_aln_generator.py -i $log_file_from_advntr -o $aln_output_file -db $vntr_db -vid $vid &>> $command_output
    done
done
