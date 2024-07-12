#!/bin/bash


genotype_log_dir="/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_10k_short"
aln_log_dir="/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_pairwise_align_10k_short"
command_output="$aln_log_dir/run_aln_command_output_batch_0.txt"

# We work with two databases for now. Use one at each run.
#all_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/combined_trf_hg38/hg38_VNTRs_by_TRF.db"
#long_vntrs="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.db"
short_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.db"
vntr_db=$short_vntr_db

# Path to run the pairwise align script from.
adnvtr_path="/nucleus/projects/saraj/vntr/sources/COH_analysis/scripts/advntr-1.5.0/adVNTR/advntr"

target_genotype="/nucleus/projects/saraj/vntr/sources/pangenome_project/target_vntrs/10k_short.txt"
# Create a directory for the logs and outputs.
mkdir -p $aln_log_dir
# Remove previous output file.
#rm -f $command_output; touch $command_output
for log_file_from_advntr in $(ls -t ${genotype_log_dir}/*.log | head -n 10); do
    for vid in $(cat $target_genotype); do
        echo "----------------" >> $command_output
        echo "working with input file $log_file_from_advntr and vntr id $vid" >> $command_output
        log_file_base_name=$(basename ${log_file_from_advntr})
        aln_output_file="${aln_log_dir}/${log_file_base_name}_vid_${vid}.txt"
        /usr/bin/time -v python ${adnvtr_path}/pairwise_aln_generator.py -i $log_file_from_advntr -o $aln_output_file -db $vntr_db -vid $vid &>> $command_output
    done
done
