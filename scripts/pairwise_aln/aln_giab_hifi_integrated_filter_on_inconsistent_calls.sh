#!/bin/bash


genotype_log_dir_base="/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_giab_hifi_integrated_inconsistent_vntrs"
aln_log_dir="/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_pairwise_aln_giab_inconsistent"
command_output="$aln_log_dir/run_aln_command_output.txt"

# We work with two databases for now. Use one at each run.
#all_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/combined_trf_hg38/hg38_VNTRs_by_TRF.db"
#long_vntrs="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.db"
short_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.db"
vntr_db=$short_vntr_db

# Path to run the pairwise align script from.
advntr_path="/nucleus/projects/saraj/vntr/sources/advntr/cloned_advntr_aug_3_merge_pipelines/adVNTR/advntr"

inconsistent_vntrs_file="/nucleus/projects/saraj/vntr/sources/pangenome_project/scripts/mendelian_consistency/inconsistent_mc_giab_vids_short.txt"
target_genotype=$inconsistent_vntrs_file

mkdir -p $aln_log_dir


#export REF_PATH="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
#export REF_CACHE="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/ref/cache/%2s/%2s/%s"

#export REF_PATH=$reference_data
#export REF_CACHE=$reference_data

# To avoid the error about too many open files, increase ulimit from the default value of 1024.
ulimit -n 4096

# For illumina reads:
for sample_name in HG002 HG003 HG004; do

    genotype_log_dir="${genotype_log_dir_base}/${sample_name}/short_vntr_db"
    for log_file_from_advntr in $(ls -t ${genotype_log_dir}/*.log); do
        for vid in $(cat $target_genotype); do
            echo "----------------" >> $command_output
            echo "working with input file $log_file_from_advntr and vntr id $vid" >> $command_output
            log_file_base_name=$(basename ${log_file_from_advntr})
            aln_output_file="${aln_log_dir}/${log_file_base_name}_vid_${vid}.txt"
            /usr/bin/time -v python ${advntr_path}/pairwise_aln_generator.py -i $log_file_from_advntr -o $aln_output_file -db $vntr_db -vid $vid &>> $command_output
        done
    done




done
