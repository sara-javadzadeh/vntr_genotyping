#!/bin/bash

genotype_log_dir="/nucleus/projects/saraj/vntr/sources/pangenome_project/scripts/mendelian_consistency/logs_illumina_wgs_5k_short_inconsistent"
hifi_samples_dir="/ribosome/projects/saraj/vntr/data/pangenome_project/pacbio_hifi_mapped_reads"
illumina_samples_dir="/ribosome/projects/saraj/vntr/data/pangenome_project/illumina_mapped_reads"

vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.db"

reference_data="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
advntr_dir="/nucleus/projects/saraj/vntr/sources/advntr/cloned_advntr_aug_3_merge_pipelines/adVNTR/advntr"

export REF_PATH="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/ref/cache/%2s/%2s/%s"

#export REF_PATH=$reference_data
#export REF_CACHE=$reference_data

# To avoid the error about too many open files, increase ulimit from the default value of 1024.
ulimit -n 20000

# VIDs that were inconsistent.
#pedigree_file="pedigree_file.txt"
#echo "FamilyID SampleID FatherID MotherID" > $pedigree_file
# For illumina reads:
#for sample_name in $(ls ${illumina_samples_dir} | grep HG | head -n 10); do
while read line ; do
    sample_name=$(echo $line | awk '{print $1}')
    vids=$(echo $line | awk '{print $2}')
    echo "working with $sample_name vids $vids"
	echo "working with sample $sample_name"
	log_dir="$genotype_log_dir/$sample_name"
    echo "log_dir $log_dir"
	if [ ! -d $log_dir ]; then
		mkdir -p $log_dir
	fi
    command_output="$log_dir/command_output.txt"

    pushd $advntr_dir
        rm $command_output
        # Genotype Illumina reads
        input_bam_file="${illumina_samples_dir}/${sample_name}/${sample_name}.final.cram"
        output_file="$log_dir/${sample_name}_illumina_short_vntrs_genotype.txt"
        #if [ ! -e $output_file ]; then
            echo "output_file $output_file"
            /usr/bin/time -v python __main__.py genotype \
            --alignment_file $input_bam_file \
            --working_directory $log_dir \
            -m $vntr_db \
            --vntr_id $vids \
            --accuracy_filter \
            > $output_file
        #fi

        # Genotype Hifi reads
        input_bam_file="${hifi_samples_dir}/${sample_name}/${sample_name}_aligned_GRCh38_winnowmap.sorted.bam"
        output_file="$log_dir/${sample_name}_hifi_short_vntrs_genotype.txt"
        #if [ ! -e $output_file ]; then
            echo "output_file $output_file"
            /usr/bin/time -v python __main__.py genotype \
            --alignment_file $input_bam_file \
            --working_directory $log_dir \
            -m $vntr_db \
            --vntr_id $vids \
            --accuracy_filter \
            --pacbio \
            --log_pacbio_reads \
            > $output_file
        #fi
    popd
    aln_log_dir="$log_dir"
    mkdir -p $aln_log_dir

    #cp logs_illumina_wgs_5k_short_sr_3_txt/$sample_name/log_${sample_name}.final.cram.log $log_dir
    #cp ../../logs/logs_pairwise_align_10k_short/log_${sample_name}_aligned_GRCh38_winnowmap.sorted.bam.log_vid_${vid}.txt $aln_log_dir
    log_file_from_advntr_illumina="$log_dir/log_${sample_name}.final.cram.log"
    log_file_from_advntr_hifi="$log_dir/log_${sample_name}_aligned_GRCh38_winnowmap.sorted.bam.log"
    echo "----------------" >> $command_output
    echo "working with input file $log_file_from_advntr_illumina and vntr id $vid" >> $command_output
    for vid in $(echo $vids | tr ',' ' '); do
        aln_output_file="${aln_log_dir}/illumina_${sample_name}_vid_${vid}.txt"
        #if [ ! -e $aln_output_file ]; then
            /usr/bin/time -v python ${advntr_dir}/pairwise_aln_generator.py \
                -i $log_file_from_advntr_illumina \
                -o $aln_output_file \
                -db $vntr_db -vid $vid &>> $command_output
        #fi
        echo "working with input file $log_file_from_advntr_hifi and vntr id $vid" >> $command_output
        aln_output_file="${aln_log_dir}/hifi_${sample_name}_vid_${vid}.txt"
        #if [ ! -e $aln_output_file ]; then
            /usr/bin/time -v python ${advntr_dir}/pairwise_aln_generator.py \
            -i $log_file_from_advntr_hifi \
            -o $aln_output_file \
            -db $vntr_db -vid $vid &>> $command_output
        #fi
    done
done < outputs/inconsistent_hprc_no_str_sr_3_per_sample.txt

# Hifi logs are stored in ../../logs/logs_pairwise_align_10k_short/log_HG00438_aligned_GRCh38_winnowmap.sorted.bam.log_vid_67173.txt
