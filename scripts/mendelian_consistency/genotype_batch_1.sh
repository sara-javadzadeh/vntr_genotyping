#!/bin/bash

# Run for Illumina reads
parent_log_dir="$(pwd)/logs_illumina_wgs_5k_short_sr_3_txt"
illumina_samples_dir="/ribosome/projects/saraj/vntr/data/pangenome_project/illumina_mapped_reads"
illumina_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.db"


reference_data="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
advntr_dir="/nucleus/projects/saraj/vntr/sources/advntr/cloned_advntr_aug_3_merge_pipelines/adVNTR/advntr"

export REF_PATH="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/ref/cache/%2s/%2s/%s"

#export REF_PATH=$reference_data
#export REF_CACHE=$reference_data

# To avoid the error about too many open files, increase ulimit from the default value of 1024.
ulimit -n 10000

#pedigree_file="pedigree_file.txt"
#echo "FamilyID SampleID FatherID MotherID" > $pedigree_file
# For illumina reads:
for sample_name in $(ls ${illumina_samples_dir} | grep HG | head -n 10); do
#for sample_name in HG01123 HG02486 HG02559; do
	pushd $advntr_dir
	echo "working with sample $sample_name"
	log_dir="$parent_log_dir/$sample_name"
	if [ ! -d $log_dir ]; then
		mkdir -p $log_dir
	fi

    ## Genotype the child
    input_bam_file="${illumina_samples_dir}/${sample_name}/${sample_name}.final.cram"
    pedigree_line="$sample_name $sample_name"
    # Remove previous output file.
    output_vcf="$log_dir/${sample_name}_short_vntr_dataset.txt"
    #if [ ! -e $output_vcf ]; then
        # To avoid the error about too many open files, increase ulimit from the default value of 1024.
        /usr/bin/time -v python __main__.py genotype \
        --alignment_file $input_bam_file \
        --working_directory $log_dir \
        -m $illumina_db \
        --accuracy_filter \
        --disable_logging \
        > $output_vcf
        #--outfmt vcf \
    #fi

	## Genotype the parents
	for parent_sample_name in $(ls ${illumina_samples_dir}/${sample_name}/parents); do
        pedigree_line="$pedigree_line $parent_sample_name"
		input_bam_file="${illumina_samples_dir}/${sample_name}/parents/$parent_sample_name/${parent_sample_name}.final.cram"
		# Remove previous output file.
		output_vcf="$log_dir/${parent_sample_name}_short_vntr_dataset.txt"
        #if [ ! -e $output_vcf ]; then
            /usr/bin/time -v python __main__.py genotype \
            --alignment_file $input_bam_file \
            --working_directory $log_dir \
            -m $illumina_db \
            --accuracy_filter \
            --disable_logging \
            > $output_vcf
            #--outfmt vcf \
        #fi
	done

    #echo $pedigree_line > $pedigree_file

	popd
done
