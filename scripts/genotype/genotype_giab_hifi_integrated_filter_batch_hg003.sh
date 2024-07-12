#!/bin/bash

# Run for Illumina reads
parent_log_dir="/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_giab_hifi_integrated_pacbcloud"
samples_dir="/nucleus/projects/saraj/vntr/data/GIAB"
short_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.db"
long_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.db"

#reference_data="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
advntr_dir="/nucleus/projects/saraj/vntr/sources/advntr/cloned_advntr_aug_3_merge_pipelines/adVNTR/advntr"

#export REF_PATH="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
#export REF_CACHE="/nucleus/projects/saraj/vntr/data/pangenome_project/hpgp_trio_ref_data/ref/cache/%2s/%2s/%s"

#export REF_PATH=$reference_data
#export REF_CACHE=$reference_data

# To avoid the error about too many open files, increase ulimit from the default value of 1024.
ulimit -n 4096

# For illumina reads:
for sample_name in HG003; do
	pushd $advntr_dir
	echo "working with sample $sample_name"
	log_dir="$parent_log_dir/$sample_name"
	if [ ! -d $log_dir ]; then
		mkdir -p $log_dir
		mkdir -p $log_dir/short_vntr_db
		mkdir -p $log_dir/long_vntr_db
	fi

    ## Genotype the child
    for input_bam_file in $(ls ${samples_dir}/${sample_name}_pacbcloud/*.bam); do
        # Short VNTR db
        bam_basename=$(basename $input_bam_file | head -c 12)
        output_vcf="$log_dir/short_vntr_db/${bam_basename}_short_vntr_db.vcf"
        echo "Genotyping $input_bam_file to output $output_vcf"
        /usr/bin/time -v python -u __main__.py genotype \
        --alignment_file $input_bam_file \
        --working_directory $log_dir/short_vntr_db \
        -m $short_vntr_db \
        --accuracy_filter \
        --outfmt vcf \
        --pacbio \
        --disable_logging \
        > $output_vcf

        # Long VNTR db
        bam_basename=$(basename $input_bam_file | head -c 12)
        output_vcf="$log_dir/long_vntr_db/${bam_basename}_long_vntr_db.vcf"
        /usr/bin/time -v python -u __main__.py genotype \
        --alignment_file $input_bam_file \
        --working_directory $log_dir/long_vntr_db \
        -m $long_vntr_db \
        --accuracy_filter \
        --outfmt vcf \
        --pacbio \
        --disable_logging \
        > $output_vcf
    done

	popd
done
