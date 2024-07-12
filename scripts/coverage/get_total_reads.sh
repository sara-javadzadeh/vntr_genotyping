#!/bin/bash

long_vntr_bed_file="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.bed"
short_vntr_bed_file="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.bed"
all_10k_vntrs_bed_file="10k_vntrs.bed"

if [ ! -e $all_10k_vntrs_bed_file ]; then
    echo "Creating $all_10k_vntrs_bed_file"
    cat $long_vntr_bed_file > $all_10k_vntrs_bed_file
    cat $short_vntr_bed_file >> $all_10k_vntrs_bed_file
fi

samples_dir="/ribosome/projects/saraj/vntr/data/pangenome_project/pacbio_hifi_mapped_reads"
for sample_dir in $(ls $samples_dir | grep HG); do
    for file in $(ls $samples_dir/$sample_dir/*.bam); do
        # Compute all reads
        #echo $sample_dir >> num_reads_in_bamfiles.txt
        #samtools view --threads 10 -c $file  >> num_reads_in_bamfiles.txt
        # Compute reads overlapping target regions
        #echo $sample_dir >> num_reads_in_bamfiles_overlapping_db.txt
        #/usr/bin/time -v samtools view -L $all_10k_vntrs_bed_file --threads 10 $file | samtools depth | awk '{sum+=$3} END { print "Average = ",sum/NR}'  >> coverage_overlapping_vntr_db.txt
        #exit 0
        # Compute coverage
        #echo $sample_dir >> sequencing_coverage.txt
        depth_file=${sample_dir}_depth.depth
        if [ ! -e $depth_file ]; then
            echo "creating depth file"
            #time samtools depth $file >> $depth_file
        fi
        output_file="num_bases.txt"
        echo "$sample_dir, $(wc -l $depth_file)" >> $output_file
        #echo $sample_dir >> $output_file
        #cat $depth_file | awk 'BEGIN {count=0} {if (int($3) >=15) count +=1} END {print count}' >> $output_file
        #time cat $depth_file | awk '{sum+=$3} END { print "Average = ",sum/NR}' >> sequencing_coverage.txt
        #echo $sample_dir >> base_coverage_histogram.txt
        #/usr/bin/time -v samtools depth -a --threads 10 $file >> sequencing_coverage.txt
        #/usr/bin/time -v cat $depth_file | awk '{print($3)}' | sort | uniq -c >> base_coverage_histogram.txt
        # Removing the depth file as it takes a lot of space
        #rm $depth_file
    done
done

