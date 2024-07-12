#!/bin/bash

# Run filtering script.
sample_ids=$(ls logs_illumina_wgs_5k_short_txt | grep HG | grep -v error)
suffix="v4_illumina_pacbio_comparison_no_str_sr_3_inc_def"
for sample_id in $sample_ids; do
    echo "Working with $sample_id"
    working_dir="$suffix/comparison_${sample_id}"
    mkdir -p $working_dir
    touch $working_dir/${sample_id}_comparison_report.txt
    output_log=$working_dir/"sample_${sample_id}_illumina_pacbio_compare_${suffix}.txt"
    echo "--------------" >> $output_log
    # Not applying filter as the filter is integrated in the advntr code now.
    python Illumina_Pacbio_comparison.py ${sample_id} $suffix > $output_log
done
