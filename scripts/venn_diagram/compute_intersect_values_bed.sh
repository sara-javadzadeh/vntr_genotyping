#!/bin/bash


function print {
    # prints unmerged and merged stats of a bed file
    filename=$1
    file_basename=$(basename $filename)
    echo "$file_basename $(wc -l $filename | awk '{print $1}') \
      after merging $(cat $filename | bedtools merge -i stdin | wc -l | awk '{print $1}')"
}

bed_files="bedfiles"

# From bedtools manual
## intersect
# -u  Write original A entry once if any overlaps found in B. In other words, just report the fact at least one overlap was found in B.
# -v Only report those entries in A that have no overlap in B.

## window
#Summary: Examines a "window" around each feature in A and
#         reports all features in B that overlap the window. For each
#         overlap the entire entry in A and B are reported.
#-l      Base pairs added upstream (left of) of each entry
#		in A when searching for overlaps in B.
#		- Allows one to define assymterical "windows".
#		- Default is 1000 bp.
#		- (INTEGER)
#
#-r      Base pairs added downstream (right of) of each entry
#		in A when searching for overlaps in B.
#		- Allows one to define assymterical "windows".
#		- Default is 1000 bp.
#		- (INTEGER)

regions="/nucleus/projects/saraj/vntr/sources/vntr_database_comparison/hg38_vntr_no_complicated_region.tab"
gene_proximal_vntrs_unmerged="$bed_files/gene_proximal_vntrs_unmerged.bed"
gene_proximal_vntrs="$bed_files/gene_proximal_vntrs.bed"
all_400k_vntrs_unsorted="/nucleus/projects/saraj/vntr/sources/vntr_database_comparison/advntr_database_no_comp_regions.bed"
all_400k_vntrs_w_comp_region="$bed_files/advntr_database_w_comp_regions_sorted.bed"
all_400k_vntrs="$bed_files/advntr_database_no_comp_regions_sorted.bed"
disease_vntrs_unsorted="/nucleus/projects/saraj/vntr/sources/vntr_database_comparison/disease_associated_sorted.bed"
disease_vntrs="$bed_files/disease_associated_vntrs.bed"
line_sine_centr_regions="/nucleus/projects/saraj/vntr/sources/vntr_database_comparison/line_sine_centr_regions_hg38_sorted.bed"
segdup_regions="/nucleus/projects/saraj/vntr/sources/vntr_database_comparison/segdup_regions_hg38_sorted.bed"
genic_low_coverage_vids_unsorted="/nucleus/projects/saraj/vntr/sources/COH_analysis/scripts/compare_short_long_vntrs/coverage_analysis/low_coverage_logs/raw_vntr_ids_with_low_coverage_across_all_7_samples.txt"
genic_low_coverage_vids="$bed_files/genic_low_coverage_vids.txt"
genic_low_coverage_low_gc_vids_unsorted="/nucleus/projects/saraj/vntr/sources/COH_analysis/scripts/compare_short_long_vntrs/coverage_analysis/low_coverage_logs/low_coverage_low_gc_vids.txt"
genic_low_coverage_low_gc_vids="$bed_files/low_coverage_low_gc_vids.txt"


# Require 50% of intervals in a to overlap with intervals in b to count as overlap
parameters=""
#parameters=" -f 0.5 "

#long_gene_vntrs="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.bed"
long_gene_vntrs="$bed_files/Pacbio_probed_long_vntrs_sorted.bed"
#short_gene_vntrs="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.bed"
short_gene_vntrs="$bed_files/illumina_probed_short_vntrs_sorted.bed"
strs="../vntrs_to_discard/strs.bed"

targeted_vntrs_unsorted="$bed_files/9999_vntr_ids.bed"
targeted_vntrs="$bed_files/9999_vntrs_sorted.bed"
bedtools sort -i $targeted_vntrs_unsorted | uniq > $targeted_vntrs


echo "Temporarily not merging before intersect calculations"
#echo "Now merging before intersect calculations"

#bedtools sort -i $bed_files/Pacbio_probed_long_vntrs.bed | uniq | bedtools merge -i stdin > $long_gene_vntrs
bedtools sort -i $bed_files/Pacbio_probed_long_vntrs.bed | uniq > $long_gene_vntrs
bedtools sort -i $bed_files/illumina_probed_short_vntrs.bed > $bed_files/illumina_probed_short_vntrs_sorted_w_strs.bed
bedtools intersect -v -a $bed_files/illumina_probed_short_vntrs_sorted_w_strs.bed -b $strs | uniq | bedtools sort -i stdin > $short_gene_vntrs
#bedtools intersect -v -a $bed_files/illumina_probed_short_vntrs_sorted_w_strs.bed -b $strs | uniq | bedtools sort -i stdin | bedtools merge -i stdin > $short_gene_vntrs

cat $long_gene_vntrs > $bed_files/gene_proximal_vntrs_unsorted.bed
cat $short_gene_vntrs >> $bed_files/gene_proximal_vntrs_unsorted.bed

bedtools sort -i $bed_files/gene_proximal_vntrs_unsorted.bed > $bed_files/gene_proximal_vntrs_sorted.bed
bedtools intersect -v -a $bed_files/gene_proximal_vntrs_sorted.bed -b $strs | uniq > $gene_proximal_vntrs_unmerged
cat $gene_proximal_vntrs_unmerged > $gene_proximal_vntrs
#bedtools merge -i $gene_proximal_vntrs_unmerged > $gene_proximal_vntrs

print $gene_proximal_vntrs


for window_size in 0; do
    echo "How many disease or genic VNTRs are in proximity ($window_size bp) of LINE/SINE/SEGDUP regions? "
    #disease_vntrs_in_complex_regions="$bed_files/disease_vntrs_overlap_with_line_sine_segdup_region.bed"
    #bedtools window -u -w $window_size -a $disease_vntrs -b $line_sine_centr_regions $segdup_regions > $disease_vntrs_in_complex_regions
    #print $disease_vntrs_in_complex_regions

    unsorted_genic_vntrs_in_complex_regions="$bed_files/unsorted_genic_vntrs_overlap_with_line_sine_segdup_region.bed"
    genic_vntrs_in_complex_regions="$bed_files/genic_vntrs_overlap_with_line_sine_segdup_region.bed"
    bedtools window -u -w $window_size -a $gene_proximal_vntrs -b $segdup_regions > $unsorted_genic_vntrs_in_complex_regions
    bedtools window -u -w $window_size -a $gene_proximal_vntrs -b $line_sine_centr_regions >> $unsorted_genic_vntrs_in_complex_regions
    cat $unsorted_genic_vntrs_in_complex_regions | sort --unique -o $unsorted_genic_vntrs_in_complex_regions
    bedtools sort -i $unsorted_genic_vntrs_in_complex_regions > $genic_vntrs_in_complex_regions

    print $genic_vntrs_in_complex_regions
    # Coverage analysis
    # Get VIDs from bed file
    genic_vids_in_complex_regions="$bed_files/genic_vids_overlap_with_line_sine_segdup_region.txt"
    cat $genic_vntrs_in_complex_regions | awk '{print $4}' > $genic_vids_in_complex_regions

    echo "low coverage VNTRs"
    #wc -l $genic_low_coverage_vids_unsorted
    # Remove the extra character in vids for low coverage vntrs
    cat $genic_low_coverage_vids_unsorted | tr -d sl | sort --unique --numeric-sort > $genic_low_coverage_vids
    wc -l $genic_low_coverage_vids
    echo "How many low coverage VNTRs are in proximity ($window_size bp) of LINE/SINE/SEGDUP regions? "
    grep -w -f $genic_low_coverage_vids $genic_vids_in_complex_regions | wc -l


    echo "low coverage low gc VNTRs"
    #wc -l $genic_low_coverage_low_gc_vids_unsorted
    # Remove the extra character in vids for low coverage vntrs
    cat $genic_low_coverage_low_gc_vids_unsorted| tr -d sl | sort --unique --numeric-sort > $genic_low_coverage_low_gc_vids
    wc -l $genic_low_coverage_low_gc_vids
    echo "How many low coverage low gc VNTRs are in proximity ($window_size bp) of LINE/SINE/SEGDUP regions? "
    grep -w -f $genic_low_coverage_low_gc_vids $genic_vids_in_complex_regions | wc -l
done



bedtools sort -i $all_400k_vntrs_unsorted | uniq | bedtools merge -i stdin > $all_400k_vntrs_w_comp_region
bedtools intersect -u -a $all_400k_vntrs_w_comp_region -b $bed_files/hg38_vntr_no_complicated_region.tab \
    > $all_400k_vntrs
print $all_400k_vntrs

bedtools sort -i $disease_vntrs_unsorted | uniq  > $disease_vntrs

print $disease_vntrs

echo "      Gene proximal and disease"
bedtools intersect $parameters -u -a $gene_proximal_vntrs -b $disease_vntrs \
    > $bed_files/genic_vntrs_in_disease.bed
    print $bed_files/genic_vntrs_in_disease.bed
bedtools intersect $parameters -u -a $disease_vntrs -b $gene_proximal_vntrs \
    > $bed_files/disease_vntrs_in_gene_proximal.bed
    print $bed_files/disease_vntrs_in_gene_proximal.bed

echo "      Gene proximal and targeted "
bedtools intersect $parameters -u -a $targeted_vntrs -b $gene_proximal_vntrs \
    > $bed_files/targeted_vntrs_in_genic.bed
    print $bed_files/targeted_vntrs_in_genic.bed
bedtools intersect $parameters -u -a $disease_vntrs -b $targeted_vntrs \
    > $bed_files/disease_vntrs_in_targeted.bed
    print $bed_files/disease_vntrs_in_targeted.bed

echo "      400k and disease"
bedtools intersect $parameters -u -a $disease_vntrs -b $all_400k_vntrs \
    > $bed_files/disease_vntrs_in_400k.bed
    print $bed_files/disease_vntrs_in_400k.bed
bedtools intersect $parameters -u -a $all_400k_vntrs -b $disease_vntrs \
    > $bed_files/400k_in_disease_vntrs.bed
    print $bed_files/400k_in_disease_vntrs.bed

echo "      all three sets"
bedtools intersect $parameters -u -a $bed_files/disease_vntrs_in_400k.bed -b $gene_proximal_vntrs \
    > $bed_files/disease_vntrs_in_400k_and_gene_proximal.bed
    print $bed_files/disease_vntrs_in_400k_and_gene_proximal.bed
bedtools intersect $parameters -u -a $gene_proximal_vntrs -b $bed_files/disease_vntrs_in_400k.bed\
    > $bed_files/gene_proximal_in_disease_and_400k.bed
    print $bed_files/gene_proximal_in_disease_and_400k.bed

echo "      Disease, excluding 400k and disease only"
bedtools intersect $parameters -v -a $disease_vntrs -b $all_400k_vntrs \
    > $bed_files/disease_not_400k.bed
    print $bed_files/disease_not_400k.bed
bedtools intersect $parameters -v -a $bed_files/disease_not_400k.bed -b $gene_proximal_vntrs \
    > $bed_files/disease_vntrs_only.bed
    print $bed_files/disease_vntrs_only.bed

echo "      Gene proximal and 400k, and excluding disease"
bedtools intersect $parameters -u -a $gene_proximal_vntrs -b $all_400k_vntrs \
    > $bed_files/gene_proximal_and_400k.bed
    print $bed_files/gene_proximal_and_400k.bed
bedtools intersect $parameters -v -a $bed_files/gene_proximal_and_400k.bed -b $disease_vntrs \
    > $bed_files/gene_proximal_and_400k_no_disease.bed
    print $bed_files/gene_proximal_and_400k_no_disease.bed
bedtools intersect -f 0.9 -u -a $bed_files/gene_proximal_and_400k_no_disease.bed -b $short_gene_vntrs \
    > $bed_files/short_gene_proximal_and_400k_no_disease.bed
    print $bed_files/short_gene_proximal_and_400k_no_disease.bed
bedtools intersect -f 0.9 -u -a $bed_files/gene_proximal_and_400k_no_disease.bed -b $long_gene_vntrs \
    > $bed_files/long_gene_proximal_and_400k_no_disease.bed
    print $bed_files/long_gene_proximal_and_400k_no_disease.bed


echo "      Gene proximal excluding 400k and disease"
bedtools intersect $parameters -v -a $gene_proximal_vntrs -b $all_400k_vntrs \
    > $bed_files/gene_proximal_no_400k.bed
    print $bed_files/gene_proximal_no_400k.bed
bedtools intersect $parameters -v -a $bed_files/gene_proximal_no_400k.bed -b $disease_vntrs \
    > $bed_files/gene_proximal_no_400k_no_disease.bed
    print $bed_files/gene_proximal_no_400k_no_disease.bed
bedtools intersect -f 0.9 -u -a $bed_files/gene_proximal_no_400k_no_disease.bed -b $short_gene_vntrs \
    > $bed_files/short_gene_proximal_no_400k_no_disease.bed
    print $bed_files/short_gene_proximal_no_400k_no_disease.bed
bedtools intersect -f 0.9 -u -a $bed_files/gene_proximal_no_400k_no_disease.bed -b $long_gene_vntrs \
    > $bed_files/long_gene_proximal_no_400k_no_disease.bed
    print $bed_files/long_gene_proximal_no_400k_no_disease.bed
bedtools intersect $parameters -u -a $gene_proximal_vntrs -b $disease_vntrs \
    > $bed_files/gene_proximal_and_disease.bed
    print $bed_files/gene_proximal_and_disease.bed


echo "      Gene proximal excluding segdup etc regions"
bedtools intersect $parameters -u -a $gene_proximal_vntrs -b $bed_files/hg38_vntr_no_complicated_region.tab \
    > $bed_files/gene_proximal_no_comp_region.bed
    print $bed_files/gene_proximal_no_comp_region.bed
bedtools intersect $parameters -u -a $bed_files/gene_proximal_no_400k_no_disease.bed -b $bed_files/hg38_vntr_no_complicated_region.tab \
    > $bed_files/gene_proximal_no_400k_no_disease_no_comp_region.bed
    print $bed_files/gene_proximal_no_400k_no_disease_no_comp_region.bed
bedtools intersect $parameters -v -a $bed_files/gene_proximal_no_400k_no_disease.bed -b $bed_files/hg38_vntr_no_complicated_region.tab \
    > $bed_files/gene_proximal_no_400k_no_disease_with_comp_region.bed
    print $bed_files/gene_proximal_no_400k_no_disease_with_comp_region.bed
bedtools intersect $parameters -u -a $disease_vntrs -b $bed_files/hg38_vntr_no_complicated_region.tab \
    > $bed_files/disease_vntr_no_comp_region.bed
    print $bed_files/disease_vntr_no_comp_region.bed

echo "      400k only"
bedtools intersect $parameters -v -a $all_400k_vntrs -b $gene_proximal_vntrs \
    > $bed_files/all_400k_no_gene_proximal.bed
    print $bed_files/all_400k_no_gene_proximal.bed
bedtools intersect $parameters -v -a $bed_files/all_400k_no_gene_proximal.bed -b $disease_vntrs \
    > $bed_files/all_400k_only.bed
    print $bed_files/all_400k_only.bed
