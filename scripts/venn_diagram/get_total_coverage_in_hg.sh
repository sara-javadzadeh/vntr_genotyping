#!/bin/bash

disease_vntrs="/nucleus/projects/saraj/vntr/sources/vntr_database_comparison/disease_associated_sorted.bed"
gene_proximal_vntrs="gene_proximal_vntrs.bed"
disease_union_gene_proximal="disease_union_gene_proximal.bed"
cat $disease_vntrs $gene_proximal_vntrs > combined_disease_gene.bed
bedtools sort -i combined_disease_gene.bed > combined_disease_gene_sorted.bed
bedtools merge -i combined_disease_gene_sorted.bed > $disease_union_gene_proximal

bp_covered=$(cat $gene_proximal_vntrs | awk '{print($3-$2)}' | awk '{sum+=$0} END {print(sum)}')
echo "gene proximal vntrs covers $bp_covered bps"
bp_covered=$(cat $gene_proximal_vntrs | bedtools sort -i stdin | bedtools merge -i stdin | awk '{print($3-$2)}' | awk '{sum+=$0} END {print(sum)}')
echo "gene proximal vntrs (merged) covers $bp_covered bps"
bp_covered=$(cat $disease_vntrs | awk '{print($3-$2)}' | awk '{sum+=$0} END {print(sum)}')
echo "disease vntrs covers $bp_covered bps"
bp_covered=$(cat $disease_vntrs | bedtools sort -i stdin | bedtools merge -i stdin | awk '{print($3-$2)}' | awk '{sum+=$0} END {print(sum)}')
echo "disease vntrs covers (merged) $bp_covered bps"
bp_covered=$(cat $disease_union_gene_proximal | awk '{print($3-$2)}' | awk '{sum+=$0} END {print(sum)}')
echo "union (merged) covers $bp_covered bps"

