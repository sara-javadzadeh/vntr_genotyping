# Analysis of Targeted Gene-proximal and Phenotype-associated VNTRs

This repository provides the tools to analyze and plot a selection of gene-proximal (G-VNTRs) as well as a selection of phenotype-associated (P-VNTRs). G-VNTRs are VNTRs with motif length of at least six basepairs, within the gene boundaries and 500bp upstream of the transcription start site including promoters and other regulatory regions. VNTRs within the segmental duplication regions, and LINE and SINE elements are excluded from G-VNTRs. P-VNTRs are extracted collected based on the literature and might include VNTRs within segdup regions.

# Requirements
[AdVNTR-1.5.0](https://github.com/mehrdadbakhtiari/adVNTR) should be installed. Part of the data files are stored in a compressed file (scripts/barchart/results.tar.gz) and should be decompressed before running the corresponding script.
```
cd scripts/barchart
tar -xzvf results.tar.gz
cd ../..
```

# Genotyping
To make the repository lightweight, the genotyping log files are excluded. Instead, the dataframes created based on the genotype log files are stored and will be loaded by default. In order to run on a new dataset, run AdVNTR genotype as follows:
```
AdVNTR --alignment_file $alignment_file \
        --working_directory $working_directory\
        -m $database\
        --accuracy_filter \
        --disable_logging \
> $output
```
Where `alignment_file` is the input BAM file. For input BAM files including HiFi reads, use `--pacbio`. `working_directory` is the path to a directory where all the log files will be stored. `database` files are availble in `target_vntrs` directory. 


# Publication
The manuscript is under review. We will soon include a link to medRxiv and published manuscript (once available).
