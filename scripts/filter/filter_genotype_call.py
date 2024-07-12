import glob
import sys
from collections import Counter

def process_genotype_file(filename, targetted_vntrs, is_expanded=False):
    genotype_file = open(filename, "r")
    lines = genotype_file.readlines()
    genotyped_vntrs = {}
    genotyped_vntr_ids = set()
    current_vntr_id = None
    for line in lines:
        if not (line[0].isdigit() or line.startswith("None")):
            # There might be other lines including warnings or logs, skip those lines.
            continue
        if line.startswith("None") or ("/" in line):
            # It's a genotype line
            if current_vntr_id is None:
                print("Error in parsing genotype file. current VNTR id not set")
                return
            if current_vntr_id in targetted_vntrs:
                vntr_genotype = line.strip()
                if is_expanded:
                    vntr_genotype = vntr_genotype.split()[0]
                if not vntr_genotype.startswith("None"):
                    genotyped_vntrs[current_vntr_id] = vntr_genotype
                    genotyped_vntr_ids.add(current_vntr_id)
            current_vntr_id = None
        else:
            # It's a line with vntr_id
            current_vntr_id = line.strip()
    return genotyped_vntrs, genotyped_vntr_ids

def process_targetted_vntrs_file(targetted_vntrs_filename) :
	targetted_vntrs = []
	with open(targetted_vntrs_filename, "r") as targetted_vntrs_file:
		lines = targetted_vntrs_file.readlines()
		targetted_vntrs = [line.strip().split(' ')[0] for line in lines]
	return targetted_vntrs

def get_genotype_alignment_for_one_dataset(vntr_aln_files_pattern,
                                           original_genotypes,
                                           target_vntrs,
                                           frac_mismatch_in_flanking_regions_threshold=0.1,
                                           min_supporting_reads=2,
                                           min_left_flanking_len=4,
                                           min_right_flanking_len=4,
                                           verbose=False):
    processed_genotypes = {}
    processed_vntr_ids = set()
    aln_file_names = glob.glob(vntr_aln_files_pattern)
    for aln_file_name in aln_file_names:
        repeat_counts = []
        aln_file = open(aln_file_name, "r")
        vntr_id = None
        lines = aln_file.readlines()
        #assert(lines[0].startswith("#VID"))
        # Skip empty files
        if len(lines) == 0:
            continue
        # Skip VNTRs that are not targetted
        vntr_id = lines[0].split()[1].strip()
        if vntr_id not in target_vntrs:
            continue
        # Example: "#VID: 912322 chr9:124947625-124947664"
        idx = 1
        # Extract repeat counts for reads where flanking region is well matched.
        while idx < len(lines):
            if idx + 4 >= len(lines):
                # Each entry for one read should have 5 lines. If less than that is available, it's not a valid entry. Skip it.
                print("Skipping the read starting at line " + \
                       str(idx + 1) + \
                       " in file " + \
                       aln_file_name + \
                       " because the alignment information is incomplete.")
                break
            line = lines[idx]
            assert(line.startswith(">") and "_RC" in line)
            # Example: ">1_RC:6 SEQLEN:244 VID:912322 REFRC:2 REPEATS:6 SR"
            repeat_count = int(line.split()[0].split(":")[1])
            vntr_id = line.split()[2].split(":")[1]
            idx += 4 # skipping query sequence, alignment line and reference sequence
            line = lines[idx]
            assert(line.startswith("# Mismatch")) # This is the ending line for each read
            # Example of flanking region log: "# Mismatch in flanking regions: 29/127 0.23, L:0/93 0.00, R:29/34 0.85"
            frac_mismatch_in_left_flanking_region = float(line.split()[8].replace(",", ""))
            frac_mismatch_in_right_flanking_region = float(line.split()[10].replace(",", ""))
            len_left_flanking = int(line.split()[7].replace("L:", "").split("/")[1])
            len_right_flanking = int(line.split()[9].replace("R:", "").split("/")[1])
            if frac_mismatch_in_left_flanking_region <= frac_mismatch_in_flanking_regions_threshold and\
               frac_mismatch_in_right_flanking_region <= frac_mismatch_in_flanking_regions_threshold and\
               len_left_flanking >= min_left_flanking_len and len_right_flanking >= min_right_flanking_len:
                # This is a high quality read in terms of alignment in flanking regions.
                repeat_counts.append(repeat_count)
            idx += 1
        if vntr_id == None:
            print("No VNTR ID found for file " + aln_file_name + " . Skipping file")
            continue
        # Compare the filtered read genotypes with the original genotypes (extracted from adVNTR direct output).
        if not vntr_id in original_genotypes.keys():
            print("ERROR: vntr id {} not found in original genotypes".format(vntr_id))
            continue
        original_genotype_values = [int(allele) for allele in original_genotypes[vntr_id].split("/")]
        processed_genotype_values = []
        for value in original_genotype_values:
            if value not in repeat_counts:
                if verbose:
                    print("Value " + str(value) + \
                    " from adVNTR original genotype values not supported by any reads processed on the alignment file " + \
                    aln_file_name)
            elif Counter(repeat_counts).get(value) >= min_supporting_reads:
                processed_genotype_values.append(value)
        if len(processed_genotype_values) == 0:
            if verbose:
                print("Not enough supporting reads or not high quality flanking region alignment for genotype " + \
                original_genotypes[vntr_id] + \
                " for vntr " + vntr_id)
        else:
            if len(processed_genotype_values) == 1:
                # This is if one of the alleles is not well supported or doesn't have high quality flanking region alignment,
                # we change the genotype to be homozygote.
                processed_genotype_values = [processed_genotype_values[0], processed_genotype_values[0]]
            processed_genotypes[vntr_id] = "/".join([str(value) for value in processed_genotype_values])
            processed_vntr_ids.add(vntr_id)
            if len(set(processed_genotype_values)) == 1 and len(set(original_genotype_values)) > 1:
                if verbose:
                    print("Genotype updated for vntr " + vntr_id + \
                    " . Not enough supporting reads or not high quality flanking region alignment.")
    return processed_genotypes, processed_vntr_ids

def print_comparison(sample, original_genotypes, filtered_genotypes, skip_same_values):
    printed_str = ""
    vntr_ids = filtered_genotypes.keys() & original_genotypes.keys()
    for vntr_id in vntr_ids:
        if skip_same_values:
            if original_genotypes[vntr_id] == filtered_genotypes[vntr_id]:
                continue
        printed_str += "sample {} id {}: {} -> {}\n".format(
                       sample,
                       vntr_id,
                       original_genotypes[vntr_id],
                       filtered_genotypes[vntr_id])
    print(printed_str)

if __name__== "__main__" :

    if len(sys.argv) < 2:
        print("Usage: python {} sample_id".format(sys.argv[0]))
        exit(1)
    sample = sys.argv[1]
    #print("For sample id {}".format(sample))
    skip_filtere_original_same_values = True
    verbose = False

    if skip_filtere_original_same_values and verbose:
        print("Skipping the printout of vntrs if original and filtered are the same")

    # Disease associated short
    #targetted_vntrs_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/target_vntrs/disease_associated_short_vntrs.txt"
    #pairwise_aln_filename = '/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_pairwise_align_disease_short/log_{}_aligned_GRCh38_winnowmap.sorted.bam.log_vid_*.txt'.format(sample)
    #genotype_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_disease_short_w_advntr_1.5.0/output_per_sample_{}.txt".format(sample)

    # Disease associated long
    #targetted_vntrs_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/target_vntrs/disease_associated_long_vntrs.txt"
    #pairwise_aln_filename = '/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_pairwise_align_disease_long/log_{}_aligned_GRCh38_winnowmap.sorted.bam.log_vid_*.txt'.format(sample)
    #genotype_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_disease_long_w_advntr_1.5.0/output_per_sample_{}.txt".format(sample)

    # Disease associated all
    targetted_vntrs_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/target_vntrs/disease_associated_all_vntrs_v3.txt"
    pairwise_aln_filename = '/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_pairwise_align_disease_all_vntrs_v3/log_{}_aligned_GRCh38_winnowmap.sorted.bam.log_vid_*.txt'.format(sample)
    genotype_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_disease_all_vntrs_v3/output_per_sample_{}.txt".format(sample)
    # EnsembleTR validation set
    #targetted_vntrs_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/target_vntrs/ensemble_validation_set.txt"
    #pairwise_aln_filename = '/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_pairwise_align_ensemble_validation_set/log_{}_aligned_GRCh38_winnowmap.sorted.bam.log_vid_*.txt'.format(sample)
    #genotype_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_ensemble_validation_set/output_per_sample_{}.txt".format(sample)

    targetted_vntrs = process_targetted_vntrs_file(targetted_vntrs_filename)
    original_genotypes, original_vntr_ids = process_genotype_file(genotype_filename, targetted_vntrs)
    genotypes, vntr_ids = get_genotype_alignment_for_one_dataset(pairwise_aln_filename, original_genotypes, targetted_vntrs)
    if sorted(vntr_ids) != sorted(original_vntr_ids) and verbose:
        print("Warning! vntr ids do not match between pairwise_aln and genotype output." + \
              " only in genotype output: {}".format(original_vntr_ids - vntr_ids) + \
              " only in parwise_aln: {}".format(vntr_ids- original_vntr_ids))
    print_comparison(sample=sample,
                     original_genotypes=original_genotypes,
                     filtered_genotypes=genotypes,
                     skip_same_values=skip_filtere_original_same_values)
