import glob
import sys
from collections import Counter

def process_genotype_file(filename, targetted_vntrs, is_expanded=False):
    genotype_file = open(filename, "r")
    lines = genotype_file.readlines()
    genotyped_vntrs = {}
    genotyped_vntr_ids = set()
    for idx in range(0, len(lines), 2):
        vntr_id = lines[idx].strip()
        if vntr_id in targetted_vntrs:
            vntr_genotype = lines[idx+1].strip()
            if is_expanded:
                vntr_genotype = vntr_genotype.split()[0]
            if not vntr_genotype.startswith("None"):
                genotyped_vntrs[vntr_id] = vntr_genotype
                genotyped_vntr_ids.add(vntr_id)
    return genotyped_vntrs, genotyped_vntr_ids

def process_targetted_vntrs_file(targetted_vntrs_filename) :
	targetted_vntrs = []
	with open(targetted_vntrs_filename, "r") as targetted_vntrs_file:
		lines = targetted_vntrs_file.readlines()
		targetted_vntrs = [line.strip().split(' ')[0] for line in lines]
		print (targetted_vntrs)
	return targetted_vntrs

def get_genotype_alignment_for_one_dataset(vntr_aln_files_pattern,
                                           original_genotypes,
                                           target_vntrs,
                                           frac_mismatch_in_flanking_regions_threshold=0.1,
                                           min_supporting_reads=2,
                                           min_left_flanking_len=4,
                                           min_right_flanking_len=4):
    processed_genotypes = {}
    processed_vntr_ids = set()
    aln_file_names = glob.glob(vntr_aln_files_pattern)
    print("len(aln_file_names): ", len(aln_file_names))
    print("len(original_genotypes): ", len(original_genotypes))
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
                print("Value " + str(value) + \
                " from adVNTR original genotype values not supported by any reads processed on the alignment file " + \
                aln_file_name)
            elif Counter(repeat_counts).get(value) >= min_supporting_reads:
                processed_genotype_values.append(value)
        if len(processed_genotype_values) == 0:
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
                print("Genotype updated for vntr " + vntr_id + \
                " . Not enough supporting reads or not high quality flanking region alignment.")
    return processed_genotypes, processed_vntr_ids



if __name__== "__main__" :

	sample = sys.argv[1]
	apply_filter = sys.argv[2]

	short_read_genotype_filename = '/home/shockyoung/VNTR/analysis/short_reads/illumina_%s.out' %sample
	long_read_genotype_filename = '/home/shockyoung/VNTR/analysis/long_reads/output_sample_%s-merge_cells_1-4_presetCCS_N1_l2000.bam' %sample
	targetted_vntrs_filename = "./target_vntr_ids/target_vntr_ids_illumina_probed_short_vntrs.txt"

	targetted_vntrs = process_targetted_vntrs_file(targetted_vntrs_filename)

	short_read_genotypes, short_read_vntr_ids = process_genotype_file(short_read_genotype_filename,targetted_vntrs,is_expanded=False)
	long_read_genotypes, long_read_vntr_ids = process_genotype_file(long_read_genotype_filename,targetted_vntrs)

	if apply_filter == 'True' :
		folder = 'filtered'
		short_read_genotypes_original = short_read_genotypes
		long_read_genotypes_original = long_read_genotypes
		short_read_genotypes, short_read_vntr_ids = get_genotype_alignment_for_one_dataset('/home/shockyoung/VNTR/analysis/short_reads/vntr_alignment/%s/vntr_*.txt'%sample, short_read_genotypes, targetted_vntrs)
		long_read_genotypes, long_read_vntr_ids = get_genotype_alignment_for_one_dataset('/home/shockyoung/VNTR/analysis/long_reads/vntr_alignment/%s/vntr_*.txt'%sample, long_read_genotypes, targetted_vntrs)
	else :
		folder = 'unfiltered'
		short_read_genotypes_original = short_read_genotypes
		long_read_genotypes_original = long_read_genotypes

	vntr_ids_union = short_read_vntr_ids | long_read_vntr_ids
	vntr_ids_intersection = short_read_vntr_ids & long_read_vntr_ids
	vntr_ids_short_read_only = short_read_vntr_ids - long_read_vntr_ids
	vntr_ids_long_read_only = long_read_vntr_ids - short_read_vntr_ids

	consistent_calls = [int(v) for v in short_read_vntr_ids if v in long_read_vntr_ids and short_read_genotypes[v] == long_read_genotypes[v]]
	inconsistent_calls = [int(v) for v in short_read_vntr_ids if v in long_read_vntr_ids and short_read_genotypes[v] != long_read_genotypes[v]]
	inconsistent_calls.sort()

	report = open('comparison_%s/%s/%s_comparison_report.txt' %(sample,folder,sample), 'w')

	report.write('numer of VNTRs genotyped by Illumina\t' + str(len(short_read_vntr_ids)) + '\n')
	report.write('numer of VNTRs genotyped by Pacbio\t' + str(len(long_read_vntr_ids)) + '\n')
	report.write('numer of VNTRs Union\t' + str(len(vntr_ids_union)) + '\n')
	report.write('numer of VNTRs Illumina only\t' + str(len(vntr_ids_short_read_only)) + '\n')
	report.write('numer of VNTRs Pacbio only\t' + str(len(vntr_ids_long_read_only)) + '\n')
	report.write('numer of VNTRs Intersection\t' + str(len(vntr_ids_intersection)) + '\n')
	report.write('numer of consistent calls\t' + str(len(consistent_calls)) + '\n')
	report.write('numer of inconsistent calls\t' + str(len(inconsistent_calls)) + '\n')

	with open("comparison_%s/%s/%s_VNTRs_genotyped_by_illumina_or_pacbio.txt" %(sample,folder,sample), "w") as f:
		for vntr_id in sorted(vntr_ids_union):
			f.write(str(vntr_id) + "\n")
	f.close()

	with open("comparison_%s/%s/%s_VNTRs_genotyped_by_illumina_missed_by_pacbio.txt" %(sample,folder,sample), "w") as f:
		for vntr_genotype in [(v, short_read_genotypes[v]) for v in vntr_ids_short_read_only]:
			f.write("^" + str(vntr_genotype[0]) + " \n")
	f.close()

	with open("comparison_%s/%s/%s_VNTRs_genotyped_by_pacbio_missed_by_illumina.txt" %(sample,folder,sample), "w") as f:
		for vntr_genotype in [(v, long_read_genotypes[v]) for v in vntr_ids_long_read_only]:
			f.write("^" + str(vntr_genotype[0]) + " \n")
	f.close()

	print(inconsistent_calls)
	with open("comparison_%s/%s/%s_VNTRs_with_inconsistent_genotype_calls.txt" %(sample,folder,sample), "w") as f:
		for id in inconsistent_calls :
			f.write('\t'.join([str(id), short_read_genotypes_original[str(id)], long_read_genotypes_original[str(id)], short_read_genotypes[str(id)], long_read_genotypes[str(id)]])+'\n')
	f.close()
	report.close()
