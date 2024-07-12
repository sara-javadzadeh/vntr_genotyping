import glob
import pickle
import sys
from collections import Counter
#from cyvcf2 import VCF

def process_vcf_genotype_file(filename, targetted_vntrs, samples):

    vcf = VCF(filename, samples=samples)
    genotyped_vntrs = {}
    genotyped_vntr_ids = set()
    for variant in vcf:
        vntr_id = variant.INFO.get("VID")
        # Compute the repeat count based on reference and alt repeat counts.
        # Note that in the VCF file, the genotypes does not represent the repeat counts.
        # Instead, the genotypes represent which of the ref or alt alleles are present in the sample.
        # Compute the repeat count in the reference
        repeat_count = {}
        motif_len = len(variant.INFO.get('RU'))
        ref_allele_len = len(variant.REF)
        repeat_count[0] = round(ref_allele_len/motif_len)
        # Compute the repeat count in the alternative alleles
        if len(variant.ALT) > 0: # No alternative alleles detected here
            for idx, alt_allele in enumerate(variant.ALT):
                alt_allele_len = len(alt_allele)
                # Alt alleles index in the genotype call start from 1.
                repeat_count[idx+1] = round(alt_allele_len / motif_len)
        print("vid {} motif len {} repeat counts {} motif {} ref allele {} alt alleles {}".format(
               vntr_id, motif_len, repeat_count, variant.INFO.get('RU'), variant.REF, variant.ALT))
        if vntr_id in targetted_vntrs:
            vcf_genotype = variant.genotypes[0]
            # example genotype: [0, 0, True]
            if vcf_genotype[0] == -1 or vcf_genotype[1] == -1:
                # None genotype indicated as value -1 or ./. in vcf file.
                # Not storing vntrs with None genotype value
                pass
            else:
                repeat_counts = sorted([repeat_count[vcf_genotype[0]], repeat_count[vcf_genotype[1]]])
                genotype = "/".join([str(repeat_counts[0]), str(repeat_counts[1])])
                print("vid {} vcf genotype {} and genotype {}".format(vntr_id, vcf_genotype, genotype))
                genotyped_vntr_ids.add(vntr_id)
                genotyped_vntrs[vntr_id] = genotype
    return genotyped_vntrs, genotyped_vntr_ids

def process_genotype_file(filename, targetted_vntrs, skip_low_complexity_vntrs):
    str_ids = []
    with open("../vntrs_to_discard/str_ids.txt", "r") as str_ids_file:
        lines = str_ids_file.readlines()
        str_ids = [int(line) for line in lines]

    genotype_file = open(filename, "r")
    lines = genotype_file.readlines()
    genotyped_vntrs = {}
    n_str_like_vntrs = 0
    genotyped_vntr_ids = set()
    vid_to_complexity_score_map = {}
    if skip_low_complexity_vntrs:
        with open("data/vid_to_complexity_score.bin", "rb") as complexity_score_file:
            vid_to_complexity_score_map = pickle.load(complexity_score_file)
    for idx in range(0, len(lines)):
        # Skip header lines or log lines.
        if (not lines[idx][0].isnumeric()) and (not lines[idx].startswith("None")):
            continue
        if lines[idx].startswith("None") or "/" in lines[idx]:
            # This is a genotype line
            vntr_genotype = lines[idx].strip()
            if vntr_id in targetted_vntrs and \
               int(vntr_id) not in str_ids: # We have some STRs in the genotyping database. Filtering them out here.
                # Skip str-like vntrs if indicated
                if skip_low_complexity_vntrs:
                    is_valid, score = is_valid_vntr(vntr_id, 0.8, vid_to_complexity_score_map)
                    if not is_valid:
                        n_str_like_vntrs += 1
                        continue
                if not vntr_genotype.startswith("None"):
                    genotyped_vntrs[vntr_id] = vntr_genotype
                    genotyped_vntr_ids.add(vntr_id)
        else:
            # This is a vntr id line
            vntr_id = int(lines[idx].strip())
    print("n_str_like_vntrs: ", n_str_like_vntrs)
    return genotyped_vntrs, genotyped_vntr_ids

def process_targetted_vntrs_file(targetted_vntrs_filename) :
    targetted_vntrs = []
    with open(targetted_vntrs_filename, "r") as targetted_vntrs_file:
        lines = targetted_vntrs_file.readlines()
        targetted_vntrs = [int(line.strip()) for line in lines]
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


def get_consistency(short_read_vntr_ids, long_read_vntr_ids, short_read_genotypes, long_read_genotypes):
    shared_vntr_ids = short_read_vntr_ids & long_read_vntr_ids
    num_both_h_call, num_hifi_h_call, num_illumina_h_call = 0, 0, 0
    consistent_calls, inconsistent_calls, partially_consistent_calls = [], [], []
    for vid in shared_vntr_ids:
        if short_read_genotypes[vid] == long_read_genotypes[vid]:
            consistent_calls.append(vid)
            continue
        elif short_read_genotypes[vid] == "None":
            # only in long reads
            continue
        elif long_read_genotypes[vid] == "None":
            # only in long reads
            continue
        short_read_alleles = set(short_read_genotypes[vid].split("/"))
        long_read_alleles = set(long_read_genotypes[vid].split("/"))
        if len(short_read_alleles & long_read_alleles) == 0:
            if len(long_read_alleles) == 1 and len(short_read_alleles) == 1:
                # E.g. 3/3 and 4/4
                #partially_consistent_calls.append(vid)
                inconsistent_calls.append(vid)
                num_both_h_call += 1
            else:
                # E.g. 3/3 vs 4/2
                # no alleles in common
                inconsistent_calls.append(vid)
        elif len(short_read_alleles & long_read_alleles) == 1:
            # one allele is shared
            # check if one genotype call is homozygote and the other one heterozygote
            if len(long_read_alleles) == 1 or len(short_read_alleles) == 1:
                # E.g. 2/2 vs 2/3
                if len(short_read_alleles) == 1:
                    num_illumina_h_call += 1
                elif len(long_read_alleles) == 1:
                    num_hifi_h_call += 1
                partially_consistent_calls.append(vid)
            else:
                # E.g. 2/4 vs 2/3
                inconsistent_calls.append(vid)
    return consistent_calls, inconsistent_calls, partially_consistent_calls,\
           num_both_h_call, num_hifi_h_call, num_illumina_h_call

def is_valid_vntr(vntr_id, motif_complexity_threshold, vid_to_complexity_score_map):
    if vntr_id in vid_to_complexity_score_map:
        score = vid_to_complexity_score_map[vntr_id]
    else:
        print("cannot calculate motif complexity score for vid", vntr_id)
        return True, None
    if score > motif_complexity_threshold:
        # VNTR is very much STR like. Skip this VNTR.
        return False, score
    return True, score

def compare_genotypes(sample, output_directory, skip_low_complexity_vntrs):
    prefix = ""
    if skip_low_complexity_vntrs:
        prefix = "_no_str_like_vntr" #logs_illumina_wgs_5k_short_sr_3_txt
    short_read_genotype_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/scripts/mendelian_consistency/" + \
                                   "logs_illumina_wgs_5k_short_sr_3_txt/%s/%s_short_vntr_dataset.txt" % (sample, sample)
    long_read_genotype_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/logs/logs_genotype_10k_short/" + \
                                  "output_per_sample_%s.txt" % sample
    targetted_vntrs_filename = "/nucleus/projects/saraj/vntr/sources/advntr/cloned_advntr_aug_30_merge_pipelines/target_vids_5k_short.txt"
    try:
        targetted_vntrs = process_targetted_vntrs_file(targetted_vntrs_filename)
    except Exception:
        return
    short_read_genotypes, short_read_vntr_ids = process_genotype_file(
        short_read_genotype_filename,
        targetted_vntrs,
        skip_low_complexity_vntrs)
    long_read_genotypes, long_read_vntr_ids = process_genotype_file(
        long_read_genotype_filename,
        targetted_vntrs,
        skip_low_complexity_vntrs)
    #print("len short_read_vntr_ids", len(short_read_vntr_ids))
    #print("head short_read_vntr_ids", sorted(list(short_read_vntr_ids))[:10])
    #print("len short_read_vntr_genotypes", len(short_read_genotypes))
    #print("head short_read_vntr_genotypes", sorted(list(short_read_genotypes))[:10])
    #print("len long_read_vntr_ids", len(long_read_vntr_ids))
    short_read_genotypes_original = short_read_genotypes
    long_read_genotypes_original = long_read_genotypes

    working_directory = "{}/comparison_{}".format(output_directory, sample)

    vntr_ids_union = short_read_vntr_ids | long_read_vntr_ids
    vntr_ids_intersection = short_read_vntr_ids & long_read_vntr_ids
    vntr_ids_short_read_only = short_read_vntr_ids - long_read_vntr_ids
    vntr_ids_long_read_only = long_read_vntr_ids - short_read_vntr_ids

    consistent_calls, inconsistent_calls, partially_consistent_calls, \
           num_both_h_call, num_hifi_h_call, num_illumina_h_call \
        = get_consistency(
        short_read_vntr_ids, long_read_vntr_ids, short_read_genotypes, long_read_genotypes)
    inconsistent_calls.sort()

    #if len(vntr_ids_intersection) > 0:
    inconsistent_percentage = float(len(inconsistent_calls)) / len(vntr_ids_intersection)*100
    partially_consistent_percentage = float(len(partially_consistent_calls)) / len(vntr_ids_intersection) * 100
    #else:
    #    inconsistent_percentage = 0
    #    partially_consistent_percentage = 0
    report = open('%s/%s_comparison_report%s.txt' %(working_directory,sample, prefix), 'w')
    report.write('number of VNTRs genotyped by Illumina\t' + str(len(short_read_vntr_ids)) + '\n')
    report.write('number of VNTRs genotyped by Pacbio\t' + str(len(long_read_vntr_ids)) + '\n')
    report.write('number of VNTRs Union\t' + str(len(vntr_ids_union)) + '\n')
    report.write('number of VNTRs Illumina only\t' + str(len(vntr_ids_short_read_only)) + '\n')
    report.write('number of VNTRs Pacbio only\t' + str(len(vntr_ids_long_read_only)) + '\n')
    report.write('number of VNTRs Intersection\t' + str(len(vntr_ids_intersection)) + '\n')
    report.write('number of consistent calls\t' + str(len(consistent_calls)) + \
                    "\t({:.1f}% of intersection)\n".format(
                        float(len(consistent_calls)) / len(vntr_ids_intersection)*100))
    report.write('number of inconsistent calls\t' + str(len(inconsistent_calls)) + \
                    "\t({:.1f}% of intersection)\n".format(
                    inconsistent_percentage))
    report.write('number of partially consistent calls\t' + str(len(partially_consistent_calls)) + \
                    "\t({:.1f}% of intersection {:.1f}% illumina_h_call {:.1f}% hifi_h_call {:.1f}% both_h_call)\n".format(
                        partially_consistent_percentage,
                        num_illumina_h_call / len(partially_consistent_calls) *100.0,
                        num_hifi_h_call / len(partially_consistent_calls) *100.0,
                        num_both_h_call / len(partially_consistent_calls) *100.0))

    with open("%s/%s_VNTRs_genotyped_by_illumina_or_pacbio%s.txt" %(working_directory,sample, prefix), "w") as f:
        for vntr_id in sorted(vntr_ids_union):
            f.write(str(vntr_id) + "\n")
    f.close()

    with open("%s/%s_VNTRs_genotyped_by_illumina_missed_by_pacbio%s.txt" %(working_directory,sample, prefix), "w") as f:
        for vntr_genotype in [(v, short_read_genotypes[v]) for v in vntr_ids_short_read_only]:
            f.write(str(vntr_genotype[0]) + " \n")
    f.close()

    with open("%s/%s_VNTRs_genotyped_by_pacbio_missed_by_illumina%s.txt" %(working_directory,sample, prefix), "w") as f:
        for vntr_genotype in [(v, long_read_genotypes[v]) for v in vntr_ids_long_read_only]:
            f.write(str(vntr_genotype[0]) + " \n")
    f.close()

    print(inconsistent_calls)
    with open("%s/%s_VNTRs_with_inconsistent_genotype_calls%s.txt" %(working_directory,sample, prefix), "w") as f:
        f.write("vid\tillumina genotype\thifi genotype\n")
        for id in inconsistent_calls:
            f.write('\t'.join([str(id),
                short_read_genotypes[id],
                long_read_genotypes[id]])+'\n')
    f.close()
    with open("%s/%s_VNTRs_with_partially_consistent_genotype_calls.txt" %(working_directory,sample), "w") as f:
        for id in partially_consistent_calls:
            try:
                f.write('\t'.join([str(id),
                    short_read_genotypes[id],
                    long_read_genotypes[id]])+'\n')
            except Exception:
                print("Error while writing vid ", id)
    f.close()
    report.close()


if __name__== "__main__" :

    sample = sys.argv[1]
    output_directory = sys.argv[2]
    for skip_low_complexity_vntrs in [False, True]:
        compare_genotypes(sample, output_directory, skip_low_complexity_vntrs)
