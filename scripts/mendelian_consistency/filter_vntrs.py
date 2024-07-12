from cyvcf2 import VCF, Writer
import numpy as np
import sys
import copy
import math
import itertools
import pickle
from pprint import pprint


# Function used for filtering VNTRs
def get_normalized_sequence_similarity(a, b):
    # Compute Hamming distance between two short sequences.
    sequences_len = len(a)
    if len(a) != len(b):
        print("Warning: Computing sequence similarity between two sequences " +\
              "of different lengths: {} and {}".format(len(a), len(b)))
        sequences_len = min(len(a), len(b))
    similarity_score = 0
    for idx in range(sequences_len):
        # "M" character represent the masked character that matches
        # any other character for the purpose of this filter.
        if a[idx] == b[idx] or a[idx] == "M" or b[idx] == "M":
            similarity_score += 1
    return similarity_score / sequences_len

# Function used for filtering VNTRs
def get_motif_complexity_score(motif):
    self_match_score = 0
    for idx_in_mask in range(len(motif)):
        # Masking a single character at a time and computing the
        # max similarity score among all masked motifs.
        # Masked character matches any other character.
        masked_motif = motif[:idx_in_mask] + "M" + motif[idx_in_mask + 1:]
        # Creating a rolling window to compare the masked motif with itself.
        motif_window = masked_motif + masked_motif
        for idx in range(1, len(masked_motif)):
            end_idx_in_window = idx + len(masked_motif)
            # Compute the max score among all possible positions
            # when sliding the motif within the motif window to compare.
            self_match_score = max(self_match_score,
                                   get_normalized_sequence_similarity(masked_motif,
                                                         motif_window[idx:end_idx_in_window]))

    return self_match_score


# Function used for filtering VNTRs with variable number of masked characters based on motif length
def get_motif_complexity_score_flexible_masking(motif):
    self_match_score = 0
    if len(motif) > 40:
        # Too long motif. These computations will take along time.
        # Compute a single character masking instead.
        #print("Too long motif. Skipping computing the score for {}".format(motif))
        score = get_motif_complexity_score(motif)
        print("motif map for > 40 for motif {} len {} score {}".format(motif, len(motif), score))
        return score
    #TODO: temporarily return 0
    return 0
    num_masked_positions = math.ceil(len(motif)/10.0)
    all_combinations = itertools.combinations(range(len(motif)), num_masked_positions)
    #c_idx = 0
    for combination in all_combinations:
        #if c_idx % 100 == 0:
        #    print("in combination idx ", c_idx)
        #c_idx += 1
        masked_motif = motif
        #print(combination, num_masked_positions)
        for idx in combination:
            masked_motif = masked_motif[:int(idx)] + "M" + masked_motif[int(idx)+1:]
        #print(masked_motif)
        motif_window = masked_motif + masked_motif
        for idx in range(1, len(masked_motif)):
            end_idx_in_window = idx + len(masked_motif)
            # Compute the max score among all possible positions
            # when sliding the motif within the motif window to compare.
            self_match_score = max(self_match_score,
                                   get_normalized_sequence_similarity(masked_motif,
                                                         motif_window[idx:end_idx_in_window]))
            if self_match_score == 1.0:
                # This is the maximum possible score. No need to search for other combinations.
                #print("motif length {} num masked positions {} score {} {}".format(
                #       len(motif), num_masked_positions, round(self_match_score, 2), motif))
                return self_match_score

    #print("motif length {} num masked positions {} score {} {}".format(
    #       len(motif), num_masked_positions, round(self_match_score, 2), motif))
    return self_match_score


# Function used for filtering VNTRs
def isValidVNTR(motif, motif_complexity_threshold, load_motif_complexity_scores, motif_complexity_score_map):
    if load_motif_complexity_scores and motif in motif_complexity_score_map:
        score = motif_complexity_score_map[motif]
    elif motif not in motif_complexity_score_map:
        #print("Motif {} not in motif_complexity_score_map".format(motif))
        score = get_motif_complexity_score_flexible_masking(motif)
        motif_complexity_score_map[motif] = score
    else:
        score = get_motif_complexity_score_flexible_masking(motif)
        motif_complexity_score_map[motif] = score
    if score > motif_complexity_threshold:
        # VNTR is very much STR like. Skip this VNTR.
        #print(motif)
        return False, score
    return True, score

def load_spanning_reads(chromosome):
    #spanning_reads_filename = "1000Genomes-TR-Analysis/mend/intermediate_data/spanning_reads_{}.bin".format(chromosome)
    spanning_reads_vid_filename = "1000Genomes-TR-Analysis/mend/intermediate_data/spanning_reads_vid_{}.bin".format(chromosome)
    spanning_reads_vid_map = {}
    #with open(spanning_reads_filename, "rb") as spanning_reads_file:
    #    spanning_reads_map = pickle.load(spanning_reads_file)
    with open(spanning_reads_vid_filename, "rb") as spanning_reads_vid_file:
        spanning_reads_vid_map = pickle.load(spanning_reads_vid_file)
    return spanning_reads_vid_map

def bgzip_output_file(input_filename):
    import bgzip
    output_filename = input_filename + ".gz"
    with open(input_filename, 'r') as f_in, bgzip.BGZipWriter(output_filename, 'wb') as f_out:
        f_out.write(f_in.encode())

def gzip_output_file(input_filename):
    import gzip
    output_filename = input_filename + ".gz"
    with open(input_filename, 'rb') as f_in, gzip.open(output_filename, 'wb') as f_out:
        f_out.writelines(f_in)

def get_sample_ids():
    pass

def create_filtered_vcf(chromosome, write_to_output=False, load_motif_complexity_map=False):
    input_vcf_filename = "/nucleus/projects/saraj/vntr/sources/EnsembleTR/filtered_advntr_calls_no_motif_in_ref_filter/ensemble_input_{}_advntr.vcf.gz".format(chromosome)
    # This is for testing input
    #input_vcf_filename = "/nucleus/projects/saraj/vntr/sources/EnsembleTR/filtered_advntr_call_plus_str_like_filter/all_pop_filtered_{}.vcf.gz".format(chromosome)
    output_vcf_filename = "/nucleus/projects/saraj/vntr/sources/EnsembleTR/advntr_call_str_like_flexible_filter_sr_5/all_pop_filtered_plus_str_like_flexible_filter_spanning_reads_{}_.vcf".format(chromosome)
    motif_complexity_score_map = {}
    # Load pre-computed motif complexity score for each motif.
    motif_complexity_score_filename = "1000Genomes-TR-Analysis/mend/intermediate_data/motif_complexity_scores_flexible_mask_ceil_div_10_skip_40_motif.bin"
    if load_motif_complexity_map:
        with open(motif_complexity_score_filename, "rb") as motif_complexity_score_file:
            motif_complexity_score_map = pickle.load(motif_complexity_score_file)
    motif_complexity_threshold = 0.8
    # Load pre-computed number of spanning reads for each call.
    spanning_reads_vid_map = load_spanning_reads(chromosome)
    #print(list(spanning_reads_vid_map.keys())[:10])
    input_vcf = VCF(input_vcf_filename)
    samples = input_vcf.samples
    if write_to_output:
        output_vcf = Writer(output_vcf_filename, input_vcf)
        output_vcf.add_info_to_header({
                    "ID": "MotifComplexity",
                    "Description": "motif complexity score in [0, 1]. Higher value indicate more STR-like.",
                    "Type": "Float",
                    "Number": "1"})
    input_vcf.add_info_to_header({
                "ID": "MotifComplexity",
                "Description": "motif complexity score in [0, 1]. Higher value indicate more STR-like.",
                "Type": "Float",
                "Number": "1"})
    if write_to_output:
        output_vcf.add_format_to_header(
                    {"ID": "FilteredOut",
                    "Description": "Reason for being filtered out, if applicable. -1 if filtered out because of not enough spanning reads. 0 if not filtered out.",
                    "Type": "Integer",
                    "Number": "1"}
        )
    input_vcf.add_format_to_header(
                {"ID": "FilteredOut",
                "Description": "Reason for being filtered out, if applicable. -1 if filtered out because of not enough spanning reads. 0 if not filtered out.",
                "Type": "Integer",
                "Number": "1"}
    )
    num_surviving_vntr = 0
    num_surviving_calls = 0
    spanning_reads_threshold = 5
    no_spanning_read_info = set()
    for row in input_vcf:
        motif = row.INFO["RU"]
        updated_row_sr = copy.deepcopy(row.format("SR"))
        #updated_row_gt = copy.deepcopy(row.format("GT"))
        #updated_row_dp = copy.deepcopy(row.format("DP"))
        #updated_row_fr = copy.deepcopy(row.format("FR"))
        #updated_row_ml = copy.deepcopy(row.format("ML"))
        updated_row_filter = np.zeros(shape=np.shape(updated_row_sr), dtype=int)
        #pprint(dir(row))
        #print(type(row.format("SR")))
        #print(type(row.genotype()))
        #print(type(row.set_format()))
        is_valid, score = isValidVNTR(motif, motif_complexity_threshold, True, motif_complexity_score_map)
        if is_valid:
            # Valid VNTR passing the motif complexity filter
            # Write this row to output_vcf
            #print(motif, motif_complexity_score, " passing filter")
            row.INFO["MotifComplexity"] = round(score, 2)
            vid = row.INFO["VID"]
            for idx, sr_in_call in enumerate(row.format("SR")):
                sample_id = samples[idx]
                if sr_in_call == -1:
                    # pull the information from the pre-loaded map
                    call_pair = (vid, sample_id)
                    if call_pair not in spanning_reads_vid_map:
                        no_spanning_read_info.add(call_pair)
                    else:
                        row.format("SR")[idx] = spanning_reads_vid_map[call_pair]
                        updated_row_sr[idx] = spanning_reads_vid_map[call_pair]
                else:
                    # spanning read information is there
                    pass
                if updated_row_sr[idx] < spanning_reads_threshold:
                    # Print a no call if failing spanning reads threshold threshold
                    #GT:DP:SR:FR:ML:FilteredOut
                    #updated_row_gt[idx] = -1
                    #updated_row_dp[idx] = -1
                    #updated_row_sr[idx] = -1
                    #updated_row_fr[idx] = -1
                    #updated_row_ml[idx] = -1
                    # Setting FilteredOut format field below

                    if updated_row_sr[idx] == -1:
                        print("Error! Spanning read info not available for vntr {} sample idx {}".format(vid, idx))
                    updated_row_filter[idx] = -1
                else:
                    num_surviving_calls += 1

            num_surviving_vntr += 1
            #print(updated_row_sr[:10])
            #print(np.shape(updated_row_sr))
            # This is for writing filteredOut and SR without a no call.
            updated_row_sr = updated_row_sr.astype(int)
            row.set_format("SR", updated_row_sr)
            #row.set_format("GT", updated_row_gt)
            #row.set_format("DP", updated_row_dp)
            #row.set_format("FR", updated_row_fr)
            #row.set_format("ML", updated_row_ml)
            row.set_format("FilteredOut", updated_row_filter)
            if write_to_output:
                output_vcf.write_record(row)
    print("len call pairs not in spanning reads map {} num distinct vntrs {} num samples starting with HG or NA {}".format(
        len(no_spanning_read_info),
        len(set([pair[0] for pair in no_spanning_read_info])),
        len(set([pair[0] for pair in no_spanning_read_info if (pair[1].startswith("HG") or pair[1].startswith("NA"))]))))
    input_vcf.close()
    if write_to_output:
        output_vcf.close()
    # bgzip and index in a bash script
    #bgzip_output_file(input_filename=output_vcf_filename)
    print("For {} {} vntrs and {} calls survived for motif complexity threshold {}".format(
           chromosome, num_surviving_vntr, num_surviving_calls, motif_complexity_threshold))
    return num_surviving_vntr, num_surviving_calls

def creat_no_call_vcf(chromosome, write_to_output=False):
    input_vcf_filename = "/nucleus/projects/saraj/vntr/sources/EnsembleTR/advntr_call_str_like_flexible_filter_sr_10/all_pop_filtered_plus_str_like_flexible_filter_spanning_reads_{}.vcf".format(chromosome)
    output_vcf_filename = "/nucleus/projects/saraj/vntr/sources/EnsembleTR/advntr_call_str_like_flexible_filter_sr_10_no_call/all_pop_filtered_plus_str_like_flexible_filter_spanning_reads_{}_.vcf".format(chromosome)
    input_file = open(input_vcf_filename, "r")
    if write_to_output:
        output_file = open(output_vcf_filename, "w+")
    for line in input_file.readlines():
        output_line = ""
        if line.startswith("#"):
            # Copy the header lines as is.
            output_line = line
        else:
            # For other lines, check out the filteredOut format field and change the call to no-call if set to -1.
            words = line.split()
            for idx, word in enumerate(words):
                if ":" in word and (not word.startswith("GT")):
                    # Make sure it's a call and not other info
                    filtered_out = word.split(":")[5]
                    if int(filtered_out) == -1:
                        output_line += "\t.:.:.:.:-1"
                    else:
                        output_line += "\t" + word
                elif idx == 0:
                    output_line += word
                else:
                    output_line += "\t" + word
                if idx == len(words) - 1:
                    # Add a new line after the last call
                    output_line += "\n"
        if write_to_output:
            output_file.write(output_line)
    input_file.close()
    if write_to_output:
        output_file.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python {} <chromosome e.g.: chr1>".format(sys.argv[0]))
        exit(1)
    if sys.argv[1] == "all":
        num_surviving_vntrs_all_chr = 0
        num_surviving_calls_all_chr = 0
        for chromosome in range(1, 23):
            num_surviving_vntrs, num_surviving_calls = create_filtered_vcf(chromosome="chr"+str(chromosome))
            num_surviving_vntrs_all_chr += num_surviving_vntrs
            num_surviving_calls_all_chr += num_surviving_vntrs
            #creat_no_call_vcf(chromosome="chr"+str(chromosome))
        print("For all chromosomes {} vntrs and {} calls survived".format(num_surviving_vntrs_all_chr, num_surviving_calls_all_chr))
    else:
        #create_filtered_vcf(chromosome=sys.argv[1])
        #creat_no_call_vcf(chromosome=sys.argv[1])
        pass

if __name__ == "__main__":
    main()
