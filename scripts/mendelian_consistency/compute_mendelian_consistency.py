import pandas as pd
from cyvcf2 import VCF
import numpy as np
import itertools
from collections import defaultdict
from matplotlib import pyplot as plt
import sys
import math
import pickle


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
def is_valid_vntr(motif, motif_complexity_threshold, load_motif_complexity_scores, motif_complexity_score_map):
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


def check_mc(sample_gt, mother_gt, father_gt):
    if sample_gt[0] in mother_gt[0:2] and sample_gt[1] in father_gt[0:2]:
        return True
    if sample_gt[1] in mother_gt[0:2] and sample_gt[0] in father_gt[0:2]:
        return True
    return False


def print_stats_map(stats_map, num_trios):
    num_total_vntrs = stats_map["num_vntrs"]
    print("Total {} VNTRs pre-filters on the script. {} VNTRs survive all filters".format(num_total_vntrs, stats_map["num_chosen_vntrs"]))
    for key in stats_map:
        if key == "num_vntrs":
            continue
        if 'vntr' in key:
            percent_vntrs_filtered_out = round(stats_map[key]/float(num_total_vntrs)*100, 0)
            print("{} vntrs {} ({}%)".format(stats_map[key], key, percent_vntrs_filtered_out))
        elif "call" in key:
            percent_calls_filtered_out = round(stats_map[key]/float(num_total_vntrs * num_trios)*100, 0)
            print("{} calls {} ({}%)".format(stats_map[key], key, percent_calls_filtered_out))


def return_mc_stats(mc_list, stats_map, num_trios, verbose=False):
    # Print number of VNTRs filtered at each step.
    num_total_vntrs = stats_map["num_vntrs"]
    if verbose:
        print_stats_map(stats_map, num_trios)
    # Print MC stats.
    num_mendelian_consistent = mc_list.count(True)
    num_mendelian_inconsistent = mc_list.count(False)
    percentage = num_mendelian_consistent/float(num_mendelian_consistent + num_mendelian_inconsistent) * 100
    percentage = round(percentage, 5)
    if verbose:
        print("num_mendelian_consistent: {} num_mendelian_consistent: {}, MC {}%".format(
               num_mendelian_consistent, num_mendelian_inconsistent, percentage))
    return percentage


def histogram(mc_per_variant, suffix):
    plt.clf()
    plt.hist(mc_per_variant, bins=50, cumulative=True)
    ax = plt.gca()
    ax.set_xlabel("MC percentage among trios")
    ax.set_ylabel("Count variants")
    plt.savefig("MC_per_variant_hist_{}.pdf".format(suffix), bbox_inches="tight", format="pdf")

def get_mc_percentage(mcs, verbose=False):
    num_mendelian_consistent = mcs.count(True)
    num_mendelian_inconsistent = mcs.count(False)
    percentage = num_mendelian_consistent/float(num_mendelian_consistent + num_mendelian_inconsistent) * 100
    percentage = round(percentage, 2)
    if verbose and (percentage < 50):
        print("vid {} num_MC {} num_MI {} percentage {} ".format(
              num_mendelian_consistent,
              num_mendelian_inconsistent,
              percentage))
    return percentage

def is_non_polymorphic(sample_genotype, mother_genotype, father_genotype):
    return int(sample_genotype[0]) == 0 and\
           int(sample_genotype[1]) == 0 and \
           int(mother_genotype[0]) == 0 and \
           int(mother_genotype[1]) == 0 and \
           int(father_genotype[0]) == 0 and \
           int(father_genotype[1]) == 0

def compute_mc(vcf_filename, pedigree_filename, load_motif_complexity_scores=True, verbose=False):
    # Getting pedigree information
    pedigree = pd.read_csv(pedigree_filename, delim_whitespace=True)
    trios = pedigree[(pedigree['FatherID'] != "0") & (pedigree['MotherID'] != "0")]
    child_in_trios = set(trios['SampleID'].to_list())
    mother_in_trios = set(trios['MotherID'].to_list())
    father_in_trios = set(trios['FatherID'].to_list())
    all_ids = (child_in_trios.union(mother_in_trios)).union(father_in_trios)
    trio_ids = []
    for index, row in pedigree.iterrows():
        trio_ids.append((row['SampleID'], row['MotherID'], row['FatherID']))

    print("Working with {} trios".format(len(trios)))


    # Load pre-computed motif complexity score for each motif.
    motif_complexity_score_map = {}
    vid_to_complexity_score = {}
    motif_complexity_score_filename = "data/motif_complexity_scores.bin"
    if load_motif_complexity_scores:
        with open(motif_complexity_score_filename, "rb") as motif_complexity_score_file:
            motif_complexity_score_map = pickle.load(motif_complexity_score_file)

    str_ids = []
    with open("../vntrs_to_discard/str_ids.txt", "r") as str_ids_file:
        lines = str_ids_file.readlines()
        str_ids = [int(line) for line in lines]

    for skip_str_like_vntrs in [True, False]:
        for discard_non_polymorphic_vntrs in [True, False]:
            print("discard_non_polymorphic_vntrs is {} skip_str_like_vntrs {}".format(
                   discard_non_polymorphic_vntrs, skip_str_like_vntrs))
            # These thresholds can vary based on the experiment
            # Threshold for the score on complexity of the motif
            vntr_filter_threshold=0.8
            # Threshold on the minimum number of spanning reads.
            # This is not necessary anymore as we do the filtering in adVNTR.
            #spanning_reads_threshold=5
            # Min number of trios should be present before we compute the MC/MI for the VNTR.
            # For example of there is only one trio for a VNTR, and the MC = 0.0,
            # it's not representative of the MC of VNTR. So we don't include that VNTR in overall MC.
            min_required_trios_for_mc = 5
            overall_chrs_stats_map = defaultdict(int)
            all_mcs_summary = []

            vcf = VCF(vcf_filename, samples = list(all_ids))
            samples = vcf.samples

            stats_map = defaultdict(int)
            all_mcs = []
            vids = set()
            mc_summary_per_row = []
            for variant in vcf:
                mc_in_row = []
                motif = variant.INFO.get('RU')
                vntr_id = variant.INFO.get("VID")
                position = "chr{}:{}".format(variant.CHROM, variant.start)
                if vntr_id is None:
                    print("VNTR ID is None ", variant)
                    continue
                if int(vntr_id) in str_ids:
                    # This is an STR and is not included in the analysis.
                    stats_map["num_strs"] += 1
                    continue
                # Skipping the STR-like VNTR filter because the filter is applied to the input
                if skip_str_like_vntrs:
                    is_valid, score = is_valid_vntr(motif,
                                   motif_complexity_threshold=vntr_filter_threshold,
                                   load_motif_complexity_scores=load_motif_complexity_scores,
                                   motif_complexity_score_map=motif_complexity_score_map)
                    vid_to_complexity_score[vntr_id] = score
                    if not is_valid:
                        # VNTR is filtered out due to being str-like.
                        stats_map["num_str_like_vntrs"] += 1
                        continue

                # VNTR is not discarded.
                stats_map["num_vntrs"] += 1
                num_no_call_for_vntr = 0
                for family in trio_ids:
                    sample_index = samples.index(family[0])
                    mother_index = samples.index(family[1])
                    father_index = samples.index(family[2])
                    fam_indices = [sample_index, mother_index, father_index]
                    sample_genotype = variant.genotypes[sample_index]
                    mother_genotype = variant.genotypes[mother_index]
                    father_genotype = variant.genotypes[father_index]
                    if sample_genotype[0] == -1 or \
                       mother_genotype[0] == -1 or \
                       father_genotype[0] == -1: # No call
                        stats_map["num_no_call"] += 1
                        num_no_call_for_vntr += 1
                        continue
                    if discard_non_polymorphic_vntrs and \
                        is_non_polymorphic(sample_genotype, mother_genotype, father_genotype):
                        stats_map["num_non_polymorphic"] += 1
                        continue
                    stats_map["num_chosen_call"] += 1
                    vids.add(vntr_id)
                    mc_val = check_mc(sample_genotype, mother_genotype, father_genotype)
                    if verbose and mc_val == 0:
                        print(vntr_id)
                        #print("vid {} genotype {} parents {} {}".format(
                        #       vntr_id, sample_genotype, mother_genotype, father_genotype))
                    all_mcs.append(mc_val)
                    mc_in_row.append(mc_val)
                if len(mc_in_row) >= min_required_trios_for_mc:
                    mc_percentage_per_vntr = get_mc_percentage(mc_in_row)
                    if mc_percentage_per_vntr < 60 and False:
                        print("MC {} for {} num calls {} motif {}".format(
                            mc_percentage_per_vntr,
                            vntr_id,
                            len(mc_in_row),
                            motif))
                    mc_summary_per_row.append(mc_percentage_per_vntr)
                    stats_map["num_vntrs_passing_min_required_trios_threshold"] += 1

            mc_percentage = return_mc_stats(all_mcs, stats_map, len(trios))
            #all_mcs_summary.append((mc_percentage, np.mean(mc_summary_per_row)))
            print("MC: mean_per row {}% median_per row {}% mean_all {}% and number of calls {} for {} vntrs".format(
                    round(np.mean(mc_summary_per_row), 5),
                    round(np.median(mc_summary_per_row), 5),
                    round(mc_percentage, 5),
                    len(all_mcs),
                    len(list(vids))))
            #print("len(all_mcs_summary): ", len(all_mcs_summary))
            print("len(mc_summary_per_row): ", len(mc_summary_per_row))

    # Save motif complexity scores.
    with open("data/motif_complexity_scores_genic.bin", "wb") as motif_complexity_file:
        pickle.dump(motif_complexity_score_map, motif_complexity_file)
    with open("data/vid_to_complexity_score_genic.bin", "wb") as vid_to_complexity_file:
        pickle.dump(vid_to_complexity_score, vid_to_complexity_file)
    # plot histogram for MCs
    #histogram(mc_summary_per_row, suffix="")
    #low_mc_variants = [mc for mc in mc_summary_per_row if mc < 90]
    #histogram(low_mc_variants, suffix="low_mc_lt_90")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python {} <path to merged vcf file>".format(sys.argv[0]))
        exit(1)
    vcf_filename = sys.argv[1]
    #pedigree_filename = "pedigree_giab.txt"
    pedigree_filename = "pedigree_hprc.txt"

    compute_mc(vcf_filename, pedigree_filename, verbose=False)
