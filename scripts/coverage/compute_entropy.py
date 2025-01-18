# Need to run conda activate bioinfo
import numpy as np
import argparse
from scipy.stats import differential_entropy, entropy
import random


def compute_entropy(samples, spanning_reads, min_sr, max_sr, is_continuous=True):
    selected_srs = []
    for sample in samples:
        selected_srs.append(spanning_reads[sample])
    if is_continuous:
        return differential_entropy(selected_srs, base=2, method="auto")
    else:
        selected_srs = [int(sr) for sr in selected_srs]
        value, counts = np.unique(selected_srs, return_counts=True)
        return entropy(counts, base=2)

def read_samples(filename):
    ids = []
    with open(filename) as samples_ids_file:
        for line in samples_ids_file.readlines():
            ids.append(line.strip())
    return ids

def compute_entropy_for_vntr(words, max_sr, targeted_samples_ids_file,
                             is_wgs,
                             downsample=True, same_scale=True):
    spanning_reads = {}
    max_scale_thr = 50
    targeted_samples = read_samples(targeted_samples_ids_file)
    for word in words:
        if is_wgs and word.startswith("HG"):
            sample_name = word
        elif not is_wgs and word.strip() in targeted_samples:
            sample_name = word
        else: # It's a spanning read number
            spanning_reads[sample_name] = float(word)
    if same_scale and max(spanning_reads.values()) > max_scale_thr:
        return None, None
    samples = spanning_reads.keys()

    entropy = compute_entropy(samples, spanning_reads, max_sr=max_sr, min_sr=0)
    if downsample:
        random.seed(0)
        num_random_experiments = 10
        downsample_size = 7
        ds_entropy = []
        for idx in range(num_random_experiments):
            subset = random.sample(samples, downsample_size)
            ds_entropy.append(compute_entropy(subset, spanning_reads, max_sr=max_sr, min_sr=0))
        return entropy, np.median(ds_entropy)
    else:
        return entropy, entropy

def compute_entropy_for_dataset(is_wgs, same_scale, targeted_samples_ids_file):
    if is_wgs:
        # WGS
        sr_file = open("data/sr_for_entropy_wgs.txt", "r")
        downsample = True
        max_sr = 50
    else:
        # Targeted
        sr_file = open("data/sr_for_entropy_targeted.txt", "r")
        downsample = False
        max_sr = 250

    entropies = []
    down_sample_entropies = []
    for line in sr_file.readlines():
        words = line.split()
        entropy, down_sample_entropy = compute_entropy_for_vntr(words[1:], downsample,
                                            targeted_samples_ids_file, is_wgs,
                                            max_sr, same_scale=same_scale)
        if entropy and down_sample_entropy:
            entropies.append(entropy)
            down_sample_entropies.append(down_sample_entropy)
    print("Num VNTRs ", len(entropies))
    print("Median entropy across different VNTRs: {}".format(np.median(entropies)))
    print("Median down_sample_entropy across different VNTRs: {}".format(np.median(down_sample_entropies)))
    sr_file.close()


def parse_arguments():
    parser = argparse.ArgumentParser(
                prog="compute_entropy",
                description="Compute entropy of the spanning reads.")

    parser.add_argument("--targeted-samples-id-file", help="Path to a txt file with sample ids in one line.",
                        type=str,
                        default="data/targeted_samples_ids.txt")
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    for same_scale in [False, True]:
        print("Same scale of max 50 SR: ", same_scale)
        print("For WGS")
        compute_entropy_for_dataset(targeted_samples_ids_file=args.targeted_samples_id_file,
                                    is_wgs=True, same_scale=same_scale)

if __name__ == "__main__":
    main()
