import numpy as np
from matplotlib import pyplot as plt
from collections import Counter

def histplot(genotype_counter, suffix):
    plt.clf()
    genotype_counter_list = [(key, genotype_counter[key]) for key in genotype_counter.keys()]
    genotype_counter_list = sorted(genotype_counter_list, key=lambda x: x[1])
    x, y = zip(*genotype_counter_list)
    # Changing x to show an index instead of vid
    x = range(1, len(x)+1)
    plt.plot(x,
             y)
    ax = plt.gca()
    ax.set_xlabel("VNTR id")
    ax.tick_params(axis="x", rotation=90)
    ax.set_ylabel("Count")
    plt.savefig("{}_genotypes_across_samples.pdf".format(suffix),
                format="pdf",
                bbox_inches="tight")

def aggregate_comparison_reports():
    sample_names = ["HG00438", "HG00621", "HG00673", "HG00735", "HG00741", "HG01106", "HG01175", "HG01258", "HG01358",
                    "HG01361", "HG01891", "HG01928", "HG01952", "HG01978", "HG02148", "HG02257", "HG02572",
                    "HG02622", "HG02630", "HG02717", "HG02886", "HG03453", "HG03516", "HG03540", "HG03579"]
    # Skipping HG01123 HG02486 HG02559, no parent data
    # Get median consistency value over all samples
    skip_str_like_vntrs = True
    if skip_str_like_vntrs:
        filename_suffix = "_no_str_like_vntr"
    else:
        filename_suffix = ""
    n_consistent, percent_consistent = [], []
    n_partially_consistent, percent_partially_consistent = [], []
    percent_illumina_h_call, percent_hifi_h_call, percent_both_h_call = [], [], []
    n_inconsistent, percent_inconsistent = [], []
    n_genotyped_illumina, n_genotyped_hifi = [], []
    #input_dir = "v2_illumina_pacbio_comparison/"
    #input_dir = "v3_illumina_pacbio_comparison_no_str_sr_3/"
    input_dir = "v4_illumina_pacbio_comparison_no_str_sr_3_inc_def/"

    for sample_name in sample_names:
       input_filename = input_dir + \
                "comparison_{}/{}_comparison_report{}.txt".format(
                 sample_name, sample_name, filename_suffix)
       try:
           with open(input_filename, "r") as input_file:
                for line in input_file.readlines():
                    if line.startswith("number of VNTRs genotyped by Illumina"):
                        n_genotyped_illumina.append(int(line.split()[6]))
                    elif line.startswith("number of VNTRs genotyped by Pacbio"):
                        n_genotyped_hifi.append(int(line.split()[6]))
                    elif line.startswith("number of consistent calls"):
                        n_consistent.append(int(line.split()[4]))
                        percent_consistent.append(float(line.split()[5].replace("%", "").replace("(", "")))
                    elif line.startswith("number of inconsistent calls"):
                        n_inconsistent.append(int(line.split()[4]))
                        percent_inconsistent.append(float(line.split()[5].replace("%", "").replace("(", "")))
                    elif line.startswith("number of partially consistent calls"):
                        # number of partially consistent calls 34
                        # (2.0% of intersection 82.4% illumina_h_call 17.6% hifi_h_call 0.0% both_h_call)
                        n_partially_consistent.append(int(line.split()[5]))
                        percent_partially_consistent.append(float(line.split()[6].replace("%", "").replace("(", "")))
                        percent_illumina_h_call.append(float(line.split()[9].replace("%", "")))
                        percent_hifi_h_call.append(float(line.split()[11].replace("%", "")))
                        percent_both_h_call.append(float(line.split()[13].replace("%", "")))
       except Exception:
           pass
    median_consistent = np.median(percent_consistent)
    median_n_consistent = np.median(n_consistent)
    median_partially_consistent = np.median(percent_partially_consistent)
    median_n_partially_consistent = np.median(n_partially_consistent)
    median_inconsistent = np.median(percent_inconsistent)
    median_n_inconsistent = np.median(n_inconsistent)
    median_n_illumina = np.median(n_genotyped_illumina)
    median_n_hifi = np.median(n_genotyped_hifi)
    median_illumina_h_call = np.median(percent_illumina_h_call)
    median_hifi_h_call = np.median(percent_hifi_h_call)
    median_both_h_call = np.median(percent_both_h_call)
    print("Median num vntrs genotyped in illumina {} hifi {}".format(
            median_n_illumina, median_n_hifi))
    print("Median over WGS samples, {}% consistent, {}% partially consistent and {}% inconsistent".format(
            median_consistent, median_partially_consistent, median_inconsistent))
    print("Among partially consistent calls, {}% illumina_h_call, {}% hifi_h_call and {}% both_h_call".format(
            median_illumina_h_call, median_hifi_h_call, median_both_h_call))
    print("Median over WGS samples, {} n_consistent, {} n_partially consistent and {} n_inconsistent".format(
            median_n_consistent, median_n_partially_consistent, median_n_inconsistent))

    # Get VNTRs with inconsistency
    for inconsistency_type in ["inconsistent", "partially_consistent"]:
        vids_across_samples = []
        for sample_name in sample_names:
           input_filename = input_dir + \
                    "comparison_{}/{}_VNTRs_with_{}_genotype_calls.txt".format(
                     sample_name, sample_name, inconsistency_type)
           try:
               with open(input_filename, "r") as input_file:
                   for line in input_file.readlines():
                        vntr_id = line.split()[0]
                        if vntr_id == "vid":
                            # header line
                            continue
                        vids_across_samples.append(vntr_id)
                        if inconsistency_type == "inconsistent":
                            #print("{} {} {}".format(
                            #      inconsistency_type, sample_name, line.strip()))
                            pass
           except Exception:
               pass
        vids_counter = Counter(vids_across_samples)
        unique_vntrs = len(list(set(vids_across_samples)))
        singleton_vntrs = Counter(vids_counter.values())[1]
        print("{} {} vntrs {} ({}%) singleton vntrs (only appear in one sample)".format(
                unique_vntrs,
                inconsistency_type,
                singleton_vntrs,
                round(singleton_vntrs/unique_vntrs*100)))
        #histplot(vids_counter, suffix=inconsistency_type)

if __name__ == "__main__":
    aggregate_comparison_reports()
