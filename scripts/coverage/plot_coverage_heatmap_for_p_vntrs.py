import os
import glob

import seaborn as sns
from matplotlib import pyplot as plt
from collections import defaultdict
import pandas as pd

def plot(dataframe, plot_suffix, num_vntrs, num_samples):
    sns.set(font_scale=2)
    font_size = 60
    fig, ax = plt.subplots(figsize=(num_samples, num_vntrs))
    heatmap = sns.heatmap(dataframe, ax=ax, #annot=True,
                          vmin=0, vmax=30)
    ticks = list(range(0, 31, 5))
    ticklabels = [str(tick) for tick in ticks]
    ticklabels[-1] = ">={}".format(ticklabels[-1])
    print(len(ticks), len(ticklabels))
    colorbar_axis = heatmap.figure.axes[-1]
    colorbar_axis.set_yticks(ticks)
    colorbar_axis.set_yticklabels(ticklabels)
    colorbar_axis.tick_params(labelsize = font_size)
    plt.xlabel("Sample", fontsize=font_size)
    plt.ylabel("VNTR", fontsize=font_size)

    plt.savefig("figures/spanning_reads_heatmap_no_annot_" + plot_suffix + ".pdf", format="pdf", bbox_inches='tight')
    #plt.savefig("figures/spanning_reads_heatmap_no_annot_" + plot_suffix + ".jpg", format="jpg", bbox_inches='tight')

def get_sample_ids():
    samples_dir = "/ribosome/projects/saraj/vntr/data/pangenome_project/pacbio_hifi_mapped_reads/HG*"
    sample_ids = []
    for file_name in glob.glob(samples_dir):
        sample_ids.append(os.path.basename(file_name))
    return sample_ids

def get_vntr_ids(file_name):
    vntr_ids = []
    with open(file_name, "r") as target_file:
        for line in target_file.readlines():
            vntr_ids.append(line.strip())
    return vntr_ids

def is_high_gc_vntr(gene_name):
    high_gc_vntrs = [
        "PRNP",
        "DRD4",
        "ABCA7",
        "C9orf72",
        "MUC1",
        "NACA",
        "SLC6A3(DAT1)",
        "SLC6A3",
        "CEL",
        "NOP56",
        "HIC1",
        "INS",
        "VWA1",
        "CSTB",
        "SLC6A4",
        "TCHH",
        "EIF4A3",
        "FZD8",
        "MAOA",
        "TAF1",
        ]
    return gene_name in high_gc_vntrs

def get_vntr_name_map(target_set):
    vntr_map = {}
    vntr_names_df = pd.read_csv("../../target_vntrs/phenotype_associated_vntrs.tsv", sep="\t")
    for idx, row in vntr_names_df.iterrows():
        vid = row["VNTR ID (TRF, hg38) - at /nucleus/projects/saraj/vntr/sources/COH_analysis/databases/combined_trf_hg38/hg38_VNTRs_by_TRF.db"]
        vntr_name = row["Gene, loci"]
        vntr_map[str(vid)] = vntr_name

    return vntr_map

    vntr_map = None
    if target_set == "short":
        vntr_map = {
                "826396": "NOS3",
                "493003": "PRNP",
                "886134": "C9orf72",
                "492236": "NOP56",
                "90143": "FZD8-1",
                "90147": "FZD8-2",
                "936687": "MAOA",
            }
    elif target_set == "long":
        vntr_map = {
                "122595": "DRD4",
                "385941": "ABCA7",
                "45940": "MUC1",
                "666226": "SLC6A3",
                "915594": "CEL",
                "123860": "INS",
                "358459": "EIF4A3",
                #DUX4 665101
                "731268": "MUC22",
                "331665": "GP1BA",
                "290964": "ACAN",
                "4237": "PER3",
                "674126": "PRDM9",
                "628764": "MUC7",
                "731247": "MUC21",
            }
    elif target_set == "all_vntrs":
        vntr_map = {
                "290964": "ACAN",
                "4237": "PER3",
                "674126": "PRDM9",
                "628764": "MUC7",
                "376172": "WDR7",
                "983809": "PLIN4",
                "450626": "IL1RN",
                "49244": "TMCO1",
                "239917": "CUL4A",
                #"239917": "CUL4A-1",
                #"239918": "CUL4A-2",
                "705182": "IL4",
                "339765": "SLC6A4",
                #"339765": "SLC6A4-1",
                #"339766": "SLC6A4-2",
                #"339767": "SLC6A4-3",
                "983810": "MUC6",
                "165038": "CACNA1C",
                "75781": "NLRP3",
                "826396": "NOS3",
                "665101": "DUX4",
                "731268": "MUC22",
                "331665": "GP1BA",
                "731247": "MUC21",
                "583186": "CNBP",
                "493003": "PRNP",
                "122595": "DRD4",
                "385941": "ABCA7",
                "886134": "C9orf72",
                "45940": "MUC1",
                "181577": "NACA",
                "666226": "SLC6A3-1",
                "666245": "SLC6A3-2",
                "915594": "CEL",
                "492236": "NOP56",
                #"492236": "NOP56-1",
                #"492237": "NOP56-2",
                "123860": "INS",
                "527656": "CSTB",
                "358459": "EIF4A3", # There are more than 2 alternatives here.
                #"358459": "EIF4A3-1",
                #"358460": "EIF4A3-2",
                "90147": "FZD8",
                #"90143": "FZD8-1", # Not processed in the heatmap
                #"90147": "FZD8-2",
                #"90148": "FZD8-3",
                #"90149": "FZD8-4",
                #"90150": "FZD8-5",
                #"90151": "FZD8-6",
                "936686": "MAOA-1",
                #"936687": "MAOA-2", # Not processed in the heatmap
                "936688": "MAOA-3",
                "983819": "TAF1",
            }


    for vntr_id in vntr_map:
        if is_high_gc_vntr(vntr_map[vntr_id]):
            vntr_map[vntr_id] = vntr_map[vntr_id] + "(high-GC)"
    return vntr_map


def get_total_spanning_reads(filename):
    num_spanning_reads = 0
    with open(filename, "r") as log_file:
        for line in log_file.readlines():
            if line.startswith(">"):
                num_spanning_reads += 1
    return num_spanning_reads

def analyze_disease_vntrs(vntr_target_set):
    vntr_target_set_name_map = {"short": "short_vntrs",
                                "long": "long_vntrs",
                                "all_vntrs": "all_vntrs_v5_unique"}
    target_vntr_filename = "/nucleus/projects/saraj/vntr/sources/pangenome_project/target_vntrs/disease_associated_" +\
                            vntr_target_set_name_map[vntr_target_set] + ".txt"
    vntr_ids = get_vntr_ids(file_name=target_vntr_filename)
    sample_ids = get_sample_ids()
    print("Num samples: ", len(sample_ids))
    vntr_name_map = get_vntr_name_map(vntr_target_set)

    removed_from_vntr_name_map = set()
    spanning_reads_df = pd.DataFrame(columns=vntr_name_map.values(), index=sample_ids)
    for vntr_id in vntr_ids:
        for sample_id in sample_ids:
            log_file_name = "/nucleus/projects/saraj/vntr/sources/pangenome_project" +\
                "/logs/logs_pairwise_align_disease_all_vntrs_v5/" + \
                "log_{}_aligned_GRCh38_winnowmap.sorted.bam.log_vid_{}.txt".format(sample_id, vntr_id)
            spanning_reads = 0
            if os.path.exists(log_file_name):
                spanning_reads = get_total_spanning_reads(log_file_name)
            if vntr_id not in vntr_name_map:
                removed_from_vntr_name_map.add(vntr_id)
                continue
            vntr_gene_name = vntr_name_map[vntr_id]
            spanning_reads_df.at[sample_id, vntr_gene_name] = spanning_reads

    print("vntrs removed from vntr_name_map {}".format(len(removed_from_vntr_name_map)))
    # Sort and prepare the dataframe for plotting.
    spanning_reads_df = spanning_reads_df.fillna(0)
    spanning_reads_df["median over VNTRs"] = spanning_reads_df.median(axis=1)
    spanning_reads_df.loc["median over\nsamples"] = spanning_reads_df.median(axis=0)
    #spanning_reads_df.loc["std over\nsamples"] = spanning_reads_df.std(axis=0)
    #print("median std over samples: ", spanning_reads_df.loc["std over\nsamples"].median())
    spanning_reads_df.at["median over\nsamples", "median over VNTRs"] = 0
    #spanning_reads_df.at["std over\nsamples", "median over VNTRs"] = 0
    spanning_reads_df = spanning_reads_df.sort_values(by="median over VNTRs", ascending=False)
    spanning_reads_df = spanning_reads_df.sort_values(by="median over\nsamples", ascending=False, axis=1)
    print(spanning_reads_df.loc["median over\nsamples"])
    all_median_over_samples = list(spanning_reads_df.loc["median over\nsamples"].astype(int))

    print("num vntrs with >=15 SR {} and len all_median_over_samples {}".format(len(
        [sr for sr in all_median_over_samples if sr >=15]),
        len(all_median_over_samples)))
    # Print basic dataframe information.
    #print(spanning_reads_df.head())
    print("Dataframe shape: ", spanning_reads_df.shape)

    # Store the dataframe and plot.
    spanning_reads_df.to_csv("spanning_reads_disease_{}.csv".format(vntr_target_set))
    plot(dataframe=spanning_reads_df.transpose(),
         plot_suffix=vntr_target_set,
         num_samples=spanning_reads_df.shape[0],
         num_vntrs=spanning_reads_df.shape[1])

def main():
    #analyze_disease_vntrs("long")
    #analyze_disease_vntrs("short")
    analyze_disease_vntrs("all_vntrs")

if __name__ == "__main__":
    main()
