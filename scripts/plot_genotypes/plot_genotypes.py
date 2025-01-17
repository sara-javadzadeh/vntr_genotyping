import os
import numpy as np
import argparse
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
from pandas import DataFrame
from collections import defaultdict
from scipy.stats import entropy

from Bio.SeqUtils import GC

from advntr.models import load_unique_vntrs_data

def discard_strs(vids, filename):
    str_ids = []
    with open(filename, "r") as str_file:
        str_ids = str_file.readlines()
    str_ids = [int(single_id) for single_id in str_ids]
    vids = [vid for vid in vids if vid not in str_ids]
    return vids

def get_gc_content(ref_vntr, flanking_threshold, round_to_int=True, verbose=False):
    """ return GC contents of vntr motif and flanking region"""
    # flanking regions
    left_flanking_region = ref_vntr.left_flanking_region[-1*flanking_threshold:]
    right_flanking_region = ref_vntr.right_flanking_region[:flanking_threshold]
    # VNTR sequence in reference model
    vntr_seq = "".join(ref_vntr.repeat_segments)
    whole_region = "".join([left_flanking_region,vntr_seq,right_flanking_region])

    whole_region_GC = GC(whole_region)

    if verbose:
        print("whole region", whole_region, whole_region_GC)

    if round_to_int:
        return int(round(whole_region_GC))
    else:
        return whole_region_GC

def process_genotype_file(filename, targetted_vntrs, is_expanded=False):
    genotype_file = open(filename, "r")
    lines = genotype_file.readlines()
    genotyped_vntrs = {}
    genotyped_vntr_ids = set()
    current_vntr_id = None
    none_genotypes = 0
    for line in lines:
        if not (line[0].isdigit() or line.startswith("None")):
            # There might be other lines including warnings or logs, skip those lines.
            continue
        if line.startswith("None") or ("/" in line):
            # It's a genotype line
            if current_vntr_id is None:
                print("Error in parsing genotype file. current VNTR id not set")
                return
            if current_vntr_id in targetted_vntrs or int(current_vntr_id) in targetted_vntrs:
                vntr_genotype = line.strip()
                if is_expanded:
                    vntr_genotype = vntr_genotype.split()[0]
                if not vntr_genotype.startswith("None"):
                    genotyped_vntrs[current_vntr_id] = vntr_genotype
                    genotyped_vntr_ids.add(current_vntr_id)
                    none_genotypes += 1
            current_vntr_id = None
        else:
            # It's a line with vntr_id
            current_vntr_id = line.strip()
    return genotyped_vntrs, genotyped_vntr_ids

def process_targetted_vntrs_file(targetted_vntrs_filename):
	targetted_vntrs = []
	with open(targetted_vntrs_filename, "r") as targetted_vntrs_file:
		lines = targetted_vntrs_file.readlines()
		targetted_vntrs = [line.strip().split(' ')[0] for line in lines]
	return targetted_vntrs

def get_vntr_name_map():
    vntr_name_map = {
            "290964": "ACAN",
            "4237": "PER3",
            "674126": "PRDM9",
            "628764": "MUC7",
            "376172": "WDR7",
            "983809": "PLIN4",
            "450626": "IL1RN",
            "49244": "TMCO1",
            "239917": "CUL4A",
            "705182": "IL4",
            "983816": "SLC6A4-1",
            "339765": "SLC6A4-2",
            "983810": "MUC6",
            "983808": "AVPR1A",
            "983811": "FXN",
            "75781": "NLRP3",
            "983812": "EIF3H-1",
            "983813": "EIF3H-2",
            "826396": "NOS3",
            "331665": "GP1BA",
            "731278": "HCG22",
            "731268": "MUC22",
            "665101": "DUX4",
            "583186": "CNBP",
            "493003": "PRNP",
            "122595": "DRD4",
            "385941": "ABCA7",
            "886134": "C9orf72",
            "45940": "MUC1",
            "181577": "NACA",
            "666226": "SLC6A3-1",
            "666245": "SLC6A3-2",
            "666244": "SLC6A3-3",
            "915594": "CEL",
            "492236": "NOP56",
            "983814": "HIC1",
            "123860": "INS",
            "983815": "VWA1",
            "527656": "CSTB",
            "983817": "TCHH",
            "358459": "EIF4A3",
            "936686": "MAOA-1",
            "936688": "MAOA-2",
            "944729": "TAF1",
            "983820": "FRA16B",
            "983821": "PCSK6",
            "746692": "TENT5A",
            "983808": "AVPR1A",
        }
    return vntr_name_map

def apply_filtered_genotypes():
    filter_log_filename = "../filter/output_disease_all.txt"
    # Example in the filter output file:
    filtered_gt_map = {}
    with open(filter_log_filename, "r") as filter_file:
        for line in filter_file.readlines():
            if not line.startswith("sample"):
                # Skip empty lines.
                continue
            _, sample, _, vid, _, _, genotype = line.strip().split()
            vid = vid.replace(":", "")
            filtered_gt_map[(sample, vid)] = genotype
    return filtered_gt_map

def get_samples(filename):
    samples = []
    with open(filename) as samples_ids_file:
        for line in samples_ids_file.readlines():
            for sample_id in line.split(","):
                samples.append(sample_id.strip())
    print("Working with {} samples".format(len(samples)))
    return samples

def get_id_to_genic_vntr_map(long_vntr_db, short_vntr_db):
    id_to_vntr_map = {}
    vntrs = load_unique_vntrs_data(long_vntr_db)
    for vntr in vntrs:
        id_to_vntr_map[vntr.id] = vntr
    #print("len of map {} after adding long vntrs".format(len(id_to_vntr_map)))
    vntrs = load_unique_vntrs_data(short_vntr_db)
    for vntr in vntrs:
        id_to_vntr_map[vntr.id] = vntr
    #print("len of map {} after adding short vntrs".format(len(id_to_vntr_map)))

    # Finding the number of VNTRs per chromosome
    chr_based_num_vntrs = defaultdict(float)
    num_vntrs = len(id_to_vntr_map)
    for vid in id_to_vntr_map.keys():
        chr_based_num_vntrs[id_to_vntr_map[vid].chromosome] += 1
    for chromosome in chr_based_num_vntrs.keys():
        chr_based_num_vntrs[chromosome] = round(chr_based_num_vntrs[chromosome]/num_vntrs* 100, 2)
    #print("chr_based_num_vntrs: ", sorted(chr_based_num_vntrs.items(), key=lambda x: x[1], reverse=True))
    return id_to_vntr_map

def get_variation_category_map():
    import re
    variation_category_map = {}
    p_vntrs_file = "../../target_vntrs/phenotype_associated_vntrs.tsv"
    p_vntrs_df = pd.read_csv(p_vntrs_file, sep="\t")
    for idx, row in p_vntrs_df.iterrows():
        vntr_name = row["Gene, loci"]
        vntr_id = row["VNTR ID (TRF, hg38) - at /nucleus/projects/saraj/vntr/sources/COH_analysis/databases/combined_trf_hg38/hg38_VNTRs_by_TRF.db"]
        variation_type = row["Associated variant category (major length variation is >=10 RC different in the RC allele range reported)"]
        variation_type = str(variation_type).lower()
        variation_type = re.sub(r"\(.*", "", variation_type).strip()
        vntr_name = str(vntr_name).strip()
        vntr_name = re.sub(r"\(.*", "", vntr_name)
        vntr_name = re.sub(r"/.*", "", vntr_name).strip()
        variation_category_map[str(vntr_id)] = variation_type
    return variation_category_map

def get_genic_genotype_distributions(long_vntr_db,
                                     short_vntr_db,
                                     strs_filename,
                                     flanking_threshold,
                                     g_vntrs_log_dir,
                                     verbose=False):
    id_to_vntr_map = get_id_to_genic_vntr_map(long_vntr_db=long_vntr_db,
                                              short_vntr_db=short_vntr_db)
    vntr_ids = id_to_vntr_map.keys()
    vntr_ids = discard_strs(vntr_ids, strs_filename)
    print("len vntr_ids: ", len(vntr_ids))
    samples = get_samples(os.path.join(g_vntrs_log_dir, "samples_ids.txt"))
    genotypes_df = DataFrame(columns=samples, index=vntr_ids)
    haplotypes = [sample + "-1" for sample in samples] + \
                     [sample + "-2" for sample in samples]
    allele_columns = haplotypes + \
                     ["vntr_id", "motif_len", "gc_content"] + \
                     ["median_rc", "var_rc", "cov_rc"] + \
                     ["median_est_len", "start_coordinate"]
    alleles_df = DataFrame(columns=allele_columns, index=vntr_ids)
    skipped_vntrs = set()
    for sample in samples:
        genotype_filenames = [g_vntrs_log_dir +\
                              "output_per_sample_{}.txt".format(sample)]
        for genotype_filename in genotype_filenames:
            genotyped_vntrs, genotyped_vntr_ids = process_genotype_file(
                    filename=genotype_filename,
                    targetted_vntrs=vntr_ids)
            for vntr_id in vntr_ids:
                if str(vntr_id) not in genotyped_vntr_ids:
                    skipped_vntrs.add(vntr_id)
                    continue
                genotype = genotyped_vntrs[str(vntr_id)]
                genotypes_df.at[vntr_id, sample] = genotype
                alleles_df.at[vntr_id, sample + "-1"] = int(genotype.split("/")[0])
                alleles_df.at[vntr_id, sample + "-2"] = int(genotype.split("/")[1])
                alleles_df.at[vntr_id, "vntr_id"] = vntr_id
                alleles_df.at[vntr_id, "start_coordinate"] = "{}:{}".format(
                        id_to_vntr_map[vntr_id].chromosome,
                        id_to_vntr_map[vntr_id].start_point)
                alleles_df.at[vntr_id, "motif_len"] = len(id_to_vntr_map[vntr_id].pattern)
                alleles_df.at[vntr_id, "gc_content"] = get_gc_content(id_to_vntr_map[vntr_id], flanking_threshold)


    #print("Len skipped vntrs {}".format(
    #    len(list(skipped_vntrs))))

    # Sort VNTRs in the alleles_df based on median allele value
    print("Computing median_rc")
    alleles_df["median_rc"] = alleles_df[haplotypes].median(axis=1)
    print("Computing median est len")
    alleles_df["median_est_len"] = alleles_df["median_rc"] * alleles_df["motif_len"]
    print("Computing cov")
    alleles_df["cov_rc"] = alleles_df[haplotypes].std(axis=1) / alleles_df[haplotypes].mean(axis=1)
    print("Computing variance")
    alleles_df["var_rc"] = alleles_df[haplotypes].var(axis=1)
    alleles_df = alleles_df.sort_values(by=["median_rc"])

    # print the vntrs with motif len between 80 and 90. Why are they so many?
    if verbose:
        print_columns = ["vntr_id", "motif_len", "gc_content"] + \
                         ["median_rc", "var_rc", "cov_rc"] + \
                         ["median_est_len", "start_coordinate"]
        print(alleles_df[(alleles_df["motif_len"] >= 83) & (alleles_df["motif_len"] <= 84)][print_columns].to_string())
        print("median GC in motif len limit: ",
              alleles_df[(alleles_df["motif_len"] >= 83) & (alleles_df["motif_len"] <= 84)]["gc_content"].median())

    sorted_vntrs = alleles_df.index
    return genotypes_df, alleles_df, sorted_vntrs, haplotypes, samples


def get_disease_genotype_distributions(flanking_threshold, vntr_db, p_vntrs_ids_file, p_vntrs_log_dir):
    vntr_name_map = get_vntr_name_map()
    vntr_name_to_vntr = get_vntr_name_to_vntr(vntr_name_map.values(), vntr_db)
    variation_category_map = get_variation_category_map()
    samples = get_samples(os.path.join(p_vntrs_log_dir, "samples_ids.txt"))
    targetted_vntrs_filename = p_vntrs_ids_file
    filtered_gt_map = apply_filtered_genotypes()
    vntr_ids = process_targetted_vntrs_file(targetted_vntrs_filename)
    genotypes_df = DataFrame(columns=samples, index=vntr_name_map.values())
    haplotypes = [sample + "-1" for sample in samples] + \
                     [sample + "-2" for sample in samples]
    allele_columns = haplotypes + \
                     ["vntr_gene_name", "motif_len", "gc_content"] + \
                     ["median_rc", "var_rc", "cov_rc", "median_est_len", "variation category"]
    alleles_df = DataFrame(columns=allele_columns, index=vntr_name_map.values())
    skipped_vntrs = set()
    for sample in samples:
        genotype_filename = "{}/output_per_sample_{}.txt".format(p_vntrs_log_dir, sample)
        genotyped_vntrs, genotyped_vntr_ids = process_genotype_file(
                filename=genotype_filename,
                targetted_vntrs=vntr_ids)
        for vntr_id in vntr_ids:
            if vntr_id in genotyped_vntr_ids:
                if vntr_id not in vntr_name_map:
                    skipped_vntrs.add(vntr_id)
                    continue
                vntr_name = vntr_name_map[vntr_id]
                genotype = genotyped_vntrs[vntr_id]
                if (sample, vntr_id) in filtered_gt_map:
                    # Update the genotype based on filtered results.
                    genotype = filtered_gt_map[(sample, vntr_id)]
                genotypes_df.at[vntr_name, sample] = genotype
                alleles_df.at[vntr_name, sample + "-1"] = int(genotype.split("/")[0])
                alleles_df.at[vntr_name, sample + "-2"] = int(genotype.split("/")[1])
                alleles_df.at[vntr_name, "vntr_gene_name"] = vntr_name
                alleles_df.at[vntr_name, "motif_len"] = len(vntr_name_to_vntr[vntr_name].pattern)
                alleles_df.at[vntr_name, "gc_content"] = get_gc_content(
                            vntr_name_to_vntr[vntr_name], flanking_threshold)
                vntr_name_for_var = vntr_name.replace("-1", "").replace("-2", "").replace("-3", "")
                alleles_df.at[vntr_name, "variation category"] = variation_category_map[vntr_id]

    genotypes_df = genotypes_df.drop_duplicates()
    alleles_df = alleles_df.drop_duplicates()
    # Sort VNTRs in the alleles_df based on median allele value
    alleles_df["median_rc"] = alleles_df[haplotypes].median(axis=1)
    alleles_df["median_est_len"] = alleles_df["median_rc"] * alleles_df["motif_len"]
    alleles_df["cov_rc"] = alleles_df[haplotypes].std(axis=1) / alleles_df[haplotypes].mean(axis=1)
    alleles_df["var_rc"] = alleles_df[haplotypes].var(axis=1)
    alleles_df = alleles_df.sort_values(by=["median_rc"])

    sorted_vntrs = alleles_df.index
    return genotypes_df, alleles_df, sorted_vntrs, haplotypes


def get_vntr_name_to_vntr(vntr_names, vntr_db):
    vntr_name_map = get_vntr_name_map()
    vntr_name_to_vntr = {}
    vntrs = load_unique_vntrs_data(vntr_db)
    for vntr in vntrs:
        if str(vntr.id) in vntr_name_map:
            vntr_name_to_vntr[vntr_name_map[str(vntr.id)]] = vntr
    return vntr_name_to_vntr

def multi_box_plot(alleles_df, haplotypes, output_postfix,
             palette="RdBu", title=None,
             font_size=8, rotation=90,
             ):
    plt.clf()
    sns.set_style("darkgrid")
    filename = "{}_multi_box_plot".format(output_postfix)
    fig, ax = plt.subplots(nrows=1,
                   ncols=3,
                   sharey=True,
                   gridspec_kw={'width_ratios': [2.5, 2.5, 1]})

    # Remove axis lines in the middle
    ax[1].spines["left"].set_visible(False)
    ax[2].spines["left"].set_visible(False)
    ax[0].spines["right"].set_visible(False)
    ax[1].spines["right"].set_visible(False)
    ax[1].tick_params(axis="y", left=False)
    ax[2].tick_params(axis="y", left=False)

    # Plot the three subplots
    subcategory_alleles_df = alleles_df[alleles_df["variation category"] == "major length variation"][haplotypes]
    sns.boxplot(data=subcategory_alleles_df.transpose(),
                ax=ax[0])
    ax[0].set(xlabel="Major RC\nvariation")
    ax[0].set_ylabel("Alleles in repeat counts")
    ax[0].tick_params(axis="x", rotation=rotation, labelsize=font_size)
    subcategory_alleles_df = alleles_df[alleles_df["variation category"] == "minor length variation"][haplotypes]
    sns.boxplot(data=subcategory_alleles_df.transpose(),
                ax=ax[1])
    ax[1].set(xlabel="Minor RC\nvariation")
    ax[1].tick_params(axis="x", rotation=rotation, labelsize=font_size)
    subcategory_alleles_df = alleles_df[alleles_df["variation category"] == "specific allele sequence"][haplotypes]
    sns.boxplot(data=subcategory_alleles_df.transpose(),
                ax=ax[2])
    ax[2].set(xlabel="Specific\nallele\nsequence")
    ax[2].tick_params(axis="x", rotation=rotation, labelsize=font_size)
    # Save figure
    plt.savefig(filename + ".pdf", bbox_inches="tight", format="pdf")

def box_plot(alleles_df, plot_swarm, output_postfix,
             palette="RdBu", title=None,
             font_size=10, rotation=90,
             violins=False):
    plt.clf()
    filename = "{}_box_plot".format(output_postfix)
    if plot_swarm:
        filename += "_swarm"
        plot = sns.swarmplot(data=alleles_df.transpose(),
                         color="grey",
                         size=2)
    if violins:
        plot = sns.violinplot(data=alleles_df.transpose(),
                       palette=palette,
                       fill=False)
    else:
        plot = sns.boxplot(data=alleles_df.transpose(),
                   palette=palette)
    if output_postfix.startswith("all"):
        font_size = 10
        rotation = 90
    plt.xticks(rotation=rotation, size=font_size)
    plot.set(xlabel="VNTRs",
             ylabel="Alleles in repeat counts",
             title=title)

    plot.figure.savefig(filename + ".pdf", bbox_inches="tight", format="pdf")

def histogram(x, xlabel, ylabel, num_bins, filename, cumulative=False, xlim=None):
    # Apply xlim as upperbound
    if xlim is not None:
        x = [min(value, xlim) for value in x]
    plt.clf()
    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 16})
    ax = plt.hist(x,
                      cumulative=cumulative,
                      bins=num_bins)

    xticks = range(0, int(max(x))+1, 10)
    xticklabels = [str(xtick) for xtick in xticks]
    xticklabels[-1] = ">={}".format(xticklabels[-1])
    ax = plt.gca()
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=45)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # Set y tick labels as percentage
    ax = plt.gca()
    num_ticks = 10
    num_items = len(x)
    steps = int(num_items / num_ticks)
    if cumulative is False:
        # Only pick yticks up to 15%
        yticks = list(range(0, int(num_items), steps))[:4]
    else:
        yticks = list(range(0, int(num_items), steps))
    ylabels = [str(round(tick/num_items*100)) + "%" for tick in yticks]
    ax.set(yticks=yticks, yticklabels=ylabels)

    plt.savefig(filename, format="pdf", bbox_inches="tight")

    # Change fonts back to 10 for other figures
    plt.rcParams.update({'font.size': 10})


def plot_side_histogram(dataframe, ax, axis, column):
    color = sns.color_palette("Blues", 6)[3]
    if axis == "x":
        data = dataframe[column]
        orientation = "vertical"
        ax.set_ylabel("Count VNTRs", fontsize=8)
        ax.yaxis.label.set_color(color)
        ax.tick_params(axis="y", colors=color)
    if axis == "y":
        data = dataframe[column]
        orientation = "horizontal"
        ax.set_xlabel("Count VNTRs", fontsize=8)
        ax.yaxis.label.set_color(color)
        ax.tick_params(axis="y", colors=color)
    ax.hist(data,
             bins=40,
             color=color,
             orientation=orientation
             )


def scatterplot(raw_dataframe, x, y, xlabel, ylabel, filename,
                hue, xlim=None, ylim=None,
                add_noise_line=False, add_quantiles=False,
                add_median=False, is_histplot=False,
                legend=False, place_curve=False):
    plt.clf()
    alpha=1
    linewidth=0
    dataframe = raw_dataframe.copy()
    if add_quantiles:
        # Make points lighter, so the quantile likes are more clear.
        alpha = 0.5
    if add_median:
        alpha = 0.8
        linewidth=0.5
    if xlim is not None:
        dataframe[x] = raw_dataframe[x].clip(0, xlim)
    if ylim is not None:
        dataframe[y] = raw_dataframe[y].clip(0, ylim)
    if is_histplot:
        # Create figure base and adjust ratios
        fig = plt.figure(figsize=(6, 6))
        gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                              left=0.1, right=0.9, bottom=0.1, top=0.9,
                              wspace=0.05, hspace=0.05)
        # Adjust spacing between the two subplots
        plt.subplots_adjust(wspace=0.2)

        color = sns.color_palette("Blues", 6)[5]
        # Create the axes
        ax = fig.add_subplot(gs[1, 0])
        ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
        ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

        # Remove double labels
        ax_histx.tick_params(axis="x", labelbottom=False)
        ax_histy.tick_params(axis="y", labelleft=False)
        # Plot side histograms
        plot_side_histogram(dataframe=dataframe,
                       column=x,
                       ax=ax_histx,
                       axis="x")
        plot_side_histogram(dataframe=dataframe,
                       column=y,
                       ax=ax_histy,
                       axis="y")
        # Plot main histplot
        plot = sns.histplot(data=dataframe,
                            x=x,
                            y=y,
                            ax=ax,
                            color=color,
                            bins=60,
                           )
        plot = ax
    else:
        plot = sns.scatterplot(data=dataframe,
                            x=x,
                            y=y,
                            hue=hue,
                            linewidth=linewidth,
                            alpha=alpha,
                            )

    if xlim is not None:
        xticks, xticklabels = plt.xticks()
        xticklabels[-1] = ">={}".format(xticklabels[-1])
        plot.set_xticks(xticks)
        plot.set_xticklabels(xticklabels)
    if ylim is not None:
        yticks, yticklabels = plt.yticks()
        yticklabels = [str(int(ytick)) for ytick in yticks]
        yticklabels[-2] = ">={}".format(yticklabels[-2])
        plot.set_yticks(yticks)
        plot.set_yticklabels(yticklabels)
    # Add a line plot with a noise driven from a normal distribution,
    # As a way to compare.
    if add_noise_line:
        line_xs, line_ys = [], []
        dataframe_xs = sorted(list(dataframe[x]))
        for line_x in dataframe_xs:
            line_xs.append(line_x)
            line_ys.append(np.random.normal(line_x,1,1))
        plt.plot(line_xs,
                 line_ys)

    # Add lines showing quantiles of y per each x.
    # This works with x with integer values.
    if add_quantiles or add_median:
        # Select more color as I need, as we're descarding light most colors.
        if add_quantiles:
            colors = sns.color_palette("Blues_r", 8)
        else:
            colors = ["black"]
        unique_xs = sorted(list(set(dataframe[x].dropna())))
        x_bins = np.linspace(start = min(unique_xs),
                            stop = max(unique_xs),
                            num = 70)
        bin_width = x_bins[1] - x_bins[0]
        quantile_xs = []
        if add_quantiles:
            quantiles = [10, 25, 50, 75, 90]
        else:
            quantiles = [50]
        quantile_values = defaultdict(list)
        for x_bin in x_bins:
            x_bin_min, x_bin_max = x_bin, x_bin + bin_width
            ys_for_current_x = dataframe[(dataframe[x] >= x_bin_min) & (dataframe[x] < x_bin_max)][y]
            # Skip empty or small bins
            if len(ys_for_current_x) > 10:
                quantile_xs.append(x_bin_min + bin_width/2.0)
                for q in quantiles:
                    quantile_values[q].append(np.quantile(ys_for_current_x, q/100.0))
        for idx, q in enumerate(quantiles[::-1]):
            if add_quantiles:
                label="{}th percentile".format(q)
            else:
                label="Median"
            plot.plot(quantile_xs,
                    quantile_values[q],
                    label=label,
                    color = colors[idx])
            if add_median and False:
                print(quantile_values[q])
                print(quantile_xs)
        plot.legend()
    # Add a diagonal curve indicating VNTRs with length equal to a certain value. y = n/x.
    if place_curve:
        for idx, mult_value in enumerate([150, 1000, 2000, 4000]):
            # Minimum X such that y is less than 150, so figure is centered on points.
            min_x = mult_value/150.0
            med_x = max(min_x, 50)
            # med_x is a value of x to have logarithmically many points before.
            # This is to avoid gaps in the curve for small x values.
            curve_x = list(np.logspace(np.log10(min_x), np.log10(med_x), 1000, base=10))
            for curr_x in np.linspace(med_x, 500, num=1000):
                curve_x.append(curr_x)
            curve_y = [mult_value/curr_x for curr_x in curve_x]
            sns.scatterplot(x=curve_x,
                            y=curve_y,
                            color="grey",
                            s=0.5)
            # Print the value of each mult_value on the plot
            if "disease" in filename:
                # Don't print for genic VNTRs, too many items to print
                plot.text(s=mult_value, x= 200 + idx*10, y=mult_value/200.0, size=5)
        filename = filename.replace(".pdf", "_w_curve.pdf")

    if legend:
        # Put a legend with 3 columns to fit all the names
        plot.legend(ncol=3)
        filename = filename.replace(".pdf", "_legend.pdf")
    else:
        mult_value_for_text = 150
        # Write the gene name for each point next to it.
        plt.legend([],[], frameon=False)
        max_x = dataframe[x].max()
        max_y = dataframe[y].max()
        delta_x = max_x/100.0
        delta_y = max_y/100.0
        for idx, row in dataframe.iterrows():
            text_x = row[x] + delta_x
            text_y = row[y] + delta_y
            # For MUC22 VNTR, change the X and Y to move the text inside the plot
            if text_x >= 500:
                text_x = min(row[x], max_x) - 5 * (delta_x)
                text_y = row[y] + (2 * delta_y)
            # Only write the text for spaced points
            if "disease" in filename:
                # Don't print for genic VNTRs, too many items to print
                plot.text(s=row[hue], x=text_x, y=text_y, size=5)

    plot.set(xlabel=xlabel,
             ylabel=ylabel)
    plot.figure.savefig(filename, bbox_inches="tight", format="pdf")

def group_vntrs_median_repeat_count(alleles_df, vntrs):
    low_median_vntrs, medium_median_vntrs, high_median_vntrs = [], [], []
    for vntr in vntrs:
        if vntr not in alleles_df.index:
            print("{} not in alleles_df".format(vntr))
            continue
        if len(alleles_df.loc[vntr].dropna()) == 0:
            print("No alleles for {} in alleles_df".format(vntr))
            continue
        if alleles_df.loc[vntr].median() < 10:
            low_median_vntrs.append(vntr)
        elif alleles_df.loc[vntr].median() < 20:
            medium_median_vntrs.append(vntr)
        else:
            high_median_vntrs.append(vntr)
    return low_median_vntrs, medium_median_vntrs, high_median_vntrs

def group_vntrs_variance(alleles_df, vntrs):
    high_var, low_var = [], []
    for vntr in vntrs:
        if vntr not in alleles_df.index:
            print("{} not in alleles_df".format(vntr))
            continue
        if len(alleles_df.loc[vntr].dropna()) == 0:
            print("No alleles for {} in alleles_df".format(vntr))
            continue
        if alleles_df.loc[vntr].var() < 20:
            low_var.append(vntr)
        else:
            high_var.append(vntr)
    return low_var, high_var

def group_vntrs_entropy(alleles_df, vntrs):
    high_ent, low_ent = [], []
    entropies = []
    for vntr in vntrs:
        if vntr not in alleles_df.index:
            print("{} not in alleles_df".format(vntr))
            continue
        if len(alleles_df.loc[vntr].dropna()) == 0:
            print("No alleles for {} in alleles_df".format(vntr))
            continue
        entropy_val = entropy(list(alleles_df.loc[vntr].value_counts()))
        print(alleles_df.loc[vntr].value_counts())
        entropies.append(entropy_val)
        if entropy_val < 20:
            low_ent.append(vntr)
        else:
            high_ent.append(vntr)
    print("entropies: ", [round(item, 2) for item in entropies])
    return low_ent, high_ent

def group_vntrs_cv(alleles_df, vntrs):
    #Coefficient of Variation
    high_cv, low_cv = [], []
    cv_vals = []
    for vntr in vntrs:
        if vntr not in alleles_df.index:
            print("{} not in alleles_df".format(vntr))
            continue
        if len(alleles_df.loc[vntr].dropna()) == 0:
            print("No alleles for {} in alleles_df".format(vntr))
            continue
        cv_val = alleles_df.loc[vntr].std() / alleles_df.loc[vntr].mean()
        cv_vals.append(cv_val)
        #print("{} std {} cv value {} for distribution {}".format(
        #       vntr,alleles_df.loc[vntr].std(), cv_val, alleles_df.loc[vntr]))
        if cv_val < 0.3:
            low_cv.append(vntr)
        else:
            high_cv.append(vntr)
    print("max CV {} min CV {}".format(max(cv_vals), min(cv_vals)))
    return low_cv, high_cv


def plot_genotypes(genotypes_df, alleles_df, vntrs, haplotypes,
                    filename_prefix, title_prefix,
                    vntr_id_column,
                    plot_median_based_plots=False, plot_swarm=False):
    ### Box plots
    print("Plotting boxplots")
    if "disease" in filename_prefix:
        # Plot all disease VNTRs
        box_plot(alleles_df[haplotypes],
                    plot_swarm=plot_swarm,
                    title=title_prefix,
                    rotation=90,
                    output_postfix=filename_prefix+"_all_sorted")
        # Plot all disease VNTRs separated by polymorphism type
        multi_box_plot(alleles_df,
                    haplotypes=haplotypes,
                    rotation=90,
                    output_postfix=filename_prefix+"_all_sorted_polymorphism")
        # Plot figues based on low, medium or high median of the distribution
        if plot_median_based_plots:
            low_median = alleles_df[alleles_df["median_rc"] < 10]
            medium_median = alleles_df[(alleles_df["median_rc"] >= 10 & alleles_df["median_rc"] < 20)]
            high_median = alleles_df[alleles_df["median_rc"] >= 20]
            box_plot(low_median[haplotypes],
                        plot_swarm=plot_swarm,
                        title=title_prefix+" with Median Repeats < 10",
                        output_postfix=filename_prefix+"_low_median_sorted")
            box_plot(medium_median[haplotypes],
                        plot_swarm=plot_swarm,
                        title=title_prefix+" with 10 <= Median Repeats < 20",
                        output_postfix=filename_prefix+"_medium_median_sorted")
            box_plot(high_median[haplotypes],
                        plot_swarm=plot_swarm,
                        title=title_prefix+" with Median Repeats >= 20",
                        output_postfix=filename_prefix+"_high_median_sorted")

    # Plot figure based on high variance in the VNTR genotype distribution
    if "disease" in filename_prefix:
        var_threshold = 10
        var_rc_lim = 300
    else: # genic VNTRs
        var_threshold = 150
        var_rc_lim = 1000
    high_var = alleles_df[alleles_df["var_rc"] >= var_threshold]
    low_var = alleles_df[alleles_df["var_rc"] < var_threshold]
    #print("len of alleles_df {} len of low var {} len of high var {}".format(
    #        len(alleles_df), len(low_var), len(high_var)))
    palette="RdBu"
    if "disease" in filename_prefix:
        box_plot(high_var[haplotypes],
                plot_swarm=plot_swarm,
                palette=palette,
                title=title_prefix+" with variance >= {} in repeats".format(var_threshold),
                output_postfix=filename_prefix+"_high_var_sorted")
        box_plot(low_var[haplotypes],
                palette=palette,
                plot_swarm=plot_swarm,
                title=title_prefix+" with variance < {} in repeats".format(var_threshold),
                font_size=5,
                rotation=90,
                output_postfix=filename_prefix+"_low_var_sorted")
    if "disease" in filename_prefix:
        cv_threshold = 0.3
    else: # genic VNTRs
        cv_threshold = 0.6
    high_cv = alleles_df[alleles_df["cov_rc"] >= cv_threshold]
    low_cv = alleles_df[alleles_df["cov_rc"] < cv_threshold]
    box_plot(high_cv[haplotypes],
                plot_swarm=plot_swarm,
                title=title_prefix+" with High CoV in Repeats",
                output_postfix=filename_prefix+"_high_cov_sorted")



    ### Scatter plots
    print("Plotting scatterplots")

    if "genic" in filename_prefix:
        scatterplot(alleles_df,
                    x="median_rc",
                    y="var_rc",
                    xlim=500,
                    xlabel="Repeat counts",
                    ylabel="variance in repeat counts in samples",
                    hue=vntr_id_column,
                    add_quantiles=True,
                    filename=filename_prefix+"_median_rc_var_rc_w_quantiles.pdf",
                    )

    scatterplot(alleles_df,
                x="motif_len",
                y="median_rc",
                xlabel="VNTR motif length (bp)",
                ylabel="VNTR repeat count",
                filename=filename_prefix+"_displot_motif_len_rc_w_median.pdf",
                add_median=True,
                is_histplot=True,
                hue=None,
                legend=True,
                ylim=100,
                )

def parse_arguments():
    parser = argparse.ArgumentParser(
                prog="vntr_genotypes_plot",
                description="Plot the vntr genotypes.")

    parser.add_argument("--vntr-db", help="Path to the db file for vntrs.",
                        type=str,
                        default="../../target_vntrs/p_vntrs_g_vntrs.db")
    parser.add_argument("--g-vntrs-short-db", help="Path to the db file for short g-vntrs.",
                        type=str,
                        default="../../target_vntrs/illumina_probed_short_vntrs.db")
    parser.add_argument("--g-vntrs-long-db", help="Path to the db file for long g-vntrs.",
                        type=str,
                        default="../../target_vntrs/pacbio_probed_long_vntrs.db")
    parser.add_argument("--p-vntrs-ids", help="Path to the text file including all p-vntr ids.",
                        type=str,
                        default="../../target_vntrs/disease_associated_all_vntrs_v5_unique.txt")
    parser.add_argument("--g-vntrs-log-dir", help="Path to the directory including g_vntrs genotypes.",
                        type=str, required=True)
    parser.add_argument("--p-vntrs-log-dir", help="Path to the directory including p_vntrs genotypes.",
                        type=str, required=True)
    parser.add_argument("--ids-to-discard-file", help="Path to the text file including ids to be discarded (e.g. test strs).",
                        type=str,
                        default="../../target_vntrs/str_ids_to_discard.txt")
    parser.add_argument("--gc-flanking-threshold", help="Threshold for how many nucleotides" + \
                            " into each of the flanking regions is considered for GC content.",
                        type=int,
                        default=100)
    args = parser.parse_args()
    return args

def main():

    args = parse_arguments()

    # Set flags on which figures to plot on this run
    plot_swarm = False
    plot_median_based_plots = False
    plot_disease_vntrs = True
    plot_genic_vntrs = True

    plt.rcParams.update({'font.size': 10})
    if plot_disease_vntrs:
        # Graphs for disease associated VNTRs
        print("Working with phenotype associated VNTRs")
        genotypes_df, alleles_df, vntrs, haplotypes = get_disease_genotype_distributions(
                        p_vntrs_ids_file=args.p_vntrs_ids,
                        vntr_db=args.vntr_db,
                        p_vntrs_log_dir=args.p_vntrs_log_dir,
                        flanking_threshold=args.gc_flanking_threshold)
        plot_genotypes(genotypes_df, alleles_df, vntrs, haplotypes,
                        filename_prefix="disease_vntr",
                        title_prefix="Disease associated VNTRs",
                        vntr_id_column="vntr_gene_name",
                        plot_median_based_plots=plot_median_based_plots,
                        plot_swarm=plot_swarm)
    if plot_genic_vntrs:
        # Graphs for Genic VNTRs
        print("Working with gene proximal VNTRs")
        genotypes_df, alleles_df, vntrs, haplotypes, samples = get_genic_genotype_distributions(
                    long_vntr_db=args.g_vntrs_long_db,
                    short_vntr_db=args.g_vntrs_short_db,
                    g_vntrs_log_dir=args.g_vntrs_log_dir,
                    strs_filename=args.ids_to_discard_file,
                    flanking_threshold=args.gc_flanking_threshold,
                    verbose=False)
        vntrs_genotyped = alleles_df.dropna(subset=["median_rc"])
        print("len alleles df {} len vntrs genotyped {}".format(
            len(alleles_df), len(vntrs_genotyped)))
        print("Num VNTRs genotyped {}. Num VNTRs with RC <= 40 {} RC<=20 {} with RC<=10 {}".format(
            len(vntrs_genotyped),
            len(alleles_df[alleles_df["median_rc"] <= 40]),
            len(alleles_df[alleles_df["median_rc"] <= 20]),
            len(alleles_df[alleles_df["median_rc"] <= 10]),
            ))
        #print(alleles_df["median_rc"].describe())

        plot_genotypes(genotypes_df, alleles_df, vntrs, haplotypes,
                        filename_prefix="genic_vntr",
                        title_prefix="Gene proximal VNTRs",
                        vntr_id_column="vntr_id",
                        plot_median_based_plots=plot_median_based_plots,
                        plot_swarm=plot_swarm)

if __name__ == "__main__":
    main()
