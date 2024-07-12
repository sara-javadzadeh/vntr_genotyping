import os
from os.path import join
import subprocess
from functools import reduce

from pandas import DataFrame
import pandas as pd
from scipy.stats import entropy
from collections import Counter
from math import log
import numpy as np
import random
import math
import seaborn as sns
from matplotlib import pyplot as plt

import pysam
from Bio.SeqUtils import GC
from advntr.models import load_unique_vntrs_data



def get_gc_contents(ref_vntr, round_to_int=True, verbose=False):
    """ return GC contents of vntr motif and flanking region"""
    # flanking 100 bps
    left_flanking_region = ref_vntr.left_flanking_region[-100:]
    right_flanking_region = ref_vntr.right_flanking_region[:100]
    # VNTR sequence in reference model
    vntr_seq = "".join(ref_vntr.repeat_segments)
    whole_region = "".join([left_flanking_region,vntr_seq,right_flanking_region])

    left_flanking_GC = GC(left_flanking_region)
    right_flanking_GC = GC(right_flanking_region)
    vntr_motif_GC = GC(vntr_seq)
    whole_region_GC = GC(whole_region)

    if verbose:
        print("left flanking", left_flanking_region, left_flanking_GC)
        print("right flanking", right_flanking_region, right_flanking_GC)
        print("vntr seq", vntr_seq, vntr_motif_GC)
        print("whole region", whole_region, whole_region_GC)

    if round_to_int:
        return int(round(left_flanking_GC)), int(round(right_flanking_GC)), int(round(vntr_motif_GC)), int(round(whole_region_GC))
    else:
        return left_flanking_GC, right_flanking_GC, vntr_motif_GC, whole_region_GC

def count_spanning_reads(chromosome, start, end, alignment_file):
    spanning_read_count = 0
    min_spanned_flanking_region = 10
    for aligned_read in alignment_file.fetch(chromosome, start, end):
        if aligned_read.reference_start <= start - min_spanned_flanking_region and \
           end + min_spanned_flanking_region <= aligned_read.reference_end:
            spanning_read_count += 1
    return spanning_read_count

def get_id_to_vntr_map_for_combined():
    # This is a much larger DB. Not recommended for normal use.
    id_to_vntr_map = {}
    vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/combined_trf_hg38/hg38_VNTRs_by_TRF.db"
    vntrs = load_unique_vntrs_data(vntr_db)
    for vntr in vntrs:
        id_to_vntr_map[vntr.id] = vntr
    print("len of map {} after adding long vntrs".format(len(id_to_vntr_map)))
    return id_to_vntr_map

def get_id_to_vntr_map():
    id_to_vntr_map = {}
    long_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.db"
    vntrs = load_unique_vntrs_data(long_vntr_db)
    for vntr in vntrs:
        id_to_vntr_map[vntr.id] = vntr
    print("len of map {} after adding long vntrs".format(len(id_to_vntr_map)))
    short_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.db"
    vntrs = load_unique_vntrs_data(short_vntr_db)
    for vntr in vntrs:
        id_to_vntr_map[vntr.id] = vntr
    print("len of map {} after adding short vntrs".format(len(id_to_vntr_map)))
    return id_to_vntr_map

def get_sample_bam_files():
    samples_dir="/ribosome/projects/saraj/vntr/data/pangenome_project/pacbio_hifi_mapped_reads"
    files = os.listdir(samples_dir)
    sample_names = [sample for sample in files if sample.startswith("HG")]
    print("sample_names len {}".format(len(sample_names)))
    sample_name_to_bam_file = {}
    for sample_name in sample_names:
        sample_name_to_bam_file[sample_name] = join(samples_dir, sample_name, sample_name + "_aligned_GRCh38_winnowmap.sorted.bam")
    return sample_names, sample_name_to_bam_file


def get_y_value_for_preset_threshold(dataframe, spanning_reads_threshold, column):
    # Assumes the dataframe is sorted in descending order.
    # Finds the first row that is greater than or equal to the threshold.
    row_counter = 0.0
    for idx, row in dataframe.iterrows():
        if row[column] < spanning_reads_threshold:
            return row_counter
        row_counter += 1

def plot_gc_quadrant_scatter_deprecated(dataframe, is_log):
    plt.clf()
    fig, ax = plt.subplots(nrows=2,
                           ncols=2,
                           sharey='row',
                           sharex='col',
                           gridspec_kw={'width_ratios': [5, 1], 'height_ratios': [5, 1]})
    # Adjust spacing between the two subplots
    plt.subplots_adjust(wspace=0.2)

    # Add title
    #fig.suptitle("Spanning reads for variable GC content of VNTR region")

    # Plot the main density plot
    plot_displot(dataframe=dataframe,
                 ax=ax[0][0],
                 is_log=is_log)

    # Plot side histograms
    plot_histogram(dataframe=dataframe,
                   ax=ax[1][0],
                   axis="x")
    plot_histogram(dataframe=dataframe,
                   ax=ax[0][1],
                   axis="y")
    plt.txt
    # Save figure
    plt.savefig("figures/gc_sr_quadrant_w_histograms.pdf", format="pdf", bbox_inches="tight")

def plot_gc_quadrant_scatter(dataframe, reads_ratio, is_log):
    plt.clf()

    # Create figure base and adjust ratios
    fig = plt.figure(figsize=(6, 6))
    gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    # Add title
    #fig.suptitle("Spanning reads and GC content of VNTR region")

    # Create the axes
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

	# Remove double labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # Plot the main density plot
    plot_displot(dataframe=dataframe,
                 reads_ratio=reads_ratio,
                 ax=ax,
                 is_log=is_log)

    # Plot side histograms
    plot_histogram(dataframe=dataframe,
                   ax=ax_histx,
                   axis="x")
    plot_histogram(dataframe=dataframe,
                   ax=ax_histy,
                   axis="y")
    # Save figure
    plt.savefig("figures/gc_sr_quadrant_w_histograms.pdf", format="pdf", bbox_inches="tight")

def plot_histogram(dataframe, ax, axis):
    color = sns.color_palette("Blues", 6)[3]
    if axis == "x":
        data = dataframe["median_spanning_reads"]
        orientation = "vertical"
        ax.set_ylabel("Count VNTRs", fontsize=8)
        ax.yaxis.label.set_color(color)
        ax.tick_params(axis="y", colors=color)
    if axis == "y":
        data = dataframe["gc_content"]
        orientation = "horizontal"
        ax.set_xlabel("Count VNTRs", fontsize=8)
        ax.xaxis.label.set_color(color)
        ax.tick_params(axis="x", colors=color)
        #ax_cumulative = ax.twinx()
        #ax_cumulative.set_ylabel("Cumulative", fontsize=8)
        #ax_cumulative.yaxis.label.set_color("grey")
        #ax_cumulative.tick_params(axis="y", colors="grey")
    # Use the color palette rocket as a discrete pallete
    # Pick the second to last color (orange) when having 6 discrete levels
    ax.hist(data,
             bins=40,
             color=color,
             orientation=orientation
             )
    if axis == "x":
        # Add a cumulative histogram as a step line (not filled) on the opposite axis.
        color = sns.color_palette("Blues", 6)[5]
        ax_cumulative = ax.twinx()
        ax_cumulative.set_ylabel("Cumulative\ncount VNTRs", fontsize=8)
        ax_cumulative.yaxis.label.set_color(color)
        ax_cumulative.tick_params(axis="y", colors=color)
        ax_cumulative.hist(data,
                 bins=40,
                 cumulative=-1,
                 histtype="step",
                 color=color,
                 orientation=orientation
                 )



def plot_displot(dataframe, ax, reads_ratio, is_log, gc_content_threshold=60):
    spanning_reads_threshold = 4
    x = "median_spanning_reads"
    y = "gc_content"
    close_to_lines = False
    on_the_corners = True
    color = sns.color_palette("Blues", 6)[5]
    # Note: Displot is a figure level functiona nd does not take ax as input
    sns.histplot(data=dataframe,
                    x=x,
                    y=y,
                    color=color,
                    ax=ax,
                    bins=60)
    # Set the quadrant lines and texts
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.vlines(x=spanning_reads_threshold,
              ymin=ymin,
              ymax=ymax,
              color='black',
              linewidth=1)
    ax.vlines(x=spanning_reads_threshold + 2,
              ymin=ymin,
              ymax=ymax,
              linestyles="dashed",
              color='grey',
              linewidth=1)
    ax.hlines(y=gc_content_threshold,
              xmin=xmin,
              xmax=xmax,
              color='black',
              linewidth=1)
    # Find values for each quadrant and write it.
    hi_gc_hi_coverage = len(dataframe[\
        (dataframe[y] >= gc_content_threshold) &\
        (dataframe[x] >= spanning_reads_threshold)])
    low_gc_hi_coverage = len(dataframe[\
        (dataframe[y] < gc_content_threshold) &\
        (dataframe[x] >= spanning_reads_threshold)])
    hi_gc_low_coverage = len(dataframe[\
        (dataframe[y] >= gc_content_threshold) &\
        (dataframe[x] < spanning_reads_threshold)])
    low_gc_low_coverage = len(dataframe[\
        (dataframe[y] < gc_content_threshold) &\
        (dataframe[x] < spanning_reads_threshold)])

    # Two possible setting of in-figure texts: either close to the lines or on the furthest corners
    if close_to_lines:
        # Write the values for each quadrant on the plot
        ax.text(s=hi_gc_hi_coverage, x=spanning_reads_threshold + 3, y=gc_content_threshold + 2)
        ax.text(s=low_gc_hi_coverage, x=spanning_reads_threshold + 3, y=gc_content_threshold - 5)
        ax.text(s=hi_gc_low_coverage, x=spanning_reads_threshold - 1.5, y=gc_content_threshold + 2)
        ax.text(s=low_gc_low_coverage, x=spanning_reads_threshold - 1.5, y=gc_content_threshold - 5)
        # Write the description for each quadrant on the plot
        ax.text(s="Low  GC", x=0.25, y=gc_content_threshold - 5)
        ax.text(s="High GC", x=0.25, y=gc_content_threshold + 2)
        ax.text(s="Low\nspanning\nreads", x=spanning_reads_threshold - 2.5, y=10)
        ax.text(s="High\nspanning\nreads", x=spanning_reads_threshold + 0.5, y=10)
    if on_the_corners:
        # Write the description and values for each quadrant on the plot
        ax.text(s="{:,} VNTRs\nlow GC\nlow SR".format(low_gc_low_coverage), x=0.25, y=5)
        ax.text(s="{:,} VNTRs\nhigh GC\nlow SR".format(hi_gc_low_coverage), x=0.25, y=75)
        ax.text(s="{:,} VNTRs\nlow GC\nhigh SR".format(low_gc_hi_coverage), x=7.5, y=5)
        ax.text(s="{:,} VNTRs\nhigh GC\nhigh SR".format(hi_gc_hi_coverage), x=7.5, y=75)
        ax.text(s="{:.2f}%\ncovering\nreads".format(reads_ratio), x=7.5, y=40)

    # Translate xtick labels to show actual value of spanning reads, instead of log
    xticks = list(range(0, 12, 2))
    xticklabels = [str((2**xtick)-1) for xtick in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    # Set X and Y labels
    if is_log:
        xlabel = "Spanning reads"
    else:
        xlable = "Spanning reads"
    ax.set_xlabel(xlabel)
    ax.set_ylabel("GC content")

    print("Num total VNTRs ", len(dataframe))
    #print("Num hi GC hi coverage ", hi_gc_hi_coverage)
    #print("Num low GC hi coverage ", low_gc_hi_coverage)
    #print("Num hi GC low coverage ", hi_gc_low_coverage)
    #print("Num low GC low coverage ", low_gc_low_coverage)


def plot_spanning_reads(dataframe, sample_ids, normalization_factor, spanning_reads_threshold=4, plot_cov=False):
    plt.clf()
    fig, ax = plt.subplots(nrows=1,
                           ncols=2,
                           sharey=True,
                           #width_ratios=[4,1,1])
                           gridspec_kw={'width_ratios': [3, 1]})
    # Adjust spacing between the two subplots
    plt.subplots_adjust(wspace=0.2)

    # Add title
    #fig.suptitle("Spanning reads")

    # Remove the outliers of spanning reads
    if False:
        print("Dataframe shape before removing quantile: ", dataframe.shape)
        dataframe_no_outlier = dataframe.copy()
        top_quantile = np.quantile(dataframe_no_outlier["median_spanning_reads"], q=0.9)
        dataframe_no_outlier = dataframe_no_outlier[dataframe_no_outlier["median_spanning_reads"] < top_quantile]
        num_vntrs_removed = len(dataframe) - len(dataframe_no_outlier)
        print("Dataframe shape after removing quantile: {} with {} vntrs removed".format(dataframe_no_outlier.shape, num_vntrs_removed))

    plot_heatmap(dataframe=dataframe, #.astype(float).transform(lambda x: log(x, 2)),
                 sample_ids=sample_ids,
                 ax=ax[0],
                 spanning_reads_threshold=spanning_reads_threshold)
    if plot_cov:
        plot_scatter(dataframe=dataframe,
                     column="CoV across samples",
                     xlabel="CoV",
                     ax=ax[1])
    plot_scatter(dataframe=dataframe,
                 column="variance across samples clip 30",
                 xlabel="Variance",
                 ax=ax[1])
    # Save the figure.
    plt.savefig("figures/spanning_reads_heatmap_10k_clip_30.pdf", format="pdf", bbox_inches='tight')

def plot_scatter(dataframe, column, xlabel, ax):
    #print("dataframe.index: ", dataframe.index)
    #print("column: ", dataframe[column])
    # For each VNTR, the row number is the y value in this plot.
    # Have to reverse the order because the dataframe is plotted from the top first,
    # But the sns plot, plots the first value at the bottom.
    y = list(range(len(dataframe.index)))
    # Shift the y axis so each value falls in the middle of each VNTR bin.
    y = [y_value + 0.5 for y_value in y]
    #print("Y for scatter plot: ", y)
    sns.scatterplot(y=y,
                 s=2,
                 legend=None,
                 x=dataframe[column],
                 linewidth=0, # Removes the edge color
                 ax=ax)
    # Put a vertical line for the median value
    median_of_column = dataframe[column].median()
    ymin, ymax = ax.get_ylim()
    print("Median of column {} is {}".format(column, median_of_column))
    high_sr = dataframe[dataframe["median_spanning_reads"] >= 4]
    print("Median of column {} is {:.4f} mean is {:.4f} for >=15 SR count is {} median is {:.4f} mean is {:.4f}".format(
        column, median_of_column,
        dataframe[column].mean(),
        len(high_sr),
        high_sr[column].median(),
        high_sr[column].mean()
        ))
    print("Number of values in column {} < 0.2 is {} and len values {}".format(column,
        len(dataframe[dataframe[column] < 0.2]),
        len(dataframe)))
    ax.vlines(x=median_of_column,
              ymin=ymin,
              ymax=ymax,
              color="grey")
    # Adjust tick labels
    xticks = np.linspace(min(dataframe[column]),
                              max(dataframe[column]),
                              4)
    xticks = list(xticks) + [median_of_column]
    xticks = [round(xtick, 2) for xtick in xticks]
    #print("xticks: ", xticks)
    #ax.set_xticks(xticks)
    ax.tick_params(axis='x', rotation=90, labelsize=6)
    ax.tick_params(left=False)
    #ax.set(yticklabels=[])
    # Set xlabel
    ax.set_xlabel(xlabel)


def plot_heatmap(dataframe, sample_ids, ax, spanning_reads_threshold):
    # Adjust colorbar ticks
    cbar_ticks = [0,1] + list(np.arange(2, 5.5, 0.5))
    cbar_tick_labels = [str(int(2**tick) - 1) for tick in cbar_ticks]
    cbar_label = "Spanning reads (capped at 30)"

    # Plot the heatmap
    heatmap = sns.heatmap(dataframe[sample_ids],
                          ax=ax,
                          yticklabels=False,
                          xticklabels=True,
                          cbar_kws = {"use_gridspec":False,
                                      "location": "left",
                                      'pad': 0.15,
                                      'label': cbar_label,
                                      'ticks': cbar_ticks},
                          vmax=5
                          )
    # Adjust colorbar tick labels
    cbar_axis = heatmap.collections[0].colorbar
    cbar_axis.set_ticklabels(cbar_tick_labels)
    cbar_axis.set_label(cbar_label, fontsize=10)

    # Set labels and titles
    ax.set_xlabel("Samples")
    ax.set_ylabel("{:,} VNTRs".format(len(dataframe)))
    # Adjust fonts
    ax.tick_params(axis='x', labelsize=4)

    # Add hline for median.
    hline_y_value = get_y_value_for_preset_threshold(dataframe, 3,#spanning_reads_threshold,
                                                     column="median_spanning_reads")
    print("for median spanning_reads_threshold: {} y_value {}".format(spanning_reads_threshold, hline_y_value))
    print("SR for first row {} middle row {} last row {}".format(
            dataframe["median_spanning_reads"].iloc[0],
            dataframe["median_spanning_reads"].iloc[5000],
            dataframe["median_spanning_reads"].iloc[-1]))
    xmin, xmax = heatmap.get_xlim()
    print("xmin, xmax, hline_y_value: ", xmin, xmax, hline_y_value)
    ax.hlines(y=hline_y_value, xmin=xmin, xmax=xmax, color="black", linewidth=1)
    ax.text(xmin, hline_y_value - 40, 'High (>=15) SRs', color="black", size=10)
    # Not enough space for displaying this text
    #ax.text(xmin, hline_y_value + 2, 'low (<15) mean spanning reads', color="white")

    # Add hline for 90th quantile.
    hline_y_value = get_y_value_for_preset_threshold(dataframe, spanning_reads_threshold,
                                                     column="90th_quantile_spanning_reads")
    print("for 10th quantile spanning_reads_threshold: {} y_value {}".format(spanning_reads_threshold, hline_y_value))

    # Adjust colorbar
    cbar_ax = ax.collections[0].colorbar.ax
    cbar_ax.yaxis.label.set_size(6)
    cbar_ax.tick_params(axis="y", labelsize=6)

def get_sample_id_to_overlapping_reads_ratio():
    # Return number of reads overlapping any target VNTR divided by total number of reads.
    overlap_filename = "num_reads_in_bamfiles_overlapping_db.txt"
    all_reads_filename = "num_reads_in_bamfiles.txt"
    sample_id_to_reads_map = {}
    # First compute the number of overlapping reads
    overlap_reads_file = open(overlap_filename, "r")
    lines = overlap_reads_file.readlines()
    sample_id = None
    for idx, line in enumerate(lines):
        if idx % 2 == 0:
            # This should be a line with sample id
            assert(line.startswith("HG"))
            sample_id = line.strip()
        elif idx % 2 == 1:
            # This should be a line with the number of reads
            assert(not line.startswith("HG"))
            sample_id_to_reads_map[sample_id] = float(line.strip())
            #print("{}: {}".format(sample_id, sample_id_to_num_reads_map[sample_id]))
    overlap_reads_file.close()
    # Next, divide this number by total reads
    all_reads_file = open(all_reads_filename, "r")
    lines = all_reads_file.readlines()
    sample_id = None
    for idx, line in enumerate(lines):
        if idx % 2 == 0:
            # This should be a line with sample id
            assert(line.startswith("HG"))
            sample_id = line.strip()
        elif idx % 2 == 1:
            # This should be a line with the number of reads
            assert(not line.startswith("HG"))
            sample_id_to_reads_map[sample_id] = sample_id_to_reads_map[sample_id] / float(line.strip())
            #print("{}: {}".format(sample_id, sample_id_to_num_reads_map[sample_id]))
    all_reads_file.close()
    return sample_id_to_reads_map

def get_sample_id_to_million_reads():
    # return total number of reads overlapping VNTR DB
    sample_id_to_num_reads_map = {}
    filename = "num_reads_in_bamfiles_overlapping_db.txt"
    with open(filename, "r") as num_reads_file:
        lines = num_reads_file.readlines()
        sample_id = None
        for idx, line in enumerate(lines):
            if idx % 2 == 0:
                # This should be a line with sample id
                assert(line.startswith("HG"))
                sample_id = line.strip()
            elif idx % 2 == 1:
                # This should be a line with the number of reads
                assert(not line.startswith("HG"))
                sample_id_to_num_reads_map[sample_id] = float(line.strip()) / (100.0 * 1000.0)
                # Temporarily removing the normalization factor
                #sample_id_to_num_reads_map[sample_id] = 1
                #print("{}: {}".format(sample_id, sample_id_to_num_reads_map[sample_id]))
    return sample_id_to_num_reads_map

def discard_strs(vids, data):
    filename = "../vntrs_to_discard/str_ids.txt"
    str_ids = []
    with open(filename, "r") as str_file:
        str_ids = str_file.readlines()
    # VNTR with id 983766 is not present in the csv file to begin with.
    str_ids = [float(single_id) for single_id in str_ids if str(single_id).strip() != "983766"]
    print("vids len before discarding STRs ", len(vids))
    vids = [vid for vid in vids if vid not in str_ids]
    data = data.drop(index=str_ids)
    print("vids len after discarding STRs ", len(vids))
    return vids, data

def load_and_plot_spanning_reads(num_samples_included,
                                num_vntrs_included,
                                save_df,
                                load_df,
                                is_normalized):
    #alignment_file
    #vids = [4237, #PER3
    #        674126, #PRDM9
    #        290964, #ACAN
    #        886134, #C9orf72
    #        ]

    # Set suffix for the csv files to load the dataframe from
    if is_normalized:
        normalized_suffix = "normalized"
        df_filename = "logs/spanning_reads_df_{}_samples_{}_vntrs_{}_no_normalization.csv".format(
                    num_samples_included, num_vntrs_included, normalized_suffix)
    else:
        normalized_suffix = "raw"
        df_filename = "logs/spanning_reads_df_{}_samples_{}_vntrs_{}.csv".format(
                    num_samples_included, num_vntrs_included, normalized_suffix)
    print("Loading dataframe from file ", df_filename)
    reference_filename = "/nucleus/projects/saraj/reference_data/grch38/hg38full.fa"

    # For test only - print GC content of P-VNTRs
    #id_to_vntr_map = get_id_to_vntr_map_for_combined()
    #with open("../../target_vntrs/disease_associated_all_vntrs_v5.txt", "r") as p_vntrs_file:
    #    p_vntr_ids = p_vntrs_file.readlines()
    #p_vntr_ids = [int(p_vntr_id.strip()) for p_vntr_id in p_vntr_ids]
    #for p_vntr_id in p_vntr_ids:
    #    if p_vntr_id in id_to_vntr_map:
    #        print("{} gc: {}".format(p_vntr_id, get_gc_contents(id_to_vntr_map[p_vntr_id])))


    id_to_vntr_map = get_id_to_vntr_map()
    vids = sorted(list(id_to_vntr_map.keys()))
    sample_ids, sample_id_to_bam_file = get_sample_bam_files()
    sample_id_to_num_reads = get_sample_id_to_million_reads()
    # Choose num_samples_included random samples
    random.seed(0)
    sample_ids = random.sample(sample_ids, num_samples_included)
    num_samples = len(sample_ids)
    vids = vids[:num_vntrs_included]

    # Populate the dataframe
    if load_df:
        spanning_reads_df = pd.read_csv(df_filename, index_col="vid")
        spanning_reads_df.drop("Unnamed: 0", axis=1, inplace=True)
        #with open("temp_vids.txt", "w") as vids_file:
        #    for vid in sorted(list(spanning_reads_df.index)):
        #        vids_file.write(str(int(vid)) + '\n')
    else:
        spanning_reads_df = DataFrame(columns=['vid'] + sample_ids, index=vids, dtype=float)
        for sample_id in sample_ids:
            print("Working with sample ", sample_id)
            for vntr_id in vids:
                vntr = id_to_vntr_map[vntr_id]
                sample_bam_file = sample_id_to_bam_file[sample_id]
                alignment_file = pysam.AlignmentFile(sample_bam_file, 'rb', reference_filename=reference_filename)
                start = vntr.start_point
                end = vntr.start_point + vntr.get_length()
                sr = count_spanning_reads(vntr.chromosome, start, end, alignment_file)
                # Temporarily disabling normalization
                #assert(sample_id_to_num_reads[sample_id] == 1)
                normalized_sr = float(sr) / sample_id_to_num_reads[sample_id]
                spanning_reads_df.at[vntr_id, 'vid'] = vntr_id
                if is_normalized:
                    spanning_reads_df.at[vntr_id, sample_id] = normalized_sr
                else:
                    spanning_reads_df.at[vntr_id, sample_id] = sr

        if save_df:
            # Save the dataframe
            spanning_reads_df.to_csv(df_filename)

    # Discard STRs
    vids, spanning_reads_df = discard_strs(vids=vids, data=spanning_reads_df)
    # Store spanning reads in log scale
    spanning_reads_df = spanning_reads_df[sample_ids].applymap(lambda x: log(float(x) + 1, 2))
    # Set GC content
    for vntr_id in vids:
        vntr = id_to_vntr_map[vntr_id]
        left_flanking_GC, right_flanking_GC, vntr_motif_GC, whole_region_GC = get_gc_contents(vntr)
        spanning_reads_df.at[vntr_id, 'gc_content'] = whole_region_GC

    #print("GC content: ", spanning_reads_df['gc_content'])
    print("Dataframe shape before adding helper columns: ", spanning_reads_df.shape)
    # Add helper columns
    #print(spanning_reads_df)
    log2_30 = math.log(30, 2)
    log2_50 = math.log(50, 2)
    spanning_reads_df["mean_spanning_reads"] = spanning_reads_df[sample_ids].mean(axis=1)
    spanning_reads_df["median_spanning_reads"] = spanning_reads_df[sample_ids].median(axis=1)
    spanning_reads_df["90th_quantile_spanning_reads"] = spanning_reads_df[sample_ids].quantile(axis=1, q=0.9)
    spanning_reads_df["median across samples"] = spanning_reads_df[sample_ids].median(axis=1)
    spanning_reads_df["std across samples"] = spanning_reads_df[sample_ids].std(axis=1)
    spanning_reads_df["variance across samples clip 30"] = spanning_reads_df[sample_ids].clip(upper=log2_30).sample(axis=1 ,n=7, random_state=0).var(axis=1).clip(upper=5)
    spanning_reads_df["variance across samples clip 50"] = spanning_reads_df[sample_ids].clip(upper=log2_50).var(axis=1)
    spanning_reads_df["CoV across samples"] = spanning_reads_df[sample_ids].std(axis=1)/spanning_reads_df[sample_ids].median(axis=1)
    entropy_file = open("logs/sr_for_entropy.txt", "w")

    # Adding information for entropy column
    for idx, row in spanning_reads_df.iterrows():
        binned_row, bins = pd.cut(row[sample_ids], right=True, bins=10, retbins=True)
        p_x = list(binned_row.value_counts(normalize=True))
        #print("{}, binned {} , bins {}, px {}".format(
        #    [round(item, 2) for item in list(row[sample_ids])],
        #    list(binned_row),
        #    [round(item, 2) for item in bins],
        #    p_x))
        spanning_reads_df.at[idx, "entropy across samples"] = float(entropy(p_x, base=2))
        # Writing down entropy to read it later on on another script
        entropy_file.write(str(idx) + " " + row[sample_ids].to_string().replace("\n", " ") + "\n")
        continue
        #spanning_reads_df.at[idx, "entropy across samples"] = float(differential_entropy(row[sample_ids]))
        #if idx % 1000 == 0:
        #    random.seed(0)
        #    down_sample = random.sample(row[sample_ids], 7)
        #    print("entropy {} entropy with downsampling {}".format(
        #           spanning_reads_df.at[idx, "entropy across samples"],
        #           differential_entropy(down_sample)))
    entropy_file.close()
    print("Dataframe shape after adding helper columns: ", spanning_reads_df.shape)
    # Sort the dataframe
    spanning_reads_df = spanning_reads_df.sort_values(by="mean_spanning_reads", ascending=False)

    # Plot GC quadrant figure
    if not is_normalized:
        sample_id_to_reads_map = get_sample_id_to_overlapping_reads_ratio()
        #print("sample_id_to_reads_map ", sample_id_to_reads_map)
        reads_ratio =  np.median(list(sample_id_to_reads_map.values())) * 100
        print("Median value of overlapping reads ratio ", reads_ratio)
        plot_gc_quadrant_scatter(spanning_reads_df, reads_ratio=reads_ratio, is_log=True)

    # Plot the heatmap
    if is_normalized:
        #plot_spanning_reads(spanning_reads_df, sample_ids, normalization_factor="100k overlapping", spanning_reads_threshold=4)
        plot_spanning_reads(spanning_reads_df, sample_ids, normalization_factor="100k overlapping", spanning_reads_threshold=4)

    print("Median of median spanning reads: ", spanning_reads_df["median across samples"].median())
    spanning_reads_df_high_coverage = spanning_reads_df[spanning_reads_df["median across samples"]>=4]
    print("num vntrs with spanning_reads >= 15: ", len(spanning_reads_df_high_coverage))
    spanning_reads_df_higher_coverage = spanning_reads_df[spanning_reads_df["median across samples"]>6]
    print("num vntrs with spanning_reads > 63: ", len(spanning_reads_df_higher_coverage))
    spanning_reads_df_low_coverage = spanning_reads_df[spanning_reads_df["median across samples"]<4]
    spanning_reads_df.to_csv("logs/spanning_reads_df.csv")
    print("num vntrs with spanning_reads < 15: ", len(spanning_reads_df_low_coverage))
    #print("len VNTRs with <15 SRs, low GC, distance_to probe <1000: ",
    #        spanning_reads_df[(spanning_reads_df["median across samples"]<4) & \
    #                          (spanning_reads_df["gc_content"]<60) & \
    #                          (spanning_reads_df["probes"]<1000)])

if __name__ == "__main__":
    # Set flags and values for this run.
    num_samples_included = 28
    num_vntrs_included = -1
    save_df = False
    load_df = True
    plot_heatmap = False

    # Plot heatmap
    if plot_heatmap:
        load_and_plot_spanning_reads(
            num_samples_included=num_samples_included,
            num_vntrs_included=num_vntrs_included,
            save_df=save_df,
            load_df=load_df,
            is_normalized=True)
    # Plot GC quadrant figure
    load_and_plot_spanning_reads(
        num_samples_included=num_samples_included,
        num_vntrs_included=num_vntrs_included,
        save_df=save_df,
        load_df=load_df,
        is_normalized=False)

