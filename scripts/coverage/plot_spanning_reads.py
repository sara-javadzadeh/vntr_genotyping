import os
from os.path import join
import subprocess
import argparse
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



def get_gc_contents(ref_vntr, gc_flanking_threshold, round_to_int=True, verbose=False):
    """ return GC contents of vntr motif and flanking region"""
    left_flanking_region = ref_vntr.left_flanking_region[-1*gc_flanking_threshold:]
    right_flanking_region = ref_vntr.right_flanking_region[:gc_flanking_threshold]
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

def get_id_to_vntr_map(long_vntr_db, short_vntr_db):
    id_to_vntr_map = {}
    vntrs = load_unique_vntrs_data(long_vntr_db)
    for vntr in vntrs:
        id_to_vntr_map[vntr.id] = vntr
    #print("len of map {} after adding long vntrs".format(len(id_to_vntr_map)))
    vntrs = load_unique_vntrs_data(short_vntr_db)
    for vntr in vntrs:
        id_to_vntr_map[vntr.id] = vntr
    #print("len of map {} after adding short vntrs".format(len(id_to_vntr_map)))
    return id_to_vntr_map

def get_sample_bam_files(samples_dir):
    files = os.listdir(samples_dir)
    sample_names = [sample for sample in files if sample.startswith("HG")]
    if len(sample_names) == 0:
        print("Error in loading the input bam files")
        print("For HPRC, --samples-dir should be the parent directory with one directory per sample.")
        print("For targeted samples, --samples-dir should point to the directory where all bam files are located.")
        exit(1)
    print("Number of samples {}".format(len(sample_names)))
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

def plot_gc_quadrant_scatter(dataframe, reads_ratio, is_log):
    plt.clf()

    # Create figure base and adjust ratios
    fig = plt.figure(figsize=(6, 6))
    gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
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



def plot_displot(dataframe, ax, reads_ratio, is_log,
                 gc_content_threshold=60,
                 spanning_reads_threshold = 4):
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
    print("Num hi GC hi coverage ", hi_gc_hi_coverage)
    print("Num low GC hi coverage ", low_gc_hi_coverage)
    print("Num hi GC low coverage ", hi_gc_low_coverage)
    print("Num low GC low coverage ", low_gc_low_coverage)


def plot_spanning_reads(dataframe, sample_ids, spanning_reads_threshold=4, plot_cov=False):
    plt.clf()
    fig, ax = plt.subplots(nrows=1,
                           ncols=2,
                           sharey=True,
                           gridspec_kw={'width_ratios': [3, 1]})
    # Adjust spacing between the two subplots
    plt.subplots_adjust(wspace=0.2)


    plot_heatmap(dataframe=dataframe,
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

def plot_scatter(dataframe, column, xlabel, ax, verbose=False):
    # For each VNTR, the row number is the y value in this plot.
    # Have to reverse the order because the dataframe is plotted from the top first,
    # But the sns plot, plots the first value at the bottom.
    y = list(range(len(dataframe.index)))
    # Shift the y axis so each value falls in the middle of each VNTR bin.
    y = [y_value + 0.5 for y_value in y]
    sns.scatterplot(y=y,
                 s=2,
                 legend=None,
                 x=dataframe[column],
                 linewidth=0, # Removes the edge color
                 ax=ax)
    # Put a vertical line for the median value
    median_of_column = dataframe[column].median()
    ymin, ymax = ax.get_ylim()
    if verbose:
        print("Median of column {} is {}".format(column, median_of_column))
    high_sr = dataframe[dataframe["median_spanning_reads"] >= 4]
    if verbose:
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
    ax.tick_params(axis='x', rotation=90, labelsize=6)
    ax.tick_params(left=False)
    # Set xlabel
    ax.set_xlabel(xlabel)


def plot_heatmap(dataframe, sample_ids, ax, spanning_reads_threshold, verbose=False):
    # Adjust colorbar ticks
    cbar_ticks = list(np.arange(0, 5, 1))
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
    hline_y_value = get_y_value_for_preset_threshold(dataframe, spanning_reads_threshold,
                                                     column="median_spanning_reads")
    if verbose:
        print("for median spanning_reads_threshold: {} y_value {}".format(spanning_reads_threshold, hline_y_value))
        print("SR for first row {} middle row {} last row {}".format(
                dataframe["median_spanning_reads"].iloc[0],
                dataframe["median_spanning_reads"].iloc[5000],
                dataframe["median_spanning_reads"].iloc[-1]))
    xmin, xmax = heatmap.get_xlim()
    if verbose:
        print("xmin, xmax, hline_y_value: ", xmin, xmax, hline_y_value)
    ax.hlines(y=hline_y_value, xmin=xmin, xmax=xmax, color="black", linewidth=1)
    ax.text(xmin, hline_y_value - 40, 'High (>=15) SRs', color="black", size=10)
    # Not enough space for displaying this text
    #ax.text(xmin, hline_y_value + 2, 'low (<15) mean spanning reads', color="white")

    # Add hline for 90th quantile.
    hline_y_value = get_y_value_for_preset_threshold(dataframe, spanning_reads_threshold,
                                                     column="90th_quantile_spanning_reads")
    if verbose:
        print("for 10th quantile spanning_reads_threshold: {} y_value {}".format(
                spanning_reads_threshold, hline_y_value))

    # Adjust colorbar
    cbar_ax = ax.collections[0].colorbar.ax
    cbar_ax.yaxis.label.set_size(6)
    cbar_ax.tick_params(axis="y", labelsize=6)

def get_sample_id_to_overlapping_reads_ratio(overlap_filename, all_reads_filename):
    # Return number of reads overlapping any target VNTR divided by total number of reads.
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
    all_reads_file.close()
    return sample_id_to_reads_map

def discard_strs(vids, data, filename):
    str_ids = []
    with open(filename, "r") as str_file:
        str_ids = str_file.readlines()
    # VNTR with id 983766 is not present in the csv file to begin with.
    str_ids = [float(single_id) for single_id in str_ids if str(single_id).strip() != "983766"]
    #print("vids len before discarding STRs ", len(vids))
    vids = [vid for vid in vids if vid not in str_ids]
    data = data.drop(index=str_ids)
    #print("vids len after discarding STRs ", len(vids))
    return vids, data

def load_and_plot_spanning_reads(num_samples_included,
                                num_vntrs_included,
                                df_filename,
                                samples_dir,
                                long_vntr_db,
                                short_vntr_db,
                                save_df,
                                load_df,
                                overlap_filename,
                                gc_flanking_threshold,
                                all_reads_filename,
                                plotting_heatmap,
                                reference_filename,
                                strs_filename,
                                verbose=False):

    print("Loading dataframe from file ", df_filename)

    id_to_vntr_map = get_id_to_vntr_map(long_vntr_db, short_vntr_db)
    vids = sorted(list(id_to_vntr_map.keys()))
    sample_ids, sample_id_to_bam_file = get_sample_bam_files(samples_dir)
    # Choose num_samples_included random samples
    random.seed(0)
    sample_ids = random.sample(sample_ids, num_samples_included)
    num_samples = len(sample_ids)
    vids = vids[:num_vntrs_included]

    # Populate the dataframe
    if load_df:
        spanning_reads_df = pd.read_csv(df_filename, index_col="vid")
        spanning_reads_df.drop("Unnamed: 0", axis=1, inplace=True)
    else:
        spanning_reads_df = DataFrame(columns=['vid'] + sample_ids, index=vids, dtype=float)
        for sample_id in sample_ids:
            print("Working with sample ", sample_id)
            for vntr_id in vids:
                vntr = id_to_vntr_map[vntr_id]
                sample_bam_file = sample_id_to_bam_file[sample_id]
                alignment_file = pysam.AlignmentFile(sample_bam_file, 'rb',
                                    reference_filename=reference_filename)
                start = vntr.start_point
                end = vntr.start_point + vntr.get_length()
                sr = count_spanning_reads(vntr.chromosome, start, end, alignment_file)
                spanning_reads_df.at[vntr_id, 'vid'] = vntr_id
                spanning_reads_df.at[vntr_id, sample_id] = sr

        if save_df:
            # Save the dataframe
            spanning_reads_df.to_csv(df_filename)

    # Discard STRs
    vids, spanning_reads_df = discard_strs(vids=vids, data=spanning_reads_df, filename=strs_filename)
    # Store spanning reads in log scale
    spanning_reads_df = spanning_reads_df[sample_ids].applymap(lambda x: log(float(x) + 1, 2))
    # Set GC content
    for vntr_id in vids:
        vntr = id_to_vntr_map[vntr_id]
        left_flanking_GC, right_flanking_GC, vntr_motif_GC, whole_region_GC = \
            get_gc_contents(vntr, gc_flanking_threshold)
        spanning_reads_df.at[vntr_id, 'gc_content'] = whole_region_GC

    print("Dataframe shape before adding helper columns: ", spanning_reads_df.shape)
    # Add helper columns
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
    entropy_file = open("data/sr_for_entropy_wgs.txt", "w")

    # Adding information for entropy column
    for idx, row in spanning_reads_df.iterrows():
        binned_row, bins = pd.cut(row[sample_ids], right=True, bins=10, retbins=True)
        p_x = list(binned_row.value_counts(normalize=True))
        spanning_reads_df.at[idx, "entropy across samples"] = float(entropy(p_x, base=2))
        # Writing down entropy to read it later on on another script
        entropy_file.write(str(idx) + " " + row[sample_ids].to_string().replace("\n", " ") + "\n")
        continue
    entropy_file.close()
    print("Dataframe shape after adding helper columns: ", spanning_reads_df.shape)
    # Sort the dataframe
    spanning_reads_df = spanning_reads_df.sort_values(by="mean_spanning_reads", ascending=False)

    # Plot GC quadrant figure
    if not plotting_heatmap:
        sample_id_to_reads_map = get_sample_id_to_overlapping_reads_ratio(overlap_filename, all_reads_filename)
        #print("sample_id_to_reads_map ", sample_id_to_reads_map)
        reads_ratio =  np.median(list(sample_id_to_reads_map.values())) * 100
        #print("Median value of overlapping reads ratio ", reads_ratio)
        plot_gc_quadrant_scatter(spanning_reads_df, reads_ratio=reads_ratio, is_log=True)

    # Plot the heatmap
    if plotting_heatmap:
        plot_spanning_reads(spanning_reads_df, sample_ids,
                    spanning_reads_threshold=4)

    print("Median of median spanning reads: ", spanning_reads_df["median across samples"].median())
    spanning_reads_df_high_coverage = spanning_reads_df[spanning_reads_df["median across samples"]>=4]
    print("num vntrs with spanning_reads >= 15: ", len(spanning_reads_df_high_coverage))
    spanning_reads_df_higher_coverage = spanning_reads_df[spanning_reads_df["median across samples"]>6]
    print("num vntrs with spanning_reads > 63: ", len(spanning_reads_df_higher_coverage))
    spanning_reads_df_low_coverage = spanning_reads_df[spanning_reads_df["median across samples"]<4]
    print("num vntrs with spanning_reads < 15: ", len(spanning_reads_df_low_coverage))


def parse_arguments():
    parser = argparse.ArgumentParser(
                prog="plot_spanning_reads",
                description="Plot the heatmap and the GC quadrant figures.")

    parser.add_argument("--samples-dir", help="Path to the directory where sample BAM files are located.",
                        type=str,
                        required=True)
    parser.add_argument("--p-vntrs-db", help="Path to the db file for p-vntrs.",
                        type=str,
                        default="../../target_vntrs/p_vntrs_g_vntrs.db")
    parser.add_argument("--g-vntrs-short-db", help="Path to the db file for short g-vntrs.",
                        type=str,
                        default="../../target_vntrs/illumina_probed_short_vntrs.db")
    parser.add_argument("--g-vntrs-long-db", help="Path to the db file for long g-vntrs.",
                        type=str,
                        default="../../target_vntrs/pacbio_probed_long_vntrs.db")
    parser.add_argument("--df-filename", help="Path to the csv file where the precomputed dataframe is stored.",
                        type=str,
                        default="data/spanning_reads_df_wgs.csv")
    parser.add_argument("--ids-to-discard-file", help="Path to the text file including ids to be discarded (e.g. test strs).",
                        type=str,
                        default="../../target_vntrs/str_ids_to_discard.txt")
    parser.add_argument("--num-samples", help="Number of samples to be included.",
                        type=int, default=28)
    parser.add_argument("--num-vntrs", help="Number of VNTRs to be included. " + \
                        "This is helpful to run a test on a subset of VNTRs to get a speedup in running. " + \
                        "Defaults to all VNTRs.",
                        type=int, default=-1)
    parser.add_argument("--gc-flanking-threshold", help="How many bps in each flanking side " + \
                        "should be included for GC computation? ",
                        type=int, default=100)
    parser.add_argument("--save-df", help="Saves the dataframe to be used later by overriding the current csv file. ",
                        action="store_true", default=False)
    parser.add_argument("--load-df", help="Loads the dataframe from a precomputed csv file. True by default. " + \
                        "Warning: Setting to False and recomputing the dataframe takes a few minutes.",
                        action="store_false", default=True)
    parser.add_argument("--reference-filename", help="Path to the GRCh38 reference fasta file.",
                        type=str,
                        required=True)
    parser.add_argument("--overlap-filename", help="Path to a file with the number of reads " + \
                        "in bamfiles overlapping the VNTR db.",
                        type=str,
                        default="data/num_reads_in_bamfiles_overlapping_db.txt")
    parser.add_argument("--all-reads-filename", help="Path to a file with the number of reads in bam files.",
                        type=str,
                        default="data/num_reads_in_bamfiles.txt")

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    # Set flags and values for this run.
    num_samples_included = args.num_samples
    num_vntrs_included = args.num_vntrs
    save_df = args.save_df
    load_df = args.load_df
    plot_heatmap = True

    if not os.path.exists("figures"):
        os.mkdir("figures")

    # Plot heatmap
    if plot_heatmap:
        load_and_plot_spanning_reads(
            num_samples_included=num_samples_included,
            num_vntrs_included=num_vntrs_included,
            df_filename=args.df_filename,
            save_df=save_df,
            load_df=load_df,
            gc_flanking_threshold=args.gc_flanking_threshold,
            overlap_filename=args.overlap_filename,
            all_reads_filename=args.all_reads_filename,
            samples_dir=args.samples_dir,
            long_vntr_db=args.g_vntrs_long_db,
            short_vntr_db=args.g_vntrs_short_db,
            plotting_heatmap=True,
            reference_filename=args.reference_filename,
            strs_filename=args.ids_to_discard_file)
    # Plot GC quadrant figure
    load_and_plot_spanning_reads(
        num_samples_included=num_samples_included,
        num_vntrs_included=num_vntrs_included,
        df_filename=args.df_filename,
        save_df=save_df,
        load_df=load_df,
        gc_flanking_threshold=args.gc_flanking_threshold,
        overlap_filename=args.overlap_filename,
        all_reads_filename=args.all_reads_filename,
        samples_dir=args.samples_dir,
        long_vntr_db=args.g_vntrs_long_db,
        short_vntr_db=args.g_vntrs_short_db,
        plotting_heatmap=False,
        reference_filename=args.reference_filename,
        strs_filename=args.ids_to_discard_file)


if __name__ == "__main__":
    main()
