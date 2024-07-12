import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from collections import defaultdict
from matplotlib.ticker import PercentFormatter

from advntr.models import load_unique_vntrs_data

def read_str_ids():
    filename="../vntrs_to_discard/str_ids.txt"
    with open(filename, "r") as str_file:
        str_ids = [int(line) for line in str_file.readlines()]
    return str_ids

def get_id_to_vntr_map():
    str_ids = read_str_ids()
    id_to_vntr_map = {}
    long_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.db"
    vntrs = load_unique_vntrs_data(long_vntr_db)
    for vntr in vntrs:
        if int(vntr.id) in str_ids:
            continue
        id_to_vntr_map[vntr.id] = vntr
    print("len of map {} after adding long vntrs".format(len(id_to_vntr_map)))
    short_vntr_db="/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.db"
    vntrs = load_unique_vntrs_data(short_vntr_db)
    for vntr in vntrs:
        if int(vntr.id) in str_ids:
            continue
        id_to_vntr_map[vntr.id] = vntr
    print("len of map {} after adding short vntrs".format(len(id_to_vntr_map)))
    combined_trf_vntr_db = "/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/combined_trf_hg38/hg38_VNTRs_by_TRF.db"
    vntrs = load_unique_vntrs_data(combined_trf_vntr_db)
    for vntr in vntrs:
        if int(vntr.id) in str_ids:
            continue
        id_to_vntr_map[vntr.id] = vntr
    print("len of map {} after adding combined vntrs (with overlap with both previous VNTRs)".format(len(id_to_vntr_map)))
    return id_to_vntr_map


# Target set can be either phenotype-associated or 10k-genic
def get_target_vntrs(id_to_vntr_map, target_set):
    str_ids = read_str_ids()
    target_vntrs_filenames = {
    "phenotype associated VNTRs": \
        "/nucleus/projects/saraj/vntr/sources/pangenome_project/target_vntrs/disease_associated_all_vntrs_v5_unique.txt",
    "10k short genic VNTRs":\
        "/nucleus/projects/saraj/vntr/sources/pangenome_project/target_vntrs/10k_short.txt",
    "10k long genic VNTRs":\
        "/nucleus/projects/saraj/vntr/sources/pangenome_project/target_vntrs/10k_long.txt"
        }
    target_vids = []
    if target_set == "phenotype-associated":
        with open(target_vntrs_filenames["phenotype associated VNTRs"], "r") as target_file:
            vids = target_file.readlines()
            vids = [int(vid) for vid in vids]
        target_vids += vids
        target_vids = list(set(target_vids))
        print("Returning target_vntrs of len ", len(target_vids))
        return target_vids
    elif target_set == "10k-genic":
        with open(target_vntrs_filenames["10k short genic VNTRs"], "r") as target_file:
            s_vids = target_file.readlines()
            s_vids = [int(vid) for vid in s_vids if int(vid) not in str_ids]
        target_vids += s_vids
        with open(target_vntrs_filenames["10k long genic VNTRs"], "r") as target_file:
            l_vids = target_file.readlines()
            l_vids = [int(vid) for vid in l_vids if int(vid) not in str_ids]
        target_vids += l_vids
        target_vids = list(set(target_vids))
        print("Returning target_vntrs of len ", len(target_vids))
        return target_vids

'''def add_trendline_plot_deprecated(dataframe, x, y):
    # Add a trend line for the scatter plot.
    bin_width_trendline = 1
    binned_x_label = x + "_binned"
    bins = np.arange(0, dataframe[x].max() + bin_width_trendline, bin_width_trendline)
    labels = [bin_value + bin_width_trendline/2.0 for bin_value in bins][:-1]
    dataframe[binned_x_label] = pd.cut(dataframe[x], bins=bins, labels=labels, include_lowest=True)
    q = dataframe[y].quantile(0.9)
    dataframe = dataframe[dataframe[y] < q]
    # Cut functions uses the interval (min, max] i.e, left exclusive and right inclusive.
    # Except for the first bin if include_lowest is True
    filtered_labels = []
    filtered_bin_mins = []
    y_trend_line = []
    y_trend_line_interval_based = []
    for label_value in labels:
        current_ys = dataframe[dataframe[binned_x_label] == label_value][y]
        if len(current_ys) > 0:
            current_y = current_ys.median()
            y_trend_line.append(current_y)
            filtered_labels.append(label_value)
            filtered_bin_mins.append(label_value - bin_width_trendline/2.0)
        else:
            pass
    for bin_start in bins[:-1]:
        bin_end = bin_start + bin_width_trendline
        y_trend_line_interval_based.append(
            dataframe[(dataframe[x] > bin_start) & \
                      (dataframe[x] <= bin_end)][y].median())
    rolling_window_size = 5
    weights = np.ones(rolling_window_size)/float(rolling_window_size)
    y_rolling_avg = np.convolve(y_trend_line, weights, mode="same")
    sns.lineplot(x=filtered_bin_mins,
                 y=y_rolling_avg,
                 linewidth=1,
                 color='black'
    )

def set_hue_deprecated(literature_df, highlighted_genes):
    # Set specific hue colors
    literature_df["Gene"] = "Other"
    for highlighted_gene in highlighted_genes:
        mask = literature_df["Gene, loci"].str.startswith(highlighted_gene)
        literature_df.loc[mask, "Gene"] = highlighted_gene

def plot_scatter_deprecated(add_trendline=False):
    plt.clf()
    sns.set_style("whitegrid")
    # ToDo: Decide to go with VNTR len by AdVNTR or TRF and set length accordingly.
    # Right now for histograms AdVNTR VNTR len is used and for scatter plot TRF is used.
    # Remove non-completed entries
    literature_df = pd.read_csv("Resources - Disease_phenotype associated TRs.tsv", sep="\t").dropna(
        subset=["Total Length", "Reference motif count"])
    # Remove STR entries
    literature_df = literature_df[literature_df["Class"] != "Short Tandem repeat - skipping"]
    # Clean up the total length column
    literature_df["Total Length"] = literature_df["Total Length"].apply(
                lambda x: int(x.replace("~", "").replace("42-48", "42")))
    # Set specific hue colors
    highlighted_genes = ["INS", "DUX4", "G4", "CEL", "MUC1", "ACAN"]
    set_hue(literature_df, highlighted_genes)
    sns.scatterplot(data=literature_df,
                    x="Year (have to check again)",
                    y="Total Length",
                    hue="Gene",
                    hue_order=["Other"] + highlighted_genes
                    )
    if add_trendline:
        add_trendline_plot(literature_df,
                      x="Year (have to check again)",
                      y="Total Length")
    if False:
        # Add trendline
        rolling_window_size = 7
        weights = np.ones(rolling_window_size)/float(rolling_window_size)
        q = literature_df["Total Length"].quantile(0.9)
        x_no_outliers = literature_df[literature_df["Total Length"] < q]["Year (have to check again)"]
        y_no_outliers = literature_df[literature_df["Total Length"] < q]["Total Length"]
        y_rolling_avg = np.convolve(y_no_outliers, weights, mode="same")
        sns.lineplot(x=x_no_outliers,
                     y=y_rolling_avg,
                     linewidth=1,
                     color='black'
        )

    plt.savefig("vntr_year_len.pdf", format="pdf", bbox_inches="tight")'''

def plot_year():
    plt.clf()
    color = sns.color_palette("colorblind")[0]
    sns.set_style("whitegrid")
    # Remove non-completed entries
    literature_df = pd.read_csv("Resources - Disease_phenotype associated TRs.tsv", sep="\t").dropna(
        subset=["Total Length", "Reference motif count"])
    # Remove STR entries
    literature_df = literature_df[literature_df["Class"] != "Short Tandem repeat - skipping"]

    # Plot the histogram
    x = literature_df["Year (have to check again)"]
    plt.hist(x,
                cumulative=True,
                color=color,
                bins=20)
    # Set y tick labels as percentage
    ax = plt.gca()
    num_ticks = 10
    num_items = len(x)
    steps = int(num_items / num_ticks)
    yticks = list(range(0, num_items, steps))
    ylabels = [str(round(tick/num_items*100)) + "%" for tick in yticks]
    ax.set(yticks=yticks, yticklabels=ylabels)

    # Set xticks with orientation
    plt.xticks(rotation=45)

    # Set axis labels and figure title
    ax.set_xlabel("Year association reported")
    ax.set_ylabel("Cumulative counts (%)")
    # Save figure
    plt.savefig("vntr_year.pdf", format="pdf", bbox_inches="tight")

def plot_non_cumulative_histogram(data, x, xlabel):
    print("len data before any filtering", len(data))
    data = data[data["motif_len"] == 84]
    print("len data with motif len 84", len(data))
    data = data[data["chromosome"] == "chr19"]
    print("len data in chr 19", len(data))
    plt.hist(data[x],
             bins=100,
             label=xlabel)
    plt.xticks(fontsize=10, rotation=90)
    #for idx, row in data.iterrows():
    #    print(row["motif"])

def plot_histogram(x, xlabel, vline_x, is_capped, is_gene_proximal):
    # Color the histogram based on the target vntrs.
    if is_gene_proximal:
        legend_label = "Gene proximal\nVNTRs (n={:,})".format(len(x))
        color = sns.color_palette("colorblind")[-1]
    else:
        legend_label = "Phenotype associated\nVNTRs (n={})".format(len(x))
        color = sns.color_palette("colorblind")[0]
    if "motif" in xlabel:
        bins = 40
    else:
        bins = len(x)
    # Compute the histogram in reverse cumulative (percentage)
    hist, bin_edges = np.histogram(x, bins=bins)
    print("x: {} y: {} percent: {} out of total {}".format(
        vline_x,
        len([value for value in x if value >= vline_x]),
        len([value for value in x if value >= vline_x]) / float(len(x)) * 100.0,
        len(x)))
    # Adjusting bin values to be at the center of each bin
    bins = []
    for idx in range(len(bin_edges)-1):
        bins.append((bin_edges[idx] + bin_edges[idx+1])/2.0)

    # Compute the reverse cumulative hist
    total_values = sum(hist)
    assert(total_values == len(x))
    reverse_cumulative_hist = []
    for idx in range(len(hist)):
        reverse_cumulative_hist.append(sum(hist[idx:]) / float(total_values) * 100)

    plt.plot(bins,
                 reverse_cumulative_hist,
                 color=color,
                 label=legend_label)
    '''plt.hist(x,
                 bins=40,
                 cumulative=-1,
                 histtype="step",
                 color=color,
                 label=legend_label,
                 #stat="density",
                 )'''
    '''
    # Set y tick labels as percentage
    ax = plt.gca()
    num_ticks = 10
    num_items = len(x)
    #print("num_items", num_items)
    steps = int(num_items / num_ticks)
    #print("steps", steps)
    yticks = list(range(0, num_items + 1, steps))
    #print("yticks", yticks)
    ylabels = [str(round(tick/num_items*100)) + "%" for tick in yticks]
    ax.set(yticks=yticks, yticklabels=ylabels)
    #ax.set_yticklabels(ylabels)
    '''
    # Put vertical indicating line
    ax = plt.gca()
    ymin, ymax = ax.get_ylim()
    if not is_gene_proximal:
        ax.vlines(x=vline_x,
                  ymin=0,
                  ymax=100,
                  colors='black',
                  linestyles="dashed")
        xticks, xticklabels = plt.xticks()
        xticks = [int(tick) for tick in xticks if tick >= x.min() and tick <= x.max()] + [0]
        xticks = sorted([xtick for xtick in xticks if abs(xtick-vline_x) > 0.2*vline_x] + [vline_x])
        xticklabels = [str(tick) for tick in xticks]
        if is_capped:
            xticklabels[-1] = ">={}".format(xticklabels[-1])

        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, rotation=45)

    """# Make sure the xtick label is present for the vertical line
    xticks, xticklabels = plt.xticks()
    # Remove extra ticks
    xticks = [0] + [tick for tick in xticks if tick >= x.min() and tick <= x.max()]
    xticklabels = [str(int(xtick)) for xtick in xticks]
    # Change xtick label to represent >=2000bp in the last bin of histogram.
    if vline_x not in xticks:
        # Remove all labels very close to the threshold label, to avoid labell overlapping.
        xticks = sorted([xtick for xtick in xticks if abs(xtick-vline_x) > 0.2*vline_x] + [vline_x])
        #xticks = sorted(list(xticks) + [vline_x])[1:-1]
        xticklabels = [str(int(xtick)) for xtick in xticks]
    if is_capped:
        xticklabels[-1] = ">={}".format(xticklabels[-1])

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=45)"""

    # Set axis labels and figure title
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Reverse cumulative counts (%)")


def main():
    # Increase the font in all plots
    plt.rcParams.update({'font.size': 16})
    sns.set_style("whitegrid")

    # Plot year associated.
    #plot_year()
    # Get VNTR DB info.
    id_to_vntr_map = get_id_to_vntr_map()

    target_set_to_vids = defaultdict(list)
    all_vids = []
    # Get list of target VNTRs.
    for target_set in ["10k-genic", "phenotype-associated"]:
        target_vids = get_target_vntrs(id_to_vntr_map=id_to_vntr_map,
                                        target_set=target_set)
        target_vids = list(set(target_vids))
        target_set_to_vids[target_set] = target_vids
        all_vids.extend(target_vids)
    all_vids=list(set(all_vids))
    print("len(all_vids): ", len(all_vids))
    vntrs_df = pd.DataFrame(columns=["motif_len", "vntr_len", "vntr_len_max_2000"], index=all_vids, dtype=int)
    for vid in all_vids:
        vntrs_df.at[vid, "motif_len"] = len(id_to_vntr_map[vid].pattern)
        vntrs_df.at[vid, "motif_len_max_100"] = min(100, len(id_to_vntr_map[vid].pattern))
        vntrs_df.at[vid, "vntr_len"] = id_to_vntr_map[vid].get_length()
        vntrs_df.at[vid, "vntr_len_max_2000"] = min(2000, id_to_vntr_map[vid].get_length())
        vntrs_df.at[vid, "chromosome"] = id_to_vntr_map[vid].chromosome
        vntrs_df.at[vid, "start_position"] = id_to_vntr_map[vid].start_point
        vntrs_df.at[vid, "motif"] = id_to_vntr_map[vid].pattern
        # Double check the motif len and vntr len are computed correctly.
        #for test_vid in [290964, 826396, 731247]: # ACAN, NOS3 (ecNOS), MUC21
        #    print("VNTR with vid {} has motif len {} and vntr len {}".format(
        #           test_vid,
        #           vntrs_df.at[test_vid, "motif_len"],
        #           vntrs_df.at[test_vid, "vntr_len"]))
    print("vntr length describe: ", vntrs_df["vntr_len"].describe())
    # Plot VNTR lengh histogram.
    plt.clf()
    for target_set in ["10k-genic", "phenotype-associated"]:
        vids = target_set_to_vids[target_set]
        is_gene_proximal = (target_set == "10k-genic")
        print("Plotting vntr_len_max_2000")
        plot_histogram(x=vntrs_df.loc[vids]["vntr_len_max_2000"],
                       xlabel="VNTR length (bp)",
                       vline_x=150,
                       is_capped=True,
                       is_gene_proximal=is_gene_proximal,
                       #filename=target_set + "_vntrs_total_len_max_2000.pdf"
                       )
    plt.legend(loc="upper right")
    # Save figure
    plt.savefig("vntrs_total_len_max_2000.pdf", format="pdf", bbox_inches="tight")
    # Plot VNTR motif length histogram.
    plt.clf()
    vids = target_set_to_vids["phenotype-associated"]

    for target_set in ["10k-genic", "phenotype-associated"]:
        vids = target_set_to_vids[target_set]
        is_gene_proximal = (target_set == "10k-genic")
        print("Plotting motif_len_max_100")
        plot_histogram(x=vntrs_df.loc[vids]["motif_len_max_100"],
                       vline_x=20,
                       xlabel="VNTR motif length (bp)",
                       is_capped=False,
                       is_gene_proximal=is_gene_proximal,
                       #filename=target_set + "_vntrs_motif_len.pdf"
                       )
    plt.legend(loc="upper right")
    # Save figure
    plt.savefig("vntrs_motif_len.pdf", format="pdf", bbox_inches="tight")
    plt.clf()
    # Study the strange case of G-VNTRs with motif len = 84
    vids = target_set_to_vids["10k-genic"]
    plot_non_cumulative_histogram(data=vntrs_df.loc[vids],
                    x="start_position",
                   xlabel="VNTR start position",
                   )
    plt.savefig("vntr_positions.pdf", format="pdf", bbox_inches="tight")
    plt.clf()


if __name__ == "__main__":
    main()
