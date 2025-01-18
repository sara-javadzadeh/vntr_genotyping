import pandas as pd
import numpy as np
import seaborn as sns
import argparse
from matplotlib import pyplot as plt
from collections import defaultdict
from matplotlib.ticker import PercentFormatter

from advntr.models import load_unique_vntrs_data

def read_str_ids(filename):
    with open(filename, "r") as str_file:
        str_ids = [int(line) for line in str_file.readlines()]
    return str_ids

def get_id_to_vntr_map(long_vntr_db, short_vntr_db, p_vntr_db, str_ids_file):
    str_ids = read_str_ids(str_ids_file)
    id_to_vntr_map = {}
    vntrs = load_unique_vntrs_data(long_vntr_db)
    for vntr in vntrs:
        if int(vntr.id) in str_ids:
            continue
        id_to_vntr_map[vntr.id] = vntr
    print("len of map {} after adding long vntrs".format(len(id_to_vntr_map)))
    vntrs = load_unique_vntrs_data(short_vntr_db)
    for vntr in vntrs:
        if int(vntr.id) in str_ids:
            continue
        id_to_vntr_map[vntr.id] = vntr
    print("len of map {} after adding short vntrs".format(len(id_to_vntr_map)))
    combined_trf_vntr_db = p_vntr_db
    vntrs = load_unique_vntrs_data(combined_trf_vntr_db)
    for vntr in vntrs:
        if int(vntr.id) in str_ids:
            continue
        id_to_vntr_map[vntr.id] = vntr
    print("len of map {} after adding combined vntrs (with overlap with both previous VNTRs)".format(len(id_to_vntr_map)))
    return id_to_vntr_map


# Target set can be either phenotype-associated or 10k-genic
def get_target_vntrs(id_to_vntr_map, target_set, p_vntr, g_vntr_short, g_vntr_long, str_ids_file):
    str_ids = read_str_ids(str_ids_file)
    target_vntrs_filenames = {
    "phenotype associated VNTRs": \
        p_vntr,
    "10k short genic VNTRs":\
        g_vntr_short,
    "10k long genic VNTRs":\
        g_vntr_long,
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

def plot_non_cumulative_histogram(data, x, xlabel):
    print("len data before any filtering", len(data))
    data = data[data["motif_len"] == 84]
    print("len data with motif len 84", len(data))
    plt.hist(data[x],
             bins=100,
             label=xlabel)
    plt.xticks(fontsize=10, rotation=90)

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

    # Set axis labels and figure title
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Reverse cumulative counts (%)")


def parse_arguments():
    parser = argparse.ArgumentParser(
                prog="vntr_length_plot",
                description="Plot the vntr lengths.")

    parser.add_argument("--p-vntrs-db", help="Path to the db file for p-vntrs.",
                        type=str,
                        default="../../target_vntrs/p_vntrs_g_vntrs.db")
    parser.add_argument("--p-vntrs-ids", help="Path to the text file including all p-vntr ids.",
                        type=str,
                        default="../../target_vntrs/disease_associated_all_vntrs_v5_unique.txt")
    parser.add_argument("--g-vntrs-short-db", help="Path to the db file for short g-vntrs.",
                        type=str,
                        default="../../target_vntrs/illumina_probed_short_vntrs.db")
    parser.add_argument("--g-vntrs-short-ids", help="Path to the text file including short g-vntr ids.",
                        type=str,
                        default="../../target_vntrs/g_vntr_short.txt")
    parser.add_argument("--g-vntrs-long-db", help="Path to the db file for long g-vntrs.",
                        type=str,
                        default="../../target_vntrs/pacbio_probed_long_vntrs.db")
    parser.add_argument("--g-vntrs-long-ids", help="Path to the text file including long g-vntr ids.",
                        type=str,
                        default="../../target_vntrs/g_vntr_long.txt")
    parser.add_argument("--ids-to-discard-file", help="Path to the text file including ids to be discarded (e.g. test strs).",
                        type=str,
                        default="../../target_vntrs/str_ids_to_discard.txt")
    args = parser.parse_args()
    return args

def main():

    args = parse_arguments()

    # Increase the font in all plots
    plt.rcParams.update({'font.size': 16})
    sns.set_style("whitegrid")

    # Get VNTR DB info.
    id_to_vntr_map = get_id_to_vntr_map(long_vntr_db=args.g_vntrs_long_db,
                                       short_vntr_db=args.g_vntrs_short_db,
                                       p_vntr_db=args.p_vntrs_db,
                                       str_ids_file=args.ids_to_discard_file)

    target_set_to_vids = defaultdict(list)
    all_vids = []
    # Get list of target VNTRs.
    for target_set in ["10k-genic", "phenotype-associated"]:
        target_vids = get_target_vntrs(id_to_vntr_map=id_to_vntr_map,
                                        target_set=target_set,
                                        p_vntr=args.p_vntrs_ids,
                                        g_vntr_short=args.g_vntrs_short_ids,
                                        g_vntr_long=args.g_vntrs_long_ids,
                                        str_ids_file=args.ids_to_discard_file)
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


if __name__ == "__main__":
    main()
