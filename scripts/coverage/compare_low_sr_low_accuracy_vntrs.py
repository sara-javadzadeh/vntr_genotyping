import os
import pandas as pd

targeted_samples_filename="../../../sample_ids.txt"

def print_median(df, column):
    print("Median of {} is {:.2f}".format(
            column, df[column].median()))

def compare(wgs_sr, targeted_spanning_reads):
    print_columns = False
    comparison_df = wgs_sr.merge(
                        targeted_spanning_reads,
                        how="left",
                        on="vid",
                        )
    print("num vntrs in WGS: ", len(comparison_df))
    print_median(comparison_df, "vntr_length")
    comparison_df = comparison_df[comparison_df["wgs_log_sr"] < comparison_df["targeted_log_sr"]]
    print("num vntrs in WGS with higher coverage in targeted: ", len(comparison_df))
    if print_columns:
        print(comparison_df.columns)
    print(comparison_df[["vid", "wgs_log_sr", "targeted_log_sr",
                         "vntr_length", "gc_content_x"]])
    print_median(comparison_df, "wgs_log_sr")
    print_median(comparison_df, "targeted_log_sr")
    print_median(comparison_df, "vntr_length")
    print_median(comparison_df, "gc_content_x")

    return comparison_df

def compare_sr_in_vntrs():
    targeted_low_coverage_path = "/nucleus/projects/saraj/vntr/sources/COH_analysis/scripts/compare_short_long_vntrs/coverage_analysis/low_coverage_logs/low_coverage_vntrs.csv"
    targeted_full_df_path = "/nucleus/projects/saraj/vntr/sources/COH_analysis/scripts/compare_short_long_vntrs/coverage_analysis/logs/spanning_reads_df.csv"
    wgs_low_coverage_path = "logs/vntrs_low_coverage.csv"
    wgs_full_df_path = "logs/spanning_reads_df.csv"

    wgs_sr = pd.read_csv(wgs_full_df_path, index_col="vid")
    wgs_sr.rename({"median across samples": "wgs_log_sr"},
                            axis=1,
                            inplace=True)
    print(wgs_sr.shape)
    targeted_spanning_reads = pd.read_csv(targeted_full_df_path)
    targeted_spanning_reads.rename(
            {"median_spanning_reads_across_samples": "targeted_sr",
             "median_log_spanning_reads_across_samples": "targeted_log_sr"},
             axis=1,
             inplace=True)
    targeted_spanning_reads["vid"] = targeted_spanning_reads["vntr_id"].str.replace("s|l", "").astype(int)
    print(targeted_spanning_reads.shape)


    compare(wgs_sr=wgs_sr,
            targeted_spanning_reads=targeted_spanning_reads)

    merged_dataframe =  wgs_sr.merge(
                        targeted_spanning_reads,
                        how="outer",
                        on="vid",
                        )
    return merged_dataframe

def get_wgs_samples():
    samples_dir="/ribosome/projects/saraj/vntr/data/pangenome_project/pacbio_hifi_mapped_reads"
    files = os.listdir(samples_dir)
    sample_names = [sample for sample in files if sample.startswith("HG")]
    return sample_names


def read_targeted_samples(filename):
    ids = []
    with open(filename) as samples_ids_file:
        for line in samples_ids_file.readlines():
            ids.append(line.strip())
    return ids

# Get inconsistent VNTRs in at least one sample
def get_inconsistent_vntrs(filename, sample, verbose):
    vids = []
    if not os.path.exists(filename):
        if verbose:
            print("File for sample {} does not exist".format(sample))
        return []
    with open(filename, "r") as inconsistent_vntrs_file:
        for line in inconsistent_vntrs_file.readlines():
            if line.startswith("vid"):
                # Skip header line
                continue
            vids.append(int(line.split()[0]))
    return vids

def get_inconsistent_wgs_vids():
    samples = get_wgs_samples()
    # This is including str_like vntrs
    sample_template = "/nucleus/projects/saraj/vntr/sources/pangenome_project/scripts/mendelian_consistency/v3_illumina_pacbio_comparison_no_str_sr_3/comparison_#/#_VNTRs_with_inconsistent_genotype_calls.txt"
    vids = []
    for sample in samples:
        inconsistent_vntrs_filename = sample_template.replace("#", sample)
        vids.extend(get_inconsistent_vntrs(inconsistent_vntrs_filename, sample, verbose=False))
    # Get unique vids
    vids = list(set(vids))
    print("Vids with at least one inconsitent call in WGS samples ", len(vids))
    return vids

def get_inconsistent_targeted_vids():
    samples = read_targeted_samples(targeted_samples_filename)
    sample_template = "/nucleus/projects/saraj/vntr/sources/COH_analysis/scripts/compare_short_long_vntrs/deprecated/v14_integrated_updates_no_str_sr_3/comparison_#/unfiltered/#_VNTRs_with_inconsistent_genotype_calls.txt"
    vids = []
    for sample in samples:
        inconsistent_vntrs_filename = sample_template.replace("#", sample)
        vids.extend(get_inconsistent_vntrs(inconsistent_vntrs_filename, sample, verbose=False))
    # Get unique vids
    vids = list(set(vids))
    print("Vids with at least one inconsitent call in targeted samples ", len(vids))
    return vids

def print_consistency_comparison(dataframe, name_for_printing):
    print("vntrs with {}: \n{}".format(
        name_for_printing,
        len(dataframe[["vid", "wgs_log_sr", "targeted_log_sr",
                             "vntr_length", "gc_content_x"]])))
    print(dataframe[["vid", "wgs_log_sr", "targeted_log_sr",
                                 "vntr_length", "gc_content_x"]])
    dataframe = dataframe[dataframe["gc_content_x"] < 60]
    print("low gc vntrs with {}: \n{}".format(
        name_for_printing,
        len(dataframe[["vid", "wgs_log_sr", "targeted_log_sr",
                             "vntr_length", "gc_content_x"]])))
    dataframe = dataframe[dataframe["targeted_log_sr"]\
                            >= 4]
    print("low gc vntrs with {} and >= 15 srs on targeted: \n{}".format(
        name_for_printing,
        len(dataframe[["vid", "wgs_log_sr", "targeted_log_sr",
                             "vntr_length", "gc_content_x"]])))

def compare_low_accuracy_vntrs(coverage_vntr_df, verbose=True):
    inconsistent_wgs_vids = get_inconsistent_wgs_vids()
    inconsistent_targeted_vids = get_inconsistent_targeted_vids()

    # Get intersection and subtracted sets
    print("vntrs inconsistent in both wgs and targeted ", len(
        set(inconsistent_wgs_vids) & set(inconsistent_targeted_vids)))
    inc_wgs_con_targeted_vids = [vid for vid in inconsistent_wgs_vids \
        if vid not in inconsistent_targeted_vids]
    print("vntrs inconsistent in wgs and consistent in targeted",
        len(inc_wgs_con_targeted_vids))
    inc_targeted_con_wgs_vids = [vid for vid in inconsistent_targeted_vids \
        if vid not in inconsistent_wgs_vids]
    print("vntrs inconsistent in targeted and consistent in wgs",
        len(inc_targeted_con_wgs_vids))

    # Select dataframes for printing
    inconsistent_wgs = coverage_vntr_df[coverage_vntr_df["vid"]\
                            .isin(inconsistent_wgs_vids)]
    if verbose:
        print("vntrs with inconsistent calls on wgs: \n{}".format(
        inconsistent_wgs[["vid", "wgs_log_sr", "targeted_log_sr",
                             "vntr_length", "gc_content_x"]]))

    inconsistent_targeted = coverage_vntr_df[coverage_vntr_df["vid"]\
                            .isin(inconsistent_targeted_vids)]
    if verbose:
        print("vntrs with inconsistent calls on targeted: \n{}".format(
        inconsistent_targeted[["vid", "wgs_log_sr", "targeted_log_sr",
                             "vntr_length", "gc_content_x"]]))
    if verbose:
        # Print detailed comparison with low gc, etc.
        inc_wgs_con_targeted = \
            coverage_vntr_df[coverage_vntr_df["vid"]\
                                .isin(inc_wgs_con_targeted_vids)]
        print_consistency_comparison(inc_wgs_con_targeted, "inconsistent calls on wgs and consistent on targeted")
        inc_targeted_con_wgs = \
            coverage_vntr_df[coverage_vntr_df["vid"]\
                                .isin(inc_targeted_con_wgs_vids)]
        print_consistency_comparison(inc_targeted_con_wgs, "inconsistent calls on targeted and consistent on wgs")

if __name__ == "__main__":
    coverage_vntr_df = compare_sr_in_vntrs()
    compare_low_accuracy_vntrs(coverage_vntr_df)
