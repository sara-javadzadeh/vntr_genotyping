from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from advntr.models import update_gene_name_and_annotation_in_database, load_unique_vntrs_data
from advntr.vntr_annotation import intersect, include, get_genes_info, \
                        get_refseq_id_to_gene_name_map, get_gene_name_from_refseq_id
from advntr.vntr_annotation import get_exons_info as get_region_info
#from advntr.vntr_annotation import get_gene_name_and_annotation_of_vntr


def read_g_database():
    long_vntr_database_filename = "/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.txt"
    short_vntr_database_filename = "/nucleus/projects/saraj/vntr/sources/COH_analysis/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.txt"
    short_vntrs_df = read_single_g_database(short_vntr_database_filename)
    print("len short vntrs df ", len(short_vntrs_df))
    long_vntrs_df = read_single_g_database(long_vntr_database_filename)
    print("len long vntrs df ", len(long_vntrs_df))
    concat_df = pd.concat([short_vntrs_df, long_vntrs_df], ignore_index=True, sort=True)
    concat_df = concat_df.drop_duplicates(subset=["vid"])
    print("len concat_df before reset index ", len(concat_df))
    concat_df = concat_df.reset_index(drop=True)
    print("len concat_df after reset index ", len(concat_df))
    #with open("temp_vids.txt", "w") as vids_file:
    #    for vid in list(concat_df["vid"]):
    #        vids_file.write(str(vid) + '\n')
    vids = concat_df["vid"]
    vids, concat_df = discard_strs(vids, concat_df)
    print("len concat_df after dropping STRs ", len(concat_df))
    #print(concat_df.head())
    #print(concat_df.columns)
    return concat_df

def read_p_dataframe(verbose=False, specify_coding=False):
    filename = "phenotype_associated_vntrs.tsv"
    dataframe = pd.read_csv(filename, sep="\t", nrows=55)
    # Replace any character coming after these characters which are meant to provide more information.
    dataframe = dataframe.replace("\(.*", "", regex=True)
    dataframe = dataframe.replace("/.*", "", regex=True)
    # Remove heading or trailing spaces
    dataframe = dataframe.replace("^ +| +$", "", regex=True)
    dataframe = dataframe.replace("-", "", regex=True)
    dataframe = dataframe.dropna(how="all")
    # Get vids
    vid_colname="VNTR ID (TRF, hg38) - at /nucleus/projects/saraj/vntr/sources/COH_analysis/databases/combined_trf_hg38/hg38_VNTRs_by_TRF.db"
    vids = dataframe[vid_colname].astype(int)
    print(vids)
    print("len vids ", len(vids))
    # No need to write v5 data. Upon checking, v4 and v5 are equal
    with open("disease_associated_all_vntrs_v5.txt", "w+") as p_vntr_file:
        for vid in vids:
            p_vntr_file.write(str(vid)+"\n")
    return dataframe

if __name__ == "__main__":
    # Phenotype associated dataframe
    p_vntr_dataframe = read_p_dataframe()
