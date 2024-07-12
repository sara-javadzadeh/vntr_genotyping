from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
from Bio.SeqUtils import GC
from advntr.models import update_gene_name_and_annotation_in_database, load_unique_vntrs_data
from advntr.vntr_annotation import intersect, include, get_genes_info, \
                        get_refseq_id_to_gene_name_map, get_gene_name_from_refseq_id
from advntr.vntr_annotation import get_exons_info as get_region_info
#from advntr.vntr_annotation import get_gene_name_and_annotation_of_vntr


polymorphism_type_colname = "Associated variant category (major length variation is >=10 RC different in the RC allele range reported)"

ANNOTATION_DIR = 'results/hg38_annotation/'
EXONS = ANNOTATION_DIR + '%s_gene_coding_exons.bed'
INTRONS = ANNOTATION_DIR + '%s_gene_introns.bed'
UTR5 = ANNOTATION_DIR + '%s_gene_5utr.bed'
UTR3 = ANNOTATION_DIR + '%s_gene_3utr.bed'
GENES = ANNOTATION_DIR + '%s_genes.bed'
ENSEMBL_TO_GENE = ANNOTATION_DIR + 'ensemblToGeneName.txt'
UCSC_TO_ENSMBL = ANNOTATION_DIR + 'knownToEnsembl.txt'
REFSEQ_TO_GENE = ANNOTATION_DIR + 'Refseq2Gene.txt'
PROMOTER_RANGE = 500


def get_gc_content(ref_vntr, round_to_int=True, verbose=False):
    """ return GC contents of vntr motif and flanking region"""
    # flanking 100 bps
    left_flanking_region = ref_vntr.left_flanking_region[-100:]
    right_flanking_region = ref_vntr.right_flanking_region[:100]
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

def discard_strs(vids, data):
    filename = "../vntrs_to_discard/str_ids.txt"
    str_ids = []
    with open(filename, "r") as str_file:
        str_ids = str_file.readlines()
    str_ids = [int(single_id) for single_id in str_ids]
    print("vids len before discarding STRs ", len(vids))
    vids = [vid for vid in vids if vid not in str_ids]
    #print(str_ids[:10], type(str_ids[0]))
    data["vid"] = data["vid"].astype(int)
    #print(data["vid"].head())
    data = data[~data["vid"].isin(str_ids)]
    print("vids len after discarding STRs ", len(vids))
    return vids, data

def print_map(input_map):
    unique_values = set(input_map.values())
    for value in unique_values:
        num_vntrs = len([vntr for vntr in input_map if input_map[vntr] == value])
        print("num vntrs within {} are {}".format(value, num_vntrs))


def get_annotations(db_file):
    genes_info = get_genes_info()
    exons_info, _ = get_region_info(EXONS)
    introns_info, _ = get_region_info(INTRONS)
    utr5_info, _ = get_region_info(UTR5)
    utr3_info, _ = get_region_info(UTR3)
    name_mapping = get_refseq_id_to_gene_name_map()
    # translate_ranges = get_translate_ranges(exons_info)
    vntr_to_annotation_map = {}
    reference_vntrs = load_unique_vntrs_data(db_file)
    for ref_vntr in reference_vntrs:
        end = ref_vntr.start_point + ref_vntr.get_length()
        new_gene, new_annotation = get_gene_name_and_annotation_of_vntr(
            ref_vntr.chromosome,
            ref_vntr.start_point,
            end, genes_info, exons_info,
            introns_info, utr3_info,
            utr5_info, name_mapping)
        vntr_to_annotation_map[ref_vntr.id] = new_annotation
    #print_map(vntr_to_annotation_map)
    return vntr_to_annotation_map

def get_gene_name_and_annotation_of_vntr(vntr_chromosome, vntr_start, vntr_end, genes, exons, introns, utr3, utr5, name_mapping=None, gene_reference='refseq'):
    if name_mapping is None:
        name_mapping = get_refseq_id_to_gene_name_map()

    def get_annotation(vntr_start, vntr_end, regions, region_name):
        gene_name, annotation = 'None', 'None'
        for start, end, identifier, direction, seg_number in regions:
            if intersect(start, end, vntr_start, vntr_end):
                if gene_reference == 'ucsc':
                    gene_name = get_gene_name_from_ucsc_id(identifier.split('_')[0])
                else:
                    gene_name = get_gene_name_from_refseq_id(identifier.split('.')[0], name_mapping)
                annotation = region_name
                break
            if start > vntr_end:
                break
        return gene_name, annotation

    gene_names = []
    annotations = []
    gene_name, annotation = "None", "None"
    # Try out all the regions
    for regions, region_txt in [(exons[vntr_chromosome], 'exon'),
                                (introns[vntr_chromosome], 'intron'),
                                (utr5[vntr_chromosome], 'utr'),
                                (utr3[vntr_chromosome], 'utr'),
                                ]:
        i_gene_name, i_annotation = get_annotation(vntr_start, vntr_end, regions, region_txt)
        if i_gene_name != 'None':
            gene_names.append(i_gene_name)
            annotations.append(i_annotation)
    if len(set(gene_names)) > 1:
        # Multiple gene names found for a VNTR
        #print("Gene names for vntr: {}".format(gene_names))
        # Pick the shorter gene name. E.g. 'CACNA1C', 'CACNA1C-AS4'
        # or 'PRPH', 'LOC101927267'
        gene_name = min(gene_names, key=len)
    elif len(set(gene_names)) > 0:
        gene_name = gene_names[0]
    if len(set(annotations)) > 1:
        # VNTR spans multiple regions
        #print("annotations for vntr: {}".format(annotations))
        if "exon" in annotations and "intron" in annotations:
            annotation = "multiple"
        elif "exon" in annotations:
            annotation = "exon"
        elif "intron" in annotations:
            annotation = "intron"
        else:
            print("annotations for vntr outside exon and introns: {}, genes {}".format(annotations, gene_names))
    elif len(set(annotations)) == 1:
        # Otherwise, there is a single region to use
        annotation = annotations[0]

    # promoter region
    if len(gene_names) == 0:
        for start, end, identifier, direction in genes[vntr_chromosome]:
            if direction == '-':
                start = end
                end += PROMOTER_RANGE
            else:
                end = start
                start -= PROMOTER_RANGE
            if intersect(start, end, vntr_start, vntr_end):
                if gene_reference == 'ucsc':
                    gene_name = get_gene_name_from_ucsc_id(identifier.split('_')[0])
                else:
                    gene_name = get_gene_name_from_refseq_id(identifier.split('.')[0], name_mapping)
                annotation = 'promoter'
                break
            if start - PROMOTER_RANGE > vntr_end:
                break
    #if annotation == "None":
    #    print("gene {} annotation {}  {}:{}-{}".format(
    #           gene_name, annotation, vntr_chromosome, start, end))
    return gene_name, annotation

def read_single_g_database(g_vntr_database_filename):
    columns = "vid, nonoverlapping, chromosome, ref_start, gene_name, "\
              "annotation, pattern, left_flanking, right_flanking, repeats"\
              .replace(",", "").split()
    dataframe = pd.read_csv(g_vntr_database_filename,
                            sep=' ',
                            header=None,
                            names=columns,
                            dtype=str)
    # Update the annotations based on the updated map
    db_filename = g_vntr_database_filename.replace(".txt", ".db")
    vntr_to_annotation_map = get_annotations(db_filename)
    #vntr_to_annotation_map = {}
    num_vntrs_with_missing_anotation = 0
    for idx, row in dataframe.iterrows():
        vid = int(row['vid'])
        if vid in vntr_to_annotation_map:
            if vntr_to_annotation_map[vid] != 'None':
                row.at['annotation'] = vntr_to_annotation_map[vid]
        else:
            num_vntrs_with_missing_anotation += 1
            #print("no updated annotation for vid {} {}".format(vid, row['annotation']))
    print("num_vntrs_with_missing_updated_anotation: ", num_vntrs_with_missing_anotation)
    # Rename annotations to match the annotations for the figure
    for value, replacement in [("Coding", "exon"),
                               ("Exon", "exon"),
                               ("Intron", "intron"),
                               ("Promoter", "promoter"),
                               ("UTR", "utr"),
                               ("None", "none")
                                ]:
        dataframe = dataframe.replace(value, replacement)
    # Print count of each annotation
    unique_annotations = dataframe["annotation"].unique()
    print("Unique values in annotation: ", unique_annotations)
    for annotation in unique_annotations:
        print("Num vntrs for {}: {}".format(
                annotation,
                len(dataframe[dataframe["annotation"] == annotation])))
    # Remove VNTRs with None annotation
    #dataframe = dataframe[dataframe["annotation"] != "None"]
    # Add a column with all values the same which helps to create a single (stacked) column in barchart
    dataframe["group"] = " "
    # Rename the column annotation to location
    dataframe["Location"] = dataframe["annotation"]
    #print(dataframe.columns)
    return dataframe

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
    dataframe = pd.read_csv(filename, sep="\t")
    # Replace any character coming after these characters which are meant to provide more information.
    dataframe = dataframe.replace("\(.*", "", regex=True)
    dataframe = dataframe.replace("/.*", "", regex=True)
    # Remove heading or trailing spaces
    dataframe = dataframe.replace("^ +| +$", "", regex=True)
    # Fix upper case and lower case discrepencies
    dataframe["Pathogenic"] = dataframe["Pathogenic"].str.lower()
    dataframe["Location"] = dataframe["Location"].str.lower()
    dataframe["Coding or non-coding exon"] = dataframe["Coding or non-coding exon"].str.lower()
    dataframe[polymorphism_type_colname] = dataframe[polymorphism_type_colname].str.lower()
    vid_column="VNTR ID (TRF, hg38) - at /nucleus/projects/saraj/vntr/sources/COH_analysis/databases/combined_trf_hg38/hg38_VNTRs_by_TRF.db"
    print("len p_vntr dataframe before dropping duplicates", len(dataframe))
    dataframe = dataframe.drop_duplicates(subset=[vid_column])

    print("len p_vntr dataframe ", len(dataframe))
    # Delete details about which intron, which exon etc for smoother aggregation
    dataframe = dataframe.replace("intron.*", "intron", regex=True)
    dataframe = dataframe.replace("exon.*", "exon", regex=True)
    dataframe = dataframe.replace("5' utr", "5'utr", regex=True)
    dataframe = dataframe.replace(".* upstream.*", "5'utr", regex=True)
    dataframe = dataframe.replace(".* downstream.*", "3'utr", regex=True)
    # Combine the 3'utr and 5'utr to match the g-vntrs
    dataframe = dataframe.replace("3'utr", "utr")
    dataframe = dataframe.replace("5'utr", "utr")
    # Rename the classes to match the literature
    dataframe = dataframe.replace("major length variation", "high\n(>=10)\nRC\nvariation", regex=True)
    dataframe = dataframe.replace("minor length variation", "low\n(<10)\nRC\nvariation", regex=True)

    # Print the count of sub categories
    high_rc_category = "high\n(>=10)\nRC\nvariation"
    df_no_high_rc = dataframe[~(dataframe[polymorphism_type_colname] == high_rc_category)]
    df_no_high_rc["Total Length"] = df_no_high_rc["Total Length"].astype(int)
    print("len vntrs in the minor length and specific allele frequencies {} of those {} are >150bp".format(
                len(df_no_high_rc),
                len(df_no_high_rc[df_no_high_rc["Total Length"] > 150])))
    #print(dataframe["GC content (%)"])
    #print("len high gc vntrs ", len(dataframe[dataframe["GC content (%)"].astype(int) >= 60]))

    print(dataframe[dataframe["Gene, loci"] == "TAF1"])
    # Specify coding or non-coding exon in the same column as location
    if specify_coding:
        # Looks like all exon are coding. So I'm skipping this for now.
        for idx, row in dataframe.iterrows():
            if row["Coding or non-coding exon"] == "coding":
                if dataframe.at[idx, "Location"] == "exon":
                    dataframe.at[idx, "Location"] = "coding exon"
            elif row["Coding or non-coding exon"] == "non-coding":
                if dataframe.at[idx, "Location"] == "exon":
                    row.at["Location"] = "non-coding exon"
            elif row["Coding or non-coding exon"] == "-" or isinstance(row["Coding or non-coding exon"], float):
                # outside exon regions
                pass
            else:
                print("unclassified exon region: -{}-".format(row["Coding or non-coding exon"]))
    if verbose:
        print(dataframe[["Gene, loci", "Pathogenic", polymorphism_type_colname,
                         "Location", "Coding or non-coding exon"]])
    print("Printing TAF1 VNTR: ", dataframe[dataframe["Gene, loci"] == "TAF1"])
    return dataframe


def individual_barchart(dataframe, ax, x, y, ylim,
                        hide_legend,
                        ylabel="count",
                        is_multi_category=True):
    colors = sns.color_palette("colorblind")

    # Change the order of labels in x and hue
    categorical_dataframe = dataframe.copy()
    #hue_order = ["multiple", "3'utr","5'utr", 'promoter', 'intron', 'exon'][::-1]
    hue_order = ["multiple", "utr", 'promoter', 'intron', 'exon']
    new_colors_order = sns.color_palette("colorblind", n_colors=len(hue_order))
    new_colors_order = new_colors_order[::-1]
    if is_multi_category:
        categorical_dataframe[x] = pd.Categorical(dataframe[x], ["high\n(>=10)\nRC\nvariation",
                                                  "low\n(<10)\nRC\nvariation",
                                                  "specific\nallele\nsequence"])
        plot = sns.histplot(data=categorical_dataframe,
                ax=ax,
                x=x,
                stat="count",
                multiple="stack",
                hue=y,
                binwidth=10,
                hue_order=hue_order,
                palette=new_colors_order,
                shrink=0.8,
                alpha=0.9)
    else:
        plot = sns.histplot(data=dataframe,
                ax=ax,
                x=x,
                stat="count",
                multiple="stack",
                hue=y,
                hue_order=hue_order,
                #bins=[0.5, 0.7],
                binwidth=2,
                palette=new_colors_order,
                alpha=0.9)

    ax.set(ylabel=ylabel)
    if hide_legend:
        ax.legend([],[], frameon=False)
    else:
        # Re-order handles and labels to see the desired order
        ax.legend(ax.legend_.legendHandles[::-1], hue_order[::-1], title="Region", fontsize=15)
        #print(dataframe[(dataframe[x] == "major\nlength\nvariation") & (dataframe["Location"] == "exon")])
    # Set tick label size
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    if is_multi_category:
        ax.set_yticks(range(0, ylim, 2))
    else:
        ax.yaxis.set_label_position("right")
        #ax.tick_params(axis='x',labelbottom='off')
        ax.set_xlim(-0.75, 0.75)

def plot_p_barchart_only(dataframe, prefix):
    plt.clf()
    sns.set(font_scale=1.7)
    # Create figure base and adjust ratios
    fig = plt.figure(figsize=(8, 6))
    fig.suptitle("Phenotype or disease associated VNTRs", y=0.95)
    # Two axis for plots and one axis to have continous grid line
    gs = fig.add_gridspec(1, 3,  width_ratios=(4, 1, 4),
                          #left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0, hspace=0)
    ylim=18

    # Adjust grid
    sns.set_style("whitegrid")
    #print(plt.rcParams.keys())
    # Create the axes
    ax = fig.add_subplot(gs[0, 0])
    ax_middle = fig.add_subplot(gs[0, 1], sharey=ax)
    ax_secondary = fig.add_subplot(gs[0, 2], sharey=ax)


    # Remove axis lines in the middle
    ax_secondary.spines["left"].set_visible(False)
    ax_middle.spines["left"].set_visible(False)
    ax_middle.spines["right"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Extend y axis so that legend fits without overlapping with the bars
    plt.rcParams["axes.grid"] = True
    ax_middle.grid(axis='x')
    ax.grid(axis='x')
    ax_secondary.grid(axis='x')

    # Set labels
    ax.set_xlabel("Pathogenic (39 VNTRs)", fontdict={'weight': 'bold'})
    ax_secondary.set_xlabel("Non-pathogenic (9 VNTRs)", fontdict={'weight': 'bold'})

    # Remove double ticks
    ax_secondary.tick_params(axis="y", labelleft=False)
    ax_middle.tick_params(axis="y", labelleft=False)
    ax_middle.tick_params(axis="x", labelbottom=False)
    #ax_secondary.get_yaxis().set_visible(False)

    # Rename the values for the figure purposes.
    plotted_dataframe = dataframe.replace(" ", "\n", regex=True)
    # Get sub data for each smaller plot
    pathogenic_df = plotted_dataframe[plotted_dataframe["Pathogenic"] == "pathogenic"]
    non_pathogenic_df = plotted_dataframe[plotted_dataframe["Pathogenic"] == "non-pathogenic"]

    print("len of pathogenic_df: ", len(pathogenic_df))
    print("len of non_pathogenic_df: ", len(non_pathogenic_df))

    # Plot the left chart
    individual_barchart(dataframe=pathogenic_df,
                        ax=ax,
                        x=polymorphism_type_colname,
                        y="Location",
                        ylim=ylim,
                        hide_legend=True)
    # Plot the right chart
    individual_barchart(dataframe=non_pathogenic_df,
                        ax=ax_secondary,
                        x=polymorphism_type_colname,
                        y="Location",
                        ylim=ylim,
                        hide_legend=True)

    plt.savefig("{}_category_chart.pdf".format(prefix), format="pdf", bbox_inches="tight")

def plot_barchart(p_vntr_dataframe, g_vntr_dataframe, prefix, save_fig=True):
    plt.clf()
    sns.set(font_scale=1.2)
    # Create figure base and adjust ratios
    fig = plt.figure(figsize=(10, 6))
    #fig.suptitle("Phenotype or disease associated VNTRs", y=0.95)
    # Two axis for plots and one axis to have continous grid line
    gs = fig.add_gridspec(1, 5,  width_ratios=(3, 0.5, 2, 0.5, 1),
                          #left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0, hspace=0)
    ylim=18

    # Adjust grid
    sns.set_style("whitegrid")
    #print(plt.rcParams.keys())
    # Create the axes for p_vntrs
    ax_p_vntr_1 = fig.add_subplot(gs[0, 0])
    ax_middle_1 = fig.add_subplot(gs[0, 1], sharey=ax_p_vntr_1)
    ax_p_vntr_2 = fig.add_subplot(gs[0, 2], sharey=ax_p_vntr_1)
    # Create the axes for g_vntrs
    ax_g_vntr = fig.add_subplot(gs[0, 4])
    ax_middle_2 = fig.add_subplot(gs[0, 3])


    # Remove axis lines in the middle
    ax_p_vntr_2.spines["left"].set_visible(False)
    ax_p_vntr_2.spines["right"].set_visible(False)
    ax_middle_1.spines["left"].set_visible(False)
    ax_middle_1.spines["right"].set_visible(False)
    ax_p_vntr_1.spines["right"].set_visible(False)

    ax_middle_2.spines["left"].set_visible(False)
    ax_middle_2.spines["right"].set_visible(False)
    ax_g_vntr.spines["left"].set_visible(False)
    ax_g_vntr.spines["right"].set_visible(True)

    # Extend y axis so that legend fits without overlapping with the bars
    # Have only horizontal grid, no vertical grid
    plt.rcParams["axes.grid"] = True
    ax_middle_1.grid(axis='x')
    ax_p_vntr_1.grid(axis='x')
    ax_p_vntr_2.grid(axis='x')
    ax_middle_2.grid(False)
    ax_g_vntr.grid(axis='x')


    # Remove double ticks
    ax_p_vntr_2.tick_params(axis="y", labelleft=False)
    ax_middle_1.tick_params(axis="y", labelleft=False)
    ax_middle_1.tick_params(axis="x", labelbottom=False)
    ax_middle_2.tick_params(axis="y", labelleft=False)
    ax_middle_2.tick_params(axis="x", labelbottom=False)
    ax_g_vntr.tick_params(axis="y", labelright=True, labelleft=False)
    #ax_g_vntr.tick_params(axis="x", labelbottom=False)

    # Rename the values for the figure purposes.
    plotted_dataframe = p_vntr_dataframe.replace(" ", "\n", regex=True)
    g_vntr_dataframe = g_vntr_dataframe.replace(" ", "\n", regex=True)

    # Get sub data for each smaller plot
    pathogenic_df = plotted_dataframe[plotted_dataframe["Pathogenic"] == "pathogenic"]
    non_pathogenic_df = plotted_dataframe[plotted_dataframe["Pathogenic"] == "non-pathogenic"]

    # Set labels
    ax_p_vntr_1.set_xlabel("Pathogenic\n({} P-VNTRs)".format(len(pathogenic_df)), fontdict={'weight': 'bold'})
    ax_p_vntr_2.set_xlabel("Non-pathogenic\n({} P-VNTRs)".format(len(non_pathogenic_df)), fontdict={'weight': 'bold'})
    ax_g_vntr.set_xlabel("{:,} G-VNTRs".format(len(g_vntr_dataframe)), fontdict={'weight': 'bold'})

    # Plot the left p-vntrs
    individual_barchart(dataframe=pathogenic_df,
                        ax=ax_p_vntr_1,
                        x=polymorphism_type_colname,
                        y="Location",
                        ylim=ylim,
                        ylabel="count P-VNTR",
                        is_multi_category=True,
                        hide_legend=True)
    # Plot the right p-vntrs
    individual_barchart(dataframe=non_pathogenic_df,
                        ax=ax_p_vntr_2,
                        x=polymorphism_type_colname,
                        y="Location",
                        ylim=ylim,
                        ylabel="count P-VNTR",
                        is_multi_category=True,
                        hide_legend=False)
    # Plot g-vntrs
    individual_barchart(dataframe=g_vntr_dataframe,
                        ax=ax_g_vntr,
                        x="group",
                        y="Location",
                        ylim=None,
                        ylabel="count G-VNTR",
                        hide_legend=True,
                        is_multi_category=False)

    if save_fig:
        plt.savefig("{}_category_chart.pdf".format(prefix), format="pdf", bbox_inches="tight")

def plot_g_barchart_only(dataframe, prefix):
    plt.clf()
    sns.set(font_scale=1.5)
    colors = sns.color_palette("colorblind")
    # Create figure base and adjust ratios
    fig = plt.figure(figsize=(8, 6))
    fig.suptitle("Gene proximal VNTRs", y=0.95)

    # Adjust grid
    sns.set_style("whitegrid")
    #print(plt.rcParams.keys())
    # Create the axes
    ax = plt.gca()

    # Change the order of labels in x and hue
    hue_order = ["3'utr","5'utr", 'promoter', 'intron', 'exon']
    print("unique values for annotation in df", dataframe["annotation"].unique())
    plot = sns.histplot(data=dataframe,
                ax=ax,
                x='annotation',
                stat="count",
                #hue="annotation",
                binwidth=10,
                hue_order=hue_order,
                palette="colorblind",
                shrink=0.8,
                alpha=0.9)

    plt.savefig("{}_category_chart.pdf".format(prefix), format="pdf", bbox_inches="tight")

if __name__ == "__main__":
    # Phenotype associated dataframe
    p_vntr_dataframe = read_p_dataframe()
    # Gene proximal database
    g_vntr_dataframe = read_g_database()
    print("len of g_vntr_dataframe: ", len(g_vntr_dataframe))
    print("len unique vids in g_vntr_dataframe: ", len(g_vntr_dataframe["vid"].unique()))
    plot_barchart(g_vntr_dataframe=g_vntr_dataframe,
                    p_vntr_dataframe=p_vntr_dataframe,
                    prefix="p_and_g_vntr")
