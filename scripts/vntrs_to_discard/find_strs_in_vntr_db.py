from advntr.models import load_unique_vntrs_data

def find_strs_in_vntr_db(db_path, write=False):
    trs = load_unique_vntrs_data(db_path)
    strs, str_ids = [], []
    for tr in trs:
        if tr.id == 331663:
            print("for vid 331663 {} {} {}".format(
                tr.chromosome,
                tr.start_point,
                tr.start_point + tr.get_length()))
        if len(tr.pattern) < 6:
            strs.append(tr)
            str_ids.append(tr.id)
    print("In total {} TRs in the database were STRs".format(len(str_ids)))
    if len(str_ids) > 0 and write:
        with open("str_ids.txt", "w") as non_target_vntrs:
            for str_id in str_ids:
                non_target_vntrs.write(str(str_id) + "\n")
        with open("strs.bed", "w") as non_target_vntrs:
            for str_info in strs:
                non_target_vntrs.write("{}\t{}\t{}\t{}\n".format(
                    str_info.chromosome,
                    str_info.start_point,
                    str_info.start_point + str_info.get_length(),
                    str_info.id))

if __name__ == "__main__":
    short_vntr_db = "/nucleus/projects/saraj/vntr/sources/COH_analysis" +\
        "/databases/illumina_vntr_db_used_for_probe_design/illumina_probed_short_vntrs.db"
    print("Working on short vntr db")
    find_strs_in_vntr_db(short_vntr_db, write=True)
    long_vntr_db = "/nucleus/projects/saraj/vntr/sources/COH_analysis" +\
        "/databases/pacbio_vntr_db_used_for_probe_design_exact_match/Pacbio_probed_long_vntrs.db"
    print("Working on long vntr db")
    find_strs_in_vntr_db(long_vntr_db)
