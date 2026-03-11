#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
import sys
import numpy as np
from colorama import Fore
from guarnatee_lib.window_srna import WindowSRNA
from guarnatee_lib.rna_classifier import RNAClassifier
from guarnatee_lib.differential_classifier import DifferentialClassifier
from guarnatee_lib.candidates_merger import CandidatesMerger
from guarnatee_lib.helpers import Helpers
from guarnatee_lib.wiggle import Wiggle
from guarnatee_lib.gff import GFF
import argparse
import argcomplete
import pandas as pd
import os
from guarnatee_lib.fasta import Fasta
import pybedtools as pybed
import configparser
import matplotlib.pyplot as plt
__author__ = ("Muhammad Elhossary <elhossary@zbmed.de> "
              "Konrad Förstner <konrad@foerstner.org> ")
__copyright__ = "2021 by Muhammad Elhossary <elhossary@zbmed.de>"
__license__ = "ISC license"
__email__ = "elhossary@zbmed.de"
__version__ = "0.1.0"
__maintainer__ = "Muhammad Elhossary"

def main():
    welcome_text = \
    r"""
===============================================================================
||    ____                    ____    _   _      _       _                   ||
||   / ___|  _   _    __ _   |  _ \  | \ | |    / \     | |_    ___    ___   ||
||  | |  _  | | | |  / _` |  | |_) | |  \| |   / _ \    | __|  / _ \  / _ \  ||
||  | |_| | | |_| | | (_| |  |  _ <  | |\  |  / ___ \   | |_  |  __/ |  __/  ||
||   \____|  \__,_|  \__,_|  |_| \_\ |_| \_| /_/   \_\   \__|  \___|  \___|  ||
===============================================================================""" + "\n\n"
    print(Fore.RED + welcome_text)
    print(Fore.WHITE + "")
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.set_defaults(func=predict)
    parser.add_argument(
        "--gffs", required=True, type=str, nargs="+", help="GFF files (space separated)"
    )
    parser.add_argument(
        "--fastas",
        required=True,
        type=str,
        nargs="+",
        help="Fasta files (space separated)",
    )
    parser.add_argument(
        "--wigs",
        required=True,
        type=str,
        action="append",
        help="Annotated .wig files (space separated)",
    )
    parser.add_argument("--threshold", type=int, default=1)
    parser.add_argument(
        "--config_file", default=f"{os.path.dirname(__file__)}/config.cfg", type=str, help="Configuration file"
    )
    parser.add_argument("--out_dir", required=True, type=str, help="")

    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    if "func" in dir(args):
        args.func(args)
    else:
        parser.print_help()


def predict(args):
    # load files
    conf_dict = {}
    conf_file_path = os.path.abspath(args.config_file)
    conf_parse = configparser.ConfigParser(strict=True)
    print(f"Loading config file: {conf_file_path}")
    conf_parse.read(conf_file_path)
    for sect in conf_parse.sections():
        conf_dict.update(dict(conf_parse.items(sect)), )
    gff_obj = GFF(gff_paths=args.gffs)
    fastas = Fasta(fasta_paths=args.fastas)
    seqid_groups = fastas.organism_seqid_groups

    # TODO: assert correctness of wig files
    wig_anno_cols = ["file_path", "condition", "replicate", "strand", "lib_mode"]
    wigs_df = pd.DataFrame([x.split(":") for x in args.wigs], columns=wig_anno_cols)
    assert wigs_df["lib_mode"].isin(["FL", "P1", "P2", "5E", "3E", "d5E"]).all()

    fl_wigs_df = wigs_df[wigs_df["lib_mode"] == "FL"].copy()
    pairs_wigs_df = wigs_df[wigs_df["lib_mode"].isin(["P1", "P2"])].copy()
    ends_wigs_df = wigs_df[wigs_df["lib_mode"].isin(["3E", "5E", "d5E"])].copy()


    export_df = pd.DataFrame()
    stats_df = pd.DataFrame()


    if not fl_wigs_df.empty:
        fl_wigs_df["file_path_3e"] = fl_wigs_df["file_path"]
        fl_wigs_df.rename({"file_path": "file_path_5e"}, axis=1, inplace=True)
        fl_wigs_df.drop(["lib_mode"], axis=1, inplace=True)
        dispatch_pair(gff_obj, fastas, conf_dict, args, fl_wigs_df, "full_length")

    if not pairs_wigs_df.empty:
        pairs_wigs_df = pd.merge(
            pairs_wigs_df.loc[pairs_wigs_df["lib_mode"] == "P1", wig_anno_cols[:-1]],
            pairs_wigs_df.loc[pairs_wigs_df["lib_mode"] == "P2", wig_anno_cols[:-1]],
            on=["condition", "replicate", "strand"], how="inner", suffixes=('_5e', '_3e'))
        dispatch_pair(gff_obj, fastas, conf_dict, args, pairs_wigs_df, "Paired")

    if not ends_wigs_df.empty:
        diff_ends_wigs_df = ends_wigs_df.loc[ends_wigs_df["lib_mode"] == "d5E", wig_anno_cols[:-1]].copy()

        ends_wigs_df = pd.merge(
            left=ends_wigs_df.loc[ends_wigs_df["lib_mode"] == "5E", wig_anno_cols[:-1]],
            right=ends_wigs_df.loc[ends_wigs_df["lib_mode"] == "3E", wig_anno_cols[:-1]],
            on=["condition", "replicate", "strand"], how="inner", suffixes=('_5e', '_3e'))

        ends_wigs_df = pd.merge(
            left=diff_ends_wigs_df,
            right=ends_wigs_df,
            on=["condition", "replicate", "strand"], how="outer", indicator=True)

        ends_wigs_df = ends_wigs_df.loc[ends_wigs_df["_merge"] != "both", wig_anno_cols]
        diff_ends_wigs_df = ends_wigs_df.loc[ends_wigs_df["_merge"] == "both", wig_anno_cols]
        diff_ends_wigs_df.rename({"file_path": "file_path_d5e"}, axis=1, inplace=True)

        if diff_ends_wigs_df.empty:
            dispatch_pair(gff_obj, fastas, conf_dict, args, ends_wigs_df, "dual_lib")
        else:
            dispatch_diff(gff_obj, fastas, conf_dict, args, diff_ends_wigs_df, "differential")




    return None

def dispatch_pair(gff_obj, fastas, conf_dict, args, wigs_df, lib_mode):
    seqid_groups = fastas.organism_seqid_groups
    export_df = pd.DataFrame()
    stats_df = pd.DataFrame()

    for i in wigs_df.index:
        desc = f"{wigs_df.at[i, "condition"]}_rep_{wigs_df.at[i, "replicate"]}"
        wigs_df.at[i, "file_desc"] = desc

        strand = wigs_df.at[i, "strand"]
        strand_sign = "+" if strand == "f" else "-"
        tmp_exp_df, tmp_stat_df = _call_srnas(
            wigs_df.at[i, "file_path_5e"],
            wigs_df.at[i, "file_path_3e"],
            conf_dict, args.threshold)
        tmp_exp_df["strand"] = strand_sign
        tmp_exp_df["TSS_lib_type"] = "pair"
        tmp_exp_df["condition"] = desc
        tmp_exp_df = Helpers.get_gff_df(tmp_exp_df, anno_source="GuaRNAtee", anno_type="candidate", new_id=True)
        tmp_exp_df = Helpers.warp_non_gff_columns(RNAClassifier(gff_obj, tmp_exp_df, fastas, conf_dict).classes)

        tmp_stat_df["strand"] = strand
        tmp_stat_df["TSS_lib_type"] = "pair"
        tmp_stat_df["file_desc"] = desc

        export_df = pd.concat([export_df, tmp_exp_df], ignore_index=True)
        stats_df = pd.concat([stats_df, tmp_stat_df], ignore_index=True)


    export_df = Helpers.warp_non_gff_columns(export_df)
    export_pb = pybed.BedTool.from_dataframe(export_df).sort().cluster(s=True, d=0)
    export_df = export_pb.to_dataframe(names=gff_obj.column_names + ["cluster_id"])
    clusters = export_df["cluster_id"].unique().size
    export_df = Helpers.warp_non_gff_columns(export_df)
    #generate_orf_stats(Helpers.expand_attributes_to_columns(export_df), seqid_groups, args)
    #if conf_dict["detailed_output"] == "False":
    #    export_df = CandidatesMerger(export_df, float(conf_dict["merge_similarity_ratio"])).merge()
    export_df.sort_values(["seqid", "start", "end"], inplace=True)
    export_df.reset_index(inplace=True, drop=True)
    export_df["source"] = "GuaRNAtee"

    print(f"Total {export_df.shape[0]} candidates in {clusters} unique regions to be exported")
    # Exports
    # GFFs
    for seqid_group, seqids in seqid_groups.items():
        export_df[export_df["seqid"].isin(seqids)].to_csv(
            os.path.abspath(f"{args.out_dir}/{seqid_group}_candidates.gff"),
            index=False,
            sep="\t",
            header=False)
    # Stats
    for i in stats_df.index:
        for k, v in seqid_groups.items():
            if stats_df.at[i, "seqid"] in v:
                stats_df.at[i, "Organism"] = k
    """
    all_thres = []
    for i in stats_df.index:
        all_thres.extend(stats_df.at[i, "TSS_peaks_thresholds"])
        all_thres.extend(stats_df.at[i, "TTS_peaks_thresholds"])
    fig = plt.figure(figsize=(16, 9))
    plt.grid()
    plt.plot(all_thres)
    plt.title("Variable threshold factors of IQR")
    fig.savefig(os.path.abspath(f"{args.out_dir}/dist_thres_iqr.jpg"))
    """
    stats_df = stats_df.groupby(["Organism", "file_desc", "TSS_lib_type"], as_index=False).agg(
        {"TSS_lib_windows_count": "sum",
         "TSS_lib_peaks_count": "sum",
         "TTS_lib_windows_count": "sum",
         "TTS_lib_peaks_count": "sum",
         "peaks_connections_count": "sum"})
    stats_df.to_csv(os.path.abspath(f"{args.out_dir}/stats.tsv"), sep='\t', index=False)

    # Excel
    drop_cols = ["source", "type", "score", "phase", "id",
                 "ss_id", "ss_diff_height", "ss_height", "ss_upstream", "ss_downstream",
                 "ts_id", "ts_diff_height", "ts_height", "ts_upstream", "ts_downstream"]

    export_df = Helpers.expand_attributes_to_columns(export_df)
    for col in drop_cols:
        if col in export_df.columns:
            export_df.drop(columns=[col], inplace=True)
    rename_cols = {c: c.replace("_", " ") for c in export_df.columns}

    classes_groups = {"intergenic": [], "ORF_int": [], "others": []}
    for cls in export_df["annotation_class"].unique():
        if ("cross" in cls and "ncRNA" in cls) or "antisense_to" in cls:
            classes_groups["others"] += [cls]
        elif "ncRNA" in cls or "intergenic" in cls:
            classes_groups["intergenic"] += [cls]
        elif "ORF_int" in cls:
            classes_groups["ORF_int"] += [cls]
        else:
            classes_groups["others"] += [cls]

    ### Tables
    for seqid_group, seqids in seqid_groups.items():
        with pd.ExcelWriter(os.path.abspath(f"{args.out_dir}/{seqid_group}_candidates.xlsx"), engine="openpyxl") \
                as writer:
            for classes_group, classes in classes_groups.items():
                seqid_type_df = export_df[(export_df["seqid"].isin(seqids)) &
                                          (export_df["annotation_class"].isin(classes))].copy()
                if seqid_type_df.empty:
                    continue
                seqid_type_df.reset_index(inplace=True, drop=True)
                seqid_type_df = seqid_type_df.replace("", np.nan).infer_objects()

                seqid_type_df.dropna(how='all', axis=1, inplace=True)
                seqid_type_df.rename(columns=rename_cols, inplace=True)
                seqid_type_df.to_excel(
                    excel_writer=writer,
                    sheet_name=f"{classes_group}",
                    index=True,
                    header=True,
                    na_rep="",
                    index_label="index")

    return export_df, stats_df


def dispatch_diff(gff_obj, fastas, conf_dict, args, wigs_df, lib_mode):
    pass

def _call_srnas(five_end_path, three_end_path, conf_dict, threshold):
    if five_end_path == three_end_path:
        wiggle = Wiggle(five_end_path)
        srnas = WindowSRNA(wiggle, wiggle)
    else:
        srnas = WindowSRNA(Wiggle(five_end_path), Wiggle(three_end_path))
    srnas.call_window_srna(conf_dict, thres_factor=threshold)

    return srnas.srna_candidates, srnas.log_df

def generate_orf_stats(export_df, seqid_groups, args):
    ### ORF int Stats
    orf_int_stats_df = export_df[export_df["sub_class"] != ""].copy()
    orf_int_stats_df["Specie"] = ""
    for seqid_group, seqids in seqid_groups.items():
        orf_int_stats_df.loc[orf_int_stats_df["seqid"].isin(seqids), "Specie"] = seqid_group
    orf_int_stats_df = orf_int_stats_df.value_counts(
        subset=["Specie", "lib_type", "condition", "sub_class"]).reset_index()
    orf_int_stats_df.columns = ["Specie", "lib_type", "condition", "sub_class", 'count']
    with pd.ExcelWriter(os.path.abspath(f"{args.out_dir}/ORF_int_stats.xlsx"), engine="openpyxl") as writer:
        for sp in orf_int_stats_df["Specie"].unique():
            tmp_stat_df = orf_int_stats_df[orf_int_stats_df["Specie"] == sp].copy()
            # tmp_stat_df.drop(columns=[sp], inplace=True)
            tmp_stat_df.to_excel(excel_writer=writer, sheet_name=sp, index=False, header=True, na_rep="", verbose=True)

if __name__ == "__main__":
    main()