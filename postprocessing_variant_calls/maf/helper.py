#!/usr/bin/env python
# imports
import os
import sys
import csv
import re

from pathlib import Path
from contextlib import ExitStack
from typing import List, Optional
import typer
import pandas as pd
import numpy as np
from .resources import tsg_genes

from postprocessing_variant_calls.maf.tag.tag_constants import (
    MAF_DUMMY_COLUMNS2,
    MAF_COLUMNS_SELECT,
    MAF_DUMMY_COLUMNS,
    MAF_TSV_COL_MAP,
    EXONIC_FILTERED,
    SILENT_FILTERED,
    NONPANEL_EXONIC_FILTERED,
    NONPANEL_SILENT_FILTERED,
    DROPPED,
    ALLOWED_EXONIC_VARIANT_CLASS,
)


def process_paths(paths):
    file = open(paths, "r")
    files = []
    for line in file.readlines():
        files.append(line.rstrip("\n"))
    file.close
    return files


def IS_EXONIC_CLASS(Gene, Variant_Classification, Coordinate):
    """
    Determine whether a variant can be considered as exonic
    based on user-defined conditions. Multiple user-defined
    conditions can be added to the conditional block.
    """
    if any(
        [
            Variant_Classification in ALLOWED_EXONIC_VARIANT_CLASS,
            Gene == "TERT" and Variant_Classification == "5'Flank",
        ]
    ):
        return (Gene, Variant_Classification, Coordinate)
    elif any(
        [
            Gene == "MET"
            and Variant_Classification == "Intron"
            and Coordinate >= 116411708
            and Coordinate <= 116414935
        ]
    ):
        return (Gene, "Splice_Site", Coordinate)
    else:
        return None


def check_maf(files: List[Path]):
    acceptable_extensions = [".maf", ".txt", ".csv", "tsv"]
    # return non if argument is empty
    if files is None:
        return None
    # check that we have a list of mafs after reading off the cli
    extensions = [os.path.splitext(f)[1] for f in files]
    for ext in extensions:
        if ext not in acceptable_extensions:
            typer.secho(
                f"If using files argument, all files must be mafs using the same extension.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()
    return files


# NOTE: move these helper functions over to tag_helpers and import from there
def add_dummy_columns(maf, columns):
    """
    Temporary function to add dummy columns
    to meet DMP requirements
    """
    for col in columns:
        if not col in maf.columns:
            maf[col] = ""
    return maf


def customize_cosmic(cosmic_id, cosmic_occurrence):
    """
    helper function to customize cosmic annotation.
    If cosmic_id is defined, but not occurrence, then
    a generic value of "1(unknown)" will be used.
    """
    if cosmic_id is not np.nan and cosmic_id != "":
        # OCCURENCE spelled incorrectly by design
        return (
            "ID="
            + cosmic_id
            + ";OCCURENCE="
            + (cosmic_occurrence if cosmic_occurrence is not np.nan else "1(unknown)")
        )
    else:
        return ""


def get_exon(maf_exon, maf_intron):
    """ "
    helper function to determine the exonic or
    intronic location of a variant.
    """
    try:
        exon, total_exon = str(maf_exon).split("/")
        return "exon" + str(exon)
    except ValueError:  # Not exonic
        try:
            intron, total_intron = str(maf_intron).split("/")
            return "intron" + str(intron)
        except ValueError:  # Not intronic
            return ""


def filter_by_annotations(maf, ref_tx_file, project_name):

    return 0


def cleanup_post_filter(df_post_filter):
    # Change duplex columns to have D_
    df_post_filter.rename(
        columns={
            "t_alt_count_fragment": "D_t_alt_count_fragment",
            "t_ref_count_fragment": "D_t_ref_count_fragment",
            "t_vaf_fragment": "D_t_vaf_fragment",
        },
        inplace=True,
    )
    # Move Status column next to Hotspots
    col = list(df_post_filter)
    col.insert(col.index("D_t_alt_count_fragment"), col.pop(col.index("Status")))
    df_post_filter = df_post_filter[col]
    # Add Match Normal columns even when sample is unmatched
    if "Matched_Norm_Sample_Barcode" not in col:
        df_post_filter.insert(
            col.index("SD_t_vaf_fragment") + 1,
            "Matched_Norm_Sample_Barcode",
            "Unmatched",
        )
        df_post_filter.insert(
            col.index("SD_t_vaf_fragment") + 2, "n_alt_count_fragment", "NA"
        )
        df_post_filter.insert(
            col.index("SD_t_vaf_fragment") + 3, "n_ref_count_fragment", "NA"
        )
        df_post_filter.insert(
            col.index("SD_t_vaf_fragment") + 4, "n_vaf_fragment", "NA"
        )
    # Always add a column with Matched_Norm_Bamfile
    # TODO: ADD filename here
    col = list(df_post_filter)
    df_post_filter.insert(
        col.index("Matched_Norm_Sample_Barcode") + 1, "Matched_Norm_Bamfile", "NA"
    )
    return df_post_filter


def make_condensed_post_filter(df_post_filter):

    # creating the "condensed" MAF -- can be customized in the future
    df_condensed = df_post_filter.loc[:, :"n_vaf_fragment"]
    return df_condensed


# FindVAFandSummary
def _find_VAFandsummary(df, sample_group):  # add category as third argumnet
    # add a line of code here to rename the simplex, duplex and simplex_duplex columns with a prefix of the category they belong to.
    df = df.copy()
    # find the VAF from the fillout (the comma separated string values that the summary will later be calculated from)
    # NOTE: col [t_vaf_fragment] already calculated by traceback, no need to create column again

    if (~df["fillout_type"].isin(["MATCHED_NORMAL", "UNMATCHED_NORMAL"])).any():
        df[f"summary_fragment_{str(sample_group)}"] = (
            "DP="
            + (
                df[f"t_alt_count_fragment_{str(sample_group)}"].astype(int)
                + df[f"t_ref_count_fragment_{str(sample_group)}"].astype(int)
            ).astype(str)
            + ";RD="
            + df[f"t_ref_count_fragment_{str(sample_group)}"].astype(str)
            + ";AD="
            + df[f"t_alt_count_fragment_{str(sample_group)}"].astype(str)
            + ";VF="
            + df[f"t_vaf_fragment_{str(sample_group)}"].fillna(0).astype(str)
        )
    else:
        df["t_vaf_fragment_standard"] = (
            df["t_alt_count_fragment_standard"]
            / (
                df["t_alt_count_fragment_standard"].astype(int)
                + df["t_ref_count_fragment_standard"].astype(int)
            )
        ).round(4)
        df[f"summary_fragment_standard"] = (
            "DP="
            + (
                df[f"t_alt_count_fragment_standard"].astype(int)
                + df[f"t_ref_count_fragment_standard"].astype(int)
            ).astype(str)
            + ";RD="
            + df[f"t_ref_count_fragment_standard"].astype(str)
            + ";AD="
            + df[f"t_alt_count_fragment_standard"].astype(str)
            + ";VF="
            + df[f"t_vaf_fragment_standard"].fillna(0).astype(str)
        )

    return df


def maf_duplicates(data_frame):
    de_duplication_columns = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Variant_Classification",
        "Variant_Type",
        "HGVSc",
        "HGVSp",
        "HGVSp_Short",
    ]
    return data_frame.drop_duplicates(subset=de_duplication_columns)


def check_txt(paths: Path):
    # return None if argument is empty, kind of unfortunate we have to handle this case
    if paths is None:
        return None
    # check that we have a text file after reading off the cli
    extension = os.path.splitext(paths)[1]
    if extension != ".txt":
        typer.secho(
            f"If using paths argument, must provided an input txt file.",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    return paths


def check_separator(separator: str):
    separator_dict = {"tsv": "\t", "csv": ","}
    if separator in separator_dict.keys():
        sep = separator_dict[separator]
    else:
        typer.secho(
            f"Separator for delimited file must be 'tsv' or 'csv', not '{separator}'",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    return sep


def read_tsv(tsv, separator, canonical_tx_ref_flag=False):
    """Read a tsv file

    Args:
        maf (File): Input MAF/tsv like format file
        canonical tx ref flag: flag that if set to True, processes a TSV meant to be the
        Reference canonical transcript file, meant to be used in the SNVs/indels nextflow pipeline.

    Returns:
        data_frame: Output a data frame containing the MAF/tsv
        OR
        tx_isoform_lst: If a Reference canonical transcript file is given, a list of extracted transcript isoforms are returned.
    """
    typer.echo("Read Delimited file...")
    if canonical_tx_ref_flag != False:
        skip = get_row(tsv)
        try:
            tx_tsv = pd.read_csv(
                tsv,
                sep=separator,
                skiprows=skip,
                low_memory=False,
                usecols=["isoform", "gene_name", "refseq_id"],
            )
            tx_isoform_list = tx_tsv.isoform.values.tolist()
            return tx_tsv, tx_isoform_list
        except [IOError, KeyError] as e:
            print(
                "Transcript file requires the following columns {}.".format(
                    ",".join(["isoform", "gene_name", "refseq_id"])
                )
            )
            raise
    else:
        skip = get_row(tsv)
        return pd.read_csv(tsv, sep=separator, skiprows=skip, low_memory=False)


def get_row(tsv_file):
    """Function to skip rows

    Args:
        tsv_file (file): file to be read

    Returns:
        list: lines to be skipped
    """
    skipped = []
    with open(tsv_file, "r") as FH:
        skipped.extend(i for i, line in enumerate(FH) if line.startswith("#"))
    return skipped


def gen_id_tsv(df):
    cols = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
    ]
    if set(cols).issubset(set(df.columns.tolist())):
        df["id"] = df[cols].apply(
            lambda x: "_".join(x.replace("-", "").astype(str)), axis=1
        )
    else:
        typer.secho(
            f"tsv file must include {cols} columns to generate an id for annotating the input maf.",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    return df


class MAFFile:
    def __init__(self, file_path, separator, header=None):
        self.file_path = file_path
        self.separator = separator
        self.cols = {
            "general": [
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ],
            "blocklist": [
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele",
                "Annotation",
            ],
            "germline_status": ["t_alt_count", "t_depth"],
            "common_variant": ["gnomAD_AF"],
            "prevalence_in_cosmicDB": ["CNT"],
            "truncating_mut_in_TSG": [
                "Consequence",
                "Variant_Classification",
                "Hugo_Symbol",
            ],
            "hotspot": ["t_alt_count", "hotspot"],
            "non_hotspot": ["t_alt_count", "hotspot"],
            "not_complex": ["complexity"],
            "mappable": ["mappability"],
            "non_common_variant": ["common_variant"],
            "cmo_ch_filter": [
                "t_alt_count",
                "t_depth",
                "gnomAD_AF",
                "Consequence",
                "Variant_Classification",
                "Hugo_Symbol",
                "t_alt_count",
                "hotspot",
                "t_alt_count",
                "complexity",
                "mappability",
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ],
            "traceback": {
                "standard": [
                    "t_ref_count_standard",
                    "t_alt_count_standard",
                    "t_total_count_standard",
                ],
                "access": [
                    "t_ref_count_fragment_simplex_duplex",
                    "t_alt_count_fragment_simplex_duplex",
                    "t_total_count_fragment_simplex_duplex",
                ],
            },
        }
        self.header = self.__process_header(header) if header is not None else None
        self.data_frame = self.__read_tsv()
        self.__gen_id()
        self.tsg_genes = tsg_genes

    def convert_annomaf_to_df(self):
        if self.data_frame.empty == False:
            self.data_frame["Chromosome"] = self.data_frame["Chromosome"].astype(str)
            self.data_frame.set_index(self.cols["general"], drop=False, inplace=True)
            self.data_frame.rename(
                columns={
                    "Matched_Norm_Sample_Barcode": "caller_Norm_Sample_Barcode",
                    "t_depth": "caller_t_depth",
                    "t_ref_count": "caller_t_ref_count",
                    "t_alt_count": "caller_t_alt_count",
                    "n_depth": "caller_n_depth",
                    "n_ref_count": "caller_n_ref_count",
                    "n_alt_count": "caller_n_alt_count",
                    "set": "CallMethod",
                },
                inplace=True,
            )
            # quick cleanup of mutect columns (can be made into mini function in the MAF class)

            # marking the mutect and vardict combo columns if the condition in the line below is met.
            self.data_frame.loc[
                (self.data_frame["MUTECT"] == 1)
                & (self.data_frame["CallMethod"] != "MuTect"),
                "CallMethod",
            ] = "VarDict,MuTect"
            self.data_frame.drop(
                ["TYPE", "FAILURE_REASON", "MUTECT"], inplace=True, axis=1
            )
            return self.data_frame
        else:
            typer.secho(
                f"failed to open path to the annotation MAF file {self.file_path}.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()

    def _convert_fillout_to_df(self):
        if self.data_frame.empty == False:
            self.data_frame["Chromosome"] = self.data_frame["Chromosome"].astype(str)
            # self.data_frame.set_index(self.cols['general'], drop=False, inplace=True)
            return self.data_frame
        else:
            typer.secho(
                f"failed to open path to the fillout MAF file {self.file_path}.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()

    def __read_tsv(self):
        """Read the tsv file and store it in the instance variable 'data_frame'.

        Args:
            self

        Returns:
            pd.DataFrame: Output a data frame containing the MAF/tsv
        """
        if Path(self.file_path).is_file():
            typer.secho(
                f"Reading Delimited file: {self.file_path}",
                fg=typer.colors.BRIGHT_GREEN,
            )
            skip = self.get_row()
            df = pd.read_csv(
                self.file_path, sep=self.separator, skiprows=skip, low_memory=False
            )
            if self.header:
                df = df[df.columns.intersection(self.header)]
            return df
        else:
            typer.secho(f"failed to open {self.file_path}", fg=typer.colors.RED)
            raise typer.Abort()

    def get_row(self):
        """Function to skip rows

        Returns:
            list: lines to be skipped
        """
        skipped = []
        with open(self.file_path, "r") as FH:
            skipped.extend(i for i, line in enumerate(FH) if line.startswith("#"))
        return skipped

    def merge(self, maf, id, how):
        maf_df = self.data_frame.merge(maf, on=id, how=how)
        return maf_df

    def __gen_id(self):
        cols = self.cols["general"]
        if set(cols).issubset(set(self.data_frame.columns.tolist())):
            self.data_frame["id"] = self.data_frame[cols].apply(
                lambda x: "_".join(x.replace("-", "").astype(str)), axis=1
            )
            first_column = self.data_frame.pop("id")
            self.data_frame.insert(0, "id", first_column)
        else:
            typer.secho(
                f"maf file must include {cols} columns to generate an id for annotating the input maf.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()

    def annotate_maf_maf(self, maf_df_a, cname, values):
        self.data_frame[cname] = np.where(
            self.data_frame["id"].isin(maf_df_a["id"]), values[0], values[1]
        )
        return self.data_frame

    def extract_fillout_type(self):

        # make a call to the _convert_fillout_to_df() function since it is also within the MAF class
        df_full_fillout = self._convert_fillout_to_df()

        # run extract fillout type function on the fillout df (will result in many mini dfs)
        # extract the VAF and summary values for the curated samples
        df_curated = df_full_fillout[df_full_fillout["fillout_type"].isin(["CURATED"])]
        df_plasma = df_full_fillout[df_full_fillout["fillout_type"].isin(["PLASMA"])]
        df_tumor = df_full_fillout[df_full_fillout["fillout_type"].isin(["CASE"])]
        df_control = df_full_fillout[df_full_fillout["fillout_type"].isin(["CONTROL"])]

        df_matched_normal = df_full_fillout[
            df_full_fillout["fillout_type"].isin(["MATCHED_NORMAL"])
        ]

        # make a call to the findVAFandSummary function for each of the subgroups within curated (simplex,duplex)

        df_curated_simplex_summary_added = _find_VAFandsummary(df_curated, "simplex")
        df_curated_simplex_duplex_summary_added = _find_VAFandsummary(
            df_curated_simplex_summary_added, "duplex"
        )
        df_all_curated_SD = _find_VAFandsummary(
            df_curated_simplex_duplex_summary_added, "simplex_duplex"
        )

        df_plasma_simplex_summary_added = _find_VAFandsummary(df_plasma, "simplex")
        df_plasma_simplex_duplex_summary_added = _find_VAFandsummary(
            df_plasma_simplex_summary_added, "duplex"
        )
        df_all_plasma_SD = _find_VAFandsummary(
            df_plasma_simplex_duplex_summary_added, "simplex_duplex"
        )

        df_tumor_simplex_summary_added = _find_VAFandsummary(df_tumor, "simplex")
        df_tumor_simplex_duplex_summary_added = _find_VAFandsummary(
            df_tumor_simplex_summary_added, "duplex"
        )
        df_all_tumor_SD = _find_VAFandsummary(
            df_tumor_simplex_duplex_summary_added, "simplex_duplex"
        )

        df_control_simplex_summary_added = _find_VAFandsummary(df_control, "simplex")
        df_control_simplex_duplex_summary_added = _find_VAFandsummary(
            df_control_simplex_summary_added, "duplex"
        )
        df_all_control_SD = _find_VAFandsummary(
            df_control_simplex_duplex_summary_added, "simplex_duplex"
        )

        df_matched_normal = _find_VAFandsummary(df_matched_normal, "standard")

        df_normals = df_full_fillout[
            df_full_fillout["fillout_type"].isin(["MATCHED_NORMAL", "UNMATCHED_NORMAL"])
        ]

        df_all_normals = _find_VAFandsummary(df_normals, "standard")

        return (
            df_all_curated_SD,
            df_all_plasma_SD,
            df_all_tumor_SD,
            df_matched_normal,
            df_all_normals,
            df_all_control_SD,
        )

    def tag(self, tagging):
        cols = self.cols[tagging]
        if isinstance(cols, dict):
            dictionary = True
        if set(cols).issubset(set(self.data_frame.columns.tolist())) or dictionary:
            if tagging == "germline_status":
                self.data_frame["t_alt_freq"] = pd.to_numeric(
                    (self.data_frame["t_alt_count"])
                ) / pd.to_numeric(self.data_frame["t_depth"])
                self.data_frame["germline_status"] = np.where(
                    (self.data_frame["t_alt_freq"] > 0.35), "likely_germline", ""
                )
            if tagging == "common_variant":
                self.data_frame["gnomAD_AF"] = pd.to_numeric(
                    self.data_frame["gnomAD_AF"]
                )
                self.data_frame["common_variant"] = np.where(
                    (self.data_frame["gnomAD_AF"] > 0.05), "yes", "no"
                )
            if tagging == "prevalence_in_cosmicDB":
                self.data_frame["prevalence_in_cosmicDB"] = self.data_frame[
                    "CNT"
                ].apply(lambda x: int(x.split(",")[0]) if pd.notnull(x) else x)
                self.data_frame.drop(["CNT"], axis=1, inplace=True)
            if tagging == "truncating_mut_in_TSG":
                self.data_frame["truncating_mutation"] = np.where(
                    (self.data_frame["Consequence"].str.contains("stop_gained"))
                    | (self.data_frame["Variant_Classification"] == "Frame_Shift_Ins")
                    | (self.data_frame["Variant_Classification"] == "Nonsense_Mutation")
                    | (self.data_frame["Variant_Classification"] == "Splice_Site")
                    | (self.data_frame["Variant_Classification"] == "Frame_Shift_Del")
                    | (
                        self.data_frame["Variant_Classification"]
                        == "Translation_Start_Site"
                    ),
                    "yes",
                    "no",
                )
                self.data_frame["tumor_suppressor_gene"] = np.where(
                    self.data_frame["Hugo_Symbol"].isin(self.tsg_genes), "yes", "no"
                )
                self.data_frame["truncating_mut_in_TSG"] = np.where(
                    (
                        (self.data_frame["tumor_suppressor_gene"] == "yes")
                        & (self.data_frame["truncating_mutation"] == "yes")
                    ),
                    "yes",
                    "no",
                )
            if tagging == "traceback":
                self.tag_traceback(cols, tagging)

            if tagging == "by_variant_classification":
                self.tag_by_variant_classification(cols, tagging)

        else:
            typer.secho(
                f"missing columns expected for {tagging} tagging expects: {set(cols).difference(set(self.data_frame.columns.tolist()))}, which was missing from the input",
                fg=typer.colors.RED,
            )
            raise typer.Abort()
        return self.data_frame

    # NOTE: replicate this for access_filters
    def tag_all(self, tagging):
        if tagging == "cmo_ch_tag":
            self.tag("germline_status")
            self.tag("common_variant")
            self.tag("prevalence_in_cosmicDB")
            self.tag("truncating_mut_in_TSG")
        return self.data_frame

    def tag_traceback(self, cols, tagging):
        if set(cols["standard"] + cols["access"]).issubset(
            set(self.data_frame.columns.tolist())
        ):
            self.data_frame["t_ref_count"] = self.data_frame[
                "t_ref_count_standard"
            ].combine_first(self.data_frame["t_ref_count_fragment_simplex_duplex"])
            self.data_frame["t_alt_count"] = self.data_frame[
                "t_alt_count_standard"
            ].combine_first(self.data_frame["t_alt_count_fragment_simplex_duplex"])
            self.data_frame["t_total_count"] = self.data_frame[
                "t_total_count_standard"
            ].combine_first(self.data_frame["t_total_count_fragment_simplex_duplex"])
        elif set(cols["standard"]).issubset(set(self.data_frame.columns.tolist())):
            self.data_frame["t_ref_count"] = self.data_frame["t_ref_count_standard"]
            self.data_frame["t_alt_count"] = self.data_frame["t_alt_count_standard"]
            self.data_frame["t_total_count"] = self.data_frame["t_total_count_standard"]
        elif set(cols["access"]).issubset(set(self.data_frame.columns.tolist())):
            self.data_frame["t_ref_count"] = self.data_frame[
                "t_ref_count_fragment_simplex_duplex"
            ]
            self.data_frame["t_alt_count"] = self.data_frame[
                "t_alt_count_fragment_simplex_duplex"
            ]
            self.data_frame["t_total_count"] = self.data_frame[
                "t_total_count_fragment_simplex_duplex"
            ]
        else:
            typer.secho(
                f"missing columns expected for {tagging} tagging expects: {set(cols['standard']).difference(set(self.data_frame.columns.tolist()))} or {set(cols['access']).difference(set(self.data_frame.columns.tolist()))}, which were missing from the input",
                fg=typer.colors.RED,
            )
            raise typer.Abort()

    def tag_by_variant_classification(self, output_dir, ref_lst):
        maf = add_dummy_columns(self.data_frame, MAF_DUMMY_COLUMNS)

        # start tagging the input MAF with the 5 categories

        def is_exonic_or_silent(row, ref_lst):
            exonic_result = IS_EXONIC_CLASS(
                row.Hugo_Symbol, row.Variant_Classification, row.vcf_pos
            )
            if exonic_result:
                if row.Transcript_ID in ref_lst:
                    return "exonic"
                else:
                    return "nonpanel_exonic"
            else:
                if row.Transcript_ID in ref_lst:
                    return "silent"
                else:
                    return "nonpanel_silent"

        def dropped(row):
            return pd.notna(row.Status) and row.Status != ""

        def tag_row(row, canonical):
            tags = []

            # Evaluate "dropped"
            if dropped(row):
                tags.append("dropped")

            # Evaluate "exonic" or "silent" independently of "dropped"
            is_exonic_or_silent_tag = is_exonic_or_silent(row, canonical)
            if is_exonic_or_silent_tag:
                tags.append(is_exonic_or_silent_tag)

            return ", ".join(tags)

        tags_list = []
        for row in maf.itertuples(index=False):
            tags_list.append(tag_row(row, ref_lst))

        maf["classification"] = tags_list

        return maf

        # file_names = [
        #     EXONIC_FILTERED,
        #     SILENT_FILTERED,
        #     NONPANEL_EXONIC_FILTERED,
        #     NONPANEL_SILENT_FILTERED,
        #     DROPPED
        #     ]

        # Create exonic, silent, and nonpanel files.
        # file_paths = [f"{output_dir}/{name}" for name in file_names]

        # headers = "\t".join(MAF_TSV_COL_MAP.values()) + "\n"

        # with ExitStack() as stack:
        #     files = [stack.enter_context(open(path, "w")) for path in file_paths]

        #     for file in files:
        #         file.write(headers)

        # with open(output_dir + "/" + EXONIC_FILTERED, "w") as exonic, open(
        #     output_dir + "/" + SILENT_FILTERED, "w") as silent, open(
        #     output_dir + "/" + NONPANEL_EXONIC_FILTERED, "w") as nonpanel_exonic, open(
        #     output_dir + "/" + NONPANEL_SILENT_FILTERED, "w") as nonpanel_silent, open(
        #     output_dir + "/" + DROPPED, "w") as dropped:

        #     headers = "\t".join(MAF_TSV_COL_MAP.values()) + "\n"
        #     for f in [exonic, silent, nonpanel_exonic, nonpanel_silent, dropped]:
        #         f.write(headers)

        return headers

        # typer.secho(
        #     f"missing columns expected for {tagging} tagging",
        #     fg=typer.colors.RED,
        # )
        # raise typer.Abort()

    def tag_by_variant_annotations(self, rules_df):
        if rules_df is not None:
            for index, row in rules_df.iterrows():
                condition = True
                is_list = lambda var: isinstance(var, list)

                (
                    Tag_Column_Name,
                    Hugo_Symbol,
                    Variant_Classification,
                    Start_Position,
                    End_Position,
                ) = row[
                    [
                        "Tag_Column_Name",
                        "Hugo_Symbol",
                        "Variant_Classification",
                        "Start_Position",
                        "End_Position",
                    ]
                ]

                if Hugo_Symbol != "none":
                    condition &= self.data_frame["Hugo_Symbol"].isin(Hugo_Symbol)
                if Variant_Classification != "none":
                    condition &= self.data_frame["Variant_Classification"].isin(
                        Variant_Classification
                    )
                if Start_Position != "none":
                    condition &= self.data_frame["Start_Position"] >= float(
                        Start_Position
                    )
                if End_Position != "none":
                    condition &= self.data_frame["End_Position"] <= float(End_Position)

                colname = " ".join(Tag_Column_Name)
                tag_column_name = f"is_{colname}_variant"
                self.data_frame[tag_column_name] = "No"
                self.data_frame.loc[condition, tag_column_name] = "Yes"
        else:
            typer.secho(
                f"MAF File is empty. Please check your inputs again.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()
        return self.data_frame

    def filter(self, filter):
        cols = self.cols[filter]
        if set(cols).issubset(set(self.data_frame.columns.tolist())):
            if filter == "hotspot":
                self.data_frame["t_alt_count"] = pd.to_numeric(
                    self.data_frame["t_alt_count"]
                )
                self.data_frame["hotspot_retain"] = np.where(
                    (
                        (self.data_frame["hotspot"] == "yes")
                        & (self.data_frame["t_alt_count"] >= 3)
                    ),
                    "yes",
                    "no",
                )
                self.data_frame = self.data_frame[
                    self.data_frame["hotspot_retain"] == "yes"
                ]
            if filter == "non_hotspot":
                self.data_frame["t_alt_count"] = pd.to_numeric(
                    self.data_frame["t_alt_count"]
                )
                self.data_frame["non_hotspot_retain"] = np.where(
                    (
                        (self.data_frame["hotspot"] == "no")
                        & (self.data_frame["t_alt_count"] >= 5)
                    ),
                    "yes",
                    "no",
                )
                self.data_frame = self.data_frame[
                    self.data_frame["non_hotspot_retain"] == "yes"
                ]
            if filter == "not_complex":
                self.data_frame = self.data_frame[self.data_frame["complexity"] == "no"]
            if filter == "mappable":
                self.data_frame = self.data_frame[
                    self.data_frame["mappability"] == "no"
                ]
            if filter == "non_common_variant":
                self.data_frame = self.data_frame[
                    self.data_frame["common_variant"] == "no"
                ]
            if filter == "cmo_ch_filter":
                self.data_frame["retain"] = np.where(
                    (
                        (
                            (self.data_frame["hotspot"] == "yes")
                            & (self.data_frame["t_alt_count"] >= 3)
                        )
                        | (
                            (self.data_frame["hotspot"] == "no")
                            & (self.data_frame["t_alt_count"] >= 5)
                        )
                    ),
                    "yes",
                    "no",
                )
                self.data_frame = self.data_frame[
                    (
                        (self.data_frame["common_variant"] == "no")
                        & (self.data_frame["mappability"] == "no")
                        & (self.data_frame["complexity"] == "no")
                        & (self.data_frame["retain"] == "yes")
                    )
                ]
        else:
            typer.secho(
                f"missing columns expected for {filter} filtering expects: {set(cols).difference(set(self.data_frame.columns.tolist()))}, which was missing from the input",
                fg=typer.colors.RED,
            )
            raise typer.Abort()
        return self.data_frame

    def __process_header(self, header):
        possible_delimiters = [",", "\t", " "]

        with open(header, "r") as file:
            first_line = file.read().rstrip("\n")
            delimiter = None
            try:
                for delim in possible_delimiters:
                    if delim in first_line:
                        delimiter = delim
                        break
                if delim is None:
                    raise ValueError(
                        "No delimiter found in the header file. Please ensure your header file is either in CSV or TSV format."
                    )
            except ValueError as e:
                raise

        file = open(header, "r")
        header = file.readline().rstrip("\n").split(delimiter)
        header = self.__check_headers(header)
        return header

    def __check_headers(self, header):
        req_columns_set = set(self.cols["general"])
        if set(req_columns_set).issubset(header):
            return header
        else:
            missing = list(req_columns_set - set(header))
            typer.secho(
                f"Header file is missing the following required column names: {missing}",
                fg=typer.colors.RED,
            )
            raise typer.Abort()


class RulesFile:
    def __init__(self, file_path):
        self.file_path = file_path
        self.data_frame = self.__read_json_to_dataframe()

    def __read_json_to_dataframe(self):
        try:
            data_frame = pd.read_json(self.file_path)
            data_frame.replace("", "none", inplace=True)

            return data_frame
        except Exception as e:
            typer.secho(
                f"Error reading input Rules File: {e}. Please check that your input file is in standard JSON form.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()
