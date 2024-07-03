#!/usr/bin/env python
# imports
import os
import sys
import csv

from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd
import numpy as np
from .resources import tsg_genes


def process_paths(paths):
    file = open(paths, "r")
    files = []
    for line in file.readlines():
        files.append(line.rstrip("\n"))
    file.close
    return files


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


def read_tsv(tsv, separator):
    """Read a tsv file

    Args:
        maf (File): Input MAF/tsv like format file

    Returns:
        data_frame: Output a data frame containing the MAF/tsv
    """
    typer.echo("Read Delimited file...")
    skip = get_row(tsv)
    return pd.read_csv(tsv, sep=separator, skiprows=skip, low_memory=True)


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

def extract_fillout_type(df_full_fillout):
    #tag fillout type in full fillout
    df_full_fillout['Fillout_Type'] = df_full_fillout['Tumor_Sample_Barcode'].apply(lambda x: "CURATED-SIMPLEX" if "-CURATED-SIMPLEX" in x else ("SIMPLEX" if "-SIMPLEX" in x else ("CURATED-DUPLEX" if "CURATED-DUPLEX" in x else "POOL")))
    #extract curated duplex and and VAf and summary column
    df_curated = (df_full_fillout.loc[df_full_fillout['Fillout_Type'] == 'CURATED-DUPLEX'])
    df_curated = __find_VAFandsummary(df_curated)
    breakpoint()
    #extract simplex and transform to simplex +duplex
    df_curatedsimplex = (df_full_fillout.loc[df_full_fillout['Fillout_Type'] == 'CURATED-SIMPLEX'])
    df_ds_curated = __create_duplexsimplex(df_curatedsimplex, df_curated)
    #extract duplex and and VAf and summary column
    df_pool = (df_full_fillout.loc[df_full_fillout['Fillout_Type'] == 'POOL'])
    df_pool = __find_VAFandsummary(df_pool)
    #extract simplex and transform to simplex +duplex
    df_simplex = (df_full_fillout.loc[df_full_fillout['Fillout_Type'] == 'SIMPLEX'])
    df_ds_tumor = __create_duplexsimplex(df_simplex, df_pool)
    #separate tumor samples from Normal samples from POOL Fillout_Type
    #This is done by assuming all tumor samples have a simplex-duplex tumor sample genotyped, while normals do not. Samples not in DS_tumor are assumed to be normals
    df_tumor, df_normal= __separate_normal_and_tumor(df_pool, df_ds_tumor)
    return df_tumor, df_normal, df_ds_tumor, df_curated, df_ds_curated

def __find_VAFandsummary(df_fillout): 
    df_fillout = df_fillout.copy()
    #find the VAF from the fillout
    df_fillout['t_vaf_fragment'] = (df_fillout['t_alt_count'] / (df_fillout['t_alt_count'].astype(int) + df_fillout['t_ref_count'].astype(int))).round(4)
    df_fillout['summary_fragment'] = 'DP='+(df_fillout['t_alt_count'].astype(int) + df_fillout['t_ref_count'].astype(int)).astype(str)+';RD='+  df_fillout['t_ref_count'].astype(str)+';AD='+ df_fillout['t_alt_count'].astype(str)+';VF='+df_fillout['t_vaf_fragment'].fillna(0).astype(str)
    return df_fillout
    
def __create_duplexsimplex(df_s, df_d):
    df_s = df_s.copy()
    df_d = df_d.copy()
    #Prep Simplex
    df_s.rename(columns = {'t_alt_count_fragment': 't_alt_count_fragment_simplex','t_ref_count_fragment':'t_ref_count_fragment_simplex'}, inplace=True)   
    df_s['Tumor_Sample_Barcode'] = df_s['Tumor_Sample_Barcode'].str.replace('-SIMPLEX','')
    df_s.set_index('Tumor_Sample_Barcode', append=True, drop=False, inplace=True)
    #Prep Duplex
    df_d.rename(columns = {'t_alt_count_fragment': 't_alt_count_fragment_duplex','t_ref_count_fragment':'t_ref_count_fragment_duplex'}, inplace=True)
    df_d['Tumor_Sample_Barcode'] = df_d['Tumor_Sample_Barcode'].str.replace('-DUPLEX', '')
    df_d.set_index('Tumor_Sample_Barcode', append=True, drop=False, inplace=True)
    #Merge
    df_ds = df_s.merge(df_d[['t_ref_count_fragment_duplex','t_alt_count_fragment_duplex']], left_index=True, right_index=True)
    ##Add
    df_ds['t_ref_count_fragment'] = df_ds['t_ref_count_fragment_simplex'] + df_ds['t_ref_count_fragment_duplex']
    df_ds['t_alt_count_fragment'] = df_ds['t_alt_count_fragment_simplex'] + df_ds['t_alt_count_fragment_duplex']
    df_ds['t_total_count_fragment'] = df_ds['t_alt_count_fragment'] + df_ds['t_ref_count_fragment']
    ##clean up
    fillout_type = df_ds['Fillout_Type']+'-DUPLEX'
    df_ds.drop(['Fillout_Type', 't_ref_count_fragment_simplex', 't_ref_count_fragment_duplex', 't_alt_count_fragment_simplex','t_alt_count_fragment_duplex'], axis=1, inplace=True)
    df_ds['Fillout_Type'] = fillout_type
    df_ds['Tumor_Sample_Barcode'] = df_ds['Tumor_Sample_Barcode']+'-SIMPLEX-DUPLEX'
    df_ds.set_index(mutation_key, drop=False, inplace=True)
    df_ds = find_VAFandsummary(df_ds)
    return df_ds
    
def __separate_normal_and_tumor (df_pool, df_ds_tumor):
    tumor_samples=[f.replace('-SIMPLEX-DUPLEX', '') for f in df_ds_tumor.Tumor_Sample_Barcode.unique().tolist()]
    df_tumor=(df_pool.loc[df_pool['Tumor_Sample_Barcode'].isin(tumor_samples)])
    df_normal=(df_pool.loc[~df_pool['Tumor_Sample_Barcode'].isin(tumor_samples)])
    df_normal.Fillout_Type='NORMAL'
    return df_tumor, df_normal
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
                "Annotation"
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
        
    def _convert_annomaf_to_df(self):
        if self.data_frame.empty == False:
            self.data_frame['Chromosome'] = self.data_frame['Chromosome'].astype(str)
            self.data_frame.set_index(self.cols['general'], drop=False, inplace=True)
            self.data_frame.rename(columns ={'Matched_Norm_Sample_Barcode':'caller_Norm_Sample_Barcode','t_depth':'caller_t_depth','t_ref_count':'caller_t_ref_count', 't_alt_count':'caller_t_alt_count', 'n_depth':'caller_n_depth','n_ref_count':'caller_n_ref_count', 'n_alt_count':'caller_n_alt_count','set':'CallMethod'}, inplace=True)
        # quick cleanup of mutect columns (can be made into mini function in the MAF class)
    
        # marking the mutect and vardict combo columns if the condition in the line below is met.
            self.data_frame.loc[(self.data_frame['MUTECT'] == 1) & (self.data_frame['CallMethod'] != 'MuTect'),'CallMethod'] = 'VarDict,MuTect'
            self.data_frame.drop(['TYPE','FAILURE_REASON','MUTECT'], inplace=True, axis=1)
            return self.data_frame
        else:
            typer.secho(f"failed to open path to the annotation MAF file {self.file_path}.", fg=typer.colors.RED)
            raise typer.Abort()
    
    def _convert_fillout_to_df(self):
        if self.data_frame.empty == False:
            self.data_frame['Chromosome'] = self.data_frame['Chromosome'].astype(str)
            self.data_frame.set_index(self.cols['general'], drop=False, inplace=True)
            return self.data_frame
        else:
            typer.secho(f"failed to open path to the fillout MAF file {self.file_path}.", fg=typer.colors.RED)
            raise typer.Abort()
            
        
    def __check_delimiter(self):
        filename = self.file_path
        try:
            with open(filename, "r", newline="") as file_in:
                reader = csv.reader(file_in, delimiter=self.separator)
                headers = next(reader)

                column_count = len(headers)
                if column_count == 1:
                    raise ValueError("Column Count is 1")

        except Exception as e:
            typer.secho(f'file \"{filename}\" is probably not separated by \"{self.separator}"', fg=typer.colors.RED)
            typer.secho(f'trying \t', fg=typer.colors.RED)
            

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
            df = pd.read_csv(self.file_path, sep=None, skiprows=skip, low_memory=True)
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


        else:
            typer.secho(
                f"missing columns expected for {tagging} tagging expects: {set(cols).difference(set(self.data_frame.columns.tolist()))}, which was missing from the input",
                fg=typer.colors.RED,
            )
            raise typer.Abort()
        return self.data_frame

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
                ].combine_first(
                    self.data_frame["t_ref_count_fragment_simplex_duplex"]
                )
                self.data_frame["t_alt_count"] = self.data_frame[
                    "t_alt_count_standard"
                ].combine_first(
                    self.data_frame["t_alt_count_fragment_simplex_duplex"]
                )
                self.data_frame["t_total_count"] = self.data_frame[
                    "t_total_count_standard"
                ].combine_first(
                    self.data_frame["t_total_count_fragment_simplex_duplex"]
                )
        elif set(cols["standard"]).issubset(
            set(self.data_frame.columns.tolist())
        ):
            self.data_frame["t_ref_count"] = self.data_frame[
                "t_ref_count_standard"
            ]
            self.data_frame["t_alt_count"] = self.data_frame[
                "t_alt_count_standard"
            ]
            self.data_frame["t_total_count"] = self.data_frame[
                "t_total_count_standard"
            ]
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
                        (self.data_frame["common_variant"] == "yes")
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
        file = open(header, "r")
        header = file.readline().rstrip("\n").split(",")
        file.close
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
