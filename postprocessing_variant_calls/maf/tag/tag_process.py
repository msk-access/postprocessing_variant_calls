#!/usr/bin/env python
# imports
from __future__ import division
import os
import sys
import vcf
import time
import logging
from pathlib import Path
from typing import List, Optional
import typer
from vcf.parser import (
    _Info as VcfInfo,
    _Format as VcfFormat,
    _vcf_metadata_parser as VcfMetadataParser,
)
from postprocessing_variant_calls.maf.helper import (
    check_maf,
    check_txt,
    check_separator,
    read_tsv,
    MAFFile,
    gen_id_tsv,
)
from utils.pybed_intersect import annotater
import pandas as pd
import numpy as np
from typing import Tuple


logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

logger = logging.getLogger("tag")

app = typer.Typer(help="post-processing command tagging maf files")
# app.add_typer(app, name="annotate", help="annotate maf files based on a given input. Currently supports bed and maf files as references.")


@app.command(
    "germline_status", help="Tag a variant in a MAF file as germline based on VAF value"
)
def germline_status(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    typer.secho(
        f"Tagging Maf with germline_status columns", fg=typer.colors.BRIGHT_GREEN
    )
    mafa = mafa.tag("germline_status")
    typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "common_variant",
    help="Tag a variant in a MAF file as common variant based on GNOMAD AF",
)
def common_variant(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    typer.secho(
        f"Tagging Maf with common_variant columns", fg=typer.colors.BRIGHT_GREEN
    )
    mafa = mafa.tag("common_variant")
    typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "prevalence_in_cosmicDB",
    help="Tag a variant in a MAF file with prevalence in COSMIC DB ",
)
def prevalence_in_cosmicDB(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    typer.secho(
        f"Tagging Maf with prevalence_in_cosmicDB columns", fg=typer.colors.BRIGHT_GREEN
    )
    mafa = mafa.tag("prevalence_in_cosmicDB")
    typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "truncating_mut_in_TSG",
    help="Tag a  truncating mutating variant in a MAF file based on its presence in the Tumor Suppressor Gene ",
)
def truncating_mut_in_TSG(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    typer.secho(
        f"Tagging Maf with truncating_mut_in_TSG columns", fg=typer.colors.BRIGHT_GREEN
    )
    mafa = mafa.tag("truncating_mut_in_TSG")
    typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "cmo_ch", help="Tag a variant in MAF file based on all the parameters listed"
)
def cmo_ch(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    typer.secho(f"Tagging Maf with cmo_ch_tag columns", fg=typer.colors.BRIGHT_GREEN)
    mafa = mafa.tag_all("cmo_ch_tag")
    typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "traceback",
    help="Generate combined count columns between standard and simplex/duplex mafs",
)
def traceback(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to tag",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
    samplesheet: List[Path] = typer.Option(
        None,
        "--samplesheet",
        "-sheet",
        help="Samplesheets in nucleovar formatting. See README for more info: `https://github.com/mskcc-omics-workflows/nucleovar/blob/main/README.md`. Used to add fillout type information to maf. The `sample_id` and `type` columns must be present.",
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)

    # Tag columns for traceback
    typer.secho(f"Tagging Maf with traceback columns", fg=typer.colors.BRIGHT_GREEN)
    mafa = mafa.tag("traceback")

    pd_samplesheet = []
    if samplesheet:
        for sheet in samplesheet:
            s = pd.read_csv(sheet)
            required_columns = ["sample_id", "type"]
            missing_columns = [col for col in required_columns if col not in s.columns]
            if len(missing_columns) == 0:
                pd_samplesheet.append(s)
            else:
                typer.secho(
                    f"Samplesheet is missing required column(s): {missing_columns}",
                    fg=typer.colors.RED,
                )
                raise typer.Abort()

        # Concatenate samplesheets
        combine_samplesheet = pd.concat(pd_samplesheet, ignore_index=True, sort=False)
        combine_samplesheet.fillna("", inplace=True)
        combine_samplesheet = combine_samplesheet[["sample_id", "type"]]

        # add in sample category columns via left merge
        typer.secho(f"Adding fillout type column", fg=typer.colors.BRIGHT_GREEN)
        mafa = pd.merge(
            mafa,
            combine_samplesheet,
            how="left",
            left_on="Tumor_Sample_Barcode",
            right_on="sample_id",
        )
        mafa.drop(columns=["sample_id"], inplace=True)
        mafa.rename(columns={"type": "fillout_type"}, inplace=True)

    # write out to csv file
    typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


# @app.command(
#     "access_filters_tag",
#     help="Tag a variant in a MAF file based on the following criteria in the access_filters python script",
# )
# def access_filters_tag(
#     df_pre_filter: Path = typer.Option(
#     ..., "--df_pre_filter",help="Pre Filtered MAF dataframe."
#     )
# ):

#     typer.secho(
#         f"Tagging Maf with germline columns", fg=typer.colors.BRIGHT_GREEN
#     )
# code for tagging for germline (if there is a matched normal and it has sufficient coverage)
# return df_pre_filter
# typer.secho(
#     f"Tagging Maf with below alt threshold values", fg=typer.colors.BRIGHT_GREEN
# )
# code for tagging for below alt threshold
# mafa = mafa.tag("below_alt_thres")
# typer.secho(
#     f"Tagging Maf for occurrence in curated samples", fg=typer.colors.BRIGHT_GREEN
# )
# code for tagging occurrence in curated
# mafa = mafa.tag("curated_occurrence")
# typer.secho(
#     f"Tagging Maf for occurrence in normal samples", fg=typer.colors.BRIGHT_GREEN
# )
# code for tagging occurrence in normal
# mafa = mafa.tag("normal_occurrence")
# typer.secho(
#     f"Tagging Maf for occurrence in blacklist", fg=typer.colors.BRIGHT_GREEN
# )
# code for tagging occurrence in blacklist
# mafa = mafa.tag("in_blacklist")


# typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
# mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
# return 0

if __name__ == "__main__":
    app()
