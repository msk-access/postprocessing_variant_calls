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
    RulesFile,
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
    sample_type: str = typer.Option(
        "--sample_type", "-t", help="String value of type of sample. Taken from meta map in traceback subworkflow."
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
    typer.secho(f"Tagging Maf with traceback columns", fg=typer.colors.BRIGHT_GREEN)
    mafa = mafa.tag("traceback")
    typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
    mafa_id_dropped = mafa.drop('id', axis=1)
    # adding in sample type as a separate column to differentiate between groups of samples
    #test_lst = [['DONOR36-TP', 'CURATED'], ['DONOR19-T', 'CURATED'], ['DONOR3-T', 'CURATED'], ['C-1V5P0L-N001-d01', 'UNMATCHED_NORMAL'], ['C-30N2P4-L004-d03', 'PLASMA'], ['C-H77MCH-L008-d08', 'PLASMA'], ['DONOR10-T', 'CURATED'], ['DONOR17-T', 'CURATED'], ['DONOR23-T', 'CURATED'], ['C-30N2P4-L003-d02', 'PLASMA'], ['C-PR83CF-L004-d04_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex', null], ['C-217F4D-N005-d05', 'UNMATCHED_NORMAL'], ['DONOR6-T', 'CURATED'], ['DONOR16-T', 'CURATED'], ['DONOR39-TP', 'CURATED'], ['C-13TFTU-N003-d02', 'UNMATCHED_NORMAL'], ['DONOR46-T', 'CURATED'], ['DONOR13-T', 'CURATED'], ['DONOR5-T', 'CURATED'], ['DONOR28-TP', 'CURATED'], ['C-001421-N001-d01', 'UNMATCHED_NORMAL'], ['C-2UW6JP-N001-d01', 'UNMATCHED_NORMAL'], ['DONOR30-TP', 'CURATED'], ['C-1DFUNN-N018-d03', 'UNMATCHED_NORMAL'], ['DONOR32-TP', 'CURATED'], ['C-2UW6JP-L008-d07', 'PLASMA'], ['DONOR38-TP', 'CURATED'], ['C-30N2P4-N001-d01', 'MATCHED_NORMAL'], ['C-09XLT2-N001-d01', 'UNMATCHED_NORMAL'], ['C-0D9K3A-N001-d01', 'UNMATCHED_NORMAL'], ['C-13TFTU-L010-d10', 'PLASMA']]
    
    temp_df = pd.DataFrame(sample_type, columns=['Tumor_Sample_Barcode', 'Type']) 
    
    mafa_final = mafa_id_dropped.merge(temp_df,on="Tumor_Sample_Barcode",how='left')
    mafa_final.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "by_variant_annotations",
    help="Tag a variant in a MAF file based on criterion provided by input rules JSON file",
)
def by_maf(
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
    rules: Path = typer.Option(
        ...,
        "--rules",
        "-r",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Intervals JSON file containing criterion to tag input MAF by",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a separator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    typer.secho(
        f"Tagging Maf with criteria from input rules JSON file",
        fg=typer.colors.BRIGHT_GREEN,
    )
    rules_file = RulesFile(rules)
    tagged_by_variant_annot_maf = mafa.tag_by_variant_annotations(rules_file.data_frame)

    typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
    tagged_by_variant_annot_maf.to_csv(
        f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t"
    )
    return 0


@app.command(
    "access",
    help="Tag a variant in a MAF file based on criterion stated by the SNV/indels ACCESS pipeline workflow",
)
def access(
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
    rules: Path = typer.Option(
        ...,
        "--rules",
        "-r",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Intervals JSON file containing criterion to tag input MAF by",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a separator for delimited data.",
        callback=check_separator,
    ),
):
    # run tag by_variant_annotations
    mafa = MAFFile(maf, separator)
    typer.secho(
        f"Tagging Maf with criteria from input rules JSON file",
        fg=typer.colors.BRIGHT_GREEN,
    )
    rules_file = RulesFile(rules)
    tagged_by_variant_annot_maf = mafa.tag_by_variant_annotations(rules_file.data_frame)

    # run tag by refseq IDs
    # run tag_by_hotspots
    # pv mag annotate mafbytsv
    # run tag by artifacts
    # xpv mag annotate mafbytsv

    # run tag by germ baseline (if applicable)

    # split into df_keep and df_drop
    # typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)


if __name__ == "__main__":
    app()
