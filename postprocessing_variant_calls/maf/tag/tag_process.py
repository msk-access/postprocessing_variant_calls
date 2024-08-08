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
from contextlib import ExitStack
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
):
    # prep maf
    mafa = MAFFile(maf, separator)
    typer.secho(f"Tagging Maf with traceback columns", fg=typer.colors.BRIGHT_GREEN)
    mafa = mafa.tag("traceback")
    typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "by_variant_classification",
    help="Tag filtered MAF file by variant classifications and subset into individual text files.",
)
def by_variant_classification(
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
        help="filtered MAF file to split by annotations with",
    ),
    canonical_tx_ref: Path = typer.Option(
        ...,
        "--canonical_tx_ref",
        "-tx_ref",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Reference canonical transcript file",
    ),
    output_dir: Path = typer.Option(
        "output_dir",
        "--output_dir",
        "-o",
        help="Output Directory to export individual text files to.",
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    """
    Parse a dataframe of annotated variants
    tag them into exonic, silent, exonic nonpanel, silent nonpanel or dropped
    write out into individual TXT/MAF output files
    """

    def format_var(variant):
        """
        Helper function to convert named tuple dervied from pandas df into tsv
        """
        try:
            columns = map(lambda x: getattr(variant, x), MAF_TSV_COL_MAP.keys())
            return "\t".join(map(str, columns)) + "\n"
        except AttributeError:
            return
            missing_columns = set(MAF_TSV_COL_MAP.keys()) - set(
                filter(lambda x: not x.startswith("_"), dir(variant))
            )
            raise Exception(
                "Missing required columns: {}".format(",".join(missing_columns))
            )

    def reformat_tx(txID, tx_df):
        """
        helper function to get reportable txID
        """
        return tx_df[tx_df.isoform == txID].refseq_id.values.tolist()

    # prep and read in maf
    typer.secho(f"Reading in input Filtered MAF file.", fg=typer.colors.BRIGHT_GREEN)
    mafa = MAFFile(maf, separator)

    # prep and read in input canonical reference transcript TSV flie
    typer.secho(
        f"Reading in input Reference canonical transcript file",
        fg=typer.colors.BRIGHT_GREEN,
    )
    tx_df, tx_isoform_lst = read_tsv(
        canonical_tx_ref, separator, canonical_tx_ref_flag=True
    )

    # start tagging by variant classification process
    final_maf = mafa.tag_by_variant_classification(output_dir, tx_isoform_lst)

    # Create exonic, silent, and nonpanel files.
    file_names = [
        EXONIC_FILTERED,
        SILENT_FILTERED,
        NONPANEL_EXONIC_FILTERED,
        NONPANEL_SILENT_FILTERED,
        DROPPED,
    ]

    file_paths = [f"{output_dir}/{name}" for name in file_names]

    headers = "\t".join(MAF_TSV_COL_MAP.values()) + "\n"

    with ExitStack() as stack:
        files = [stack.enter_context(open(path, "w")) for path in file_paths]
        for file in files:
            file.write(headers)

        # Iterate through the DataFrame rows
        for variant in final_maf.itertuples(index=False):
            tag = variant.classification

            # Reformat the row
            formatted_variant = format_var(variant)

            # Determine which file(s) to write to based on the tag
            if "exonic" in tag:
                variant_tuple = (
                    variant.Hugo_Symbol,
                    variant.Variant_Classification,
                    variant.vcf_pos,
                )
                transcript_id = reformat_tx(variant.Transcript_ID, tx_df)

                variant = variant._replace(
                    Hugo_Symbol=variant_tuple[0],  # Gene
                    Entrez_Gene_Id=variant.Entrez_Gene_Id,
                    Center=variant.Center,
                    NCBI_Build=variant.NCBI_Build,
                    Chromosome=variant.Chromosome,
                    Variant_Classification=variant_tuple[1],  # VariantClass
                    Variant_Type=variant.Variant_Type,
                    vcf_pos=variant_tuple[2],  # Start coordinate
                    Tumor_Sample_Barcode=variant.Tumor_Sample_Barcode,
                    caller_Norm_Sample_Barcode=variant.caller_Norm_Sample_Barcode,
                    Transcript_ID=transcript_id,  # TranscriptID
                )
                formatted_exonic_variant_row = format_var(variant)
                files[0].write(formatted_exonic_variant_row)
            if "silent" in tag:
                files[1].write(formatted_variant)
            if "nonpanel_exonic" in tag:
                files[2].write(formatted_variant)
            if "nonpanel_silent" in tag:
                files[3].write(formatted_variant)
            if "dropped" in tag:
                files[4].write(formatted_variant)

    return 0


if __name__ == "__main__":
    app()
