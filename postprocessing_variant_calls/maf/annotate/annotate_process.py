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
from .annotate_helpers import read_bed, read_maf
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

logger = logging.getLogger("annotate")

app = typer.Typer(help="post-processing command annotating maf file by bed a bed file.")
# app.add_typer(app, name="annotate", help="annotate maf files based on a given input. Currently supports bed and maf files as references.")


@app.command("mafbybed", help="annotate a maf column by a bed file.")
def maf_bed(
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
        help="input maf file",
    ),
    bed: Path = typer.Option(
        ...,
        "--bed",
        "-b",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="bed file to annotate maf",
    ),
    outputFile: str = typer.Option(
        "output.maf",
        "--output",
        "-o",
        help="output maf file",
    ),
    cname: str = typer.Option(
        "annotation",
        "--cname",
        "-c",
        help="name for annotation column",
    ),
):
    logger.info("starting annotation")
    typer.secho(
        f"annotating {maf} with {bed}".format(maf=maf, bed=bed), fg=typer.colors.GREEN
    )
    # read bed
    bed_df = read_bed(bed)
    # read maf
    maf_df = read_maf(maf)
    # annotate maf with processed bed file
    annotated_maf = annotater(maf_df, bed_df, cname)
    # write to csv
    typer.secho(
        f"Writing out maf file to the following location: {outputFile}.maf".format(
            outputFile=outputFile
        ),
        fg=typer.colors.GREEN,
    )
    annotated_maf.to_csv(
        f"{outputFile}".format(outputFile=outputFile), index=False, sep="\t"
    )
    return 1


@app.command("mafbytsv", help="annotate a maf column by a bed file.")
def maf_tsv(
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
    tsv: Path = typer.Option(
        ...,
        "--tsv",
        "-t",
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
    oc: str = typer.Option(
        "hotspot",
        "--outcome_column",
        "-oc",
        help="name for outcome column",
    ),
    values: Tuple[str, str] = typer.Option(
        ("yes", "no"),
        "--values",
        "-v",
        help="name for annotation column. Defaults to (Yes, No)",
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)

    # prep tsv
    tsva = read_tsv(tsv, separator)
    tsva = gen_id_tsv(tsva)

    # annotate maf with processed bed file
    annotated_maf = mafa.annotate_maf_maf(tsva, oc, values)
    # added line of code to remove ID column
    annotated_maf = annotated_maf.drop("id", axis=1)
    annotated_maf["Start_Position"] = annotated_maf["Start_Position"].astype(int)
    annotated_maf["End_Position"] = annotated_maf["End_Position"].astype(int)
    # added line of code to convert start/end coordinates from float to integer

    # write out paths
    annotated_maf.to_csv(output_maf, index=False, sep="\t")
    return 0


@app.command(
    "extract_blocklist",
    help="Extract values from an optional blocklist file if provided. Used in SNVs/indels workflow.",
)
def extract_blocklist(
    tsv: Path = typer.Option(
        ...,
        "--blocklist_file",
        "-b",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Blocklist text file to extract values from. Needs to be in TSV format",
    ),
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
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a separator for delimited data.",
        callback=check_separator,
    ),
):
    # reading in input blocklist file
    header = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele",
        "Annotation",
    ]
    mafa = MAFFile(maf, separator)
    tsva = read_tsv(tsv, separator)

    # processing the input blocklist file and extracting the blocklist values, returning a list
    if tsva.empty:
        blacklist = []
        print("Extracted blocklist values for input blocklist file.")
        # performing actual filtering using the MAF read in object
        return blacklist
    elif tsva.empty == False:
        if list(tsva.columns.values) != header:
            raise Exception(
                "Blacklist provided is in the wrong format, file should have the following in the header (in order):"
                + ", ".join(header)
            )
        else:
            tsva.drop(["Annotation"], axis=1, inplace=True)
            tsva.drop_duplicates(inplace=True)
            blacklist = [
                str(b[0])
                + "_"
                + str(b[1])
                + "_"
                + str(b[2])
                + "_"
                + str(b[3])
                + "_"
                + str(b[4])
                for b in tsva.values.tolist()
            ]
            print("Extracted blocklist values for input blocklist file.")
            # performing actual filtering using the MAF read in object
            return blacklist
    else:
        raise IOError(
            "Blacklist file provided does not exist. Please check inputs again."
        )


if __name__ == "__main__":
    app()
