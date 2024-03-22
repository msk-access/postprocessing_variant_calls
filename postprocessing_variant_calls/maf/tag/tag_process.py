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
    "by_rules",
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
    tagged_by_rules_maf = mafa.tag_by_rules(rules_file.data_frame)

    typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)
    tagged_by_rules_maf.to_csv(
        f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t"
    )
    return 0


@app.command(
    "by_access",
    help="Tag a variant in a MAF file based on criterion stated by the SNV/indels ACCESS pipeline workflow",
)
def by_access(
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
    # run tag by_rules
    mafa = MAFFile(maf, separator)
    typer.secho(
        f"Tagging Maf with criteria from input rules JSON file",
        fg=typer.colors.BRIGHT_GREEN,
    )
    rules_file = RulesFile(rules)
    tagged_by_rules_maf = mafa.tag_by_rules(rules_file.data_frame)

    # run tag_by_hotspots

    # typer.secho(f"Writing Delimited file: {output_maf}", fg=typer.colors.BRIGHT_GREEN)


if __name__ == "__main__":
    app()
