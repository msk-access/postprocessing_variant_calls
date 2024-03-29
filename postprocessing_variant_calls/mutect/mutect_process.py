#!/usr/bin/env python
# imports
import os
import sys
import pandas as pd
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
from .mutect_class import mutect_sample
import typer


app = typer.Typer(help="post-processing commands for MuTect version 1.1.5 VCFs.")
# paired sample filter
single_app = typer.Typer()
app.add_typer(
    single_app,
    name="case-control",
    help="Post-processing commands for case-control filtering of MuTect version 1.1.5 VCF input file.",
)


@single_app.command("filter")
def filter(
    inputVcf: Path = typer.Option(
        ...,
        "--inputVcf",
        "-i",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Input vcf generated by MuTect which needs to be processed",
    ),
    inputTxt: Path = typer.Option(
        ...,
        "--inputTxt",
        "-i",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Input Txt file generated by MuTect which needs to be processed",
    ),
    refFasta: Path = typer.Option(
        ...,
        "--refFasta",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Input reference fasta",
    ),
    tsampleName: str = typer.Option(
        ..., "--tsampleName", help="Name of the tumor sample."
    ),
    totalDepth: int = typer.Option(
        20,
        "--totalDepth",
        "-dp",
        min=0,
        clamp=True,
        help="Tumor total depth threshold",
    ),
    alleleDepth: int = typer.Option(
        1,
        "--alleledepth",
        "-ad",
        min=0,
        clamp=True,
    ),
    tnRatio: int = typer.Option(
        1,
        "--tnRatio",
        "-tnr",
        min=0,
        clamp=True,
        help="Tumor-Normal variant fraction ratio threshold",
    ),
    variantFraction: float = typer.Option(
        0.00005,
        "--variantFraction",
        "-vf",
        min=0,
        clamp=True,
        help="Tumor variant fraction threshold",
    ),
    outputDir: str = typer.Option(
        "", "--outDir", "-o", help="Full Path to the output dir"
    ),
):
    """
    This tool helps to filter MuTect version 1.1.5 VCFs for case-control calling
    """
    typer.secho(
        f"process_mutect1: Started the run for doing standard filter.",
        fg=typer.colors.BLACK,
    )
    # single sample case
    # create mutect object
    to_filter = mutect_sample(
        inputVcf,
        inputTxt,
        refFasta,
        outputDir,
        tsampleName,
        totalDepth,
        alleleDepth,
        variantFraction,
        tnRatio,
    )

    # check for normal
    if to_filter.has_tumor_and_normal_cols():
        vcf_out, txt_out = to_filter.filter_paired_sample()
        typer.secho(
            f"process_mutect1: Filtering for MuTect VCFs has completed. Please refer to output directory for filtered MuTect VCF and text file.",
            fg=typer.colors.BLACK,
        )
    else:
        typer.secho(
            f"Tumor and normal columns not identified in input MuTect VCF file. Please check input file again.",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    return vcf_out, txt_out


if __name__ == "__main__":
    app()
