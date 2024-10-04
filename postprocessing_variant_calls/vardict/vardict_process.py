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
from .vardict_class import var_sample

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

logger = logging.getLogger("filter")

app = typer.Typer(help="post-processing commands for VarDict version 1.4.6 VCFs.")
# single filter
single_app = typer.Typer()
app.add_typer(
    single_app,
    name="single",
    help="Post-processing commands for a single sample VarDict version 1.4.6 VCFs",
)
# case control filter
case_control_app = typer.Typer()
app.add_typer(
    case_control_app,
    name="case-control",
    help="Post-processing commands for a case-controlled VarDict version 1.4.6 VCFs",
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
        help="Input vcf generated by vardict which needs to be processed",
    ),
    sampleName: str = typer.Option(
        ...,
        "--tsampleName",
        help="Name of the tumor Sample",
    ),
    totalDepth: int = typer.Option(
        20,
        "--totalDepth",
        "-dp",
        min=20,
        help="Tumor total depth threshold",
    ),
    alleleDepth: int = typer.Option(
        "",
        "--alleledepth",
        "-ad",
        min=1,
        clamp=True,
    ),
    tnRatio: int = typer.Option(
        1,
        "--tnRatio",
        "-tnr",
        help="Tumor-Normal variant fraction ratio threshold",
    ),
    variantFraction: float = typer.Option(
        5e-05,
        "--variantFraction",
        "-vf",
        help="Tumor variant fraction threshold",
    ),
    minQual: int = typer.Option(
        0,
        "--minQual",
        "-mq",
        help="Minimum variant call quality",
    ),
    filterGermline: bool = typer.Option(
        False,
        "--filterGermline",
        "-fg",
        help="Whether to remove calls without 'somatic' status",
    ),
    outputDir: str = typer.Option(
        "", "--outDir", "-o", help="Full Path to the output dir"
    ),
):
    """
    This tool helps to filter vardict version 1.4.6 VCFs for single sample calling
    """
    logger.info("process_vardict: Started the run for doing standard filter.")
    # single sample case
    # create vardict object
    to_filter = var_sample(
        inputVcf,
        outputDir,
        sampleName,
        minQual,
        totalDepth,
        alleleDepth,
        variantFraction,
        tnRatio,
        filterGermline,
    )
    
    # check for normal
    if to_filter.has_normal():
        logger.exception(
            "`single filter` command was used, but a normal variant was detected. Did you mean to use `case_control filter`?"
        )
    else:
        # filter single
        vcf_out, vcf_complex_out, txt_out = to_filter.filter_single()
        vcf_out_sort = to_filter.sort_vcf()
        vcf_complex_out_sort = to_filter.sort_vcf_complex()
    return vcf_out_sort, vcf_complex_out_sort, txt_out


@case_control_app.command("filter")
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
        help="Input vcf generated by vardict which needs to be processed",
    ),
    sampleName: str = typer.Option(
        ...,
        "--tsampleName",
        help="Name of the tumor Sample",
    ),
    totalDepth: int = typer.Option(
        20,
        "--totalDepth",
        "-dp",
        min=20,
        help="Tumor total depth threshold",
    ),
    alleleDepth: int = typer.Option(
        "",
        "--alleledepth",
        "-ad",
        min=1,
        clamp=True,
    ),
    tnRatio: int = typer.Option(
        1,
        "--tnRatio",
        "-tnr",
        help="Tumor-Normal variant fraction ratio threshold",
    ),
    variantFraction: float = typer.Option(
        5e-05,
        "--variantFraction",
        "-vf",
        help="Tumor variant fraction threshold",
    ),
    minQual: int = typer.Option(
        0,
        "--minQual",
        "-mq",
        help="Minimum variant call quality",
    ),
    filterGermline: bool = typer.Option(
        False,
        "--filterGermline",
        "-fg",
        help="Whether to remove calls without 'somatic' status",
    ),
    outputDir: str = typer.Option(
        "", "--outDir", "-o", help="Full Path to the output dir"
    ),
):
    """
    This tool helps to filter vardict version 1.4.6 VCFs for case control calling
    """
    logger.info("process_vardict: Started the run for doing standard filter.")
    # single sample case
    # create vardict object
    to_filter = var_sample(
        inputVcf,
        outputDir,
        sampleName,
        minQual,
        totalDepth,
        alleleDepth,
        variantFraction,
        tnRatio,
        filterGermline,
    )
    # check for normal
    if to_filter.has_normal():
        # filter with normal
        vcf_out, vcf_complex_out, txt_out = to_filter.filter_case_control()
        vcf_out_sort = to_filter.sort_vcf()
        vcf_complex_out_sort = to_filter.sort_vcf_complex()
    else:
        # raise exception
        logger.exception(
            "`case_control filter` command was used, but no normal variant was detected. Did you mean to use `single filter`?"
        )
    return vcf_out_sort, vcf_complex_out_sort, txt_out


if __name__ == "__main__":
    app()
