#! python
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
from vcf.parser import _Info as VcfInfo, _Format as VcfFormat, _vcf_metadata_parser as VcfMetadataParser
from vardict_class import vardict

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

logger = logging.getLogger("process_vardict")

app = typer.Typer()

@app.command()
def process_vardict(
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
    refFasta: Path = typer.Option(
        ...,
        "--refFasta",
        "-rf",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Reference genome in fasta format",
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
        ..., "--outDir", "-o", help="Full Path to the output dir"
    ),
    ## TODO: instead boolean, maybe string where you specify single vs double, maybe true or false, required
    normalFlag: bool = typer.Option(
        False,
        "--normalFlag",
        "-n", 
        help="Indicate whether a normal sample is present in the Vardict",
    ),
):

    '''
    @Description : This tool helps to filter vardict version 1.4.6 vcf for matched calling
    @Created : 03/23/2022
    @author : Ronak H Shah


    Visual representation of how this module works:

    "Somatic" not in record['STATUS'] and filter_germline ?
    |
    yes --> DONT KEEP
    |
    no --> tumor_variant_fraction > nvfRF ?
            |
            no --> DONT KEEP
            |
            yes --> tmq >= minQual and
                    nmq >= minQual and
                    tdp >= totalDepth and
                    tad >= alleleDepth and
                    tvf >= variantFraction ?
                    |
                    no --> DONT KEEP
                    |
                    yes --> KEEP

    Note: BasicFiltering VarDict's additional filters over MuTect include:
    1. Tumor variant quality threshold
    2. Normal variant quality threshold
    3. Somatic filter (MuTect does not report germline events)
    '''
    logger.info("process_vardict: Started the run for doing standard filter.")
    if normalFlag: 
        # TODO: QC double case 
        # create vardict object 
        to_filter = vardict(
                inputVcf, outputDir, sampleName, minQual, totalDepth, 
                alleleDepth, variantFraction, tnRatio, filterGermline
        )
        # check for normal
        if to_filter.has_normal():
            # filter vardict with tumor/normal
            vcf_out, vcf_complex_out, txt_out = to_filter.filter_two()
        else: 
            logger.exception('normalFlag was set to True without a normal sample present in the vardict vcf.')
    else: 
        # TODO: contiue work on single case, check with Karthi / Ronak about single filter  
        to_filter = vardict(
                inputVcf, outputDir, sampleName, minQual, totalDepth, 
                alleleDepth, variantFraction, tnRatio, filterGermline
        )
        # check for normal
        if to_filter.has_normal():
            logger.exception('normalFlag was set to False with a normal sample present in the vardict vcf.')
        else:  
            vcf_out, vcf_complex_out, txt_out = to_filter.filter_single()
    logger.info("process_vardict: Finished the run for doing vcf processing.")
    return vcf_out, vcf_complex_out, txt_out


if __name__ == "__main__":
    start_time = time.time()
    app()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("process_vardict: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
