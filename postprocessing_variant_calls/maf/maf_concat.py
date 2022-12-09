"""
Script to merge multiple maf files produced by variants callers as part of the postprocessing process.
Tool does the following operations:
A. Read one or more maf files generated by variant callers
B. Obtain rows from multi-sample maf files based on string containing sample name
C. Merge input mafs to single maf for downstream processing
This infile is a tab-delimited MAF file. For full-spec, see: https://github.com/mskcc/vcf2maf/blob/main/docs/vep_maf_readme.txt

"""
import os
import time
import argparse


from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd

def main(
    list_of_files: Path = typer.Option(
        "--list",
        "-l",
        help="File of files, List of maf files to be concatenated, one per line, no header",
    ),
    maf: Optional[List[Path]] = typer.Option(
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Maf files to be concatenated. Can be given multiple times",
    ),
    output_maf_file_prefix: str = typer.Option(
        "concat_maf_output",
        "--prefix",
        "-p",
        help="Prefix of the output MAF",
    ),
):
    if not list_of_files:
        typer.secho(
            "File are not provided as file of files.", fg=typer.colors.BRIGHT_YELLOW
        )
        if not maf:
            typer.secho(
                "File were not provided via command line as well",
                fg=typer.colors.BRIGHT_RED,
            )
            raise typer.Abort()
    # Read file of files
    if not maf:
        maf = [line.strip() for line in open(list_of_files, "r")]
    final_df = pd.DataFrame()
    for maf_file in maf:
        if Path(maf_file).is_file():
            # Read maf file
            typer.secho(f"Reading: {maf_file}", fg=typer.colors.BRIGHT_GREEN)
            maf_df = pd.read_csv(maf_file, sep='\t', low_memory=True)
            maf_col_df = maf_df[["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode"]]
            merged_mafs = pd.concat([maf_col_df], join='inner')
        else:
            continue
    else:
        typer.secho(f"{maf_file} file does not exists", fg=typer.colors.BRIGHT_RED)
        raise typer.Abort()
    # write concatanted df to maf
    typer.secho(
        f"Done processing the CSV file writing output to {output_maf_file_prefix} in txt and excel format",
        fg=typer.colors.GREEN,
    )
    final_df.to_csv(f"{output_maf_file_prefix}.maf", index=False, sep="\t")
    # Minimum required columns for maf format
    #Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Variant_Type, Tumor_Sample_Barcode
    # Concatenate maf files
    # Export merged maf
#    merged_mafs().to_csv(merged_maf, sep="\t", index=False)

if __name__ == "__main__":
    typer.run(main)