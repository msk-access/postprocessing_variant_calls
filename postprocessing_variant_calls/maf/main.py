# imports
from postprocessing_variant_calls.maf.helper import (
    MAFFile,
    check_maf,
    check_txt,
    process_paths,
    maf_duplicates,
)
from .subset.subset_helpers import read_tsv, read_ids, filter_by_rows, check_separator
import typer
from importlib import resources
from pathlib import Path
from typing import List, Optional
import typer
import time
import os
import pandas as pd

# dir path
dir_path = os.path.dirname(os.path.realpath(__file__))

# Main app
app = typer.Typer()

# Declared apps
from .annotate import annotate_process
from .tag import tag_process
from .filter import filter_process

app.add_typer(
    annotate_process.app,
    name="annotate",
    help="annotate maf files based on a given input.",
)
app.add_typer(tag_process.app, name="tag", help="tag maf files based on a given input.")
app.add_typer(
    filter_process.app, name="filter", help="filter maf files based on a given input."
)


# Concat Features: This needs to be in main since we it doesn't have sub-commands
@app.command("concat", help="row-wise concatenation for maf files.")
def concat(
    files: List[Path] = typer.Option(
        [],
        "--files",
        "-f",
        help="MAF file to concatenate. Default assumes MAFs are tsv. MAF inputs are specified here, or using paths parameter",
        callback=check_maf,  # call back allow us to check input parameters
    ),
    paths: Path = typer.Option(
        None,
        "--paths",
        "-p",
        help="A text file containing paths of maf files to concatenate. Default assumes MAFs are tsv. MAF files are specified here, or using files parameter.",
        callback=check_txt,  # call back allow us to check input parameters
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    header: Path = typer.Option(
        None,
        "--header",
        "-h",
        help="A header file containing the columns to concatenate input mafs on. \
              It must be a subset of: \
              Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2. \
              These are also the default columns used for concatenation",
        callback=check_txt,
    ),
    deduplicate: bool = typer.Option(
        False,
        "--deduplicate",
        "-de",
        help="deduplicate outputted maf file.",
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # option to get files from text file
    if paths:
        files = process_paths(paths)
    if files:
        # create maf files
        maf_list = []
        for maf in files:
            maf = MAFFile(maf, separator, header)
            maf_list.append(maf.data_frame)
        # concat
        typer.secho(
            f"Concatenating maf files.",
            fg=typer.colors.BRIGHT_GREEN,
        )
        concat_df = pd.concat(maf_list, axis=0, ignore_index=True)
        # deduplicate
        if deduplicate:
            concat_df = maf_duplicates(concat_df)
        # write out paths
        # code to remove id column and verify it doesn't exist before writing output maf
        if "id" in list(concat_df.columns):
            try:
                concat_df = concat_df.drop(columns=["id"])
                concat_df.to_csv(output_maf, index=False, sep="\t")
            except:
                concat_df.to_csv(output_maf, index=False, sep="\t")
    else:
        typer.echo("--files or --path argument must be provided.")
        raise typer.Abort()
    return 0


# Add Subset App
@app.command("subset", help="subset maf files.")
def subset_maf(
    maf: Path = typer.Option(
        None,
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
    ids: Path = typer.Option(
        "",
        "--ids",
        "-i",
        help="List of ids to search for in the 'Tumor_Sample_Barcode' column. Header of this file is 'sample_id'",
    ),
    sid: Optional[List[str]] = typer.Option(
        None,
        "--sid",
        help="Identifiers to search for in the 'Tumor_Sample_Barcode' column. Can be given multiple times",
    ),
    output_file: str = typer.Option(
        "output_subset.maf",
        "--output",
        "-o",
        help="Name of the output file",
    ),
    col_name: str = typer.Option(
        "Tumor_Sample_Barcode",
        "--cname",
        "-c",
        help="Name of the column header to be used for sub-setting",
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
    Subset MAF/TSV file by selected sample ID/IDs

    Tool to do the following operations:
    A. Get subset of variants based on Tumor_Sample_Barcode in data_mutations_extended.txt file
    B. Mark the variants as overlapping with BED file as covered [yes/no], by appending "covered" column to the subset MAF

    Requirement:
    pandas; typing; typer; bed_lookup(https://github.com/msk-access/python_bed_lookup)

    """
    typer.echo("Subset and Annotate MAF...")
    if not ids:
        typer.echo("Identifiers were not provided in a text file")
        if not sid:
            typer.echo("Identifiers were not provided via command line as well")
            raise typer.Abort()

    maf_df = read_tsv(maf, separator)
    ids_to_subset = read_ids(sid, ids)
    subset_maf = filter_by_rows(ids_to_subset, maf_df, col_name)
    typer.echo(f"Writing to {output_file}")
    subset_maf.drop_duplicates().to_csv(output_file, sep="\t", index=False)
    typer.echo("Done!")


@app.command("mergetsv", help="merge a tsv file onto a maf by a shared id column.")
def mergetsv(
    mafa: Path = typer.Option(
        None,
        "--mafa",
        "-ma",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    mafb: Path = typer.Option(
        None,
        "--mafb",
        "-mb",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF",
    ),
    output_maf: Path = typer.Option(
        "merged.maf", "--output", "-o", help="Maf output file name."
    ),
    id: str = typer.Option("id", "--merge_id", "-id", help="id to merge mafs on."),
    how: str = typer.Option(
        "left",
        "--how",
        "-h",
        help="Type of merge to be performed on mafs. Defaults to left.",
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
    mafa = MAFFile(mafa, separator)
    mafb = read_tsv(mafb, separator)
    maf = mafa.merge(mafb, id, how)
    maf.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


if __name__ == "__main__":
    app()
