#!/usr/bin/env python
# imports
import os
import sys
import csv

from pathlib import Path
from typing import List, Optional
from .resources import acceptable_extensions, possible_maf_columns, minimal_maf_columns
import typer
import pandas as pd

def process_paths(paths): 
    file = open(paths, 'r')
    files = []
    for line in file.readlines():
        files.append(line.rstrip('\n'))
    file.close
    return files

def process_header(header):
    file = open(header, 'r')
    line = file.readline().split(',')
    file.close
    return line

def check_maf(files: List[Path]):
    # return non if argument is empty
    if files is None:
        return None
    # check that we have a list of mafs after reading off the cli
    extensions = [os.path.splitext(f)[1] for f in files]
    for ext in extensions:
        if ext not in acceptable_extensions:
            typer.secho(f"If using files argument, all files must be mafs.", fg=typer.colors.RED)
            raise typer.Abort()
    return files

def check_txt(paths: Path):
    # return None if argument is empty, kind of unfortunate we have to handle this case
    if paths is None:
        return None
    # check that we have a text file after reading off the cli
    extension = os.path.splitext(paths)[1]
    if extension != '.txt':
        typer.secho(f"If using paths argument, must provided an input txt file.", fg=typer.colors.RED)
        raise typer.Abort()
    return paths

def check_headers(maf, header):
    if set(header).issubset(set(maf.columns)):
            return maf[header]
    else:
        typer.secho(f"{maf} is not a subset of {header}. Please provide custom header file if the provided a header file or edit the current header file if you maf uses different columns names", 
                    fg=typer.colors.RED)
        raise typer.Abort()
    
def concat_mafs(files, output_maf, header):
    #TODO appending to empty data frame is slow, we sould make list of frames and flatten
    maf_list = []
    for maf in files:
        if Path(maf).is_file():
            # Read maf file
            typer.secho(f"Reading: {maf}", fg=typer.colors.BRIGHT_GREEN)
            maf_df = pd.read_csv(maf, sep="\t", low_memory=True)
            # custom header 
            if header:
                maf_col_df = check_headers(maf_df, header)
            # default header
            else:
                maf_col_df = check_headers(maf_df, minimal_maf_columns)
            maf_list.append(maf_col_df)
        else:
            typer.secho(f"failed to open {maf}", fg=typer.colors.RED)
            raise typer.Abort()
        
    # merge mafs
    merged_mafs = pd.concat(maf_list, axis=0, ignore_index=True)
    typer.secho(
        f"Done processing the concatenation of maf files writing output to: {output_maf}",
        fg=typer.colors.GREEN,
    )
    # write final df and return 0
    merged_mafs.to_csv(f"{output_maf}.maf", index=False, sep="\t")
    return 0
