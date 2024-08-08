from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd
import numpy as np
import sys

from postprocessing_variant_calls.maf.helper import (
    check_maf,
    check_txt,
    check_separator,
    read_tsv,
    MAFFile,
)

from postprocessing_variant_calls.maf.tag.tag_constants import (
    MAF_DUMMY_COLUMNS2,
    MAF_COLUMNS_SELECT,
    GNOMAD_COLUMNS,
)

app = typer.Typer()


def extract_blocklist(blocklist_file, separator):
    """
    The function `extract_blocklist` reads and processes a blocklist file, extracting specific values
    and returning them as a list.

    :param blocklist_file: The `blocklist_file` parameter is the file path to the blocklist file that
    contains information about genomic regions to be blocked. This function reads the contents of this
    file and processes it to extract the blocklist values based on the specified separator

    :param separator: The `separator` parameter in the `extract_blocklist` function is used to specify
    the delimiter that separates values in the input blocklist file. This delimiter is used when reading
    the file to properly parse the data into columns. Common separators include commas (`,`), tabs
    (`\t`), and spaces

    :return: The function `extract_blocklist` returns a list of blocklist values extracted from the
    input blocklist file after processing it.
    """
    # reading in input blocklist file
    header = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele",
        "Annotation",
    ]
    tsva = read_tsv(blocklist_file, separator)
    # processing the input blocklist file and extracting the blocklist values, returning a list
    if tsva.empty:
        blocklist = []
        # performing actual filtering using the MAF read in object
        return blocklist
    elif tsva.empty == False:
        if list(tsva.columns.values) != header:
            raise Exception(
                "Blacklist provided is in the wrong format, file should have the following in the header (in order):"
                + ", ".join(header)
            )
        else:
            tsva.drop(["Annotation"], axis=1, inplace=True)
            tsva.drop_duplicates(inplace=True)
            blocklist = [
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
            return blocklist
    else:
        raise IOError(
            "Blacklist file provided does not exist. Please check inputs again."
        )


def __generate_table_and_find_summary_stats(
    df_fillout, mutation_key, fillout_type, alt_thres, new_fillout_type
):
    """
    The function generates a summary table with statistics based on centered mutation key columns and sample
    information.

    :param df_fillout: `df_fillout` is a DataFrame containing data related to mutations and samples. It
    likely includes columns such as `mutation_key`, `Tumor_Sample_Barcode`, `summary_fragment`,
    `t_variant_frequency_standard`, `t_alt_count_standard`, `t_vaf_fragment`, and `t_alt

    :param mutation_key: The `mutation_key` parameter is a list of columns derived from the fillout MAF class object which is used as the index for the pivot table and for
    grouping the data to calculate summary statistics in the function
    `__generate_table_and_find_summary_stats`. It helps in organizing and summarizing the data based on
    the mutation key provided

    :param fillout_type: The `fillout_type` parameter in the function
    `__generate_table_and_find_summary_stats` is used to determine the type of data processing to be
    performed based on the input value. If `fillout_type` is equal to "NORMAL_", certain calculations
    and operations are carried out on the input data

    :param alt_thres: The `alt_thres` parameter is used to specify the threshold for the alternative
    count. It is used to determine the number of samples with an alternative count above this threshold
    in the dataset

    :param new_fillout_type: The `new_fillout_type` parameter is a string that represents a new type of
    fillout data that will be added to the summary table. It is used to create new column names in the
    summary table based on the type of fillout data being processed

    :return: The function `__generate_table_and_find_summary_stats` returns a summary table with various
    statistics calculated based on the input DataFrame `df_fillout`, mutation key, fillout type, alt
    threshold, and new fillout type. The summary table includes columns for median VAF, number of
    samples with alt count above the threshold, and number of samples with Total Depth > 0.
    """

    if fillout_type == "NORMAL_":
        summary_table = df_fillout.pivot_table(
            index=mutation_key,
            columns="Tumor_Sample_Barcode",
            values="summary_fragment_standard",
            aggfunc=lambda x: " ".join(x),
        )
        # Find the median VAF for the set
        summary_table[new_fillout_type + "median_VAF"] = df_fillout.groupby(
            mutation_key
        )["t_vaf_fragment_standard"].median()

        # Find the number of samples with alt count above the threshold (alt_thres)
        summary_table[
            new_fillout_type + "n_fillout_sample_alt_detect"
        ] = df_fillout.groupby(mutation_key)["t_alt_count_fragment_standard"].aggregate(
            lambda x: (x >= float(alt_thres)).sum()
        )
        # Find the number of sample with the Total Depth is >0
        # 't_vaf_fragment' column is NA for samples where mutation had no coverage, so count() will exclude it
        summary_table[new_fillout_type + "n_fillout_sample"] = df_fillout.groupby(
            mutation_key
        )["t_vaf_fragment_standard"].count()

        new_columns = {
            col: f"{col}-{new_fillout_type}"
            for col in summary_table.columns
            if not col.startswith(new_fillout_type)
        }
        summary_table.rename(columns=new_columns, inplace=True)
    else:
        summary_table = df_fillout.pivot_table(
            index=mutation_key,
            columns="Tumor_Sample_Barcode",
            values="summary_fragment",
            aggfunc=lambda x: " ".join(x),
        )
        # Find the median VAF for the set
        summary_table[new_fillout_type + "median_VAF"] = df_fillout.groupby(
            mutation_key
        )["t_vaf_fragment"].median()

        # Find the number of samples with alt count above the threshold (alt_thres)

        summary_table[new_fillout_type + "n_fillout_sample_alt_detect"] = (
            df_fillout.groupby(mutation_key)["t_alt_count_fragment"].aggregate(
                lambda x: (x >= float(alt_thres)).sum()
            )
        )
        # Find the number of sample with the Total Depth is >0
        # 't_vaf_fragment' column is NA for samples where mutation had no coverage, so count() will exclude it
        summary_table[new_fillout_type + "n_fillout_sample"] = df_fillout.groupby(
            mutation_key
        )["t_vaf_fragment"].count()

        new_columns = {
            col: f"{col}-{new_fillout_type}"
            for col in summary_table.columns
            if not col.startswith(new_fillout_type)
        }
        summary_table.rename(columns=new_columns, inplace=True)

    return summary_table


def _create_fillout_summary(df_fillout, alt_thres, mutation_key):
    """
    The function `_create_fillout_summary` processes data based on fillout type and generates summary
    statistics for duplex and simplex-duplex categories.

    :param df_fillout: The function `_create_fillout_summary` takes in three parameters: `df_fillout`,
    `alt_thres`, and `mutation_key`

    :param alt_thres: The `alt_thres` parameter in the `_create_fillout_summary` function is used to
    specify a threshold value for the alternative allele count. This threshold is likely used in the
    function to filter or process data based on the alternative allele count being above or below this
    threshold.

    :param mutation_key: The `mutation_key` parameter in the `_create_fillout_summary` function is used
    to specify the key that identifies the mutation in the data frame `df_fillout`. This key is
    essential for generating the summary statistics and tables for the mutation data.

    :return: The function `_create_fillout_summary` returns either a list containing two summary tables
    for duplex and simplex-duplex stats, or a single summary table for standard stats, depending on the
    condition of the `fillout_type`.
    """

    # make sure there is a valid fillout type value and that is suffixed with "_"
    try:
        fillout_type = df_fillout["fillout_type"].iloc[0]
        if fillout_type != "":
            fillout_type = fillout_type + "_"
    except:
        print(
            "The fillout provided to summarize was not run through extract_fillout_type"
        )
        fillout_type = ""
        raise

    # for each of the groups (exception of normals), we are going to subset to their duplex and simplex-duplex stats as well as rename the fillout type those categories
    if fillout_type not in ("NORMAL_"):
        # take the duplex stats and convert their column
        base_columns = [
            col
            for col in df_fillout.columns
            if not col.endswith(("_simplex_duplex", "_duplex", "_simplex"))
        ]
        duplex_columns = [
            col
            for col in df_fillout.columns
            if not col.endswith(("_simplex_duplex", "_simplex"))
        ]

        duplex_df = df_fillout[duplex_columns]
        duplex_df["fillout_type"] = f"{fillout_type}DUPLEX_"
        new_fillout_type_D = f"{fillout_type}DUPLEX_"

        simplex_duplex_columns = base_columns + [
            col for col in df_fillout.columns if col.endswith(("_simplex_duplex"))
        ]
        simplex_duplex_df = df_fillout[simplex_duplex_columns]
        simplex_duplex_df["fillout_type"] = f"{fillout_type}SIMPLEX_DUPLEX_"
        new_fillout_type_SD = f"{fillout_type}SIMPLEX_DUPLEX_"

        # remove the nan t alt count, t ref count and t total count and t vaf fragment columns

        duplex_df = duplex_df.rename(
            columns={
                "summary_fragment_duplex": "summary_fragment",
                "t_vaf_fragment_duplex": "t_vaf_fragment",
                "t_alt_count_fragment_duplex": "t_alt_count_fragment",
                "t_ref_count_fragment_duplex": "t_ref_count_fragment",
                "t_total_count_fragment_duplex": "t_total_count_fragment",
            }
        )

        simplex_duplex_df = simplex_duplex_df.rename(
            columns={
                "summary_fragment_simplex_duplex": "summary_fragment",
                "t_vaf_fragment_simplex_duplex": "t_vaf_fragment",
                "t_alt_count_fragment_simplex_duplex": "t_alt_count_fragment",
                "t_ref_count_fragment_simplex_duplex": "t_ref_count_fragment",
                "t_total_count_fragment_simplex_duplex": "t_total_count_fragment",
            }
        )

        for col in ["t_vaf_fragment", "t_alt_count_fragment", "t_ref_count_fragment"]:
            duplex_df[col] = duplex_df[col].astype(float)
            simplex_duplex_df[col] = simplex_duplex_df[col].astype(float)

        duplex_summary_table = __generate_table_and_find_summary_stats(
            duplex_df, mutation_key, fillout_type, alt_thres, new_fillout_type_D
        )

        simplex_duplex_summary_table = __generate_table_and_find_summary_stats(
            simplex_duplex_df,
            mutation_key,
            fillout_type,
            alt_thres,
            new_fillout_type_SD,
        )

        return [duplex_summary_table, simplex_duplex_summary_table]
    else:

        for col in [
            "t_vaf_fragment_standard",
            "t_alt_count_fragment_standard",
            "t_ref_count_fragment_standard",
        ]:
            df_fillout[col] = df_fillout[col].astype(float)
            df_normal_combined = df_fillout

        summary_table = __generate_table_and_find_summary_stats(
            df_normal_combined, mutation_key, fillout_type, alt_thres, str(fillout_type)
        )

        return summary_table


def _extract_tn_genotypes(
    df_tumor, df_matched_normal, tumor_samplename, normal_samplename
):
    """
    The function `_extract_tn_genotypes` extracts genotype information from tumor and matched normal
    samples and combines them into a single DataFrame.

    :param df_tumor: The function `_extract_tn_genotypes` takes in two dataframes `df_tumor` and
    `df_matched_normal`, as well as two sample names `tumor_samplename` and `normal_samplename`

    :param df_matched_normal: The function `_extract_tn_genotypes` takes in two dataframes `df_tumor`
    and `df_matched_normal`, as well as two sample names `tumor_samplename` and `normal_samplename`. It
    extracts specific columns related to the matched normal sample

    :param tumor_samplename: The `tumor_samplename` parameter is used to specify the sample name of the
    tumor sample for which genotypes are being extracted from the `df_tumor` DataFrame. This sample name
    is used to filter out the relevant rows from the DataFrame based on the tumor sample

    :param normal_samplename: The `normal_samplename` parameter is used to specify the sample name of
    the matched normal sample in the function `_extract_tn_genotypes`. This parameter is used to extract
    the genotype information for the matched normal sample from the `df_matched_normal` DataFrame.

    :return: The function `_extract_tn_genotypes` returns a DataFrame `df_tn_genotype_final` containing
    the genotypes of the tumor sample and, if provided, the matched normal sample.
    """

    df_tn_genotype = df_tumor[
        df_tumor["Tumor_Sample_Barcode"] == str(tumor_samplename)
    ][
        [
            "t_alt_count_fragment_simplex_duplex",
            "t_ref_count_fragment_simplex_duplex",
            "t_vaf_fragment_simplex_duplex",
            "t_alt_count_fragment_duplex",
            "t_ref_count_fragment_duplex",
            "t_vaf_fragment_duplex",
            "t_alt_count_fragment_simplex",
            "t_ref_count_fragment_simplex",
            "t_vaf_fragment_simplex",
        ]
    ]

    if df_tn_genotype.shape[0] == 0:
        raise Exception(
            "Tumor Sample ID {} not found in maf file".format(str(tumor_samplename))
        )
    df_tn_genotype.rename(
        columns={
            "t_alt_count_fragment_simplex_duplex": "SD_t_alt_count_fragment",
            "t_ref_count_fragment_simplex_duplex": "SD_t_ref_count_fragment",
            "t_vaf_fragment_simplex_duplex": "SD_t_vaf_fragment",
            "t_alt_count_fragment_duplex": "D_t_alt_count_fragment",
            "t_ref_count_fragment_duplex": "D_t_ref_count_fragment",
            "t_vaf_fragment_duplex": "D_t_vaf_fragment",
            "t_alt_count_fragment_simplex": "S_t_alt_count_fragment",
            "t_ref_count_fragment_simplex": "S_t_ref_count_fragment",
            "t_vaf_fragment_simplex": "S_t_vaf_fragment",
        },
        inplace=True,
    )

    if str(normal_samplename) != "":

        df_n_genotype = df_matched_normal[
            df_matched_normal["Tumor_Sample_Barcode"] == str(normal_samplename)
        ][
            [
                "t_alt_count_fragment_standard",
                "t_ref_count_fragment_standard",
                "t_vaf_fragment_standard",
            ]
        ]
        df_n_genotype.insert(0, "Matched_Norm_Sample_Barcode", str(normal_samplename))
        df_n_genotype.rename(
            columns={
                "t_alt_count_fragment_standard": "n_alt_count_fragment",
                "t_ref_count_fragment_standard": "n_ref_count_fragment",
                "t_vaf_fragment_standard": "n_vaf_fragment",
            },
            inplace=True,
        )
    df_tn_genotype_final = pd.concat(
        [
            df_tn_genotype.reset_index(drop=True),
            df_n_genotype.reset_index(drop=True),
        ],
        axis=1,
    )
    return df_tn_genotype_final


def make_pre_filtered_maf(
    df_annotation,
    df_all_curated,
    df_all_plasma,
    df_all_tumor,
    df_matched_normal,
    df_all_normals,
    mutation_key,
    tumor_detect_alt_thres,
    curated_detect_alt_thres,
    plasma_detect_alt_thres,
    control_detect_alt_thres,
    tumor_samplename,
    normal_samplename,
):
    """
    The function `make_pre_filtered_maf` generates a pre-filtered MAF file based on various input
    data frames and parameters.

    :param df_annotation: `df_annotation` is a DataFrame containing annotation data

    :param df_all_curated: `df_all_curated` is a DataFrame containing curated samples data.

    :param df_all_plasma: `df_all_plasma` is a DataFrame containing data related to mutations detected
    in plasma samples

    :param df_all_tumor: `df_all_tumor` is a DataFrame containing mutation data for all tumor samples

    :param df_matched_normal: The function `make_pre_filtered_maf` takes several data frames and
    parameters as input to create a pre-filtered MAF (Mutation Annotation Format) data frame. The
    `df_matched_normal` parameter is a dataframe containing information about matched normal samples
    :param df_all_normals: The `df_all_normals` parameter in the `make_pre_filtered_maf` function
    represents a DataFrame containing information about all the normal samples in the dataset. This
    DataFrame is used to generate a summary table for normal samples based on certain criteria like
    mutation detection thresholds and fillout types.

    :param mutation_key: The `mutation_key` parameter is a set of columns used as a key used to pivot the data table.

    :param tumor_detect_alt_thres: The `tumor_detect_alt_thres` parameter is used to specify the
    threshold for detecting alternate alleles in the tumor samples during the mutation analysis process.
    This threshold helps determine whether a mutation is present in the tumor sample based on the level
    of alternate alleles detected

    :param curated_detect_alt_thres: The parameter `curated_detect_alt_thres` is used to specify the
    detection threshold for alternate alleles in curated data. This threshold is likely used to filter
    out mutations based on their allele frequencies in the curated dataset.

    :param plasma_detect_alt_thres: The parameter `plasma_detect_alt_thres` is used to set the threshold
    for detecting alternate alleles in plasma samples. This threshold value determines what level of
    alternate allele frequency is considered significant for identifying mutations in the plasma
    samples.

    :param control_detect_alt_thres: The `control_detect_alt_thres` parameter likely represents the
    detection threshold for alternate alleles in the control samples. This threshold is used to
    determine if a mutation is detected in the control samples based on the alternate allele frequency.

    :param tumor_samplename: The `tumor_samplename` parameter is used to specify the name of the tumor
    sample in the dataset. It is a string that identifies the tumor sample within the dataframes
    provided as input to the `make_pre_filtered_maf` function.

    :param normal_samplename: The `normal_samplename` parameter in the `make_pre_filtered_maf` function
    is used to specify the name of the normal sample in the dataset. This parameter is important for
    identifying and processing the normal sample data within the function.

    :return: The function `make_pre_filtered_maf` returns a DataFrame `df_pre_filter` which is the
    result of merging various summary tables and genotypes extracted from input DataFrames based on
    specified conditions and thresholds.
    """

    [tumor_duplex_summary_table, tumor_simplex_duplex_summary_table] = (
        _create_fillout_summary(df_all_tumor, tumor_detect_alt_thres, mutation_key)
    )

    if df_all_normals.empty:
        df_normal_summary = pd.DataFrame(index=tumor_summary_table.index.copy())
        df_normal_summary["NORMAL_median_VAF"] = "no_normals_in_pool"
        df_normal_summary["NORMAL_n_fillout_sample_alt_detect"] = "no_normals_in_pool"
        df_normal_summary["NORMAL_n_fillout_sample"] = "no_normals_in_pool"
    else:
        df_all_normals["fillout_type"] = df_all_normals["fillout_type"].replace(
            ["UNMATCHED_NORMAL", "MATCHED_NORMAL"], "NORMAL"
        )
        normal_summary_table = _create_fillout_summary(
            df_all_normals, tumor_detect_alt_thres, mutation_key
        )

    [curated_duplex_summary_table, curated_simplex_duplex_summary_table] = (
        _create_fillout_summary(df_all_curated, curated_detect_alt_thres, mutation_key)
    )
    [plasma_duplex_summary_table, plasma_simplex_duplex_summary_table] = (
        _create_fillout_summary(df_all_plasma, curated_detect_alt_thres, mutation_key)
    )

    df_tn_genotype_final = _extract_tn_genotypes(
        df_all_tumor, df_all_normals, tumor_samplename, normal_samplename
    )

    df_anno_with_genotypes = pd.concat(
        [
            df_annotation.reset_index(drop=True),
            df_tn_genotype_final.reset_index(drop=True),
        ],
        axis=1,
    )

    df_anno_with_genotypes.index = df_annotation.index

    df_pre_filter = (
        df_anno_with_genotypes.merge(
            normal_summary_table, left_index=True, right_index=True
        )
        .merge(curated_duplex_summary_table, left_index=True, right_index=True)
        .merge(curated_simplex_duplex_summary_table, left_index=True, right_index=True)
        .merge(plasma_duplex_summary_table, left_index=True, right_index=True)
        .merge(plasma_simplex_duplex_summary_table, left_index=True, right_index=True)
    )

    return df_pre_filter


def apply_filter_maf(pre_filter_maf, **kwargs):
    """
    The `apply_filter_maf` function processes mutation data by applying various filters and updating the
    mutation status based on specified criteria.

    :param pre_filter_maf: The `pre_filter_maf` parameter in the `apply_filter_maf` function seems to
    represent a DataFrame containing mutation data. This DataFrame is used within the function to
    iterate over each mutation and update its status based on certain conditions specified in the
    tagging functions

    :return: The `apply_filter_maf` function returns a modified DataFrame `df_post_filter_final` after
    applying various tagging functions


    """

    # mini tagging functions (will need to be moved into helper.pyx)
    def tag_germline(mut, status, inner_kwargs):
        """
        The function `tag_germline` checks for germline mutations based on specified criteria and
        updates the mutation status accordingly.

        :param mut: The `mut` parameter seems to be a DataFrame containing mutation data. The function
        `tag_germline` checks certain conditions within this DataFrame and updates the `status` variable
        accordingly based on the conditions met.

        :param status: The `status` parameter is a string that represents the current status of a
        mutation. It is updated based on certain conditions in the `tag_germline` function.

        :param inner_kwargs: The `inner_kwargs` parameter is a dictionary containing inner keyword
        arguments that are used within the `tag_germline` function.

        :return: The function `tag_germline` is returning the updated `status` variable after checking
        certain conditions. If the condition for a matched normal with sufficient coverage is met, it
        appends "Germline;" to the status. If the conditions for a common variant with sufficient tumor
        coverage and VAF are met, it appends "LikelyGermline;" to the status.
        """
        # if there is a matched normal and it has sufficient coverage
        if "n_vaf_fragment" in mut.index.tolist() and mut["n_ref_count_fragment"] + mut[
            "n_alt_count_fragment"
        ] > float(inner_kwargs["normal_TD_min"]):
            status = status + "Germline;"
            return status
        elif (
            "common_variant" in mut["FILTER"]
            and (mut["SD_t_ref_count_fragment"] + mut["SD_t_alt_count_fragment"])
            > float(inner_kwargs["tumor_TD_min"])
            and mut["SD_t_vaf_fragment"]
            > float(inner_kwargs["tumor_vaf_germline_thres"])
        ):
            status = status + "LikelyGermline;"
            return status

    def tag_below_alt_threshold(mut, status, inner_kwargs):
        """
        The function `tag_below_alt_threshold` checks mutation alt counts against specified thresholds
        and updates the status accordingly.

        :param mut: The `mut` parameter seems to be a dictionary containing mutation-related
        information. The function `tag_below_alt_threshold` is checking certain conditions based on
        values in the `mut` dictionary and the `inner_kwargs` dictionary.

        :param status: The `status` parameter is a string variable that likely stores information about
        the status of a mutation based on certain conditions being met in the `tag_below_alt_threshold`
        function.

        :param inner_kwargs: The `inner_kwargs` parameter seems to be a dictionary containing the
        following keys and their corresponding values:

        :return: The function `tag_below_alt_threshold` is returning the updated `status` variable based
        on the conditions specified in the function. The `status` variable is being modified by
        appending specific strings based on the comparisons made in the function.
        """
        if mut["SD_t_alt_count_fragment"] < float(inner_kwargs["tier_one_alt_min"]) or (
            mut["hotspot_whitelist"] == False
            and mut["SD_t_alt_count_fragment"] < float(inner_kwargs["tier_two_alt_min"])
        ):
            if mut["caller_t_alt_count"] >= float(inner_kwargs["tier_two_alt_min"]) or (
                mut["hotspot_whitelist"] == True
                and mut["caller_t_alt_count"] >= float(inner_kwargs["tier_one_alt_min"])
            ):
                status = status + "BelowAltThreshold;LostbyGenotyper;"
            else:
                status = status + "BelowAltThreshold;"
        return status

    def occurrence_in_curated(mut, status, inner_kwargs):
        """
        The function checks if a specified mutation occurrence meets a minimum threshold in curated
        samples and updates the status accordingly.

        :param mut: The `mut` parameter seems to be a dictionary containing information related to a
        mutation. The function `occurrence_in_curated` checks if a specific key
        `"CURATED_DUPLEX_n_fillout_sample_alt_detect"` in the `mut` dictionary is greater than or equal
        to a minimum value specified

        :param status: The `status` parameter is a variable that holds the current status of an item or
        process. In the provided function `occurrence_in_curated`, the `status` variable is being
        updated based on a condition related to the `mut` dictionary and the `inner_kwargs` dictionary.

        :param inner_kwargs: The `inner_kwargs` parameter is a dictionary containing key-value pairs of
        additional arguments or parameters that are passed to the function `occurrence_in_curated`.

        :return: the updated `status` variable, which may have the string "InCurated;" appended to it if
        the condition `mut["CURATED_DUPLEX_n_fillout_sample_alt_detect"] >=
        float(inner_kwargs["min_n_curated_samples_alt_detected"])` is met.
        """
        if mut["CURATED_DUPLEX_n_fillout_sample_alt_detect"] >= float(
            inner_kwargs["min_n_curated_samples_alt_detected"]
        ):
            status = status + "InCurated;"
        return status

    def occurrence_in_normal(mut, status, inner_kwargs):
        """
        The function `occurrence_in_normal` checks for certain conditions related to coverage and VAF
        ratios in mutation data and updates the status string accordingly.

        :param mut: `mut` is a dictionary containing various mutation-related information such as
        counts, VAF (Variant Allele Frequency), and other metrics for a specific mutation
        :param status: The `status` parameter is a string that keeps track of the status of a mutation
        based on certain conditions being met in the `occurrence_in_normal` function. It is updated with
        specific status codes like "TNRatio-curatedmedian" or "TNRatio-matchnorm" based on the
        conditions evaluated

        :param inner_kwargs: `inner_kwargs` is a dictionary containing various parameters used in the
        function `occurrence_in_normal`. The function uses values from `inner_kwargs` to make decisions
        based on certain thresholds and ratios.

        :return: The function `occurrence_in_normal` returns the updated `status` variable after
        checking certain conditions based on the input mutation data (`mut`) and inner keyword arguments
        (`inner_kwargs`). The function modifies the `status` string by appending specific strings based
        on the conditions met within the function. The final updated `status` variable is returned as
        the output of the function.
        """
        # if normal and tumor coverage is greater than the minimal
        if mut["SD_t_ref_count_fragment"] + mut["SD_t_alt_count_fragment"] > float(
            inner_kwargs["tumor_TD_min"]
        ):
            if mut["CURATED_DUPLEX_median_VAF"] != 0:
                if mut["SD_t_vaf_fragment"] / mut["CURATED_DUPLEX_median_VAF"] < float(
                    inner_kwargs["tn_ratio_thres"]
                ):
                    status = status + "TNRatio-curatedmedian;"

            if "n_vaf_fragment" in mut.index.tolist():
                if (
                    mut["n_ref_count_fragment"] + mut["n_alt_count_fragment"]
                    > float(inner_kwargs["normal_TD_min"])
                    and mut["n_vaf_fragment"] != 0
                ):
                    if mut["SD_t_vaf_fragment"] / mut["n_vaf_fragment"] < float(
                        inner_kwargs["tn_ratio_thres"]
                    ):
                        status = status + "TNRatio-matchnorm;"
        return status

    def in_blocklist(mut, status, inner_kwargs):
        """
        The function `in_blocklist` checks if a mutation is listed in a blocklist and updates the status
        accordingly.

        :param mut: The `mut` parameter is a dictionary containing information about a mutation. It
        likely includes keys such as "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
        and "Tumor_Seq_Allele2" with corresponding values for each key

        :param status: The `status` parameter is a variable that holds the current status of a mutation.
        It is updated based on whether the mutation is listed in a blocklist

        :param inner_kwargs: `inner_kwargs` is a dictionary containing the key "blocklist_lst", which is
        a list of strings representing mutations that are in a blocklist. The function `in_blocklist`
        takes three parameters: `mut` (a dictionary representing a mutation), `status` (a string
        representing the current

        :return: the updated `status` variable, which may have the string "InBlacklist;" appended to it
        if the mutation specified in the input `mut` is found in the blocklist specified in
        `inner_kwargs["blocklist_lst"]`.
        """
        # if mutation is listed in blocklist
        if (
            str(mut["Chromosome"])
            + "_"
            + str(mut["Start_Position"])
            + "_"
            + str(mut["End_Position"])
            + "_"
            + str(mut["Reference_Allele"])
            + "_"
            + str(mut["Tumor_Seq_Allele2"])
            in inner_kwargs["blocklist_lst"]
        ):
            status = status + "InBlacklist;"
        return status

    def cleanup_post_filter(df_post_filter):
        """
        The function `cleanup_post_filter` reorganizes columns in a DataFrame and adds new columns if
        certain columns are missing.

        :param df_post_filter: The `cleanup_post_filter` function is designed to clean up and modify a
        DataFrame `df_post_filter` by moving the "Status" column next to "Hotspots" and adding
        additional columns related to matched normal samples if they do not already exist

        :return: The function `cleanup_post_filter` is returning the DataFrame `df_post_filter` after
        performing the specified operations, which include moving the "Status" column next to
        "Hotspots", adding columns for "Matched_Norm_Sample_Barcode", "n_alt_count_fragment",
        "n_ref_count_fragment", and "n_vaf_fragment" if they do not already exist, and inserting a
        "Match"
        """

        # Move Status column next to Hotspots
        col = list(df_post_filter)
        col.insert(col.index("SD_t_alt_count_fragment"), col.pop(col.index("Status")))
        df_post_filter = df_post_filter[col]
        # Add Match Normal columns even when sample is unmatched
        if "Matched_Norm_Sample_Barcode" not in col:
            df_post_filter.insert(
                col.index("SD_t_vaf_fragment") + 1,
                "Matched_Norm_Sample_Barcode",
                "Unmatched",
            )
            df_post_filter.insert(
                col.index("SD_t_vaf_fragment") + 2, "n_alt_count_fragment", "NA"
            )
            df_post_filter.insert(
                col.index("SD_t_vaf_fragment") + 3, "n_ref_count_fragment", "NA"
            )
            df_post_filter.insert(
                col.index("SD_t_vaf_fragment") + 4, "n_vaf_fragment", "NA"
            )

        col = list(df_post_filter)
        df_post_filter.insert(
            col.index("Matched_Norm_Sample_Barcode") + 1, "Matched_Norm_Bamfile", "NA"
        )
        return df_post_filter

    df_post_filter = pre_filter_maf.copy()

    df_post_filter["Status"] = ""

    for i, mut in df_post_filter.iterrows():
        status = ""
        status = tag_germline(mut, status, kwargs)
        if status is None:
            status = ""
        status = tag_below_alt_threshold(mut, status, kwargs)
        if status is None:
            status = ""
        status = occurrence_in_curated(mut, status, kwargs)
        if status is None:
            status = ""
        status = occurrence_in_normal(mut, status, kwargs)
        if status is None:
            status = ""
        status = in_blocklist(mut, status, kwargs)
        df_post_filter.loc[i, "Status"] = status

    df_post_filter_final = cleanup_post_filter(df_post_filter)

    return df_post_filter_final


def make_condensed_post_filter(df_post_filter):

    # creating the "condensed" MAF -- can be customized in the future
    df_condensed = df_post_filter.loc[:, :"n_vaf_fragment"]
    return df_condensed


def calculate_stats_for_mpath(df_post_filter_format):

    # compute various mutation depth and vaf metrics from fragment columns calculated in access_filters
    df_post_filter_format["D_t_count_fragment"] = (
        df_post_filter_format["D_t_ref_count_fragment"]
        + df_post_filter_format["D_t_alt_count_fragment"]
    )

    df_post_filter_format["SD_t_count_fragment"] = (
        df_post_filter_format["SD_t_ref_count_fragment"]
        + df_post_filter_format["SD_t_alt_count_fragment"]
    )

    df_post_filter_format["S_t_ref_count_fragment"] = (
        df_post_filter_format["SD_t_ref_count_fragment"]
        - df_post_filter_format["D_t_ref_count_fragment"]
    )

    df_post_filter_format["S_t_alt_count_fragment"] = (
        df_post_filter_format["SD_t_alt_count_fragment"]
        - df_post_filter_format["D_t_alt_count_fragment"]
    )

    df_post_filter_format["S_t_count_fragment"] = (
        df_post_filter_format["S_t_ref_count_fragment"]
        + df_post_filter_format["S_t_alt_count_fragment"]
    )

    df_post_filter_format["n_count_fragment"] = (
        df_post_filter_format["n_ref_count_fragment"]
        + df_post_filter_format["n_alt_count_fragment"]
    )

    df_post_filter_format["S_t_vaf_fragment"] = (
        df_post_filter_format["S_t_alt_count_fragment"]
        / df_post_filter_format["S_t_count_fragment"]
    ).fillna(0)

    df_post_filter_format["SD_t_vaf_fragment_over_n_vaf_fragment"] = (
        df_post_filter_format["SD_t_vaf_fragment"]
        / df_post_filter_format["n_vaf_fragment"]
    ).fillna(0)

    # convert NaN and inf computed values to 0
    computed_maf = df_post_filter_format.replace([np.inf, np.nan], 0)

    # format SNP column
    computed_maf["dbSNP_RS"] = computed_maf["dbSNP_RS"].apply(
        lambda x: x if isinstance(x, str) and x.startswith("rs") else ""
    )

    # generate occurrence stats columns
    # computed_maf["CURATED_DUPLEX_n_fillout_sample"] = (
    #         computed_maf["CURATED_DUPLEX_n_fillout_sample_alt_detect"].map(str)
    #         + ";"
    #         + computed_maf["CURATED_DUPLEX_median_VAF"].map(str)
    #     )
    # computed_maf["CURATED_SIMPLEX_DUPLEX_n_fillout_sample"] = (
    #         computed_maf["CURATED_SIMPLEX_DUPLEX_n_fillout_sample_alt_detect"].map(str)
    #         + ";"
    #         + computed_maf["CURATED_SIMPLEX_DUPLEX_median_VAF"].map(str)
    #     )
    # computed_maf["NORMAL_n_fillout_sample"] = (
    #         computed_maf["NORMAL_n_fillout_sample_alt_detect"].map(str)
    #         + ";"
    #         + computed_maf["NORMAL_median_VAF"].map(str)
    #     )

    return computed_maf


def add_dummy_columns(maf, columns):
    """
    Temporary function to add dummy columns
    to meet DMP requirements
    """
    for col in columns:
        if not col in maf.columns:
            maf[col] = ""
    return maf


def customize_cosmic(cosmic_id, cosmic_occurrence):
    """
    helper function to customize cosmic annotation.
    If cosmic_id is defined, but not occurrence, then
    a generic value of "1(unknown)" will be used.
    """
    if cosmic_id is not np.nan and cosmic_id != "":
        # OCCURENCE spelled incorrectly by design
        return (
            "ID="
            + cosmic_id
            + ";OCCURENCE="
            + (cosmic_occurrence if cosmic_occurrence is not np.nan else "1(unknown)")
        )
    else:
        return ""


def get_exon(maf_exon, maf_intron):
    """ "
    helper function to determine the exonic or
    intronic location of a variant.
    """
    try:
        exon, total_exon = str(maf_exon).split("/")
        return "exon" + str(exon)
    except ValueError:  # Not exonic
        try:
            intron, total_intron = str(maf_intron).split("/")
            return "intron" + str(intron)
        except ValueError:  # Not intronic
            return ""


def format_maf_for_mpath(df_post_filter):
    incompatible_required_column_headers = filter(
        lambda x: "-" in x and any([x.startswith("CURATED-"), x.startswith("NORMAL-")]),
        df_post_filter.columns,
    )

    renamed_maf = df_post_filter.rename(
        columns=dict(
            zip(
                incompatible_required_column_headers,
                map(
                    lambda x: x.replace("-", "_"), incompatible_required_column_headers
                ),
            )
        )
    )

    # assign potential missing expected columns in MAF
    renamed_maf_w_dummy_cols = add_dummy_columns(renamed_maf, MAF_DUMMY_COLUMNS2)

    # if a mutation does not have a flag for "Mutation_Status", classify it as Novel
    renamed_maf_w_dummy_cols["Mutation_Class"] = np.vectorize(
        lambda x: "Novel" if x is np.nan else ""
    )(renamed_maf_w_dummy_cols["Status"])

    # add modified cosmic column
    renamed_maf_w_dummy_cols["Cosmic_ID"] = np.vectorize(
        customize_cosmic, otypes=[str]
    )(
        renamed_maf_w_dummy_cols["cosmic_ID"],
        renamed_maf_w_dummy_cols["cosmic_OCCURENCE"],
    )

    # commenting this out for now since it is removing individual sample summary statistics
    # try:
    #     renamed_maf_w_dummy_cols = renamed_maf_w_dummy_cols[MAF_COLUMNS_SELECT]
    # except KeyError:
    #     missing_columns = set(MAF_COLUMNS_SELECT) - set(renamed_maf_w_dummy_cols.columns.values.tolist())
    #     missing_columns_str = ",".join(missing_columns)
    #     # TODO: add proper error handling message here for missing required columns
    #     sys.exit(f"The following required columns are missing: {missing_columns_str}")

    # compute exon and intron columns
    renamed_maf_w_dummy_cols["EXON"] = np.vectorize(get_exon, otypes=[str])(
        renamed_maf_w_dummy_cols["EXON"], renamed_maf_w_dummy_cols["INTRON"]
    )
    renamed_maf_w_dummy_cols = renamed_maf_w_dummy_cols.drop(["INTRON"], axis=1)

    # computing all the gnomad associated cols
    # "TypeError: '>=' not supported between instances of 'float' and 'str'"
    renamed_maf_w_dummy_cols[GNOMAD_COLUMNS] = renamed_maf_w_dummy_cols[
        GNOMAD_COLUMNS
    ].replace("", np.nan)

    # get max of gnomad
    renamed_maf_w_dummy_cols["gnomAD_Max_AF"] = np.max(
        renamed_maf_w_dummy_cols[GNOMAD_COLUMNS].values, axis=1
    )

    # generate occurrence stats columns
    calculated_and_formatted_maf = calculate_stats_for_mpath(renamed_maf_w_dummy_cols)

    return calculated_and_formatted_maf


if __name__ == "__main__":
    app()
