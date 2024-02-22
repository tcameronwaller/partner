"""
Supply functionality for extraction and organization of information from
reports of other tools.

This module is not directly executable.

This module within subpackage 'partner' provides executable functionality under
the management of a higher level package. Importation paths must represent this
hierarchy.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Monroe, North Carolina 28110
    United States of America

License:

    This file is part of Partner
    (https://github.com/tcameronwaller/partner/).

    Partner supports data analysis in multiple other projects.
    Copyright (C) 2024 Thomas Cameron Waller

    Partner is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    Partner is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with Partner. If not, see <http://www.gnu.org/licenses/>.
"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard
import os
import csv
import copy
import textwrap
import string
import gzip
import shutil
import textwrap
import itertools
import math

# Relevant

import pandas
import sklearn
import sklearn.preprocessing
import scipy
import numpy
import statsmodels.api

# Custom
import partner.utility as putility # this import path for subpackage

#dir()
#importlib.reload()

###############################################################################
# Functionality


# TODO: TCW; 12 December 2022; 7 August 2023
# Use the new utility function to calculate and organize 95% and 99% confidence intervals and ranges.

##########
# Extraction and collection of relevant information from the logs from
# Linkage Disequilibrium Score Regression (LDSC)


def read_extract_ldsc_heritability(
    path_file=None,
):
    """
    Reads and extracts information from log file of LDSC for estimation of SNP
    heritability from GWAS summary statistics.

    https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
    https://github.com/bulik/ldsc/blob/master/munge_sumstats.py
    https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format

    https://www.mathsisfun.com/data/confidence-interval.html

    arguments:
        path_file (str): full path to file that contains information from a
            SNP heritability analysis in LDSC

    raises:

    returns:
        (dict): information about LDSC analysis

    """

    # Initialize variables for extraction.
    variants = float("nan")
    heritability = float("nan")
    heritability_error = float("nan")
    heritability_ci95_low = float("nan")
    heritability_ci95_high = float("nan")
    heritability_ci99_low = float("nan")
    heritability_ci99_high = float("nan")
    chi_square = float("nan")
    ratio = float("nan")
    ratio_error = float("nan")

    # Read relevant lines of character strings from file.
    lines = putility.read_file_text_lines(
        path_file=path_file,
        start=22,
        stop=30,
    )

    # Define character strings that indicate relevant information.
    prefix_variants = "After merging with regression SNP LD, "
    suffix_variants = " SNPs remain."
    #prefix_heritability = "Total Observed scale h2: "
    prefix_heritability_observed = "Total Observed scale h2: "
    prefix_heritability_liability = "Total Liability scale h2: "
    prefix_chi_square = "Mean Chi^2: "
    prefix_ratio = "Ratio: "

    # Extract information from lines.
    for line in lines:
        if prefix_variants in line:
            variants = int(
                line.replace(prefix_variants, "").replace(suffix_variants, "")
            )
        elif prefix_heritability_observed in line:
            content = line.replace(prefix_heritability_observed, "")
            contents = content.split(" (")
            heritability_test = contents[0]
            if (not "NA" in heritability_test):
                heritability = float(contents[0])
                heritability_error = float(contents[1].replace(")", ""))
            pass
        elif prefix_heritability_liability in line:
            content = line.replace(prefix_heritability_liability, "")
            contents = content.split(" (")
            heritability_test = contents[0]
            if (not "NA" in heritability_test):
                heritability = float(contents[0])
                heritability_error = float(contents[1].replace(")", ""))
            pass
        elif (
            (not math.isnan(heritability)) and
            (prefix_chi_square in line)
        ):
            chi_square = float(line.replace(prefix_chi_square, ""))
            pass
        elif (
            (not math.isnan(heritability)) and
            (prefix_ratio in line)
        ):
            content = line.replace(prefix_ratio, "")
            contents = content.split(" (")
            ratio_test = contents[0]
            if (not "NA" in ratio_test):
                ratio = float(contents[0])
                ratio_error = float(
                    contents[1].replace(")", "")
                )
            pass
        pass

    # Organize information.
    if (
        (not pandas.isna(heritability)) and
        (not pandas.isna(heritability_error))
    ):
        heritability_ci95_low = (heritability - (1.960 * heritability_error))
        heritability_ci95_high = (heritability + (1.960 * heritability_error))
        heritability_ci99_low = (heritability - (2.576 * heritability_error))
        heritability_ci99_high = (heritability + (2.576 * heritability_error))
        pass
    heritability_ci95 = str(
        str(round(heritability_ci95_low, 4)) + " ... " +
        str(round(heritability_ci95_high, 4))
    )
    heritability_ci99 = str(
        str(round(heritability_ci99_low, 4)) + " ... " +
        str(round(heritability_ci99_high, 4))
    )
    summary_heritability_error = str(
        str(heritability) + " (" + str(round(heritability_error, 4)) + ")"
    )
    summary_heritability_ci95 = str(
        "(h2: " + str(heritability) + "; 95% CI: " + str(heritability_ci95) + ")"
    )
    summary_heritability_ci99 = str(
        "(h2: " + str(heritability) + "; 99% CI: " + str(heritability_ci99) + ")"
    )

    # Collect information.
    record = dict()
    record["summary_heritability_error"] = summary_heritability_error
    record["summary_heritability_ci95"] = summary_heritability_ci95
    record["summary_heritability_ci99"] = summary_heritability_ci99
    record["variants"] = variants
    record["chi_square"] = chi_square
    record["ratio"] = ratio
    record["ratio_error"] = ratio_error
    record["heritability"] = heritability
    record["heritability_error"] = heritability_error
    record["heritability_ci95_low"] = heritability_ci95_low
    record["heritability_ci95_high"] = heritability_ci95_high
    record["heritability_ci99_low"] = heritability_ci99_low
    record["heritability_ci99_high"] = heritability_ci99_high
    # Define list of variables in order.
    variables = [
        "summary_heritability_error",
        "summary_heritability_ci95",
        "summary_heritability_ci99",
        "variants",
        "chi_square",
        "ratio",
        "ratio_error",
        "heritability",
        "heritability_error",
        "heritability_ci95_low",
        "heritability_ci95_high",
        "heritability_ci99_low",
        "heritability_ci99_high",
    ]
    # Return information.
    pail = dict()
    pail["record"] = record
    pail["variables"] = variables
    return pail


def read_extract_ldsc_correlation(
    path_file=None,
):
    """
    Reads and extracts information from log file of LDSC for estimation of
    genetic correlation between two sets of GWAS summary statistics.

    https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
    https://github.com/bulik/ldsc/blob/master/munge_sumstats.py
    https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format

    https://www.mathsisfun.com/data/confidence-interval.html

    arguments:
        path_file (str): full path to file that contains information from a
            genetic correlation analysis in LDSC

    raises:

    returns:
        (dict): information about LDSC analysis

    """

    # Determine whether report has a summary table.
    # Keep index of prefix title.
    # Read relevant lines from file.
    #lines = putility.read_file_text_lines(
    #    path_file=path_file,
    #    start=50,
    #    stop=250,
    #)
    #index = 50
    #index_title = float("nan")
    #count_lines = len(lines)
    #for line in lines:
    #    if "Summary of Genetic Correlation Results" in line:
    #        summary_table = True
    #        index_title = index
    #    index += 1
    #if summary_table:
    #    # Read values from bottom table in report.
    #    table = pandas.read_csv(
    #        path_file,
    #        sep="\s+",
    #        header=(index_title + 2),
    #        #skip_blank_lines=True,
    #    )

    # Initialize variables for extraction.
    variants = float("nan")
    variants_valid = float("nan")
    covariance = float("nan")
    covariance_error = float("nan")
    correlation = float("nan")
    correlation_error = float("nan")
    correlation_ci95_low = float("nan")
    correlation_ci95_high = float("nan")
    correlation_ci99_low = float("nan")
    correlation_ci99_high = float("nan")
    correlation_ci95_not_zero = float("nan")
    correlation_ci95_not_one = float("nan")
    correlation_absolute = float("nan")
    z_score = float("nan")
    p_not_zero = float("nan")

    # Read relevant lines of character strings from file.
    lines = putility.read_file_text_lines(
        path_file=path_file,
        start=25,
        stop=57,
    )

    # Define character strings that indicate relevant information.
    prefix_variants = "After merging with summary statistics, "
    suffix_variants = " SNPs remain."
    prefix_variants_valid = ""
    suffix_variants_valid = " SNPs with valid alleles."
    prefix_covariance = "Total Observed scale gencov: "
    prefix_correlation = "Genetic Correlation: "
    prefix_z_score = "Z-score: "
    prefix_p = "P: "

    # Extract information from lines.
    for line in lines:
        if prefix_variants in line:
            variants = int(
                line.replace(prefix_variants, "").replace(suffix_variants, "")
            )
        elif suffix_variants_valid in line:
            variants_valid = int(
                line.replace(suffix_variants_valid, "")
            )
        elif prefix_correlation in line:
            content = line.replace(prefix_correlation, "")
            contents = content.split(" (")
            correlation_test = contents[0]
            if (not "nan" in correlation_test):
                correlation = float(contents[0])
                correlation_error = float(contents[1].replace(")", ""))
            pass
        elif (
            (not math.isnan(correlation)) and
            (prefix_covariance in line)
        ):
            content = line.replace(prefix_covariance, "")
            contents = content.split(" (")
            covariance_test = contents[0]
            if (not "NA" in covariance_test):
                covariance = float(contents[0])
                covariance_error = float(
                    contents[1].replace(")", "")
                )
            pass
        elif (
            (not math.isnan(correlation)) and
            (prefix_z_score in line)
        ):
            z_score = float(line.replace(prefix_z_score, ""))
            pass
        elif (
            (not math.isnan(correlation)) and
            (prefix_p in line)
        ):
            p_not_zero = float(line.replace(prefix_p, ""))
            pass
        pass

    # Organize information.
    if (
        (not pandas.isna(correlation)) and
        (not pandas.isna(correlation_error))
    ):
        correlation_absolute = math.fabs(correlation)
        correlation_ci95_low = (correlation - (1.960 * correlation_error))
        correlation_ci95_high = (correlation + (1.960 * correlation_error))
        correlation_ci99_low = (correlation - (2.576 * correlation_error))
        correlation_ci99_high = (correlation + (2.576 * correlation_error))
        # Determine whether confidence interval crosses zero or one.
        if (
            ((correlation_ci95_low > 0) and (correlation_ci95_high > 0)) or
            ((correlation_ci95_low < 0) and (correlation_ci95_high < 0))
        ):
            correlation_ci95_not_zero = 1
        else:
            correlation_ci95_not_zero = 0
            pass
        if ((correlation_ci95_low < 1.0) and (correlation_ci95_high < 1.0)):
            correlation_ci95_not_one = 1
        else:
            correlation_ci95_not_one = 0
            pass
        pass
    correlation_ci95 = str(
        str(round(correlation_ci95_low, 4)) + " ... " +
        str(round(correlation_ci95_high, 4))
    )
    correlation_ci99 = str(
        str(round(correlation_ci99_low, 4)) + " ... " +
        str(round(correlation_ci99_high, 4))
    )
    summary_correlation_error = str(
        str(correlation) + " (" + str(round(correlation_error, 4)) + ")"
    )
    summary_correlation_ci95 = str(
        "(rg: " + str(correlation) + "; 95% CI: " + str(correlation_ci95) + ")"
    )
    summary_correlation_ci99 = str(
        "(rg: " + str(correlation) + "; 99% CI: " + str(correlation_ci99) + ")"
    )

    # Collect information.
    record = dict()
    record["summary_correlation_error"] = summary_correlation_error
    record["summary_correlation_ci95"] = summary_correlation_ci95
    record["summary_correlation_ci99"] = summary_correlation_ci99
    record["variants"] = variants
    record["variants_valid"] = variants_valid
    record["covariance"] = correlation
    record["covariance_error"] = correlation_error
    record["z_score"] = z_score
    record["p_not_zero"] = p_not_zero
    record["correlation_ci95_not_zero"] = correlation_ci95_not_zero
    record["correlation_ci95_not_one"] = correlation_ci95_not_one
    record["correlation_absolute"] = correlation_absolute
    record["correlation"] = correlation
    record["correlation_error"] = correlation_error
    record["correlation_ci95_low"] = correlation_ci95_low
    record["correlation_ci95_high"] = correlation_ci95_high
    record["correlation_ci99_low"] = correlation_ci99_low
    record["correlation_ci99_high"] = correlation_ci99_high
    # Define list of variables in order.
    variables = [
        "summary_correlation_error",
        "summary_correlation_ci95",
        "summary_correlation_ci99",
        "variants",
        "variants_valid",
        "covariance",
        "covariance_error",
        "z_score",
        "p_not_zero",
        "correlation_ci95_not_zero",
        "correlation_ci95_not_one",
        "correlation_absolute",
        "correlation",
        "correlation_error",
        "correlation_ci95_low",
        "correlation_ci95_high",
        "correlation_ci99_low",
        "correlation_ci99_high",
    ]
    variables_convenience = [
        "correlation",
        "p_not_zero",
        "summary_correlation_error",
        "summary_correlation_ci95",
        "summary_correlation_ci99",
        "variants",
        "variants_valid",
        "covariance",
        "covariance_error",
        "z_score",
        "correlation_ci95_not_zero",
        "correlation_ci95_not_one",
        "correlation_absolute",
        "correlation_error",
        "correlation_ci95_low",
        "correlation_ci95_high",
        "correlation_ci99_low",
        "correlation_ci99_high",
    ]
    # Return information.
    pail = dict()
    pail["record"] = record
    pail["variables"] = variables_convenience
    return pail


def read_extract_from_all_ldsc_files_in_directory(
    path_directory=None,
    file_name_pattern=None,
    file_name_pattern_not=None,
    analysis=None,
    report=None,
):
    """
    Reads and extracts information about LDSC SNP heritability or genetic
        correlation analyses from all files within a parent directory.

        Linkage Disequilibrium Score Regression (LDSC)

    arguments:
        path_directory (str): full path to parent directory that contains
            relevant files
        file_name_pattern (str): character string in names of relevant files
        file_name_pattern_not (str): character string not in names of relevant
            files
        analysis (str): type of analysis in LDSC, either 'heritability' or
            'correlation'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table with variables across columns and
            samples or instances across rows

    """

    # Read names of relevant files within parent directory.
    names_files = putility.extract_directory_file_names_filter_by_name(
        path=path_directory,
        name=file_name_pattern,
        name_not=file_name_pattern_not,
    )

    # Collect information from analysis on each study.
    records = list()
    # Iterate on files for each study.
    counter = 0
    for name_file in names_files:
        # Define full path to file.
        path_file = os.path.join(path_directory, name_file,)
        # Extract information from LDSC analysis.
        if (str(analysis).strip() == "heritability"):
            pail = read_extract_ldsc_heritability(
                path_file=path_file,
            )
        elif (str(analysis).strip() == "correlation"):
            pail = read_extract_ldsc_correlation(
                path_file=path_file,
            )
            pass
        # Collect and organize information about study.
        # Extraction functions return the relevant variables for each type of
        # analysis.
        pail["record"]["path_directory"] = str(path_directory)
        pail["record"]["name_file"] = str(name_file)
        records.append(pail["record"])
        if (counter == 0):
            variables = pail["variables"]
        # Increment counter.
        counter += 1
        pass
    # Organize table.
    table = putility.convert_records_to_dataframe(
        records=records
    )
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # columns.insert(0, dependence)
    columns = ["path_directory", "name_file",]
    columns.extend(variables)
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    table = table[[*columns]]
    table.sort_values(
        by=["path_directory", "name_file",],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Report.
    if report:
        putility.print_terminal_partition(level=3)
        print("report: ")
        name_function = (
            "read_extract_from_all_ldsc_files_in_directory()"
        )
        print(name_function)
        putility.print_terminal_partition(level=4)
        print("Names of relevant files within parent directory:")
        print(names_files)
        putility.print_terminal_partition(level=5)
        print("Table after collection:")
        print(table)
        pass
    # Return information.
    return table


##########
# Organization of information from the extraction, collection tables for
# subsequent integration and analysis.


def read_organize_table_ldsc_correlation_single(
    path_file_table=None,
    types_columns=None,
    report=None,
):
    """
    Reads from file and organizes a table as a Pandas data frame. The original
    source table is in a tab-delimited text file.

    The source table represents information about genetic correlations between
    primary and secondary genome-wide association studies (GWAS's) after
    extraction of this information from individual text logs created by the tool
    Linkage Disequilibrium Score Regression (LDSC). This function handles a
    table from a single file that represents indentifier information about both
    primary and secondary studies within the table's column 'name_file'.

    arguments:
        path_file_table (str): path to file for original source table
        types_columns (dict<str>): types of variables in each column
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Read information from file.
    table_raw = pandas.read_csv(
        path_file_table,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=["nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",],
    )

    # Copy information in table.
    table = table_raw.copy(deep=True)
    # Extract names of primary and secondary studies.
    table["study_primary"] = table.apply(
        lambda row:
            str(row["name_file"]).strip().replace(".log", "").split("_-_")[0],
        axis="columns", # apply function to each row
    )
    table["study_secondary"] = table.apply(
        lambda row:
            str(row["name_file"]).strip().replace(".log", "").split("_-_")[1],
        axis="columns", # apply function to each row
    )
    # Remove unnecessary columns.
    table.drop(
        labels=[
            "path_directory",
            "name_file",
            "summary_correlation_error",
            "summary_correlation_ci95",
            "summary_correlation_ci99",
            "covariance",
            "covariance_error",
            "z_score",
            "correlation_ci95_not_zero",
            "correlation_ci95_not_one",
            #"correlation_absolute",
            "correlation_ci95_low",
            "correlation_ci95_high",
            "correlation_ci99_low",
            "correlation_ci99_high",
        ],
        axis="columns",
        inplace=True
    )

    # Report.
    if report:
        putility.print_terminal_partition(level=4)
        print("Source table before organization:")
        print(table_raw)
        print("Column labels:")
        labels_columns = table_raw.columns.to_list()
        print(labels_columns)
        putility.print_terminal_partition(level=4)
        print("Source table after organization:")
        print(table)
        print("Column labels:")
        labels_columns = table.columns.to_list()
        print(labels_columns)
        putility.print_terminal_partition(level=4)

    # Return information.
    return table


def read_organize_table_ldsc_correlation_multiple(
    path_directory_parent=None,
    types_columns=None,
    report=None,
):
    """
    Reads from file and organizes a table as a Pandas data frame. The original
    source tables are in a tab-delimited text file.

    The source tables represent information about genetic correlations between
    primary and secondary genome-wide association studies (GWAS's) after
    extraction of this information from individual text logs created by the tool
    Linkage Disequilibrium Score Regression (LDSC). This function handles tables
    from multiple files that represents indentifier information about the
    primary study in the name of the file for the original table and about the
    secondary study within the table's column 'name_file'.

    arguments:
        path_directory_parent (str): path to parent director for original source
            files
        types_columns (dict<str>): types of variables in each column
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Read all matching files within parent directory and organize paths to
    # these files.
    paths = putility.read_paths_match_child_files_within_parent_directory(
        path_directory_parent=path_directory_parent,
        name_file_child_prefix="table_",
        name_file_child_suffix=".tsv",
        name_file_child_not="blabbergaster",
        report=report,
    )

    # Read files as Pandas dataframe tables.
    # Iterate on names of files to read and organize tables.
    # Collect tables.
    switch = 0
    pail = dict()
    for path in paths:
        # Extract name of file and table that distinguishes it from all others.
        name_file = os.path.basename(path)
        name_table = name_file.replace(str(".tsv"), "")
        name_study_primary = name_table.replace(str("table_"), "")
        # Read information from file.
        table_raw = pandas.read_csv(
            path,
            sep="\t",
            header=0,
            dtype=types_columns,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
        )
        # Store information about primary study.
        table_raw["study_primary"] = name_study_primary
        # Extract names of primary and secondary studies.
        table_raw["study_secondary"] = table_raw.apply(
            lambda row:
                str(row["name_file"]).strip().replace(".log", ""),
            axis="columns", # apply function to each row
        )
        # Remove unnecessary columns.
        table_raw.drop(
            labels=[
                "path_directory",
                "name_file",
                "summary_correlation_error",
                "summary_correlation_ci95",
                "summary_correlation_ci99",
                "covariance",
                "covariance_error",
                "z_score",
                "correlation_ci95_not_zero",
                "correlation_ci95_not_one",
                #"correlation_absolute",
                "correlation_ci95_low",
                "correlation_ci95_high",
                "correlation_ci99_low",
                "correlation_ci99_high",
            ],
            axis="columns",
            inplace=True
        )
        # Concatenate new table with aggregation table.
        if switch == 1:
            # Concatenate new table with aggregation table.
            table = pandas.concat(
                [table, table_raw,],
                axis="index",
                join="outer",
                ignore_index=True,
                copy=True,
            )
        else:
            # Copy information in table.
            table = table_raw.copy(deep=True)
            # Change switch for all instances after first.
            switch = 1
            pass
        pass

    # Report.
    if report:
        putility.print_terminal_partition(level=4)
        print("Concatenation table after organization:")
        print(table)
        print("Column labels:")
        labels_columns = table.columns.to_list()
        print(labels_columns)
        putility.print_terminal_partition(level=4)
    # Return information.
    return table


def filter_table_ldsc_correlation_studies(
    table=None,
    studies_keep=None,
    threshold_p=None,
    report=None,
):
    """
    Filters the genetic correlations that a Pandas data-frame table represents.

    arguments:
        table (object): Pandas data-frame table
        studies_keep (list<str>): identifiers or names of primary and or
            secondary studies for which to keep information in table
        threshold_p (float): threshold on p-values in table's column
            'p_not_zero' below which to keep information in table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    table_filter = table.loc[
        (
            table["study_primary"].isin(studies_keep) &
            table["study_secondary"].isin(studies_keep)
        ), :
    ]
    table_filter = table_filter.loc[
        table_filter["p_not_zero"] < threshold_p, :
    ]

    # Report.
    if report:
        putility.print_terminal_partition(level=4)
        count_rows_table_source = (table.shape[0])
        count_rows_table_product = (table_filter.shape[0])
        print("Count of rows in source table: " + str(count_rows_table_source))
        print(
            "Count of rows in product table: " +
            str(count_rows_table_product)
        )
        putility.print_terminal_partition(level=4)
        print("Count of studies told to keep: " + str(len(studies_keep)))
        print("List of primary and or secondary studies for which to keep")
        print("genetic correlations:")
        print(studies_keep)
        putility.print_terminal_partition(level=4)
        print("Threshold on p-value: " + str(threshold_p))
        putility.print_terminal_partition(level=4)
        print("Table after filters: ")
        print(table_filter)
        putility.print_terminal_partition(level=4)
    # Return information.
    return table_filter








###############################################################################
# Procedure
# Currently, this module is not executable.

##########
