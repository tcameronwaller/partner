"""
Supply functionality for extraction and organization of information from
reports of other tools.

This module is not directly executable.

This subpackage 'promiscuity' provides executable functionality under the
management of a higher level package. Importation paths represent this
hierarchy.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of Promiscuity
    (https://github.com/tcameronwaller/promiscuity/).

    Promiscuity supports data analysis in multiple other projects.
    Copyright (C) 2022 Thomas Cameron Waller

    Promiscuity is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    Promiscuity is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with Promiscuity. If not, see <http://www.gnu.org/licenses/>.
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
import promiscuity.utility as utility # this import path for subpackage

#dir()
#importlib.reload()

###############################################################################
# Functionality


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
    lines = utility.read_file_text_lines(
        path_file=path_file,
        start=22,
        stop=30,
    )

    # Define character strings that indicate relevant information.
    prefix_variants = "After merging with regression SNP LD, "
    suffix_variants = " SNPs remain."
    prefix_heritability = "Total Observed scale h2: "
    prefix_chi_square = "Mean Chi^2: "
    prefix_ratio = "Ratio: "

    # Extract information from lines.
    for line in lines:
        if prefix_variants in line:
            variants = int(
                line.replace(prefix_variants, "").replace(suffix_variants, "")
            )
        elif prefix_heritability in line:
            content = line.replace(prefix_heritability, "")
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
            chi_square = float(line.replace(prefix_ratio, ""))
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
        "(h2: " + str(heritability) + "; 95% CI: " + str(heritability_ci99) + ")"
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



# parameters:
# 1. path to parent directory
# 2. file name suffix by which to recognize relevant files in parent directory
# 3. type of information to extract: "heritability" or "correlation"
# functionality:
# read names of all files within the parent directory
# iterate on all files, calling the extraction function for each
# collect information from the extraction record for each file
# append the file name to the extraction record
# organize the extracted information in a table
# return table


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
    names_files = utility.extract_directory_file_names_filter_by_name(
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
            print("place holder")
            pass
        # Collect and organize information about study.
        record["path_directory"] = str(path_directory)
        record["name_file"] = str(name_file)
        records.append(pail["record"])
        if (counter == 0):
            variables = pail["variables"]
        # Increment counter.
        counter += 1
        pass
    # Organize table.
    table = utility.convert_records_to_dataframe(
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
        utility.print_terminal_partition(level=3)
        print("report: ")
        name_function = (
            "read_extract_from_all_ldsc_files_in_directory()"
        )
        print(name_function)
        utility.print_terminal_partition(level=4)
        print("Names of relevant files within parent directory:")
        print(names_files)
        utility.print_terminal_partition(level=5)
        print("Table after collection:")
        print(table)
        pass
    # Return information.
    return table




###############################################################################
# Procedure
# Currently, this module is not executable.
