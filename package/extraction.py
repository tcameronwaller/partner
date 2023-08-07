"""
Supply functionality for extraction and organization of information from
reports of other tools.

This module is not directly executable.

This subpackage 'partner' provides executable functionality under the
management of a higher level package. Importation paths represent this
hierarchy.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Monroe, North Carolina 28110
    United States of America

License:

    This file is part of Partner
    (https://github.com/tcameronwaller/partner/).

    Partner supports data analysis in multiple other projects.
    Copyright (C) 2023 Thomas Cameron Waller

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
import partner.utility as utility # this import path for subpackage

#dir()
#importlib.reload()

###############################################################################
# Functionality


# TODO: TCW; 12 December 2022; 7 August 2023
# Use the new utility function to calculate and organize 95% and 99% confidence intervals and ranges.

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
    #lines = utility.read_file_text_lines(
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
    p_greater_zero = float("nan")

    # Read relevant lines of character strings from file.
    lines = utility.read_file_text_lines(
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
            p_greater_zero = float(line.replace(prefix_p, ""))
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
    record["p_greater_zero"] = p_greater_zero
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
        "p_greater_zero",
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
    # Return information.
    pail = dict()
    pail["record"] = record
    pail["variables"] = variables
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
