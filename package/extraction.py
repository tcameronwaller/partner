"""
Supply functionality for extraction of information from reports of other tools.

This module is not directly executable.

This module within subpackage 'partner' provides executable functionality under
the management of a higher level package. Importation paths must represent this
hierarchy.

Review:
On 11 March 2024, TCW checked that the extraction values matched those in the
raw report text logs from LDSC for SNP heritability and genetic correlation.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Monroe, North Carolina 28110
    United States of America

License:

    This file is part of project 'partner'
    (https://github.com/tcameronwaller/partner/).

    Project 'partner' supports data analysis in multiple other projects.
    Copyright (C) 2024 Thomas Cameron Waller

    Project 'partner' is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    Project 'partner' is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'partner'.
    If not, see <http://www.gnu.org/licenses/>.
"""

###############################################################################
# Notes

# TODO: TCW; 6 March 2024
# TODO: If the genetic correlation estimate is missing ("nan") with the explanation
# "rg out of bounds", then access the rg, se, p-value, z-score, etc from the
# tab-delimited summary information at the bottom of the report log.


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
import partner.utility as putly # this import path for subpackage
import partner.scale as pale
import partner.description as pdesc

#dir()
#importlib.reload()

###############################################################################
# Functionality


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
    lambda_gc = float("nan")
    chi_square = float("nan")
    intercept = float("nan")
    intercept_error = float("nan")
    ratio = float("nan")
    ratio_error = float("nan")
    summary_heritability_error = str("NA (NA)")
    summary_heritability_ci95 = str("(h2: NA; 95% CI: NA ... NA)")
    summary_heritability_ci99 = str("(h2: NA; 99% CI: NA ... NA)")

    # Read relevant lines of character strings from file.
    lines = putly.read_file_text_lines(
        path_file=path_file,
        start=22,
        stop=35,
    )

    # Define character strings that indicate relevant information.
    prefix_variants = "After merging with regression SNP LD, "
    suffix_variants = " SNPs remain."
    prefix_lambda_gc = "Lambda GC: "
    prefix_chi_square = "Mean Chi^2: "
    prefix_intercept = "Intercept: "
    prefix_ratio = "Ratio: "
    prefix_heritability_observed = "Total Observed scale h2: "
    prefix_heritability_liability = "Total Liability scale h2: "

    # Extract information from lines.
    for line in lines:
        if prefix_variants in line:
            variants = int(
                line.replace(prefix_variants, "").replace(
                    suffix_variants, ""
                ).strip()
            )
        pass
    for line in lines:
        if prefix_heritability_observed in line:
            content = line.replace(prefix_heritability_observed, "")
            contents = content.split(" (")
            heritability_test = contents[0]
            if (
                (not "nan" in heritability_test) and
                (not "NA" in heritability_test)
            ):
                heritability = float(contents[0].strip())
                heritability_error = float(contents[1].replace(")", "").strip())
            pass
        elif prefix_heritability_liability in line:
            content = line.replace(prefix_heritability_liability, "")
            contents = content.split(" (")
            heritability_test = contents[0]
            if (
                (not "nan" in heritability_test) and
                (not "NA" in heritability_test)
            ):
                heritability = float(contents[0].strip())
                heritability_error = float(contents[1].replace(")", "").strip())
            pass
        pass
    for line in lines:
        if (
            (not math.isnan(heritability)) and
            (prefix_lambda_gc in line)
        ):
            lambda_gc = float(line.replace(prefix_lambda_gc, "").strip())
            pass
        elif (
            (not math.isnan(heritability)) and
            (prefix_chi_square in line)
        ):
            chi_square = float(line.replace(prefix_chi_square, "").strip())
            pass
        elif (
            (not math.isnan(heritability)) and
            (prefix_intercept in line)
        ):
            content = line.replace(prefix_intercept, "")
            contents = content.split(" (")
            intercept_test = contents[0]
            if (
                (not "nan" in intercept_test) and
                (not "NA" in intercept_test)
            ):
                intercept = float(contents[0].strip())
                intercept_error = float(contents[1].replace(")", "").strip())
            pass
        elif (
            (not math.isnan(heritability)) and
            (prefix_ratio in line)
        ):
            content = line.replace(prefix_ratio, "")
            contents = content.split(" (")
            ratio_test = contents[0]
            if (
                (not "nan" in ratio_test) and
                (not "NA" in ratio_test)
            ):
                ratio = float(contents[0].strip())
                ratio_error = float(contents[1].replace(")", "").strip())
            pass
        pass

    # Organize information.
    if (
        (not pandas.isna(heritability)) and
        (not pandas.isna(heritability_error))
    ):
        # Determine confidence intervals.
        pail_ci = pdesc.determine_95_99_confidence_intervals_ranges(
            estimate=heritability,
            standard_error=heritability_error,
        )
        heritability_ci95_low = pail_ci["range_95_low"]
        heritability_ci95_high = pail_ci["range_95_high"]
        heritability_ci99_low = pail_ci["range_99_low"]
        heritability_ci99_high = pail_ci["range_99_high"]
        # Create text summaries.
        summary_heritability_error = str(
            str(heritability) + " (" + str(round(heritability_error, 4)) + ")"
        )
        summary_heritability_ci95 = str(
            "(h2: " + str(heritability) +
            "; 95% CI: " + str(pail_ci["range_95"]) + ")"
        )
        summary_heritability_ci99 = str(
            "(h2: " + str(heritability) +
            "; 99% CI: " + str(pail_ci["range_99"]) + ")"
        )
        pass

    # Collect information.
    record = dict()
    record["variants"] = variants
    record["heritability"] = heritability
    record["heritability_error"] = heritability_error
    record["heritability_ci95_low"] = heritability_ci95_low
    record["heritability_ci95_high"] = heritability_ci95_high
    record["heritability_ci99_low"] = heritability_ci99_low
    record["heritability_ci99_high"] = heritability_ci99_high
    record["lambda_gc"] = lambda_gc
    record["chi_square"] = chi_square
    record["intercept"] = intercept
    record["intercept_error"] = intercept_error
    record["ratio"] = ratio
    record["ratio_error"] = ratio_error
    record["summary_heritability_error"] = summary_heritability_error
    record["summary_heritability_ci95"] = summary_heritability_ci95
    record["summary_heritability_ci99"] = summary_heritability_ci99
    # Define list of variables in order.
    variables = [
        "variants",
        "heritability",
        "heritability_error",
        "heritability_ci95_low",
        "heritability_ci95_high",
        "heritability_ci99_low",
        "heritability_ci99_high",
        "lambda_gc",
        "chi_square",
        "intercept",
        "intercept_error",
        "ratio",
        "ratio_error",
        "summary_heritability_error",
        "summary_heritability_ci95",
        "summary_heritability_ci99",
    ]
    # Return information.
    pail = dict()
    pail["record"] = record
    pail["variables"] = variables
    return pail


def read_extract_ldsc_correlation_values_body(
    path_file=None,
):
    """
    Reads and extracts information from log file of LDSC for estimation of
    genetic correlation between two sets of GWAS summary statistics.

    This function specifically extracts values from the main lines of the text
    log report. This function recognizes the relevant lines by specific,
    consistent text prefixes.

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

    # Initialize variables for extraction.
    covariance = float("nan")
    covariance_error = float("nan")
    intercept = float("nan")
    intercept_error = float("nan")
    correlation = float("nan")
    correlation_error = float("nan")
    z_statistic_ldsc = float("nan")
    p_value_ldsc = float("nan")

    # Read relevant lines of character strings from file.
    lines_covariance = putly.read_file_text_lines(
        path_file=path_file,
        start=45,
        stop=51,
    )
    # Define character strings that indicate relevant information.
    prefix_covariance = "Total Observed scale gencov: "
    prefix_intercept = "Intercept: "
    # Extract information from lines.
    for line in lines_covariance:
        if (prefix_covariance in line):
            content = line.replace(prefix_covariance, "")
            contents = content.split(" (")
            covariance_test = contents[0]
            if (
                (not "nan" in covariance_test) and
                (not "NA" in covariance_test)
            ):
                covariance = float(contents[0].strip())
                covariance_error = float(
                    contents[1].replace(")", "").strip()
                )
            pass
        elif (prefix_intercept in line):
            content = line.replace(prefix_intercept, "")
            contents = content.split(" (")
            intercept_test = contents[0]
            if (
                (not "nan" in intercept_test) and
                (not "NA" in intercept_test)
            ):
                intercept = float(contents[0].strip())
                intercept_error = float(
                    contents[1].replace(")", "").strip()
                )
            pass
        pass

    # Read relevant lines of character strings from file.
    lines_correlation = putly.read_file_text_lines(
        path_file=path_file,
        start=51,
        stop=60,
    )
    # Define character strings that indicate relevant information.
    prefix_correlation = "Genetic Correlation: "
    prefix_z = "Z-score: "
    prefix_p = "P: "
    # Extract information from lines.
    for line in lines_correlation:
        if prefix_correlation in line:
            content = line.replace(prefix_correlation, "")
            contents = content.split(" (")
            correlation_test = contents[0]
            if (
                (not "nan" in correlation_test) and
                (not "NA" in correlation_test)
            ):
                correlation = float(contents[0].strip())
                correlation_error = float(contents[1].replace(")", "").strip())
            pass
        elif (
            (not math.isnan(correlation)) and
            (prefix_z in line)
        ):
            z_statistic_ldsc = float(line.replace(prefix_z, "").strip())
            pass
        elif (
            (not math.isnan(correlation)) and
            (prefix_p in line)
        ):
            p_value_ldsc = float(line.replace(prefix_p, "").strip())
            pass
        pass

    # Collect information.
    pail = dict()
    pail["covariance"] = covariance
    pail["covariance_error"] = covariance_error
    pail["intercept"] = intercept
    pail["intercept_error"] = intercept_error
    pail["correlation"] = correlation
    pail["correlation_error"] = correlation_error
    pail["z_statistic_ldsc"] = z_statistic_ldsc
    pail["p_value_ldsc"] = p_value_ldsc
    # Return information.
    return pail


def read_extract_ldsc_correlation(
    path_file=None,
):
    """
    Reads and extracts information from text log file that LDSC creates to
    report estimates of genetic correlation between two sets of GWAS summary
    statistics.

    Hypothesis Tests

    LDSC reports the Z-statistic and p-value corresponding to the null
    hypothesis that the genetic correlation is zero. For a reference on how to
    calculate Z-statistics and p-values for whether the genetic correlation is
    not zero or less than one respectively, refer to the Supplement and
    Supplemental Table 5 in Blokland et al, Biological Psychiatry, 2022
    (PubMed:34099189).

    Alternative hypotheses for genetic correlation (rg):
    1. rg is either less than or greater than zero (two-tailed test)
       Z = (rg) / (standard error of rg)
    2. rg is less than positive one (one-tailed test; right tail)
       Z = (1 - rg) / (standard error of rg)

    References

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
    #lines = putly.read_file_text_lines(
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
    correlation = float("nan")
    correlation_error = float("nan")
    z_statistic_ldsc = float("nan")
    p_value_ldsc = float("nan")
    z_statistic_not_zero = float("nan")
    p_value_not_zero = float("nan")
    z_statistic_less_one = float("nan")
    p_value_less_one = float("nan")
    correlation_ci95_low = float("nan")
    correlation_ci95_high = float("nan")
    correlation_ci99_low = float("nan")
    correlation_ci99_high = float("nan")

    covariance = float("nan")
    covariance_error = float("nan")
    intercept = float("nan")
    intercept_error = float("nan")
    correlation_absolute = float("nan")
    correlation_ci95_not_zero = float("nan")
    correlation_ci95_not_one = float("nan")
    summary_correlation_error = str("NA (NA)")
    summary_correlation_ci95 = str("(rg: NA; 95% CI: NA ... NA)")
    summary_correlation_ci99 = str("(rg: NA; 99% CI: NA ... NA)")

    # Extract information from main body of report log.
    pail_body = read_extract_ldsc_correlation_values_body(
        path_file=path_file,
    )

    # Read relevant lines of character strings from file.
    lines_variants = putly.read_file_text_lines(
        path_file=path_file,
        start=25,
        stop=29,
    )
    # Define character strings that indicate relevant information.
    prefix_vars = "After merging with summary statistics, "
    suffix_vars = " SNPs remain."
    prefix_variants_valid = ""
    suffix_variants_valid = " SNPs with valid alleles."
    # Extract information from lines.
    for line in lines_variants:
        if prefix_vars in line:
            variants = int(
                line.replace(prefix_vars, "").replace(suffix_vars, "").strip()
            )
        elif suffix_variants_valid in line:
            variants_valid = int(
                line.replace(suffix_variants_valid, "").strip()
            )
        pass

    # Determine whether to use values extracted from main body or summary.
    covariance = pail_body["covariance"]
    covariance_error = pail_body["covariance_error"]
    intercept = pail_body["intercept"]
    intercept_error = pail_body["intercept_error"]
    if (
        (not pandas.isna(pail_body["correlation"])) and
        (not pandas.isna(pail_body["correlation_error"]))
    ):
        correlation = pail_body["correlation"]
        correlation_error = pail_body["correlation_error"]
        z_statistic_ldsc = pail_body["z_statistic_ldsc"]
        p_value_ldsc = pail_body["p_value_ldsc"]
    else:
        # Extract information from summary of report log.
        # This block only assigns new, non-missing values if the summary exists.
        summary_check = "Summary of Genetic Correlation Results"
        with open(path_file, 'r') as file:
            line_summary_check = file.readlines()[-6]
        with open(path_file, 'r') as file:
            line_summary_values = file.readlines()[-4]
        if (summary_check in line_summary_check):
            summary_values = line_summary_values.split()
            correlation = float(summary_values[2].strip())
            correlation_error = float(summary_values[3].strip())
            z_statistic_ldsc = float(summary_values[4].strip())
            p_value_ldsc = float(summary_values[5].strip())
            pass
        pass

    # Organize information.
    if (
        (not pandas.isna(correlation)) and
        (not pandas.isna(correlation_error))
    ):
        # Determine Z-statistics and p-values for null hypotheses, respectively.
        # Do not round Z-statistics or p-values to avoid loss of information.
        z_statistic_not_zero = float(correlation / correlation_error)
        p_value_not_zero = pale.calculate_p_value_from_z_statistic(
            z_statistic=z_statistic_not_zero,
            tail="both", # two-tailed test; less than or greater than zero
        )
        z_statistic_less_one = float((1 - correlation) / correlation_error)
        p_value_less_one = pale.calculate_p_value_from_z_statistic(
            z_statistic=z_statistic_less_one,
            tail="right", # one-tailed test; less than positive one
        )
        # Determine confidence intervals.
        pail_ci = pdesc.determine_95_99_confidence_intervals_ranges(
            estimate=correlation,
            standard_error=correlation_error,
        )
        correlation_ci95_low = pail_ci["range_95_low"]
        correlation_ci95_high = pail_ci["range_95_high"]
        correlation_ci99_low = pail_ci["range_99_low"]
        correlation_ci99_high = pail_ci["range_99_high"]
        # Determine absolute value of correlation.
        correlation_absolute = math.fabs(correlation)
        # Create text summaries.
        summary_correlation_error = str(
            str(correlation) + " (" + str(round(correlation_error, 4)) + ")"
        )
        summary_correlation_ci95 = str(
            "(rg: " + str(correlation) +
            "; 95% CI: " + str(pail_ci["range_95"]) + ")"
        )
        summary_correlation_ci99 = str(
            "(rg: " + str(correlation) +
            "; 99% CI: " + str(pail_ci["range_99"]) + ")"
        )
        # Determine whether confidence interval crosses zero or one.
        if (
            (
                (pail_ci["range_95_low"] > 0) and (pail_ci["range_95_high"] > 0)
            ) or
            (
                (pail_ci["range_95_low"] < 0) and (pail_ci["range_95_high"] < 0)
            )
        ):
            correlation_ci95_not_zero = 1
        else:
            correlation_ci95_not_zero = 0
            pass
        if (
            (pail_ci["range_95_low"] < 1.0) and (pail_ci["range_95_high"] < 1.0)
        ):
            correlation_ci95_not_one = 1
        else:
            correlation_ci95_not_one = 0
            pass
        pass

    # Collect information.
    record = dict()
    record["variants"] = variants
    record["variants_valid"] = variants_valid
    record["correlation"] = correlation
    record["correlation_error"] = correlation_error
    record["z_statistic_ldsc"] = z_statistic_ldsc
    record["p_value_ldsc"] = p_value_ldsc
    record["z_statistic_not_zero"] = z_statistic_not_zero
    record["p_value_not_zero"] = p_value_not_zero
    record["z_statistic_less_one"] = z_statistic_less_one
    record["p_value_less_one"] = p_value_less_one
    record["correlation_ci95_low"] = correlation_ci95_low
    record["correlation_ci95_high"] = correlation_ci95_high
    record["correlation_ci99_low"] = correlation_ci99_low
    record["correlation_ci99_high"] = correlation_ci99_high
    record["covariance"] = covariance
    record["covariance_error"] = covariance_error
    record["intercept"] = intercept
    record["intercept_error"] = intercept_error
    record["correlation_absolute"] = correlation_absolute
    record["correlation_ci95_not_zero"] = correlation_ci95_not_zero
    record["correlation_ci95_not_one"] = correlation_ci95_not_one
    record["summary_correlation_error"] = summary_correlation_error
    record["summary_correlation_ci95"] = summary_correlation_ci95
    record["summary_correlation_ci99"] = summary_correlation_ci99
    # Define list of variables in order.
    variables = [
        "variants",
        "variants_valid",
        "correlation",
        "correlation_error",
        "z_statistic_ldsc",
        "p_value_ldsc",
        "z_statistic_not_zero",
        "p_value_not_zero",
        "z_statistic_less_one",
        "p_value_less_one",
        "correlation_ci95_low",
        "correlation_ci95_high",
        "correlation_ci99_low",
        "correlation_ci99_high",
        "covariance",
        "covariance_error",
        "intercept",
        "intercept_error",
        "correlation_absolute",
        "correlation_ci95_not_zero",
        "correlation_ci95_not_one",
        "summary_correlation_error",
        "summary_correlation_ci95",
        "summary_correlation_ci99",
    ]
    # Return information.
    pail = dict()
    pail["record"] = record
    pail["variables"] = variables
    return pail


def read_extract_from_all_ldsc_files_in_directory(
    path_directory=None,
    name_file_prefix=None,
    name_file_suffix=None,
    name_file_not=None,
    type_analysis=None,
    report=None,
):
    """
    Reads and extracts information about LDSC SNP heritability or genetic
        correlation analyses from all files within a parent directory.

        Linkage Disequilibrium Score Regression (LDSC)

    arguments:
        path_directory (str): full path to parent directory that contains
            relevant files
        name_file_prefix (str): string prefix in names of relevant child files
            within parent directory
        name_file_suffix (str): string suffix in names of relevant child files
            within parent directory
        name_file_not (str): string not in names of relevant child files within
            parent directory
        type_analysis (str): type of analysis in LDSC, either 'heritability' or
            'correlation', corresponding to information for extraction
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table with variables across columns and
            samples or instances across rows

    """

    # Extract and filter names of child files within parent directory.
    names_files = putly.extract_filter_child_file_names(
        path_directory=path_directory,
        name_file_prefix=name_file_prefix,
        name_file_suffix=name_file_suffix,
        name_file_not=name_file_not,
    )
    putly.print_terminal_partition(level=4)
    putly.print_terminal_partition(level=4)
    putly.print_terminal_partition(level=4)
    print("From extraction driver function")
    print(len(names_files))
    putly.print_terminal_partition(level=4)
    putly.print_terminal_partition(level=4)
    putly.print_terminal_partition(level=4)
    
    # Collect information from analysis on each study.
    records = list()
    # Iterate on files for each study.
    counter = 0
    for name_file in names_files:
        # Define full path to file.
        path_file = os.path.join(path_directory, name_file,)
        # Extract information from LDSC analysis.
        if (str(type_analysis).strip() == "heritability"):
            pail = read_extract_ldsc_heritability(
                path_file=path_file,
            )
        elif (str(type_analysis).strip() == "correlation"):
            pail = read_extract_ldsc_correlation(
                path_file=path_file,
            )
            pass
        # Collect and organize information about study.
        # Extraction functions return the relevant variables for each type of
        # analysis.
        pail["record"]["path_directory"] = str(path_directory).strip()
        pail["record"]["name_file"] = str(name_file).strip()
        pail["record"]["type_analysis"] = str(type_analysis).strip()
        records.append(pail["record"])
        if (counter == 0):
            variables = pail["variables"]
        # Increment counter.
        counter += 1
        pass
    # Organize table.
    table = putly.convert_records_to_dataframe(
        records=records
    )
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # columns.insert(0, dependence)
    columns = ["path_directory", "name_file", "type_analysis",]
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
        putly.print_terminal_partition(level=3)
        print("report: ")
        name_function = (
            "read_extract_from_all_ldsc_files_in_directory()"
        )
        print(name_function)
        putly.print_terminal_partition(level=4)
        #print("Names of relevant files within parent directory:")
        #print(names_files)
        putly.print_terminal_partition(level=5)
        print("Table after collection:")
        print(table)
        pass
    # Return information.
    return table


##########
# Organization of information from the extraction, collection tables for
# subsequent integration and analysis.


def define_genetic_correlation_table_column_types():
    """
    Defines the variable types of columns within table for genetic correlations.

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["path_directory"] = "string"
    types_columns["name_file"] = "string"
    types_columns["variants"] = "float"
    types_columns["variants_valid"] = "float"
    types_columns["correlation"] = "float"
    types_columns["correlation_error"] = "float"
    types_columns["z_statistic_ldsc"] = "float"
    types_columns["p_value_ldsc"] = "float"
    types_columns["z_statistic_not_zero"] = "float"
    types_columns["p_value_not_zero"] = "float"
    types_columns["z_statistic_less_one"] = "float"
    types_columns["p_value_less_one"] = "float"
    types_columns["correlation_ci95_low"] = "float"
    types_columns["correlation_ci95_high"] = "float"
    types_columns["correlation_ci99_low"] = "float"
    types_columns["correlation_ci99_high"] = "float"
    types_columns["covariance"] = "float"
    types_columns["covariance_error"] = "float"
    types_columns["intercept"] = "float"
    types_columns["intercept_error"] = "float"
    types_columns["correlation_absolute"] = "float"
    types_columns["correlation_ci95_not_zero"] = "float"
    types_columns["correlation_ci95_not_one"] = "float"
    types_columns["summary_correlation_error"] = "string"
    types_columns["summary_correlation_ci95"] = "string"
    types_columns["summary_correlation_ci99"] = "string"
    # Return information.
    return types_columns


def define_genetic_correlation_table_column_sequence():
    """
    Defines the columns in sequence within table for genetic correlations.

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify sequence of columns within table.
    columns_sequence = [
        #"path_directory",
        #"name_file",
        #"sort_primary",
        #"sort_secondary",
        "study_primary",
        "study_secondary",
        "abbreviation_primary",
        "abbreviation_secondary",
        "trait_primary",
        "trait_secondary",
        "sex_primary",
        "sex_secondary",
        #"variants",
        "variants_valid",
        "correlation",
        "correlation_error",
        "z_statistic_ldsc",
        "p_value_ldsc",
        "q_value_ldsc",
        "z_statistic_not_zero",
        "p_value_not_zero",
        "q_value_not_zero",
        "z_statistic_less_one",
        "p_value_less_one",
        "q_value_less_one",
        #"q_significance",
        "correlation_ci95_low",
        "correlation_ci95_high",
        "correlation_ci99_low",
        "correlation_ci99_high",
        #"covariance",
        #"covariance_error",
        #"intercept",
        #"intercept_error",
        #"correlation_absolute",
        #"correlation_ci95_not_zero",
        #"correlation_ci95_not_one",
        "summary_correlation_error",
        "summary_correlation_ci95",
        "summary_correlation_ci99",
    ]
    # Return information.
    return columns_sequence


def read_organize_table_ldsc_correlation_single(
    path_file_table=None,
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
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Define variable types of columns within table.
    types_columns = define_genetic_correlation_table_column_types()
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
    if False:
        table.drop(
            labels=[
                "path_directory",
                "name_file",
                "correlation_ci95_low",
                "correlation_ci95_high",
                "correlation_ci99_low",
                "correlation_ci99_high",
                "covariance",
                "covariance_error",
                "correlation_ci95_not_zero",
                "correlation_ci95_not_one",
                "correlation_absolute",
                "summary_correlation_error",
                "summary_correlation_ci95",
                "summary_correlation_ci99",
            ],
            axis="columns",
            inplace=True
        )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Source table before organization:")
        print(table_raw)
        print("Column labels:")
        labels_columns = table_raw.columns.to_list()
        print(labels_columns)
        putly.print_terminal_partition(level=4)
        print("Source table after organization:")
        print(table)
        print("Column labels:")
        labels_columns = table.columns.to_list()
        print(labels_columns)
        putly.print_terminal_partition(level=4)

    # Return information.
    return table


def read_organize_table_ldsc_correlation_multiple(
    path_directory_parent=None,
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
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Define variable types of columns within table.
    types_columns = define_genetic_correlation_table_column_types()
    # Read all matching files within parent directory and organize paths to
    # these files.
    paths = putly.read_paths_match_child_files_within_parent_directory(
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
        if False:
            table.drop(
                labels=[
                    "path_directory",
                    "name_file",
                    "correlation_ci95_low",
                    "correlation_ci95_high",
                    "correlation_ci99_low",
                    "correlation_ci99_high",
                    "covariance",
                    "covariance_error",
                    "correlation_ci95_not_zero",
                    "correlation_ci95_not_one",
                    "correlation_absolute",
                    "summary_correlation_error",
                    "summary_correlation_ci95",
                    "summary_correlation_ci99",
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
        putly.print_terminal_partition(level=4)
        print("Concatenation table after organization:")
        print(table)
        print("Column labels:")
        labels_columns = table.columns.to_list()
        print(labels_columns)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table


def filter_table_rows_ldsc_correlation(
    table=None,
    studies_primary_keep=None,
    studies_secondary_keep=None,
    name_primary=None,
    name_secondary=None,
    remove_redundancy=None,
    remove_self_pair=None,
    report=None,
):
    """
    Filters rows within a Pandas data-frame table of LDSC genetic correlations
    on the basis of pairs of primary and secondary studies, which can optionally
    be considered as interchangeable, depending on the application.

    arguments:
        table (object): Pandas data-frame table
        studies_primary_keep (list<str>): identifiers or names of primary
            studies for which to keep information in table
        studies_secondary_keep (list<str>): identifiers or names of secondary
            studies for which to keep information in table
        name_primary (str): name of column for primary identifier
        name_secondary (str): name of column for secondary identifier
        remove_redundancy (bool): whether to remove redundant pairs of primary
            and secondary studies
        remove_self_pair (bool): whether to remove all self pairs of identical
            primary and secondary studies; otherwise only remove the redundant
            occurrences if 'remove_redundancy' is True
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_filter = table.copy(deep=True)
    # Filter table's rows by identifiers of studies.
    table_filter = table_filter.loc[
        (
            table_filter[name_primary].isin(studies_primary_keep) &
            table_filter[name_secondary].isin(studies_secondary_keep)
        ), :
    ]
    # Remove redundant pairs of interchangeable primary and secondary studies.
    if remove_redundancy:
        # Organize information in table.
        # This action creates the "index" column necessary for the next filter step.
        table_filter.reset_index(
            level=None,
            inplace=True,
            drop=False, # remove index; do not move to regular columns
        )
        # Filter redundancy.
        table_filter["redundancy_studies"] = table_filter.apply(
            lambda row:
                porg.match_table_row_unique_interchangeable_identifiers_a_b(
                    table=table_filter,
                    name_index="index",
                    name_primary="study_primary",
                    name_secondary="study_secondary",
                    row_index=row["index"],
                    row_primary=row["study_primary"],
                    row_secondary=row["study_secondary"],
                    match_self_pair=remove_self_pair, # return 1 for self pairs
                    report=False,
                ),
            axis="columns", # apply function to each row
        )
        table_filter = table_filter.loc[
            (table_filter["redundancy_studies"] == 0), :
        ]
        # Remove unnecessary columns.
        table_filter.drop(
            labels=["index", "redundancy_studies",],
            axis="columns",
            inplace=True
        )
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        count_rows_table_source = (table.shape[0])
        count_rows_table_product = (table_filter.shape[0])
        print("Count of rows in source table: " + str(count_rows_table_source))
        print(
            "Count of rows in product table: " +
            str(count_rows_table_product)
        )
        putly.print_terminal_partition(level=4)
        print("Primary studies to keep: ")
        print("Count: " + str(len(studies_primary_keep)))
        print("Studies: ")
        print(studies_primary_keep)
        putly.print_terminal_partition(level=4)
        print("Secondary studies to keep: ")
        print("Count: " + str(len(studies_secondary_keep)))
        print("Studies: ")
        print(studies_secondary_keep)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_filter



# TODO: TCW; 27 March 2024
# TODO: make this more versatile and move to partner.organization
def transfer_study_attributes_table_ldsc_rg(
    table=None,
    name_primary=None,
    name_secondary=None,
    table_attribute_primary=None,
    table_attribute_secondary=None,
    attributes=None,
    report=None,
):
    """
    Transfer attributes of primary and secondary studies to a Pandas data-frame
    table of LDSC genetic correlations.

    arguments:
        table (object): Pandas data-frame table
        name_primary (str): name of column for primary identifier
        name_secondary (str): name of column for secondary identifier
        table_attribute_primary (object): table of study attributes
        table_attribute_secondary (object): table of study attributes
        attributes (list<str>): names of attributes to transfer
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_attribute_primary = table_attribute_primary.copy(deep=True)
    table_attribute_secondary = table_attribute_secondary.copy(deep=True)
    # Create dictionaries for sort sequences of studies.
    table_attribute_primary.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_attribute_secondary.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_attribute_primary.set_index(
        ["study",],
        append=False,
        drop=True,
        inplace=True,
    )
    table_attribute_secondary.set_index(
        ["study",],
        append=False,
        drop=True,
        inplace=True,
    )
    attribute_primary = table_attribute_primary.to_dict("index")
    attribute_secondary = table_attribute_secondary.to_dict("index")

    # Copy information in table.
    table_attribute = table.copy(deep=True)
    # Organize information in table.
    table_attribute.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Transfer attributes to the table of genetic correlations.
    for attribute in attributes:
        novel_primary = str(attribute + "_primary")
        novel_secondary = str(attribute + "_secondary")
        table_attribute[novel_primary] = table_attribute.apply(
            lambda row: (
                attribute_primary[row[name_primary]][attribute]
            ),
            axis="columns", # apply function to each row
        )
        table_attribute[novel_secondary] = table_attribute.apply(
            lambda row: (
                attribute_secondary[row[name_secondary]][attribute]
            ),
            axis="columns", # apply function to each row
        )
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        count_columns_source = (table.shape[1])
        count_columns_product = (table_attribute.shape[1])
        print("Count of columns in source table: " + str(count_columns_source))
        print(
            "Count of columns in product table: " +
            str(count_columns_product)
        )
        print("Table")
        print(table_attribute)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_attribute


# TODO: TCW; 26 March 2024
# TODO: function "sort_table_rows_ldsc_correlation()" is now obsolete, but I
# still need to test.
# TODO: a more versatile implementation is now in partner.organization
def sort_table_rows_ldsc_correlation(
    table=None,
    name_primary=None,
    name_secondary=None,
    table_sort_primary=None,
    table_sort_secondary=None,
    name_sort_sequence=None,
    name_sort_study=None,
    report=None,
):
    """
    Sorts rows within a Pandas data-frame table of LDSC genetic correlations
    on the basis of sort sequence of primary and secondary studies.

    arguments:
        table (object): Pandas data-frame table
        name_primary (str): name of column for primary identifier
        name_secondary (str): name of column for secondary identifier
        table_sort_primary (object): table of study sort sequence
        table_sort_secondary (object): table of study sort sequence
        name_sort_sequence (str): name of column for sort sequence
        name_sort_study (str): name of column for sort study
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_sort_primary = table_sort_primary.copy(deep=True)
    table_sort_secondary = table_sort_secondary.copy(deep=True)
    # Create dictionaries for sort sequences of studies.
    table_sort_primary.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_sort_secondary.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_sort_primary.set_index(
        [name_sort_study],
        append=False,
        drop=True,
        inplace=True,
    )
    table_sort_secondary.set_index(
        [name_sort_study],
        append=False,
        drop=True,
        inplace=True,
    )
    sort_primary = table_sort_primary.to_dict("index")
    sort_secondary = table_sort_secondary.to_dict("index")

    # Copy information in table.
    table_sort = table.copy(deep=True)
    # Organize information in table.
    table_sort.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Determine sort sequences for primary and secondary studies.
    table_sort["sort_primary"] = table_sort.apply(
        lambda row: (
            sort_primary[row[name_primary]][name_sort_sequence]
        ),
        axis="columns", # apply function to each row
    )
    table_sort["sort_secondary"] = table_sort.apply(
        lambda row: (
            sort_secondary[row[name_secondary]][name_sort_sequence]
        ),
        axis="columns", # apply function to each row
    )
    # Sort rows within table.
    table_sort.sort_values(
        by=["sort_primary", "sort_secondary",],
        axis="index",
        ascending=True,
        inplace=True,
    )
    table_sort.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Remove unnecessary columns.
    if False:
        table_sort.drop(
            labels=["sort_primary", "sort_secondary"],
            axis="columns",
            inplace=True
        )
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        count_rows_table_source = (table.shape[0])
        count_rows_table_product = (table_sort.shape[0])
        print("Count of rows in source table: " + str(count_rows_table_source))
        print(
            "Count of rows in product table: " +
            str(count_rows_table_product)
        )
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_sort


def needs_update_filter_table_ldsc_correlation_studies(
    table=None,
    studies_keep=None,
    threshold_q=None,
    order_columns=None,
    report=None,
):
    """
    Filters the genetic correlations that a Pandas data-frame table represents.

    arguments:
        table (object): Pandas data-frame table
        studies_keep (list<str>): identifiers or names of primary and or
            secondary studies for which to keep information in table
        threshold_q (float): threshold on q-values in table's new column
            'q_value' below which to keep information in table
        order_columns (list<str>): names of columns in order for sort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Filter table's rows by identifiers of studies.
    table_filter = table.loc[
        (
            table["study_primary"].isin(studies_keep) &
            table["study_secondary"].isin(studies_keep)
        ), :
    ]
    # Organize information in table.
    # This action creates the "index" column necessary for the next filter step.
    table_filter.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Remove redundant pairs of interchangeable primary and secondary studies.
    if True:
        table_filter["redundancy_studies"] = table_filter.apply(
            lambda row:
                porg.match_table_row_unique_interchangeable_identifiers_a_b(
                    table=table_filter,
                    name_index="index",
                    name_primary="study_primary",
                    name_secondary="study_secondary",
                    row_index=row["index"],
                    row_primary=row["study_primary"],
                    row_secondary=row["study_secondary"],
                    match_self_pair=True, # whether to return 1 for self pairs
                    report=False,
                ),
            axis="columns", # apply function to each row
        )
        table_filter = table_filter.loc[
            (table_filter["redundancy_studies"] == 0), :
        ]
    # Calculate Benjamini-Hochberg q-values for False Discovery Rate (FDR).
    table_filter_q = putly.calculate_table_false_discovery_rate_q_values(
        threshold=0.05,
        name_column_p_value="p_value",
        name_column_q_value="q_value",
        name_column_significance="q_significance",
        table=table_filter,
    )
    # Remove unnecessary columns.
    table_filter_q.drop(
        labels=["q_significance",],
        axis="columns",
        inplace=True
    )
    # Filter by Benjamini-Hochberg q-value for False Discovery Rate (FDR).
    table_filter_q = table_filter_q.loc[
        table_filter_q["q_not_zero"] < threshold_q, :
    ]
    # Sort table's rows.
    table_filter_q.sort_values(
        by=["study_primary", "study_secondary",],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Filter and sort table's columns.
    #table_filter = table_filter.loc[
    #    :, table_filter.columns.isin(order_columns)
    #]
    table_filter_q = table_filter_q.filter(
        items=order_columns,
        axis="columns",
    )
    table_filter_q = table_filter_q[[*order_columns]]

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        count_rows_table_source = (table.shape[0])
        count_rows_table_product = (table_filter_q.shape[0])
        print("Count of rows in source table: " + str(count_rows_table_source))
        print(
            "Count of rows in product table: " +
            str(count_rows_table_product)
        )
        putly.print_terminal_partition(level=4)
        print("Count of studies told to keep: " + str(len(studies_keep)))
        print("List of primary and or secondary studies for which to keep")
        print("genetic correlations:")
        print(studies_keep)
        putly.print_terminal_partition(level=4)
        print("Threshold on q-value: " + str(threshold_q))
        putly.print_terminal_partition(level=4)
        print("Table after filters: ")
        print(table_filter_q)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_filter_q


##########
# Queries on tables of information about genetic correlations.


def match_table_row_exclusion_identifiers_primary_secondary(
    row_primary=None,
    row_secondary=None,
    exclusions=None,
):
    """
    Dependency:
    This function is a dependency of the function below.
    partner.extraction.query_table_correlation_range_variables()

    Determines whether values of two interchangeable identifiers from a single
    row in a table match a list of exclusions.

    arguments:
        row_primary (str): current row's value of primary identifier
        row_secondary (str): current row's value of secondary identifier
        exclusions (list<dict<str>>): pairs of specific primary and secondary
            studies to exclude from the query

    raises:

    returns:
        (int): binary representation of whether the current row matches
            exclusions

    """

    # Determine whether the current combination of primary and secondary studies
    # matches any exclusions.
    indicator = 0
    for exclusion in exclusions:
        if (
            (row_primary == exclusion["primary"]) and
            (row_secondary == exclusion["secondary"])
        ):
            indicator = 1
            pass
        pass
    return indicator


def query_table_correlation_range_variables(
    table=None,
    label=None,
    studies_primary=None,
    studies_secondary=None,
    name_primary=None,
    name_secondary=None,
    exclusions=None,
    variables=None,
    remove_self_pair=None,
    report=None,
):
    """
    Queries information within a Pandas data-frame table of LDSC genetic
    correlations. Filters by a selection of primary and secondary studies and
    then describes ranges of a selection of continuous, quantitative variables.

    arguments:
        table (object): Pandas data-frame table
        label (str): character string to print at top of report
        studies_primary (list<str>): identifiers or names of primary
            studies for which to describe information in table
        studies_secondary (list<str>): identifiers or names of secondary
            studies for which to describe information in table
        name_primary (str): name of column for primary identifier
        name_secondary (str): name of column for secondary identifier
        exclusions (list<dict<str>>): pairs of specific primary and secondary
            studies to exclude from the query
        variables (list<str>): name of column for primary identifier
        remove_self_pair (bool): whether to remove all self pairs of identical
            primary and secondary studies
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_query = table.copy(deep=True)
    # Filter table's rows by identifiers of studies.
    table_query = table_query.loc[
        (
            table_query[name_primary].isin(studies_primary) &
            table_query[name_secondary].isin(studies_secondary)
        ), :
    ]
    # Remove self pairs of identical primary and secondary studies.
    if remove_self_pair:
        table_query = table_query.loc[
            (
                table_query[name_primary] != table_query[name_secondary]
            ), :
        ]
        pass
    # Exclude specific pairs of primary and secondary studies from the query.
    if (len(exclusions) > 0):
        # Filter pairs for exclusion.
        table_query["exclusion"] = table_query.apply(
            lambda row:
                match_table_row_exclusion_identifiers_primary_secondary(
                    row_primary=row["study_primary"],
                    row_secondary=row["study_secondary"],
                    exclusions=exclusions,
                ),
            axis="columns", # apply function to each row
        )
        table_query = table_query.loc[
            (table_query["exclusion"] == 0), :
        ]
        # Remove unnecessary columns.
        table_query.drop(
            labels=["exclusion",],
            axis="columns",
            inplace=True
        )
        pass
    # Organize information in table.
    table_query.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print(label)
        putly.print_terminal_partition(level=5)
        print("Primary studies in query: ")
        print("Count: " + str(len(studies_primary)))
        print("Studies: ")
        print(studies_primary)
        putly.print_terminal_partition(level=5)
        print("Secondary studies in query: ")
        print("Count: " + str(len(studies_secondary)))
        print("Studies: ")
        print(studies_secondary)
        putly.print_terminal_partition(level=5)
        check_exclusions = (len(exclusions) > 0)
        print("Exclusions: " + str(check_exclusions))
        putly.print_terminal_partition(level=5)
        count_rows_table_source = (table.shape[0])
        count_rows_table_product = (table_query.shape[0])
        print("Count of rows in source table: " + str(count_rows_table_source))
        print(
            "Count of rows in product table: " +
            str(count_rows_table_product)
        )
        #putly.print_terminal_partition(level=5)
        #print("Table after query filters: ")
        #print(table_query)

    # Describe ranges of variables within filtered, query table.
    # Report.
    if report:
        pdesc.report_table_range_variables(
            table=table_query,
            variables=variables,
        )
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_query







###############################################################################
# Procedure
# Currently, this module is not executable.

##########
