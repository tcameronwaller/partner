"""
Supply functionality for description of variables and preparation of tables and
other reports.

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

################################################################################
# Functionality


##########
# Table Attribution


def create_attribution_record(
    name_cohort=None,
    name_variable_value=None,
    variable=None,
    value=None,
    table=None,
):
    """
    Organize a record (single row in table) to describe attribution of
    categorical or discrete variable values across cohorts.

    arguments:
        name_cohort (str): name of cohort
        name_variable_value (str): name of variable's value for report
        variable (str): name of table's column for variable
        value (object): categorical or discrete value of variable
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["cohort"] = str(name_cohort)
    record["variable_value"] = str(name_variable_value)
    record["variable"] = str(variable)
    record["value"] = str(value)
    # Copy information.
    table = table.copy(deep=True)

    # Stratify table.
    # Select relevant rows of the table.
    table_variable_value = table.loc[
        (
            (~pandas.isna(table[variable])) &
            (table[variable] == value)
        ), :
    ]

    # Count records.
    count_total = table.shape[0]
    count_variable_value = table_variable_value.shape[0]

    # Calculate percentages.
    if (count_total > 0):
        percentage_variable_value = round(
            ((count_variable_value / count_total) * 100), 3
        )
    else:
        percentage_variable_value = float("nan")
        pass

    # Collect information for record.
    record["count_cohort_total_records"] = count_total
    record["count_variable_value"] = count_variable_value
    record["count_variable_value_report"] = str(
        str(count_variable_value) +
        " (" + str(percentage_variable_value) + "%)"
    )
    # Return information.
    return record


def drive_assemble_attribution_table(
    records_attribution=None,
    records_cohorts=None,
    report=None,
):
    """
    Drives the assembly of a description table from records for attribution of
    specific values of nominal, categorical, or discrete variables within
    cohorts.

    arguments:
        records_attribution (list<dict>): records with information about values
            of variables for attribution within cohorts, including entries for
            'name', 'variable', and 'value'
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Collect summary records for rows within description table.
    records_description = list()
    # Iterate on cohorts.
    for record_cohort in records_cohorts:
        # Iterate on variables.
        for record_attribution in records_attribution:
            # Organize information for description record.
            record_description = create_attribution_record(
                name_cohort=record_cohort["name"],
                name_variable_value=record_attribution["name"],
                variable=record_attribution["variable"],
                value=record_attribution["value"],
                table=record_cohort["table"],
            )
            # Collect records.
            records_description.append(record_description)
            pass
        pass
    # Organize table.
    table = pandas.DataFrame(data=records_description)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("drive_assemble_attribution_table()")
        utility.print_terminal_partition(level=3)
        print(table)
        pass
    # Return information.
    return table


##########
# Table Quantitation


def create_quantitation_record(
    name_cohort=None,
    variable=None,
    table=None,
):
    """
    Organize a record (single row in table) to describe for measures on
    quantitative variables across cohorts.

    Report percentages relative to the total count of records in the cohort.

    arguments:
        name_cohort (str): name of cohort
        variable (str): name of table's column for variable
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["cohort"] = str(name_cohort)
    record["variable"] = str(variable)
    # Copy information.
    table = table.copy(deep=True)

    # Stratify table.
    # Select relevant rows of the table.
    table_genotype = table.loc[
        (
            (~pandas.isna(table["genotype_availability"])) &
            (table["genotype_availability"] == 1)
        ), :
    ]

    # Count records.
    record["count_cohort_total_records"] = int(table.shape[0])
    record["count_cohort_total_genotypes"] = int(table_genotype.shape[0])
    count_total = record["count_cohort_total_records"]
    count_total_genotype = record["count_cohort_total_genotypes"]
    percentage_total_genotype = round(
        ((count_total_genotype / count_total) * 100), 3
    )
    record["percentage_cohort_total_genotypes"] = str(
        str(record["count_cohort_total_genotypes"]) + " (" +
        str(percentage_total_genotype) + "%)"
    )
    # Initialize missing values.
    record["count_variable_non_missing"] = float("nan")
    record["percentage_variable_non_missing"] = str("nan (nan%)")
    record["count_variable_genotypes"] = float("nan")
    record["percentage_variable_genotypes"] = str("nan (nan%)")
    record["median"] = float("nan")
    record["minimum"] = float("nan")
    record["maximum"] = float("nan")
    record["mean"] = float("nan")
    record["standard_deviation"] = float("nan")
    record["standard_error"] = float("nan")
    record["interval_95"] = float("nan")
    record["range_95_low"] = float("nan")
    record["range_95_high"] = float("nan")
    record["interval_99"] = float("nan")
    record["range_99_low"] = float("nan")
    record["range_99_high"] = float("nan")
    # Determine whether table has the column.
    if (variable in table.columns.to_list()):
        array = copy.deepcopy(table[variable].dropna().to_numpy()) # non-missing
        array_genotype = copy.deepcopy(
            table_genotype[variable].dropna().to_numpy()
        ) # non-missing variable values with genotypes
        # Determine count of valid values.
        record["count_variable_non_missing"] = int(array.size)
        count_variable = record["count_variable_non_missing"]
        percentage_variable = round(
            float((count_variable / count_total) * 100), 3
        )
        record["percentage_variable_non_missing"] = str(
            str(count_variable) + " (" +
            str(percentage_variable) + "%)"
        )
        record["count_variable_genotypes"] = int(array_genotype.size)
        count_genotype_variable = record["count_variable_genotypes"]
        percentage_genotype_variable = round(
            ((count_genotype_variable / count_total) * 100), 3
        )
        record["percentage_variable_genotypes"] = str(
            str(count_genotype_variable) + " (" +
            str(percentage_genotype_variable) + "%)"
        )
        if (count_variable > 10):
            # Determine mean, median, standard deviation, and standard error of
            # values in array.
            record["median"] = numpy.nanmedian(array)
            record["minimum"] = numpy.nanmin(array)
            record["maximum"] = numpy.nanmax(array)
            record["mean"] = numpy.nanmean(array)
            record["standard_deviation"] = numpy.nanstd(array)
            record["standard_error"] = scipy.stats.sem(array)
            # Confidence intervals and ranges.
            pail_confidence = (
                utility.determine_95_99_confidence_intervals_ranges(
                    estimate=record["mean"],
                    standard_error=record["standard_error"],
            ))
            record.update(pail_confidence)
            pass
        pass
    # Return information.
    return record


def drive_assemble_quantitation_table(
    variables=None,
    records_cohorts=None,
    report=None,
):
    """
    Drives the assembly of a description table from records of quantitative
    descriptive statistics on variables of interest.

    These descriptive statistics are most appropriate for continuous variables
    on interval, or ratio scales, but they can also be informative for discrete
    variables on ordinal scales.

    arguments:
        variables (list<str>): names of variables
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Collect summary records for rows within description table.
    records_description = list()
    # Iterate on cohorts.
    for record_cohort in records_cohorts:
        # Iterate on variables.
        for variable in variables:
            # Organize information for description record.
            record_description = create_quantitation_record(
                name_cohort=record_cohort["name"],
                variable=variable,
                table=record_cohort["table"],
            )
            # Collect records.
            records_description.append(record_description)
            pass
        pass
    # Organize table.
    table = pandas.DataFrame(data=records_description)
    # Select and sort relevant columns from table.
    columns = [
        "cohort",
        "variable",
        "count_cohort_total_records",
        "count_cohort_total_genotypes",
        "percentage_cohort_total_genotypes",
        "count_variable_non_missing",
        "percentage_variable_non_missing",
        "count_variable_genotypes",
        "percentage_variable_genotypes",
        "median",
        "minimum",
        "maximum",
        "mean",
        "standard_deviation",
        "standard_error",
        "range_95",
        "range_99",
        "interval_95",
        "range_95_low",
        "range_95_high",
        "interval_99",
        "range_99_low",
        "range_99_high",
    ]
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    table = table[[*columns]]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("drive_assemble_quantitation_table()")
        utility.print_terminal_partition(level=3)
        print(table)
        pass
    # Return information.
    return table


##########
# Confidence Intervals


def determine_range_text(
    minimum=None,
    maximum=None,
    delimiter=None,
):
    """
    Prepares a textual representation of a range.

    arguments:
        minimum (float): value of minimum
        maximum (float): value of maximum
        delimiter (str): character string to place in between values

    raises:

    returns:
        (str): textual representation of range
    """

    range = str(
        str(round((float(minimum)), 5)) +
        " ... " +
        str(round((float(maximum)), 5))
    )
    return range


def calculate_confidence_interval_range(
    confidence=None,
    standard_error=None,
    estimate=None,
):
    """
    Calculates a confidence interval from a standard error.

    This function assumes that the confidence interval is symmetrical about the
    central estimate.

    https://www.mathsisfun.com/data/confidence-interval.html
    80% Confidence Interval: (1.282 * standard_error)
    85% Confidence Interval: (1.440 * standard_error)
    90% Confidence Interval: (1.645 * standard_error)
    95% Confidence Interval: (1.960 * standard_error)
    99% Confidence Interval: (2.576 * standard_error)
    99.5% Confidence Interval: (2.807 * standard_error)
    99.9% Confidence Interval: (3.291 * standard_error)

    arguments:
        confidence (float): proportion of confidence, 0.80, 0.85, 0.90, 0.95,
            0.99, 0.995, 0.999
        standard_error (float): value of standard error
        estimate (float): value of estimate

    raises:

    returns:
        (dict<float>): values of confidence interval, minimum, and maximum
    """

    # Copy information.
    confidence = copy.deepcopy(confidence)
    standard_error = copy.deepcopy(standard_error)
    estimate = copy.deepcopy(estimate)
    # Determine factor for specific proportion of confidence.
    confidence = float(confidence)
    if (confidence == 0.80):
        factor = float(1.282)
    elif (confidence == 0.85):
        factor = float(1.440)
    elif (confidence == 0.90):
        factor = float(1.645)
    elif (confidence == 0.95):
        factor = float(1.960)
    elif (confidence == 0.99):
        factor = float(2.576)
    elif (confidence == 0.995):
        factor = float(2.807)
    elif (confidence == 0.999):
        factor = float(3.291)
    else:
        print("calculate_confidence_interval_from_standard_error")
        print(
            "Error: " +
            "Function did not recognize confidence as a float proportion."
        )
        factor = float("nan")
        pass
    # Collect information.
    pail = dict()
    # Calculate confidence interval.
    pail["interval"] = float(factor * standard_error)
    # Calculate confidence_range.
    pail["minimum"] = float(float(estimate) - float(pail["interval"]))
    pail["maximum"] = float(float(estimate) + float(pail["interval"]))
    # Return information.
    return pail


def determine_95_99_confidence_intervals_ranges(
    estimate=None,
    standard_error=None,
):
    """
    Determines information about 95% and 99% confidence intervals and their
        ranges

    arguments:
        estimate (float): value of estimate
        standard_error (float): value of standard error

    raises:

    returns:
        (dict<float, str>): information about 95% and 99% confidence intervals and
            their ranges
    """

    # Collect information.
    pail = dict()
    # Calculate confidence intervals and ranges.
    pail_95 = calculate_confidence_interval_range(
        confidence=0.95,
        standard_error=standard_error,
        estimate=estimate,
    )
    pail_99 = calculate_confidence_interval_range(
        confidence=0.99,
        standard_error=standard_error,
        estimate=estimate,
    )
    # Organize information.
    pail["interval_95"] = pail_95["interval"]
    pail["interval_99"] = pail_99["interval"]
    pail["range_95_low"] = pail_95["minimum"]
    pail["range_95_high"] = pail_95["maximum"]
    pail["range_99_low"] = pail_99["minimum"]
    pail["range_99_high"] = pail_99["maximum"]
    pail["range_95"] = determine_range_text(
        minimum=pail_95["minimum"],
        maximum=pail_95["maximum"],
        delimiter=" ... ",
    )
    pail["range_99"] = determine_range_text(
        minimum=pail_99["minimum"],
        maximum=pail_99["maximum"],
        delimiter=" ... ",
    )
    # Return information.
    return pail


##########
# Write


def write_product_table(
    name=None,
    table=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        name (str): base name for file
        table (object): table of information to write to file
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_directory, str(name + ".pickle")
    )
    path_table_text = os.path.join(
        path_directory, str(name + ".tsv")
    )
    # Write information to file.
    table.to_pickle(path_table)
    table.to_csv(
        path_or_buf=path_table_text,
        sep="\t",
        header=True,
        index=True,
    )
    pass


def write_product_tables(
    pail_write=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict<object>): collection of information to write to file,
            with keys defining names for Pandas data-frame tables as entries
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    for name in pail_write.keys():
        write_product_table(
            name=name,
            table=pail_write[name],
            path_directory=path_directory,
        )
    pass



################################################################################
# Procedure
# Currently, this module is not executable.
