"""
Supply functionality for description of variables and preparation of tables and
other reports.

This module 'description' is part of the 'partner' package.

This module is not directly executable.

This subpackage 'partner' provides executable functionality under the
management of a higher level package. Importation paths require this hierarchy.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This module file is part of the project package directory 'partner'
    (https://github.com/tcameronwaller/partner/).

    Project 'partner' supports data analysis in multiple other projects.
    Copyright (C) 2024 Thomas Cameron Waller

    The code within project 'partner' is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the GNU
    General Public License, or (at your option) any later version.

    The code within project 'partner' is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'partner'. If not, see <http://www.gnu.org/licenses/>.
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
import partner.utility as putly # this import path for subpackage
import partner.organization as porg # this import path for subpackage
import partner.scale as pscl

#dir()
#importlib.reload()

################################################################################
# Functionality


##########
# Table Attribution


def create_attribution_record(
    cohort_name=None,
    name_variable_value=None,
    variable=None,
    value=None,
    table=None,
):
    """
    Organize a record (single row in table) to describe attribution of
    categorical or discrete variable values across cohorts.

    arguments:
        cohort_name (str): name of cohort
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
    record["description_cohort_name"] = str(cohort_name)
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

    The description records and the description table that this function
    assembles preserve all information (variable names and their values) within
    each cohort record (Python dictionaries). Use variables within each cohort
    record to define details of each stratification cohort such as type of
    data records available (phenotypes, genotypes), sex (any, female, male),
    ancestry or race or ethnicity, stage of life (young, middle, old,
    premenopause, perimenopause, postmenopause), and special exclusions. At a
    minimum, the cohort record needs a variable named "cohort_name".

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
                cohort_name=record_cohort["cohort_name"],
                name_variable_value=record_attribution["name"],
                variable=record_attribution["variable"],
                value=record_attribution["value"],
                table=record_cohort["table"],
            )
            # Preserve information from stratification cohort record.
            record_description.update(record_cohort)
            del record_description["table"]
            # Collect records.
            records_description.append(record_description)
            pass
        pass
    # Organize table.
    table = pandas.DataFrame(data=records_description)
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("report: ")
        print("drive_assemble_attribution_table()")
        putly.print_terminal_partition(level=3)
        print(table)
        pass
    # Return information.
    return table


##########
# Table Quantitation


def create_quantitation_record(
    cohort_name=None,
    variable=None,
    variable_attribution=None,
    value_attribution=None,
    table=None,
):
    """
    Organize a record (single row in table) to describe for measures on
    quantitative variables across cohorts.

    Report percentages relative to the total count of records in the cohort.

    arguments:
        cohort_name (str): name of cohort
        variable (str): name of table's column for variable
        variable_attribution (str): name of table's column for a special
            variable for which to report attribution against main variable
        value_attribution (object): value of special variable for attribution
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (dict): information for summary table record on cohort

    """

    # Collect information for record.
    record = dict()
    record["description_cohort_name"] = str(cohort_name)
    record["variable"] = str(variable)
    record["variable_attribution"] = str(variable_attribution)
    record["value_attribution"] = value_attribution
    # Copy information.
    table = table.copy(deep=True)

    # Stratify table.
    # Select relevant rows of the table.
    table_attribution = table.loc[
        (
            (~pandas.isna(table[variable_attribution])) &
            (table[variable_attribution] == value_attribution)
        ), :
    ]

    # Counts.
    record["count_cohort_total_records"] = int(table.shape[0])
    record["count_cohort_total_attribution"] = int(table_attribution.shape[0])
    count_total = record["count_cohort_total_records"]
    count_total_attribution = record["count_cohort_total_attribution"]
    # Percentages.
    percentage_total_attribution = round(
        ((count_total_attribution / count_total) * 100), 3
    )
    record["percentage_cohort_total_attribution"] = str(
        str(record["count_cohort_total_attribution"]) + " (" +
        str(percentage_total_attribution) + "%)"
    )
    # Initialize missing values.
    record["count_variable_non_missing"] = float("nan")
    record["percentage_variable_non_missing"] = str("nan (nan%)")
    record["count_variable_attribution"] = float("nan")
    record["percentage_variable_attribution"] = str("nan (nan%)")
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
        array_attribution = copy.deepcopy(
            table_attribution[variable].dropna().to_numpy()
        ) # non-missing variable values with genotypes
        # Determine count of valid values.
        record["count_variable_non_missing"] = int(array.size)
        count_variable = record["count_variable_non_missing"]
        record["count_variable_attribution"] = int(array_attribution.size)
        count_attribution_variable = record["count_variable_attribution"]
        # Percentages.
        percentage_variable = round(
            float((count_variable / count_total) * 100), 3
        )
        record["percentage_variable_non_missing"] = str(
            str(count_variable) + " (" +
            str(percentage_variable) + "%)"
        )
        percentage_attribution_variable = round(
            ((count_attribution_variable / count_total) * 100), 3
        )
        record["percentage_variable_attribution"] = str(
            str(count_attribution_variable) + " (" +
            str(percentage_attribution_variable) + "%)"
        )
        if (count_variable > 5):
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
                determine_95_99_confidence_intervals_ranges(
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
    variable_attribution=None,
    value_attribution=None,
    records_cohorts=None,
    report=None,
):
    """
    Drives the assembly of a description table from records of quantitative
    descriptive statistics on variables of interest.

    These descriptive statistics are most appropriate for continuous variables
    on interval, or ratio scales, but they can also be informative for discrete
    variables on ordinal scales.

    The description records and the description table that this function
    assembles preserve all information (variable names and their values) within
    each cohort record (Python dictionaries). Use variables within each cohort
    record to define details of each stratification cohort such as type of
    data records available (phenotypes, genotypes), sex (any, female, male),
    ancestry or race or ethnicity, stage of life (young, middle, old,
    premenopause, perimenopause, postmenopause), and special exclusions. At a
    minimum, the cohort record needs a variable named "cohort_name".

    arguments:
        variables (list<str>): names of variables
        variable_attribution (str): name of table's column for a special
            variable for which to report attribution against main variable
        value_attribution (object): value of special variable for attribution
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
                cohort_name=record_cohort["cohort_name"],
                variable=variable,
                variable_attribution=variable_attribution,
                value_attribution=value_attribution,
                table=record_cohort["table"],
            )
            # Preserve information from stratification cohort record.
            record_description.update(record_cohort)
            del record_description["table"]
            # Collect records.
            records_description.append(record_description)
            pass
        pass
    # Organize table.
    table = pandas.DataFrame(data=records_description)
    # Select and sort relevant columns from table.
    columns = [
        "description_cohort_name",
        "count_cohort_total_records",
        "variable_attribution",
        "value_attribution",
        "count_cohort_total_attribution",
        "percentage_cohort_total_attribution",
        "variable",
        "count_variable_non_missing",
        "percentage_variable_non_missing",
        "count_variable_attribution",
        "percentage_variable_attribution",
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
    #    "cohort_name",
    #    "cohort_phenotypes",
    #    "cohort_genotypes",
    #    "cohort_sex",
    #    "cohort_race",
    #    "cohort_ancestry",
    #    "cohort_life_stage",
    #    "cohort_exclusions",
    record_cohort_example = copy.deepcopy(records_cohorts[0])
    del record_cohort_example["table"]
    #columns.insert(0, list(record_cohort_example.keys()))
    columns.extend(list(record_cohort_example.keys()))
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    table = table[[*columns]]
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("report: ")
        print("drive_assemble_quantitation_table()")
        putly.print_terminal_partition(level=3)
        print(table)
        pass
    # Return information.
    return table


##########
# Variance


def calculate_variance_measures(
    array=None,
):
    """
    Calculates measures of variance.

    Variance is mathematically never negative.

    A synonymous common name for the 'coefficient of variation' is the 'relative
    standard deviation (RSD)'.

    Review: TCW; 2 April 2025

    arguments:
        array (object): NumPy array of values on a quantitative, continuous
           scale of measurement, such as ratio or interval

    raises:

    returns:
        (float): value of relative variance

    """

    # Copy information.
    values_raw = copy.deepcopy(array)
    # Filter to nonmissing values.
    values_nonmissing = numpy.copy(
        values_raw[~numpy.isnan(values_raw)]
    )
    # Calculate descriptive statistical measures.
    mean = numpy.nanmean(values_nonmissing)
    mean_absolute = math.fabs(mean)
    variance = numpy.nanvar(
            values_nonmissing,
            ddof=1, # divisor is (n - 1) for sample variance
    )
    standard_deviation = numpy.nanstd(
            values_nonmissing,
            ddof=1, # divisor is (n - 1) for sample standard deviation
    )
    # ((mean > 0.0) or (mean < 0.0))
    if (mean != 0.0):
        #scipy.stats.variation()
        coefficient_variation = float(
            standard_deviation / mean
        )
    else:
        coefficient_variation = float("nan")
        pass
    # Collect information.
    pail = dict()
    pail["variance"] = variance
    pail["standard_deviation"] = standard_deviation
    pail["coefficient_variation"] = coefficient_variation
    # Return information.
    return pail



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
        str(round((float(minimum)), 4)) +
        str(delimiter) +
        str(round((float(maximum)), 4))
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
    pail["minimum"] = round(
        (float(float(estimate) - float(pail["interval"]))), 4
    )
    pail["maximum"] = round(
        (float(float(estimate) + float(pail["interval"]))), 4
    )
    # Return information.
    return pail


def determine_95_99_confidence_intervals_ranges(
    estimate=None,
    standard_error=None,
):
    """
    Determines information about 95% and 99% confidence intervals and their
        ranges

    95% Confidence Interval:
    qnorm(0.975) = 1.960

    99% Confidence Interval:
    qnorm(0.995) = 2.576

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
# Correlation


def calculate_correlations_table_columns_pair(
    table=None,
    column_primary=None,
    column_secondary=None,
    count_minimum_observations=None,
    report=None,
):
    """
    Calculates correlations between a pair of columns in a table.

    Review: 17 December 2024

    arguments:
        table (object): Pandas data-frame table of values on continuous
            interval or ratio scales of measurement corresponding to features
            across columns and observations across rows
        column_primary (str): name of column in table corresponding to first or
            primary feature for correlation
        column_secondary (str): name of second column in table corresponding to
            second feature for correlation
        count_minimum_observations (int): minimal count of observations that
            must have nonmissing values for the two features in the pair
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about correlations



    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Remove table rows with missing values in relevant columns.
    table.dropna(
        axis="index",
        how="any",
        subset=[column_primary, column_secondary],
        inplace=True,
    )

    # Collect information.
    pail = dict()
    # Define names for values of statistical measures.
    measures = list()
    measures.append("correlation_pearson")
    measures.append("probability_pearson")
    measures.append("confidence_95_low_pearson")
    measures.append("confidence_95_high_pearson")
    measures.append("correlation_spearman")
    measures.append("probability_spearman")
    #measures.append("confidence_95_low_spearman")
    #measures.append("confidence_95_high_spearman")
    # Define sequence of values.
    names = copy.deepcopy(measures)
    names.insert(0, "count_observations")
    pail["names"] = names

    # Determine count of observations for which pair of features has matching
    # values.
    pail["count_observations"] = copy.deepcopy(table.shape[0])
    # Calculate correlations.
    if (pail["count_observations"] >= count_minimum_observations):
        #if method == "pearson": ...
        # Calculate correlations.
        # Pearson.
        results_pearson = scipy.stats.pearsonr(
            table[column_primary].to_numpy(),
            table[column_secondary].to_numpy(),
            alternative="two-sided",
        )
        pail["correlation_pearson"] = results_pearson.statistic
        pail["probability_pearson"] = results_pearson.pvalue
        confidence_95_pearson = results_pearson.confidence_interval(
            confidence_level=0.95
        )
        pail["confidence_95_low_pearson"] = confidence_95_pearson.low
        pail["confidence_95_high_pearson"] = confidence_95_pearson.high
        # Spearman.
        results_spearman = scipy.stats.spearmanr(
            table[column_primary].to_numpy(),
            table[column_secondary].to_numpy(),
            nan_policy="omit",
            alternative="two-sided",
        )
        pail["correlation_spearman"] = results_spearman.statistic
        pail["probability_spearman"] = results_spearman.pvalue
        # Unfortunately the implementation of Spearman Correlation in SciPy
        # does not offer a convenient method to calculate the confidence
        # interval.
        # Kendall.
        #correlation_kendall, probability_kendall = scipy.stats.kendalltau(
        #    table[column_one].to_numpy(),
        #    table[column_two].to_numpy(),
        #    alternative="two-sided",
        #)
    else:
        for measure in measures:
            pail[measure] = float("nan")
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: description.py")
        function = str(
            "calculate_correlations_table_columns_pair()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail




##########
# Benjamini-Hochberg False-Discovery Rate (FDR) q-value


def calculate_table_false_discovery_rate_q_values(
    table=None,
    name_column_p_value=None,
    name_column_q_value=None,
    name_column_significance=None,
    threshold=None,
):
    """
    Calculates Benjamini-Hochberg q-values to indicate false discovery rates
    (FDRs) from original p-values.

    Implementations in R:
    p.adjust

    Implementations in Python:
    scipy.stats.false_discovery_control
    statsmodels.stats.multitest.multipletests
    statsmodels.stats.multitest.fdrcorrection
    statsmodels.stats.multitest.fdrcorrection_twostage

    References:
    Benjamini, Krieger, Yekutieli; Biometrika; 2006
    (https://doi.org/10.1093/biomet/93.3.491)
    Benjamini, Hochberg; Journal of the Royal Statistical Society; 1995
    site: https://www.r-bloggers.com/2023/07/
    the-benjamini-hochberg-procedure-fdr-and-p-value-adjusted-explained/
    site: https://www.statisticshowto.com/benjamini-hochberg-procedure/

    arguments:
        table (object): Pandas data frame with column of probabilities across
            observations in rows
        name_column_p_value (str): name of table's column of p-values
        name_column_q_value (str): name for table's column of q-values
        name_column_significance (str): name for table's column of FDR
            indication that null hypothesis can be rejected
        threshold (float): value of alpha, or family-wise error rate in
            calculation of false discovery rate

    raises:

    returns:
        (object): Pandas data frame of original p-values and novel q-values to
            indicate false discovery rates

    """

    # Copy information.
    table = table.copy(deep=True)

    # False name_column_q_value rate method cannot accommodate missing values.
    # Remove null values.
    table_null_boolean = pandas.isna(table[name_column_p_value])
    table_null = table.loc[table_null_boolean]
    table_valid = table.dropna(
        axis="index",
        how="any",
        subset=[name_column_p_value],
        inplace=False,
    )
    # Calculate false name_column_q_value rates from probabilities.
    p_values = copy.deepcopy(table_valid[name_column_p_value].to_numpy())
    if len(p_values) > 3:
        # Benjamini, Hochberg; 1995
        #report = statsmodels.stats.multitest.multipletests(
        #    p_values,
        #    alpha=threshold, # family-wise error rate
        #    method="fdr_bh", # fdr_bh use Benjamini-Hochberg False-Discovery Rate (FDR)
        #    is_sorted=False,
        #    #return_sorted=False, # keep q-values in original sequence
        #)
        report = statsmodels.stats.multitest.fdrcorrection(
            p_values,
            alpha=threshold, # family-wise error rate
            method="indep", # independent tests
            is_sorted=False,
        ) # Benjamini, Hochberg, 1995
        # Benjamini, Krieger, Yekutieli; 2006
        #report = statsmodels.stats.multitest.fdrcorrection_twostage(
        #    p_values,
        #    alpha=threshold, # family-wise error rate
        #    method="bky", # bky: Benjamini, Krieger, Yekuteli Definition 6
        #    maxiter=1, # 1: two-stage method
        #    is_sorted=False,
        #) # Benjamini, Krieger, Yekutieli; 2006
        significances = report[0] # whether valid to reject null hypothesis
        #significances = numpy.invert(rejects)
        q_values = report[1]
        table_valid[name_column_significance] = significances
        table_valid[name_column_q_value] = q_values
    else:
        table_valid[name_column_significance] = float("nan")
        table_valid[name_column_q_value] = float("nan")
        pass
    table_null[name_column_significance] = False
    table_null[name_column_q_value] = float("nan")

    # Combine null and valid portions of data.
    table_discoveries = pandas.concat(
        [table_valid, table_null],
        ignore_index=False,
    )
    # Return information.
    return table_discoveries


def calculate_q_values_in_table_groups(
    factors=None,
    threshold=None,
    name_column_p_value=None,
    name_column_q_value=None,
    name_column_significance=None,
    table=None,
    report=None,
):
    """
    Calculates Benjamini-Hochberg q-values to indicate false discovery rates
    (FDRs) from original p-values within stratification groups.

    arguments:
        factors (list<str>): names of columns in table by which to split groups
        threshold (float): value of alpha, or family-wise error rate in
            calculation of false discovery rate
        name_column_p_value (str): name of table's column of p-values
        name_column_q_value (str): name for table's column of q-values
        name_column_significance (str): name for table's column of FDR
            indication that null hypothesis can be rejected
        table (object): Pandas data frame with column of probabilities across
            observations in rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of original p-values and novel q-values to
            indicate false discovery rates

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # do not remove index; move to regular columns
    )
    #table.drop_duplicates(
    #    subset=None,
    #    keep="first",
    #    inplace=True,
    #)
    table.set_index(
        factors,
        append=False,
        drop=True,
        inplace=True
    )
    # Split rows within table by factor columns.
    groups = table.groupby(
        level=factors,
    )
    # Initiate table for collection of group tables after any transformations.
    table_collection = pandas.DataFrame()
    # Iterate on groups.
    for name, table_group in groups:
        # Copy information in table.
        table_group = table_group.copy(deep=True)
        # Organize table.
        table_group.reset_index(
            level=None,
            inplace=True,
            drop=False, # do not remove index; move to regular columns
        )
        # Complete further organization on each group table after split,
        # such as setting new index.
        # Report.
        if report:
            print_terminal_partition(level=4)
            print("Name of group:")
            print(name)
            print("Table for group after split:")
            print(table_group)
            print_terminal_partition(level=4)
        # Complete procedures on each group table after split.
        # For example, calculate summary statistics on each group and then
        # collect within a new summary table.

        # Transformations to group tables.
        table_q = calculate_table_false_discovery_rate_q_values(
            threshold=threshold, # alpha; family-wise error rate
            name_column_p_value=name_column_p_value,
            name_column_q_value=name_column_q_value,
            name_column_significance=name_column_significance,
            table=table_group,
        )
        # Collect group tables.
        table_collection = pandas.concat(
            [table_collection, table_q],
            ignore_index=False,
        )
        pass
    # Return information.
    return table_collection


def calculate_table_long_false_discovery_rate_q_values(
    table=None,
    name_index_type_value=None,
    type_p_value=None,
    threshold=None,
    names_indices_rows=None,
    report=None,
):
    """
    Beginning with a table of values in long format, this function calculates
    Benjamini-Hochberg q-values to control the False-Discovery Rate (FDR) for
    testing multiple hypotheses.

    Example of source table in long format:

    group_1st     group_2nd     type_value     value
    group_1_a     group_2_a     signal    -0.15
    group_1_a     group_2_a     p_value   0.001
    group_1_a     group_2_b     signal    0.25
    group_1_a     group_2_b     p_value   0.001
    group_1_a     group_2_c     signal    0.35
    group_1_a     group_2_c     p_value   0.001
    group_1_b     group_2_a     signal    0.75
    group_1_b     group_2_a     p_value   0.001
    group_1_b     group_2_b     signal    -0.35
    group_1_b     group_2_b     p_value   0.001
    group_1_b     group_2_c     signal    0.45
    group_1_b     group_2_c     p_value   0.001

    Example of product table in long format:

    group_1st     group_2nd     type_value     value
    group_1_a     group_2_a     signal    -0.15
    group_1_a     group_2_a     p_value   0.001
    group_1_a     group_2_a     q_value   0.01
    group_1_a     group_2_b     signal    0.25
    group_1_a     group_2_b     p_value   0.001
    group_1_a     group_2_b     q_value   0.01
    group_1_a     group_2_c     signal    0.35
    group_1_a     group_2_c     p_value   0.001
    group_1_a     group_2_c     q_value   0.01
    group_1_b     group_2_a     signal    0.75
    group_1_b     group_2_a     p_value   0.001
    group_1_b     group_2_a     q_value   0.01
    group_1_b     group_2_b     signal    -0.35
    group_1_b     group_2_b     p_value   0.001
    group_1_b     group_2_b     q_value   0.01
    group_1_b     group_2_c     signal    0.45
    group_1_b     group_2_c     p_value   0.001
    group_1_b     group_2_c     q_value   0.01

    Recommendations for names of indices:
    name_index_type_value="type_value"

    Recommendation for categorical name of p-values within index level
    "type_value":
    type_p_value="p_value"

    Review: TCW; 6 June 2024

    arguments:
        table (object): Pandas data-frame table in long format with definitions
            of indices across columns and rows
        name_index_type_value (str): name of index across rows that designates
            which of multiple types of values correspond to categorical indices
            across columns and rows
        type_p_value (str): name of category for p-values
        threshold (float): value of alpha, or family-wise error rate in
            calculation of false discovery rate
        names_indices_rows (list<str>): names of columns for indices across
            rows in original source table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information in table.
    table_source = table.copy(deep=True)

    ##########
    # Separate the p-values from other types of values within the table.
    #table_p = table_source.loc[
    #    (table_source[name_index_type_value] == type_p_value), :
    #].copy(deep=True)
    table_p = table_source[
        table_source.index.get_level_values(
            name_index_type_value
        ).isin([type_p_value])
    ].copy(deep=True)
    #matrix_p = numpy.copy(table_p.to_numpy())

    ##########
    # Translate names of columns.
    translations = dict()
    translations["value"] = "p_value"
    table_p.rename(
        columns=translations,
        inplace=True,
    )

    ##########
    # Calculate Benjamini-Hochberg False Discovery Rate q-values.
    table_q = calculate_table_false_discovery_rate_q_values(
        table=table_p,
        name_column_p_value="p_value",
        name_column_q_value="q_value",
        name_column_significance="q_significance",
        threshold=threshold, # alpha; family-wise error rate
    )

    ##########
    # Organize table of q-values.
    # Remove unnecessary columns.
    table_q.drop(
        labels=["p_value", "q_significance",],
        axis="columns",
        inplace=True
    )
    # Organize information in table.
    table_q.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_q[name_index_type_value] = "q_value"
    # Translate names of columns.
    translations = dict()
    translations["q_value"] = "value"
    table_q.rename(
        columns=translations,
        inplace=True,
    )

    ##########
    # Combine the novel q-values with the p-values and other information in the
    # original table.
    # Organize information in table.
    table_source.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Concatenate tables.
    table_product = pandas.concat(
        [table_source, table_q,],
        axis="index",
        join="outer",
        ignore_index=False,
        copy=True,
    )
    # Organize indices in table.
    table_product.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_product.set_index(
        names_indices_rows,
        append=False,
        drop=True,
        inplace=True,
    )

    ##########
    # Return information.
    return table_product


##########
# Quantiles


def describe_quantiles_ordinal(
    table=None,
    column_source=None,
    column_product=None,
    columns_category=None,
    report=None,
):
    """
    Describe an integer ordinal encoding of quantiles for a feature variable
    with values on a continuous interval or ratio measurement scale.

    Review: TCW; 27 August 2024

    arguments:
        table (object): Pandas data-frame table of subjects, samples, and their
            attribute features
        column_source (str): name of source column in table for a feature
            variable with values on a continuous interval or ratio measurement
            scale
        column_product (str): name of product column in table for designation
            of quantiles with encoding as integers on an ordinal scale
        columns_category (list<str>): names of columns in table for categorical
            factor variables to describe within quantiles
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Extract unique values from column of designations for quantiles.
    quantiles = table[column_product].unique()
    #sorted(quantiles, reverse=True) # error with string designators
    count = len(quantiles)
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.description.py")
        print("function: describe_quantiles_ordinal()")
        putly.print_terminal_partition(level=4)
        print(str("name of source column: " + column_source))
        putly.print_terminal_partition(level=5)
        print(str("count of quantiles: " + str(count)))
        putly.print_terminal_partition(level=5)
        print("quantile designators:")
        print(quantiles)
        putly.print_terminal_partition(level=5)
        print("count of samples in each quantile:")
        print(table[column_product].value_counts(dropna=False))
        #putly.print_terminal_partition(level=5)
        #print("description of quantiles:")
        #print(table[column_product].describe(include=columns_category))
        putly.print_terminal_partition(level=5)
        print("what follows is a description of each quantile separately")
        pass
    # Describe individual quantiles.
    for designator in quantiles:
        table_quantile = table.loc[
            (table[column_product] == designator), :
        ]
        count_rows = (table_quantile.shape[0])
        count_columns = (table_quantile.shape[1])
        minimum = table_quantile[column_source].min()
        maximum = table_quantile[column_source].max()
        # Report.
        if report:
            putly.print_terminal_partition(level=4)
            print(str("quantile: " + str(designator)))
            print(str("count rows: " + str(count_rows)))
            print(str("minimum: " + str(minimum)))
            print(str("maximum: " + str(maximum)))
            pass
        for column_category in columns_category:
            # Report.
            if report:
                putly.print_terminal_partition(level=5)
                print(str("categorical column: " + column_category))
                print(
                    table_quantile[column_category].value_counts(dropna=False)
                )
                pass
            pass
        pass
    # Return information.
    pass


def determine_describe_quantiles_ordinal(
    table=None,
    column_source=None,
    column_product=None,
    count=None,
    text_string=None,
    name_prefix=None,
    report=None,
):
    """
    Determine and describe an integer ordinal encoding of quantiles for a
    feature variable with values on a continuous interval or ratio measurement
    scale.

    See the function below as an example of applying quantiles separately by
    some factor of stratification.
    function:
    define_ordinal_stratifications_by_sex_continuous_variables
    module:
    organization.py
    package:
    sexy_age_hormones

    Review: TCW; 27 August 2024

    arguments:
        table (object): Pandas data-frame table of features across columns and
           observations across rows
        column_source (str): name of source column in table for a feature
            variable with values on a continuous interval or ratio measurement
            scale
        column_product (str): name of product column in table for designation
            of quantiles with encoding as integers on an ordinal scale
        count (int): count of quantiles to create and describe
        text_string (bool): whether to convert the integer, ordinal designators
            of quantiles to text strings
        name_prefix (str): prefix for names as designators of the tertiles
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Determine labels for quantiles.
    designators_raw = list(range(0, count, 1))
    if text_string:
        designators = list(map(
            lambda designator: str(name_prefix + str(designator)),
            designators_raw
        ))
    else:
        designators = designators_raw

    if True:
        # Determine quantiles.
        table[column_product], bins = pandas.qcut(
            table[column_source],
            q=count,
            labels=designators,
            retbins=True,
            duplicates="raise", # "raise", "drop"
        )
    else:
        # Sort rows within table.
        table.sort_values(
            by=[column_source,],
            axis="index",
            ascending=True,
            na_position="last",
            inplace=True,
        )
        # Assign ranks to the rows in the table.
        # Unambiguous ranks on the sorted rows in the data-frame table can
        # help to avoid overlapping bins for the quantiles.
        # Determine quantiles.
        table[column_product], bins = pandas.qcut(
            table.rank(method="first"),
            q=count,
            labels=designators,
            retbins=True,
        )
        pass

    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["bins"] = bins
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.description.py")
        print("function: determine_describe_quantiles_ordinal()")
        putly.print_terminal_partition(level=4)
        print(str("name of source column: " + column_source))
        putly.print_terminal_partition(level=5)
        print(str("count of quantiles: " + str(count)))
        putly.print_terminal_partition(level=5)
        print("quantile bins:")
        print(bins)
        describe_quantiles_ordinal(
            table=table,
            column_source=column_source,
            column_product=column_product,
            columns_category=[],
            report=False,
        )
        pass
    # Return information.
    return pail


##########
# Report ranges of variables from table


def report_table_range_variables(
    table=None,
    variables=None,
):
    """
    Reports the ranges of quantitative, continuous variables in table.

    arguments:
        table (object): Pandas data-frame table
        variables (list<str>): name of column for primary identifier

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Report.
    putly.print_terminal_partition(level=5)
    # Transfer attributes to the table of genetic correlations.
    for variable in variables:
        array = copy.deepcopy(table[variable].dropna().to_numpy()) # non-missing
        count_values = int(array.size)
        if (count_values > 1):
            minimum = round(float(numpy.nanmin(array)), 9)
            maximum = round(float(numpy.nanmax(array)), 9)
        elif (count_values == 1):
            minimum = round(float(array[0]), 9)
            maximum = round(float(array[0]), 9)
        # Report.
        putly.print_terminal_partition(level=5)
        print("Variable: " + str(variable))
        print("Count non-missing values: " + str(count_values))
        print("Minimum: " + str(minimum))
        print("Maximum: " + str(maximum))
        pass
    pass


##########
# Split and apply operation or transformation to groups within table


def template_split_apply_table_groups(
    factors=None,
    table=None,
    report=None,
):
    """
    Splits rows within table by groups of factor columns.
    Applies procedure to to each group of rows.

    The function serves as a sort of template.

    arguments:
        factors (list<str>): names of columns in table by which to split groups
        table (object): Pandas data-frame table of columns for feature variables
            across rows for observation records
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # do not remove index; move to regular columns
    )
    table.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    table.set_index(
        factors,
        append=False,
        drop=True,
        inplace=True
    )
    # Split rows within table by factor columns.
    groups = table.groupby(
        level=factors,
    )
    # Initiate table for collection of group tables after any transformations.
    table_collection = pandas.DataFrame()
    # Iterate on groups.
    for name, table_group in groups:
        # Copy information in table.
        table_group = table_group.copy(deep=True)
        # Organize table.
        table_group.reset_index(
            level=None,
            inplace=True,
            drop=False, # do not remove index; move to regular columns
        )
        # Complete further organization on each group table after split,
        # such as setting new index.
        # Report.
        if report:
            print_terminal_partition(level=4)
            print("Name of group:")
            print(name)
            print("Table for group after split:")
            print(table_group)
            print_terminal_partition(level=4)
        # Complete procedures on each group table after split.
        # For example, calculate summary statistics on each group and then
        # collect within a new summary table.

        # Transformations to group tables.
        # Do something...

        # Collect group tables.
        table_collection = pandas.concat(
            [table_collection, table_group],
            ignore_index=False,
        )
        pass
    # Return information.
    return table_collection


def define_sequence_columns_describe_table_features_by_groups():
    """
    Defines names of columns in sequence.

    arguments:

    raises:

    returns:
        (list<str>): names of columns in sequence by which to filter and sort
            columns in table

    """

    # Specify sequence of columns within table.
    columns_sequence = [
        #index_feature,
        #"feature_translation",
        #column_group,
        #"count_observations",
        #"count_observations_nonmissing",
        #"percentage_observations_nonmissing",
        "minimum",
        "maximum",
        "mean",
        "median",
        "standard_error",
        "standard_deviation",
        "coefficient_variation",
        "interquartile",
        "range_95_low",
        "range_95_high",
        "range_99_low",
        "range_99_high",
    ]
    # Return information.
    return columns_sequence


def describe_table_features_by_groups(
    table=None,
    column_group=None,
    columns_features=None,
    index_feature=None,
    translations_feature=None,
    threshold_observations=None,
    digits_round=None,
    report=None,
):
    """
    Describe values corresponding to columns for features within groups of rows
    for observations.

    ----------
    Format of source table (name: "table_source")
    ----------
    Format of source table is in wide format with features across columns and
    observations across rows. A special column gives identifiers corresponding
    to each observation across rows. Another special column provides names
    of categorical groups of observations, with multiple observations in each
    group. For versatility, this table does not have explicitly defined indices
    across columns or rows. Values for observations of features are on a
    quantitative, continuous, interval or ratio scale of measurement.
    ----------
    observation     group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_1   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_2   0.001     0.001     0.001     0.001     0.001
    observation_6   group_2   0.001     0.001     0.001     0.001     0.001
    ----------

    ----------
    Format of product table (name: "table_product")
    ----------
    Format of product table is in partial long format with summary statistics
    and measures across columns and features across rows. A special column
    gives identifiers corresponding to each feature across rows. Another
    special column provides names of categorical groups of observations. For
    versatility, this table does not have explicity defined indices across
    columns or rows.
    ----------
    detail    group   mean standard_error standard_deviation median interqua...
    feature
    feature_1 group_1 0.01 0.001          0.001              0.015  0.5
    feature_1 group_2 0.01 0.001          0.001              0.015  0.5
    feature_1 group_3 0.01 0.001          0.001              0.015  0.5
    feature_1 group_4 0.01 0.001          0.001              0.015  0.5
    feature_2 group_1 0.01 0.001          0.001              0.015  0.5
    feature_2 group_2 0.01 0.001          0.001              0.015  0.5
    feature_2 group_3 0.01 0.001          0.001              0.015  0.5
    feature_2 group_4 0.01 0.001          0.001              0.015  0.5
    ----------

    Review; TCW; 30 January 2025

    arguments:
        table (object): Pandas data-frame table of features across columns and
            observations across rows with values on a quantitative, continuous,
            interval or ratio scale of measurement
        column_group (str): name of column in source table for groups of
            observations
        column_features (list<str>): names of columns in original source table
            for features on a quantitative, continuous, interval, or ratio
            scale of measurement
        index_feature (str): name for column of index across rows in product
            table corresponding to names of features
        translations_feature (dict<str>): translations for names of features
        threshold_observations (int): minimal count of observations for which
            to calculate summary statistics
        digits_round (int): count of digits to right of decimal to which to
            round measures
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information.
    table_source = table.copy(deep=True)
    # Copy other information.
    columns_features = copy.deepcopy(columns_features)
    translations_feature = copy.deepcopy(translations_feature)

    # Collect records of information, which will become rows in table.
    records = list()
    # Iterate on features, apply operations, and collect information.
    for column_feature in columns_features:
        # Extract values from table for current feature across groups of
        # observations.
        pail_extract = (
            porg.extract_array_values_from_table_column_by_groups_rows(
                table=table_source,
                column_group=column_group,
                column_feature=column_feature,
                report=False,
        ))
        # Iterate on groups.
        for group in pail_extract["groups_sequence"]:
            # Collect information.
            record = dict()
            record[index_feature] = column_feature
            if (
                (translations_feature is not None) and
                (column_feature in translations_feature.keys())
            ):
                record["feature_translation"] = (
                    translations_feature[column_feature]
                )
            else:
                record["feature_translation"] = column_feature
                pass
            record[column_group] = group
            # Count.
            record["count_observations"] = int(len(
                pail_extract["groups_values"][group]
            ))
            record["count_observations_nonmissing"] = int(len(
                pail_extract["groups_values_nonmissing"][group]
            ))
            # Percentage.
            portion = record["count_observations_nonmissing"]
            total = record["count_observations"]
            if (total > 0):
                record["percentage_observations_nonmissing"] = round(float(
                    ((portion / total) * 100)), 2
                )
            else:
                record["percentage_observations_nonmissing"] = float("nan")
                pass
            # Determine whether there are adequate observations.
            count = record["count_observations_nonmissing"]
            if (count >= threshold_observations):
                # Range.
                record["minimum"] = round(numpy.nanmin(
                    pail_extract["groups_values"][group]
                ), digits_round)
                record["maximum"] = round(numpy.nanmax(
                    pail_extract["groups_values"][group]
                ), digits_round)
                # Center.
                record["mean"] = round(numpy.nanmean(
                    pail_extract["groups_values"][group]
                ), digits_round)
                record["median"] = round(numpy.nanmedian(
                    pail_extract["groups_values"][group]
                ), digits_round)
                # Spread, variance.
                record["standard_error"] = round(scipy.stats.sem(
                    pail_extract["groups_values"][group],
                    ddof=1, # divisor is (n - 1) for sample standard deviation
                    nan_policy="omit", # ignore missing values in calculation
                ), digits_round)
                record["standard_deviation"] = round(numpy.nanstd(
                    pail_extract["groups_values"][group],
                    ddof=1, # divisor is (n - 1) for sample standard deviation
                ), digits_round)
                if ((record["mean"] > 0.0) or (record["mean"] < 0.0)):
                    record["coefficient_variation"] = round(float(
                        record["standard_deviation"] / record["mean"]
                    ), digits_round)
                else:
                    record["coefficient_variation"] = float("nan")
                    pass
                record["interquartile"] = round(scipy.stats.iqr(
                    pail_extract["groups_values"][group],
                    rng=(25, 75),
                    nan_policy="omit", # ignore missing values in calculation.
                ), digits_round)
                pail_95 = calculate_confidence_interval_range(
                    confidence=0.95,
                    standard_error=record["standard_error"],
                    estimate=record["mean"],
                )
                record["range_95_low"] = round(
                    pail_95["minimum"], digits_round
                )
                record["range_95_high"] = round(
                    pail_95["maximum"], digits_round
                )
                pail_99 = calculate_confidence_interval_range(
                    confidence=0.99,
                    standard_error=record["standard_error"],
                    estimate=record["mean"],
                )
                record["range_99_low"] = round(
                    pail_99["minimum"], digits_round
                )
                record["range_99_high"] = round(
                    pail_99["maximum"], digits_round
                )
            else:
                # Fill missing values.
                columns_sequence = (
                    define_sequence_columns_describe_table_features_by_groups()
                )
                for column in columns_sequence:
                    record[column] = float("nan")
                    pass
                pass
            # Collect information.
            records.append(record)
            pass
        pass
    # Organize information in a table.
    table_product = pandas.DataFrame(data=records)
    # Specify sequence of columns within table.
    columns_sequence = [
        index_feature,
        "feature_translation",
        column_group,
        "count_observations",
        "count_observations_nonmissing",
        "percentage_observations_nonmissing",
    ]
    columns_sequence_extra = (
        define_sequence_columns_describe_table_features_by_groups()
    )
    columns_sequence.extend(columns_sequence_extra)
    # Filter and sort columns within table.
    table_product = porg.filter_sort_table_columns(
        table=table_product,
        columns_sequence=columns_sequence,
        report=False,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.description.py")
        function = str(
            "describe_table_features_by_groups" +
            "()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_product



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



#
