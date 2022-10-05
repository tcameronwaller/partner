"""
Supply functionality for transformation of distribution scales.

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


def transform_values_distribution_scale_logarithm(
    values_array=None,
    shift_minimum=None,
    base=None,
):
    """
    Shifts values within a NumPy array to positive scale and then transforms
    their distribution and scale by Logarithm at a specific base.

    Logarithm does not have a definition for negative values.

    arguments:
        values_array (object): NumPy array of original values
        shift_minimum (float): scalar value for new minimum greater than zero
            after translation
        base (float): value for logarithmic base

    raises:

    returns:
        (object): NumPy array of novel values after transformation

    """

    # Copy information.
    values_array = numpy.copy(values_array)
    # Calculate shift translation for new minimum.
    # Some contexts call this shift a pseudo count.
    minimum = numpy.nanmin(values_array)
    shift = (shift_minimum - minimum)
    # Shift values to be positive and greater than zero.
    values_array = numpy.add(values_array, shift)
    # Calculate logarithm.
    if (math.isclose(float(base), math.e)):
        values_log_array = numpy.log(values_array)
    elif (int(base) == 2):
        values_log_array = numpy.log2(values_array)
    elif (int(base) == 10):
        values_log_array = numpy.log10(values_array)
    else:
        # math.log(value, base)
        values_log_array = (numpy.log(values_array) / numpy.log(base))
    # Return information.
    return values_log_array


def drive_transform_variables_distribution_scale_logarithm(
    columns=None,
    suffix=None,
    table=None,
):
    """
    Transforms distribution and scale of variables' values by Natural Logarithm
    at Base e.

    Unlike some other transformations on the distribution and scale of a
    variable's values (Standard Z Score, Rank-Based Inverse Normal), the
    Logarithmic Transformation does not depend on the variance in the variable
    across all samples. Hence the Logarithmic Transformation is useful for
    exploratory analyses before any final filters or cohort stratifications on
    the samples.

    Table format must have samples (cases, observations) across rows and
    dependent and independent variables (features) across columns.

    This function creates new columns for values of the original variable after
    transformation.

    arguments:
        columns (list<str>): names of columns for continuous, ratio-scale
            variables to transform
        suffix (str): string suffix to append to the names of the original
            variables' columns
        table (object): Pandas data frame of variables across columns and
            samples across rows

    raises:

    returns:
        (object): Pandas data frame

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Filter columns by whether they are in the table.
    columns_relevant = list(filter(
        lambda column: (str(column) in table.columns.to_list()),
        columns
    ))
    # Apply logarithmic transformation to multiple variables.
    for column in columns_relevant:
        column_transform = str(str(column) + str(suffix))
        table[column_transform] = transform_values_distribution_scale_logarithm(
            values_array=table[column].to_numpy(),
            shift_minimum=1.0, # 0.001 or 1.0
            base=math.e, # e, 2, 10, etc
        )
    # Return information.
    return table


# Standard Z Score


def transform_values_distribution_scale_standard_z_score(
    values_array=None,
):
    """
    Transforms values within a NumPy array by Standard Z Score.

    This method ignores and propagates missing values in the original variable
    such that any missing values in the original variable do not affect the
    transformation (such as the variation in values across samples).

    An alternative would be to use sklearn.preprocessing.StandardScaler().

    arguments:
        values_array (object): NumPy array of original values

    raises:

    returns:
        (object): NumPy array of novel values after transformation

    """

    # Copy information.
    values_array = numpy.copy(values_array)
    # Calculate Standard Z Score across values in array.
    # This method inserts missing values if the standard deviation is zero.
    values_z_array = scipy.stats.zscore(
        values_array,
        axis=0,
        ddof=1, # (N - 1) Degrees of Freedom for Sample Standard Deviation.
        nan_policy="omit", # Ignore missing values in calculation.
    )
    # Return information.
    return values_z_array


def drive_transform_variables_distribution_scale_z_score(
    table=None,
    columns=None,
    report=None,
):
    """
    Transforms distribution and scale of variables' values by Standard Z Score.

    The mean of each variable's values will be zero (mean = 0).
    The standard deviation of each variable's values will be one (standard
    deviation = 1).

    arguments:
        table (object): Pandas data frame of variables across columns and
            samples across rows
        columns (list<str>): names of table's columns for which to standardize
            the scale of values
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of variables (features) across columns and
            samples (cases) across rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    table_scale = table.copy(deep=True)
    # Filter columns by whether they are in the table.
    columns_relevant = list(filter(
        lambda column: (str(column) in table.columns.to_list()),
        columns
    ))
    # Calculate standard scores by column.
    for column in columns_relevant:
        table_scale[column] = (
            transform_values_distribution_scale_standard_z_score(
                values_array=table_scale[column].to_numpy(),
        ))
        pass
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "drive_transform_variables_distribution_scale_standard_z_score()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        table_report = table.copy(deep=True)
        table_report = table_report.loc[
            :, table_report.columns.isin(columns_relevant)
        ]
        table_scale_report = table_scale.copy(deep=True)
        table_scale_report = table_scale_report.loc[
            :, table_scale_report.columns.isin(columns_relevant)
        ]
        print("Summary statistics before standardization.")
        table_mean = table_report.aggregate(
            lambda series: series.mean(),
            axis="index", # Apply function to each column of table.
        )
        print("Mean")
        print(table_mean.iloc[0:25])
        utility.print_terminal_partition(level=5)
        table_deviation = table_report.aggregate(
            lambda series: series.std(),
            axis="index", # Apply function to each column of table.
        )
        print("Standard deviation")
        print(table_deviation.iloc[0:25])
        utility.print_terminal_partition(level=4)
        print("Summary statistics after standardization.")
        table_mean = table_scale_report.aggregate(
            lambda series: series.mean(),
            axis="index", # Apply function to each column of table.
        )
        print("Mean")
        print(table_mean.iloc[0:25])
        utility.print_terminal_partition(level=5)
        table_deviation = table_scale_report.aggregate(
            lambda series: series.std(),
            axis="index", # Apply function to each column of table.
        )
        print("Standard deviation")
        print(table_deviation.iloc[0:25])
    # Return information.
    return table_scale


# Rank-Based Inverse Normalization


def transform_values_distribution_scale_rank_inverse(
    values_array=None,
):
    """
    Transforms values within a NumPy array by Rank-Based Inverse Normalization.

    This method ignores and propagates missing values in the original variable
    such that any missing values in the original variable do not affect the
    transformation (such as the variation in values across samples).

    arguments:
        values_array (object): NumPy array of original values

    raises:

    returns:
        (object): NumPy array of novel values after transformation

    """

    # Copy information.
    values_array = numpy.copy(values_array)
    # Calculate Rank-Base Inverse Normal.
    values_rank_array = numpy.squeeze(sklearn.preprocessing.quantile_transform(
        values_array,
        axis=0,
        n_quantiles=1e+6, # Use one quantile for each sample.
        output_distribution="normal",
        ignore_implicit_zeros=True,
        subsample=1e+6, # Count of Quantiles cannot exceed count of Samples.
        random_state=777,
        copy=True,
    ))
    # Return information.
    return values_rank_array


# Apply multiple transformations to a single variable.


def apply_transformations_to_variable_distribution_scale(
    column=None,
    logarithm_e=None,
    standard_z_score=None,
    rank_inverse=None,
    suffix_logarithm_e=None,
    suffix_standard_z_score=None,
    suffix_rank_inverse=None,
    table=None,
    report=None,
):
    """
    Applies transformations to the distribution and scale of a single variable's
    values.

    Transformations:
    1. Natural Logarithm at Base e
    2. Standard Z Score (mean of zero, standard deviation of one)
    3. Rank-Based Inverse Normal

    Some of these transformations (Standard Z Score, Rank-Based Inverse Normal)
    depend on the variance in the variable across all samples. Hence it is
    important to apply these transformations after any filters or cohort
    stratifications on the samples.

    The current methods for Standard Z Score and Rank-Based Inverse Normal
    ignore and propagate missing values in the original variable such that any
    missing values in the original variable do not affect the transformation
    (such as the variation in values across samples).

    Table format must have samples (cases, observations) across rows and
    dependent and independent variables (features) across columns.

    This function creates new columns for values of the original variable after
    transformation.

    arguments:
        column (str): name of table's column for the variable for transformation
        logarithm_e (bool): whether to apply the Base e Natural Logarithmic
            transformation to values of the variable
        standard_z_score (bool): whether to apply the Standard Z Score
            transformation to values of the variable
        rank_inverse (bool): whether to apply the Rank-Based Inverse Normal
            transformation to values of the variable
        suffix_logarithm_e (str): string suffix to append to the name of the
            original variable's column
        suffix_standard_z_score (str): string suffix to append to the name of
            the original variable's column
        suffix_rank_inverse (str): string suffix to append to the name of the
            original variable's column
        table (object): Pandas data frame of variables across columns and
            samples across rows
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of regression's residuals and statistics
    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Base e Natural Logarithmic Transformation.
    if ((str(column) in table.columns.to_list()) and (logarithm_e)):
        column_transform = str(str(column) + str(suffix_logarithm_e))
        table[column_transform] = transform_values_distribution_scale_logarithm(
            values_array=table[column].to_numpy(),
            shift_minimum=1.0, # 0.001 or 1.0
            base=math.e, # e, 2, 10, etc
        )

    # Standard Z-Score Transformation.
    if ((str(column) in table.columns.to_list()) and (standard_z_score)):
        column_transform = str(str(column) + str(suffix_standard_z_score))
        table[column_transform] = (
            transform_values_distribution_scale_standard_z_score(
                values_array=table[column].to_numpy(),
        ))

    # Rank-Based Inverse Normal Transformation.
    if ((str(column) in table.columns.to_list()) and (rank_inverse)):
        column_transform = str(str(column) + str(suffix_rank_inverse))
        table[column_transform] = (
            transform_values_distribution_scale_rank_inverse(
                values_array=table[column].to_numpy(),
        ))

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "apply_transformations_to_variable_distribution_scale()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table


# Drive transformations on multiple variables within multiple separate cohorts.


def drive_transformations_on_multiple_variables_in_cohorts(
    variables=None,
    records_cohorts=None,
    report=None,
):
    """
    Drives the application of Distribution Scale Transformations on multiple
    variables within separate tables for stratification cohorts.

    At a minimum, cohort records include a name of the cohort (key: "name") and
    a table (key: "table") of variables across relevant samples in the cohort.

    Table format must have samples (cases, observations) across rows and
    dependent and independent variables (features) across columns.

    arguments:
        variables (list<str>): name of columns within tables for variables on
            which to apply transformations
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Copy information.
    records_cohorts_scale = copy.deepcopy(records_cohorts)
    # Iterate across records.
    for record in records_cohorts_scale:
        # Iterate across variables.
        for variable in variables:
            record["table"] = (
                apply_transformations_to_variable_distribution_scale(
                    column=str(variable),
                    logarithm_e=True,
                    standard_z_score=True,
                    rank_inverse=True,
                    suffix_logarithm_e="_log",
                    suffix_standard_z_score="_z",
                    suffix_rank_inverse="rank",
                    table=record["table"],
                    report=False,
            ))
            pass
        pass
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "drive_transformations_on_multiple_variables_in_cohorts()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        count_cohorts = int(len(records_cohorts))
        count_cohorts_scale = int(len(records_cohorts_scale))
        print("count of original cohorts: " + str(count_cohorts))
        print("count of novel cohorts: " + str(count_cohorts_scale))
        print("variables for transformation: ")
        print(variables)
        pass
    # Return information
    return records_cohorts_scale





###############################################################################
# Procedure
# Currently, this module is not executable.
