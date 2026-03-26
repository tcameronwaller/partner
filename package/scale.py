"""
Supply functionality for transformation of distribution scales.

This module 'scale' is part of the 'partner' package.

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

# Here are a few informal (need review) notes about methods of normalization.
# robust scaling ((x - median(x) / interquartile_range(x)))
# logarithmic transformation
#   mitigate skewness, especially left
# square root transformation
#   mitigate skewness in presence of zeros
# square transformation
#   mitigate skewness, especially left
# exponential transformation
#   mitigate skewness, especially when relationship is exponential


###############################################################################
# Installation and importation

# Standard
import sys
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
import functools

# Relevant

import pandas
import sklearn
import sklearn.preprocessing
import scipy
import numpy
import statsmodels.api

# Custom
import partner.utility as putly # this import path for subpackage
import partner.organization as porg
#dir()
#importlib.reload()

###############################################################################
# Functionality


# Logarithm and exponent


# TODO: obsolete? (TCW; 13 December 2024)
def transform_values_distribution_scale_logarithm(
    values_array=None,
    shift_minimum=None,
    base=None,
):
    """
    Shifts values within a NumPy array to positive scale and then transforms
    their distribution and scale by Logarithm at a specific base.

    Logarithm does not have a definition for negative values.

    WARNING: The shift to positive values might introduce errors.

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


# TODO: obsolete? (TCW; 13 December 2024)
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
        table[column_transform] = (
            transform_values_distribution_scale_logarithm(
                values_array=table[column].to_numpy(),
                shift_minimum=1.0, # 0.001 or 1.0
                base=math.e, # e, 2, 10, etc
        ))
    # Return information.
    return table


def transform_logarithm_by_table_columns(
    table=None,
    columns=None,
    base=None,
    report=None,
):
    """
    Transform values by calculating the logarithm at a specific base to yield
    the inverse operation of the exponent.

    This function preserves the names and sequence of columns from the original
    source table.

    Review: TCW; 3 October 2025
    Review: TCW; 13 December 2024

    arguments:
        table (object): Pandas data-frame table of values on continuous ratio
            or interval scale of measurement
        columns (list<str>): names of columns in table for which to apply the
            transformation
        base (float): base for exponent
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values

    """

    # Define subordinate function.
    def calculate_base_logarithm(
        values=None,
        base=None,
    ):
        """
        Calculates the exponents by raise a specific base value to each of
        multiple values in an array.

        arguments:
            values (array): NumPy array of values to be the power in the
                exponential calculations
            base (base): value to be the base in the exponential calculations

        raises:

        returns:
            (array): novel exponents from calculating the base to the power of
                each original value

        """

        # Calculate logarithm.
        if (math.isclose(float(base), math.e, rel_tol=1e-7)):
            logarithms = numpy.log(values)
        elif (math.isclose(float(base), 2.0, rel_tol=1e-7)):
            logarithms = numpy.log2(values)
        elif (math.isclose(float(base), 10.0, rel_tol=1e-7)):
            logarithms = numpy.log10(values)
        else:
            # math.log(value, base)
            logarithms = (numpy.log(values) / numpy.log(base))
            pass
        # Handle missing values.
        logarithms = numpy.nan_to_num(
            logarithms,
            nan=numpy.nan,
            posinf=numpy.nan,
            neginf=numpy.nan,
            copy=True,
        )
        return logarithms

    # Copy information.
    table = table.copy(deep=True)
    columns = copy.deepcopy(columns)

    # Filter columns by whether they are available in the table.
    columns_relevant = list(filter(
        lambda column: (str(column) in table.columns.to_list()),
        columns
    ))

    # Calculate the standard z-score of values in each column of table.
    # This method inserts missing values if the standard deviation is zero.
    for column in columns_relevant:
        table[column] = calculate_base_logarithm(
            values=table[column].to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
            ),
            base=base,
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: scale.py")
        function = (
            "transform_logarithm_by_table_columns()"
        )
        print(str("function: " + function))
        putly.print_terminal_partition(level=5)
        print(str("base: " + str(base)))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


def transform_exponent_by_table_columns(
    table=None,
    columns=None,
    base=None,
    report=None,
):
    """
    Transform values by calculating the exponent at a specific base to yield
    the inverse operation of the logarithm.

    This function preserves the names and sequence of columns from the original
    source table.

    Review: TCW; 13 December 2024

    arguments:
        table (object): Pandas data-frame table of values on continuous ratio
            or interval scale of measurement
        columns (list<str>): names of columns in table for which to apply the
            transformation
        base (float): base for exponent
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values

    """

    # Define subordinate function.
    def calculate_base_power_exponent(
        values=None,
        base=None,
    ):
        """
        Calculates the exponents by raise a specific base value to each of
        multiple values in an array.

        arguments:
            values (array): NumPy array of values to be the power in the
                exponential calculations
            base (base): value to be the base in the exponential calculations

        raises:

        returns:
            (array): novel exponents from calculating the base to the power of
                each original value

        """

        # Calculate exponent.
        if (math.isclose(float(base), math.e, rel_tol=1e-7)):
            exponents = numpy.exp(values)
        elif (math.isclose(float(base), 2.0, rel_tol=1e-7)):
            exponents = numpy.exp2(values)
        else:
            exponents = numpy.power(base, values)
            pass
        # Handle missing values.
        exponents = numpy.nan_to_num(
            exponents,
            nan=numpy.nan,
            posinf=numpy.nan,
            neginf=numpy.nan,
            copy=True,
        )
        return exponents

    # Copy information.
    table = table.copy(deep=True)
    columns = copy.deepcopy(columns)

    # Filter columns by whether they are available in the table.
    columns_relevant = list(filter(
        lambda column: (str(column) in table.columns.to_list()),
        columns
    ))

    # Calculate the standard z-score of values in each column of table.
    # This method inserts missing values if the standard deviation is zero.
    for column in columns_relevant:
        table[column] = calculate_base_power_exponent(
            values=table[column].to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
            ),
            base=base,
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: scale.py")
        function = (
            "transform_exponent_by_table_columns()"
        )
        print(str("function: " + function))
        putly.print_terminal_partition(level=5)
        print(str("base: " + str(base)))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


def transform_product_by_table_columns(
    table=None,
    columns=None,
    factor=None,
    report=None,
):
    """
    Transform values by calculating the product with a specific factor.

    This function preserves the names and sequence of columns from the original
    source table.

    Review or revision: TCW; 7 January 2026

    arguments:
        table (object): Pandas data-frame table of values on continuous ratio
            or interval scale of measurement
        columns (list<str>): names of columns in table for which to apply the
            transformation
        base (float): base for exponent
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values

    """

    # Define subordinate functions for internal use.
    # Use this function with an "apply()" iterator.
    def calculate_product(
        factor_one=None,
        factor_two=None,
    ):
        # Determine whether there is adequate information for division.
        if (
            (factor_one is not None) and
            (factor_two is not None) and
            (not math.isnan(factor_one)) and
            (not math.isnan(factor_two))
        ):
            # Manage types of values.
            factor_one = float(factor_one)
            factor_two = float(factor_two)
            # There is adequate information for multiplication.
            product = (factor_one * factor_two)
        else:
            # There is inadequate information for division.
            product = float("nan")
            pass
        # Return information.
        return product

    # Copy information.
    table = table.copy(deep=True)
    columns = copy.deepcopy(columns)

    # Filter columns by whether they are available in the table.
    columns_relevant = list(filter(
        lambda column: (str(column) in table.columns.to_list()),
        columns
    ))

    # Calculate the standard z-score of values in each column of table.
    # This method inserts missing values if the standard deviation is zero.
    for column in columns_relevant:
        table[column] = numpy.multiply(
            table[column].to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
            ),
            factor,
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: scale.py")
        function = (
            "transform_product_by_table_columns()"
        )
        print(str("function: " + function))
        putly.print_terminal_partition(level=5)
        print(str("factor: " + str(factor)))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table




# Standard Z Score


# TODO: TCW; 13 December 2024
# Enhance this function by making it possible to apply the transformation only
# across a selection of columns within in each rows.
def transform_standard_z_score_by_table_rows(
    table=None,
    report=None,
):
    """
    Transform values to z-scores to standardize their scales and distributions
    within each separate row of a table.

    The mean of values in each row will be zero (mean = 0).
    The standard deviation of values in each row will be one (standard
    deviation = 1).

    If the table has any columns that do not correspond to values for
    transformation, then these columns must be explicitly defined as an index.

    Review: TCW; 16 October 2024

    arguments:
        table (object): Pandas data-frame table of values on continuous ratio
            or interval scale of measurement
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values

    """

    # Copy information in table.
    table = table.copy(deep=True)
    table_scale = table.copy(deep=True)

    # Calculate the standard z-score of values in each row of table.
    # This method inserts missing values if the standard deviation is zero.
    table_scale = table_scale.transform(
        lambda row: scipy.stats.zscore(
            row.to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            ),
            axis=0,
            ddof=1, # divisor is (n - 1) for sample standard deviation
            nan_policy="omit", # ignore missing values in calculation
        ),
        axis="columns", # apply function to each row
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("module: partner.scale.py")
        function = (
            "transform_standard_z_score_by_table_rows()"
        )
        print(str("function: " + function))

        putly.print_terminal_partition(level=3)
        table_report = table.copy(deep=True)
        table_scale_report = table_scale.copy(deep=True)
        print("... summary statistics before standardization ...")
        table_mean = table_report.aggregate(
            lambda series: series.mean(),
            axis="columns", # apply function to each row
        )
        print("mean:")
        print(table_mean.iloc[0:10])
        putly.print_terminal_partition(level=5)
        table_deviation = table_report.aggregate(
            lambda series: series.std(),
            axis="columns", # apply function to each row
        )
        print("standard deviation:")
        print(table_deviation.iloc[0:10])
        putly.print_terminal_partition(level=4)
        print("... summary statistics after standardization ...")
        table_mean = table_scale_report.aggregate(
            lambda series: series.mean(),
            axis="columns", # apply function to each row
        )
        print("mean:")
        print(table_mean.iloc[0:10])
        putly.print_terminal_partition(level=5)
        table_deviation = table_scale_report.aggregate(
            lambda series: series.std(),
            axis="columns", # apply function to each row
        )
        print("standard deviation:")
        print(table_deviation.iloc[0:10])
    # Return information.
    return table_scale


def transform_standard_z_score_by_table_columns(
    table=None,
    columns=None,
    report=None,
):
    """
    Transform values for features corresponding to specific columns in table
    to have mean zero (0) and standard deviation one (1) across all rows
    corresponding to observations.

    Transform values to z-scores to standardize their scales and distributions
    within each separate column of a table.

    The z score forces all values to have mean of zero (mean = 0) and standard
    deviation of one (standard deviation = 1).

    This function preserves the names and sequence of columns from the original
    source table.

    formula:
    z-score = (value - mean) / (standard deviation)

    References:

    1. Comparison of methods to adjust scale in 'Scikit Learn'
       title: 'Compare the effect of different scalers on data with outliers'
       site: https://scikit-learn.org/stable/auto_examples/preprocessing/
              plot_all_scaling.html#plot-all-scaling-minmax-scaler-section

    2. Documentation on 'StandardScaler' in 'Scikit Learn'
       title: 'StandardScaler'
       site: https://scikit-learn.org/stable/modules/generated/sklearn.
              preprocessing.StandardScaler.html#sklearn.preprocessing.
              StandardScaler

    Review: TCW; 2 April 2025

    arguments:
        table (object): Pandas data-frame table of values on continuous ratio
            or interval scale of measurement
        columns (list<str>): names of columns in table for which to apply the
            transformation
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values

    """

    # Copy information in table.
    table = table.copy(deep=True)
    table_scale = table.copy(deep=True)
    # Copy other information.
    columns = copy.deepcopy(columns)

    # Filter columns by whether they are in the table.
    columns_relevant = list(filter(
        lambda column: (str(column) in table.columns.to_list()),
        columns
    ))

    # Calculate the standard z-score of values in each column of table.
    # This method inserts missing values if the standard deviation is zero.
    for column in columns_relevant:
        table_scale[column] = scipy.stats.zscore(
            table_scale[column].to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
            ),
            axis=0,
            ddof=1, # divisor is (n - 1) for sample standard deviation
            nan_policy="omit", # ignore missing values in calculation
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: scale.py")
        name_function = (
            "transform_standard_z_score_by_table_columns()"
        )
        print(str("function: " + name_function))
        putly.print_terminal_partition(level=5)
        table_report = table.copy(deep=True)
        table_report = table_report.loc[
            :, table_report.columns.isin(columns_relevant)
        ]
        table_scale_report = table_scale.copy(deep=True)
        table_scale_report = table_scale_report.loc[
            :, table_scale_report.columns.isin(columns_relevant)
        ]
        print("... summary statistics before standardization ...")
        table_mean = table_report.aggregate(
            lambda series: series.mean(),
            axis="index", # apply function to each column
        )
        print("mean:")
        print(table_mean.iloc[0:10])
        putly.print_terminal_partition(level=5)
        table_deviation = table_report.aggregate(
            lambda series: series.std(),
            axis="index", # apply function to each column
        )
        print("standard deviation:")
        print(table_deviation.iloc[0:10])
        putly.print_terminal_partition(level=4)
        print("... summary statistics after standardization ...")
        table_mean = table_scale_report.aggregate(
            lambda series: series.mean(),
            axis="index", # apply function to each column
        )
        print("mean:")
        print(table_mean.iloc[0:10])
        putly.print_terminal_partition(level=5)
        table_deviation = table_scale_report.aggregate(
            lambda series: series.std(),
            axis="index", # apply function to each column
        )
        print("standard deviation:")
        print(table_deviation.iloc[0:10])
        pass
    # Return information.
    return table_scale


# obsolete???
def transform_values_distribution_scale_standard_z_score(
    values_array=None,
):
    """
    Transforms values within a NumPy array by Standard Z Score.

    This method ignores and propagates missing values in the original variable
    such that any missing values in the original variable do not affect the
    transformation (such as the variation in values across samples).

    An alternative would be to use sklearn.preprocessing.StandardScaler().

    Review: 24 September 2024

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
        ddof=1, # divisor is (n - 1) for sample standard deviation
        nan_policy="omit", # ignore missing values in calculation
    )
    # Return information.
    return values_z_array


def calculate_p_value_from_z_statistic(
    z_statistic=None,
    tail=None,
):
    """
    Calculate the p-value from the Z-statistic corresponding to a null
    hypothesis, assuming a normal distribution of potential values for the
    Z-statistic.

    continuous normal random distribution: "scipy.stats.norm()" in Python SciPy
    probability density function (PDF): "dnorm()" in R; "pdf()" in Python SciPy
    cumulative density function (CDF): "pnorm()" in R; "cdf()" in Python SciPy
    quantile inverse CDF: "qnorm()" in R; "ppf()" in Python SciPy

    The SciPy implementation and standard definition of the Cumulative Density
    Function (CDF) returns the probability that a random variable would take a
    value less than or equal to the given value (the left tail under the
    normal-distribution bell curve).

    Another condition of which to be aware is whether the distribution of
    potential values for the Z-statistic is symmetric about zero. If this
    distribution is symmetric about zero, then the two-tailed test p-value
    equals twice the CDF of negative one times the absolute value of the
    Z-statistic (p-value = 2 * CDF(-1 * abs(Z))).

                  .|.
                .  |  .
               .   |   .
              .    |    .
            .      |      .
         .         |         .
    .****          |          ****.

    Review: TCW; 15 March 2024

    arguments:
        z_statistic (float): value of Z-statistic for null hypothesis
        tail (str): 'left', 'right', or 'both' for tail selection

    raises:

    returns:
        (float): value of p-value assuming a normal distribution of potential
            values for the Z-statistic

    """

    # Calculate p-value from Z-statistic, assuming normal distribution.
    distribution = scipy.stats.norm() # continuous normal random distribution
    if (tail == "left"):
        # Left tail.
        p_value = float(distribution.cdf(z_statistic))
    elif (tail == "right"):
        # Right tail.
        p_value = float(1 - distribution.cdf(z_statistic))
    if (tail == "both"):
        # Both tails.
        # Only if symmetric.
        #p_value = 2 * float(distribution.cdf(-1 * abs(z_statistic)))
        # More general.
        p_value = 2 * min(
            float(distribution.cdf(z_statistic)),
            float(1 - distribution.cdf(z_statistic)),
        )
    # Return information.
    return p_value


# Unit range


def transform_unit_range_by_table_columns(
    table=None,
    columns=None,
    report=None,
):
    """
    Transform values for features corresponding to specific columns in table
    to have unit range from zero (0) to one (1) across all rows corresponding
    to observations.

    This function preserves the names and sequence of columns from the original
    source table.

    formula:
    unit-range = (value - minimum) / (maximum - minimum)

    References:

    1. Comparison of methods to adjust scale in 'Scikit Learn'
       title: 'Compare the effect of different scalers on data with outliers'
       site: https://scikit-learn.org/stable/auto_examples/preprocessing/
          plot_all_scaling.html#plot-all-scaling-minmax-scaler-section

    2. Documentation on 'MinMaxScaler' in 'Scikit Learn'
       title: 'MinMaxScaler'
       site: https://scikit-learn.org/stable/modules/generated/sklearn.
              preprocessing.MinMaxScaler.html

    Review: TCW; 2 April 2025

    arguments:
        table (object): Pandas data-frame table of values on continuous ratio
            or interval scale of measurement
        columns (list<str>): names of columns in table for which to apply the
            transformation
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values

    """

    def transform_values_unit_range(
        series=None,
    ):
        """
        Define subordinate function for internal use within dominant function.
        """
        def unit_range(value, minimum, maximum):
            return ((value - minimum) / (maximum - minimum))
        vectorized_unit_range = numpy.vectorize(unit_range)
        series = series.copy(deep=True)
        values_raw = series.to_numpy(
            dtype="float64",
            na_value=numpy.nan,
            copy=True,
        )
        minimum = numpy.nanmin(values_raw)
        maximum = numpy.nanmax(values_raw)
        values_scale = vectorized_unit_range(
            values_raw, minimum, maximum
        )
        return values_scale

    # Copy information in table.
    table = table.copy(deep=True)
    table_scale = table.copy(deep=True)
    # Copy other information.
    columns = copy.deepcopy(columns)

    # Filter columns by whether they are in the table.
    columns_relevant = list(filter(
        lambda column: (str(column) in table.columns.to_list()),
        columns
    ))

    # Calculate the standard z-score of values in each column of table.
    # This method inserts missing values if the standard deviation is zero.
    for column in columns_relevant:
        table_scale[column] = transform_values_unit_range(
            series=table_scale[column]
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: scale.py")
        name_function = (
            "transform_standard_z_score_by_table_columns()"
        )
        print(str("function: " + name_function))
        putly.print_terminal_partition(level=5)
        table_report = table.copy(deep=True)
        table_report = table_report.loc[
            :, table_report.columns.isin(columns_relevant)
        ]
        table_scale_report = table_scale.copy(deep=True)
        table_scale_report = table_scale_report.loc[
            :, table_scale_report.columns.isin(columns_relevant)
        ]
        print("... summary statistics before scale ...")
        table_minimum = table_report.aggregate(
            lambda series: series.min(),
            axis="index", # apply function to each column
        )
        print("minimum:")
        print(table_minimum.iloc[0:10])
        putly.print_terminal_partition(level=5)
        table_maximum = table_report.aggregate(
            lambda series: series.max(),
            axis="index", # apply function to each column
        )
        print("maximum:")
        print(table_maximum.iloc[0:10])
        putly.print_terminal_partition(level=4)
        print("... summary statistics after scale ...")
        table_minimum = table_scale_report.aggregate(
            lambda series: series.min(),
            axis="index", # apply function to each column
        )
        print("minimum:")
        print(table_minimum.iloc[0:10])
        putly.print_terminal_partition(level=5)
        table_maximum = table_scale_report.aggregate(
            lambda series: series.max(),
            axis="index", # apply function to each column
        )
        print("maximum:")
        print(table_maximum.iloc[0:10])
        pass
    # Return information.
    return table_scale


# Manage transformation of scale for values of features on a quantitative,
# continuous, ratio or interval scale of measurement.


def manage_transform_scale_feature_by_table_columns(
    table=None,
    features_continuity_scale=None,
    adjust_scale=None,
    method_scale=None,
    report=None,
):
    """
    Dependency:
    This function is a dependency of the functions below.
    1.
    package: partner
    module or script: organization.py
    function: prepare_table_features_observations_for_analysis()

    Manage transformation of scale for values of features on a quantitative,
    continuous, ratio or interval scale of measurement.

    ----------
    Format of source data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. For versatility, this table does not have
    explicitly defined indices across rows or columns.
    ----------
    identifiers     feature_1 feature_2 feature_3 feature_4 feature_5 ...

    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    Date, revision or review: 20 February 2026
    Date, revision or review: 8 May 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for analysis
        features_continuity_scale (list<str>): names of columns in data table
            for feature variables with values on quantitative, continuous scale
            of measurement, interval or ratio, for which to standardize the
            scale by z score
        adjust_scale (bool): whether to adjust or standardize the scale of
            values for features across observations
        method_scale (str): name of method to use to adjust the scale of values
            for features across observations, either 'z_score' or 'unit_range'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of data with features and
            observations for analysis

    """

    # Copy information.
    table = table.copy(deep=True)
    features_continuity_scale = copy.deepcopy(features_continuity_scale)

    # Standardize scale of values for observations of features.
    if (
        (adjust_scale) and
        (len(features_continuity_scale) > 0) and
        (str(features_continuity_scale).strip().lower() != "none") and
        (method_scale is not None) and
        (str(method_scale).strip().lower() != "none")
    ):
        if (
            (str(method_scale).strip().lower() == "z_score")
        ):
            table = transform_standard_z_score_by_table_columns(
                    table=table,
                    columns=features_continuity_scale,
                    report=report,
            )
        elif (
            (str(method_scale).strip().lower() == "unit_range")
        ):
            table = transform_unit_range_by_table_columns(
                    table=table,
                    columns=features_continuity_scale,
                    report=report,
            )
            pass
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: scale.py")
        name_function = (
            "manage_transform_scale_feature_by_table_columns()"
        )
        print(str("function: " + name_function))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


##########
# Scale alignment between samples by sum total signal intensity


def scale_signals_between_samples_by_sum_total_signal_intensity(
    table=None,
    name_index_columns=None,
    name_index_rows=None,
    columns_signal_scale=None,
    scale_factor_batch=None,
    report=None,
):
    """
    This function is a dependency of the function below.
    package: partner
    module or script: 1_organize_table_signal_proteomics
    function: align_scale_signal_between_samples_by_sum_total()

    Align scale of values for signal intensity between separate samples using
    the sum total signal intensity across all analytes in each sample.

    Date, revision or review: 4 February 2026

    arguments:
        table (object): Pandas data-frame table of signals for analytes
            from their measurements across samples
        name_index_columns (str): name of index across columns in table
        name_index_rows (str): name of index across rows in table
        columns_signal_scale (list<str>): names of column in table
        scale_factor_batch (float): extra factor for adjustment of signals for
            all samples in batch
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table = table.copy(deep=True)
    table_total = table.copy(deep=True)
    columns_signal_scale = copy.deepcopy(columns_signal_scale)

    # Define indices in table.
    table_total = porg.explicate_table_indices_columns_rows_single_level(
        table=table_total,
        index_columns=name_index_columns,
        index_rows=name_index_rows,
        explicate_indices=True,
        report=report,
    )
    # Calculate sum total values of signal intensity across all analytes in
    # each sample.
    intensities_sum_total = table_total[columns_signal_scale].aggregate(
        lambda column: numpy.nansum(
            column.to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            )
        ),
        axis="index", # apply function to each column
    )

    # Calculate mean of sum total values of signal intensity across all
    # analytes in all samples.
    intensity_total_mean = numpy.nanmean(
        intensities_sum_total.to_numpy(
            dtype="float64",
            na_value=numpy.nan,
            copy=True,
        )
    )
    # Convert Pandas series to data-frame table.
    table_factor = intensities_sum_total.to_frame(
        name="intensity",
    )
    # Define indices in table.
    table_factor.index.set_names("sample", inplace=True)
    table_factor.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Calculate the scale factor.
    table_factor["intensity_total"] = intensity_total_mean
    table_factor["scale_factor_1"] = table_factor.apply(
        lambda row: putly.divide_accommodate_zero_missing(
            numerator=row["intensity_total"],
            denominator=row["intensity"],
        ),
        axis="columns", # apply function to each row
    )
    table_factor["scale_factor_2"] = table_factor.apply(
        lambda row: numpy.multiply(
            row["scale_factor_1"],
            scale_factor_batch,
        ),
        axis="columns", # apply function to each row
    )
    # Define indices in table.
    table_factor = porg.explicate_table_indices_columns_rows_single_level(
        table=table_factor,
        index_columns="column",
        index_rows="sample",
        explicate_indices=True,
        report=report,
    )
    # Calculate the scale adjustments for each analyte.
    for column_signal in columns_signal_scale:
        #name_column = str(column_signal + "_scale")
        table[column_signal] = table.apply(
            lambda row: numpy.multiply(
                row[column_signal],
                table_factor["scale_factor_2"][column_signal],
            ),
            axis="columns", # apply function to each row
        )
        pass
    # Remove unnecessary columns from table.
    table.drop(
        labels=[],
        axis="columns",
        inplace=True
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("scale.py")
        print("module: " + name_module)
        name_function = str(
            "scale_signals_between_samples_by_sum_total_signal_intensity()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


##########
# Scale alignment between batches by internal standard bridge samples


def scale_signals_between_batches_by_bridge_internal_standard(
    table=None,
    name_index_columns=None,
    name_index_rows=None,
    columns_bridges_batch=None,
    columns_bridges_all=None,
    columns_observations_batch=None,
    adjust_scale_bridges_batch=None,
    report=None,
):
    """
    This function is a dependency of the function below.
    package: partner
    module or script: 1_organize_table_signal_proteomics
    function: manage_merge_scale_alignment_distribution_normalization()

    Adjust scale of values for signals between separate batches using internal
    standards called bridge samples.

    Date, revision or review: 28 January 2026

    arguments:
        table (object): Pandas data-frame table of signals for analytes
            from their measurements across samples
        name_index_columns (str): name of index across columns in table
        name_index_rows (str): name of index across rows in table
        columns_bridges_batch (list<str>): names of column in table
        columns_bridges_all (list<str>): names of column in table
        columns_observations_batch (list<str>): names of columns in table
        adjust_scale_bridges_batch (bool): whether to adjust the scale of
            signals from bridge samples in the batch
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Define subordinate functions for internal use.
    def calculate_product(
        factor_one=None,
        factor_two=None,
    ):
        # Determine whether there is adequate information for division.
        if (
            (factor_one is not None) and
            (factor_two is not None) and
            (not math.isnan(factor_one)) and
            (not math.isnan(factor_two))
        ):
            # Manage types of values.
            factor_one = float(factor_one)
            factor_two = float(factor_two)
            # There is adequate information for multiplication.
            product = (factor_one * factor_two)
        else:
            # There is inadequate information for division.
            product = float("nan")
            pass
        # Return information.
        return product

    # Copy information.
    table = table.copy(deep=True)
    columns_bridges_batch = copy.deepcopy(columns_bridges_batch)
    columns_bridges_all = copy.deepcopy(columns_bridges_all)
    columns_observations_batch = copy.deepcopy(columns_observations_batch)
    columns_signals_scale = copy.deepcopy(columns_observations_batch)

    # Calculate mean values of signals from bridge and observation samples in
    # batch and across all batches.
    #table["mean"] = ((table["column_1"] + table["column_2"]) / 2)
    #float((row["column_1"] + row["column_2"]) / 2)
    if (len(columns_bridges_batch) > 1):
        table["mean_bridges_batch"] = table.apply(
            lambda row: numpy.nanmean(row[columns_bridges_batch]),
            axis="columns", # apply function to each row
        )
    else:
        table["mean_bridges_batch"] = table[columns_bridges_batch[0]]
        pass
    table["mean_bridges_all"] = table.apply(
        lambda row: numpy.nanmean(row[columns_bridges_all]),
        axis="columns", # apply function to each row
    )
    table["mean_observations"] = table.apply(
        lambda row: numpy.nanmean(row[columns_observations_batch]),
        axis="columns", # apply function to each row
    )
    # Calculate the scale factor.
    table["scale_factor"] = table.apply(
        lambda row: putly.divide_accommodate_zero_missing(
            numerator=row["mean_bridges_all"],
            denominator=row["mean_bridges_batch"],
        ),
        axis="columns", # apply function to each row
    )
    # Determine whether to adjust the scale of signals from bridge samples in
    # the batch.
    if (adjust_scale_bridges_batch):
        columns_signals_scale.extend(columns_bridges_batch)
        pass
    # Calculate the scale adjustments for each analyte.
    for column_signal in columns_signals_scale:
        #name_column = str(column_signal + "_scale")
        #table[column_signal] = table.apply(
        #    lambda row: calculate_product(
        #        factor_one=row[column_signal],
        #        factor_two=row["scale_factor"],
        #    ),
        #    axis="columns", # apply function to each row
        #)
        table[column_signal] = table.apply(
            lambda row: numpy.multiply(
                row[column_signal],
                row["scale_factor"],
            ),
            axis="columns", # apply function to each row
        )
        pass

    # Remove unnecessary columns from table.
    table.drop(
        labels=[
            "mean_bridges_batch",
            "mean_bridges_all",
            "mean_observations",
            "scale_factor",
        ],
        axis="columns",
        inplace=True
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("scale.py")
        print("module: " + name_module)
        name_function = str(
            "scale_signals_between_batches_by_bridge_internal_standard()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


def calculate_difference_signals_between_batches_by_bridge(
    table=None,
    name_index_columns=None,
    name_index_rows=None,
    columns_bridges_batch_1=None,
    columns_bridges_batch_2=None,
    filter_sort_columns=None,
    sort_rows=None,
    report=None,
):
    """
    Calculate the difference in signals for individual analytes between mean of
    bridge samples from each of two batches.

    Date, revision or review: 3 February 2026

    arguments:
        table (object): Pandas data-frame table of signals for analytes
            from their measurements across samples
        name_index_columns (str): name of index across columns in table
        name_index_rows (str): name of index across rows in table
        columns_bridges_batch_1 (list<str>): names of column in table
        columns_bridges_batch_2 (list<str>): names of column in table
        filter_sort_columns (bool): whether to filter and sort to minimal
            columns in table that are most relevant
        sort_rows (bool): whether to sort rows in table
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table = table.copy(deep=True)
    columns_bridges_batch_1 = copy.deepcopy(columns_bridges_batch_1)
    columns_bridges_batch_2 = copy.deepcopy(columns_bridges_batch_2)
    columns_bridges_all = copy.deepcopy(columns_bridges_batch_1)
    columns_bridges_all.extend(copy.deepcopy(columns_bridges_batch_2))

    # Calculate mean values of signal intensity for individual analytes from
    # bridge samples in each batch.
    #table["mean"] = ((table["column_1"] + table["column_2"]) / 2)
    #float((row["column_1"] + row["column_2"]) / 2)
    if (len(columns_bridges_batch_1) > 1):
        table["mean_bridges_batch_1"] = table.apply(
            lambda row: numpy.nanmean(row[columns_bridges_batch_1]),
            axis="columns", # apply function to each row
        )
    else:
        table["mean_bridges_batch_1"] = table[columns_bridges_batch_1[0]]
        pass
    if (len(columns_bridges_batch_2) > 1):
        table["mean_bridges_batch_2"] = table.apply(
            lambda row: numpy.nanmean(row[columns_bridges_batch_2]),
            axis="columns", # apply function to each row
        )
    else:
        table["mean_bridges_batch_2"] = table[columns_bridges_batch_2[0]]
        pass
    table["mean_bridges_all"] = table.apply(
        lambda row: numpy.nanmean(row[columns_bridges_all]),
        axis="columns", # apply function to each row
    )

    # Calculate the difference in values of signal intensity for individual
    # analytes from bridge samples in each batch.
    table["difference_bridge_batch"] = table.apply(
        lambda row: (
            row["mean_bridges_batch_1"] - row["mean_bridges_batch_2"]
        ),
        axis="columns", # apply function to each row
    )
    # Calculate the relative difference.
    table["difference_bridge_batch_relative"] = table.apply(
        lambda row: putly.divide_accommodate_zero_missing(
            numerator=row["difference_bridge_batch"],
            denominator=row["mean_bridges_all"],
        ),
        axis="columns", # apply function to each row
    )
    table["difference_bridge_batch_percent"] = numpy.multiply(
        table["difference_bridge_batch_relative"], 100
    )

    # Filter rows in table.
    columns_sequence = [
        "difference_bridge_batch_percent",
        "difference_bridge_batch_relative",
        "difference_bridge_batch",
        "mean_bridges_all",
        "mean_bridges_batch_1",
        "mean_bridges_batch_2",
    ]
    table.dropna(
        axis="index",
        how="all",
        subset=columns_sequence,
        inplace=True,
    )
    # Sort sequence of rows in table.
    if (sort_rows):
        table.sort_values(
            by=["difference_bridge_batch_percent"],
            axis="index",
            ascending=True,
            na_position="last",
            inplace=True,
        )
        # Restore or reset indices in table to generic default.
        table.reset_index(
            level=None,
            inplace=True,
            drop=False, # remove index; do not move to regular columns
        )
        table.columns.rename(
            None,
            inplace=True,
        ) # single-dimensional index
        pass
    # Ranks of analytes in sort sequence.
    table["rank"] = range(1, len(table) + 1)

    # Specify sequence of columns within table.
    columns_sequence = [
        name_index_rows,
        "rank",
        "difference_bridge_batch_percent",
        "difference_bridge_batch_relative",
        "difference_bridge_batch",
        "mean_bridges_all",
        "mean_bridges_batch_1",
        "mean_bridges_batch_2",
    ]
    # Filter and sort sequence of columns within table.
    if (filter_sort_columns):
        table = porg.filter_sort_table_columns(
            table=table,
            columns_sequence=columns_sequence,
            report=report,
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("scale.py")
        print("module: " + name_module)
        name_function = str(
            "calculate_difference_signals_between_batches_by_bridge()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


##########
# Scale alignment between samples by the DESeq method of median-ratio scaling


def calculate_ratios_to_geometric_mean_across_feature_observations(
    table=None,
    columns=None,
):
    """
    Dependency:
    This function is a dependency of the function below.
    partner.scale.scale_feature_values_between_observations_by_deseq()

    Calculates for each feature the ratio of observation values to the
    geometric mean across all observations.

    This function creates a new column named "mean_geometric". This function
    also creates a new column with suffix "_ratio" for each column in the
    original table. This function also creates a new column with suffix
    "_ratio_log" for each column in the original table.

    It is important and necessary that there are not any values in the table
    that are missing or less than zero. Careful and thorough filter and
    imputation operations on the data ought to happen before calling this
    function.

    For this method, it is necessary to exclude any features that have a
    geometric mean of zero across samples.

    Table's format and orientation

    Table has values for each feature oriented across rows with their
    observations oriented across columns. All values are on a continous ratio
    scale of measurement and have the same type, such as corresponding to
    intensities from a particular type of measurement.

    The table must have a single-level index (potentially with name
    "observations") across columns and a single-level index (potentially with
    name "features") across rows.

    observations   observation_1 observation_2 observation_3 observation_4
    features
    feature_1      ...           ...           ...           ...
    feature_2      ...           ...           ...           ...
    feature_3      ...           ...           ...           ...
    feature_4      ...           ...           ...           ...
    feature_5      ...           ...           ...           ...
    feature_6      ...           ...           ...           ...
    feature_7      ...           ...           ...           ...
    feature_8      ...           ...           ...           ...
    feature_9      ...           ...           ...           ...
    feature_10     ...           ...           ...           ...

    This function does not modify the names of indices across columns or rows
    from the original table.

    In its implementation, this function attempts to prioritize clarity over
    efficiency. While matrix transformations would be more efficient, they
    would also be more difficult to follow and critique.

    Review: 29 April 2025
    Review: 17 September 2024

    arguments:
        table (object): Pandas data-frame table of values for observations
            across columns and for features across rows
        columns (list<str>): names of columns corresponding to values of
            features across observations

    raises:

    returns:
        (object): Pandas data-frame table of values across observations in
            columns and across features in rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Calculate geometric mean of values for each feature across all
    # observations.
    table["mean_geometric"] = table.apply(
        lambda row:
            scipy.stats.mstats.gmean(
                row[columns].to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
                ),
                nan_policy="omit",
            ),
        axis="columns", # apply function to each row
    )

    #table["mean_geometric"] = table["mean_geometric"].replace(
    #    to_replace=0.0,
    #    value=pandas.NA,
    #)
    #table.dropna(
    #    axis="index",
    #    how="all",
    #    subset=["mean_geometric"],
    #    inplace=True,
    #)

    # Calculate ratio of each observation of each feature to that feature's
    # geometric mean value across all observations.
    # For this method, it is necessary to return missing values for any
    # features that have a geometric mean of zero across samples.
    for column in columns:
        table[str(column + "_ratio")] = table.apply(
            lambda row:
                (row[column] / row["mean_geometric"])
                if (
                    (not pandas.isna(row["mean_geometric"])) and
                    (row["mean_geometric"] != 0.0)
                ) else pandas.NA,
            axis="columns", # apply function to each row
        )
        pass
    # Copy information in table.
    table = table.copy(deep=True)
    # Return information.
    return table


def scale_feature_values_between_observations_by_deseq(
    table=None,
    name_columns=None,
    name_rows=None,
    report=None,
):
    """
    Scale values of features between observations by their
    observation-specific median ratio to their geometric mean value across all
    observations.

    This scaling method is defined for comparisons of very many features,
    such as thousands to tens of thousands of genes or proteins, between
    comparatively few observations or samples. The goal of this method is to
    remove the influence of technical confounders so that comparisons between
    observations are more relevant to real differences, such as those due to
    biological regulation of expression of genes and proteins. A fundamental
    assumption of this scaling method is that between observations the vast
    majority of features, such as genes or proteins, have similar values.

    For this method, it is necessary to exclude any features that have a
    geometric mean of zero across observations.

    The popular DESeq2 tool applies this strategy of scaling or normalization.

    References:
    1. Anders et al, Genome Biology, 2010; PubMed:20979621
    2. Bernstein, "Median-ratio normalization for bulk RNA-seq data", 2023;
       <https://mbernste.github.io/posts/median_ratio_norm/>

    It is important and necessary that there are not any values in the table
    that are missing or less than zero. Careful and thorough filter and
    imputation operations on the data ought to happen before calling this
    function.

    Table's format and orientation

    Table has values for each feature oriented across rows with their
    observations oriented across columns. All values are on a continous ratio
    scale of measurement and have the same type, such as corresponding to
    signal intensities from a particular type of measurement.

    The table must have a single-level index (potentially with name
    "observations") across columns and a single-level index (potentially with
    name "features") across rows.

    observations   observation_1 observation_2 observation_3 observation_4
    features
    feature_1      ...           ...           ...           ...
    feature_2      ...           ...           ...           ...
    feature_3      ...           ...           ...           ...
    feature_4      ...           ...           ...           ...
    feature_5      ...           ...           ...           ...
    feature_6      ...           ...           ...           ...
    feature_7      ...           ...           ...           ...
    feature_8      ...           ...           ...           ...
    feature_9      ...           ...           ...           ...
    feature_10     ...           ...           ...           ...

    This function does not modify the names of indices across columns or rows
    from the original table.

    In its implementation, this function attempts to prioritize clarity over
    efficiency. While matrix transformations would be more efficient, they
    would also be more difficult to follow and critique.

    Date, revision or review: 3 February 2026
    Date, revision or review: 29 April 2025
    Date, revision or review: 13 September 2024

    arguments:
        table (object): Pandas data-frame table of values for observations
            across columns and for features across rows
        name_columns (str): name of single-level index across columns
        name_rows (str): name of single-level index across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values for observations across
            columns and for features across rows

    """

    # Define subordinate functions for internal use.
    def calculate_scale_value(value_source=None, median_ratio=None):
        if (
            (not pandas.isna(median_ratio)) and
            (median_ratio != 0.0)
        ):
            value_product = (
                value_source / median_ratio
            )
        else:
            value_product = pandas.NA
        return value_product

    # Copy information in table.
    table_scale = table.copy(deep=True)
    # Copy names of columns and rows in original table.
    names_columns = copy.deepcopy(
        table_scale.columns.get_level_values(name_columns).to_list()
    )
    names_rows = copy.deepcopy(
        table_scale.index.get_level_values(name_rows).to_list()
    )
    names_columns_ratio = list(map(
        lambda name: str(name + "_ratio"),
        names_columns,
    ))
    names_columns_scale = list(map(
        lambda name: str(name + "_scale"),
        names_columns,
    ))
    # Calculate geometric mean of values for each feature across all
    # observations.
    # Calculate ratio of each observation of each feature to that feature's
    # geometric mean value across all observations.
    table_scale = (
        calculate_ratios_to_geometric_mean_across_feature_observations(
            table=table_scale,
            columns=names_columns,
    ))

    # Note: TCW; 3 February 2026
    # The median is independent of the distribution of values. It does not
    # require normality.

    # In each observation, a few features will vary greater than the geometric
    # mean, while a few features will vary less than the geometric mean. An
    # assumption is that most features will not vary much from the true
    # mean, and any variation in these features will correspond to technical
    # artifacts. The next goal is to identify the feature in each observation
    # that varies the least: median ratio, neither greater nor less than the
    # geometric mean. This most stable feature will represent a sort of
    # internal reference standard for that observation.

    # Calculate median value of ratios for each observation across all
    # features.
    # These median ratios will be the scaling factor for each observation of
    # all features.
    median_ratios = table_scale[names_columns_ratio].aggregate(
        lambda column: numpy.nanmedian(
            column.to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            )
        ),
        axis="index", # apply function to each column
    )
    # Calculate the scale values via division by the median-ratio scale factors
    # for each observation.
    for column in names_columns:
        table_scale[str(column + "_scale")] = table_scale.apply(
            lambda row:
                calculate_scale_value(
                    value_source=row[column],
                    median_ratio=median_ratios[str(column + "_ratio")],
                ),
            axis="columns", # apply function to each row
        )
        pass

    # Organize information in table.
    # Copy information in table.
    table_scale_clean = table_scale[names_columns_scale].copy(deep=True)
    # Translate names of columns.
    translations = dict()
    names_iterator = zip(
        names_columns_scale,
        names_columns,
    )
    for (column_scale, column) in names_iterator:
        translations[column_scale] = column
        pass
    table_scale_clean.rename(
        columns=translations,
        inplace=True,
    )
    # Remove unnecessary columns.
    if False:
        table_scale.drop(
            labels=["mean_geometric",],
            axis="columns",
            inplace=True
        )
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("module: partner.scale.py")
        print(
            "function: " +
            "scale_feature_values_between_observations_by_deseq()"
        )
        putly.print_terminal_partition(level=3)
        print(
            "Original source table of values for observations across " +
            "columns and for features across rows:")
        print(table)
        putly.print_terminal_partition(level=4)
        print(
            "Series of geometric means across all observations for each " +
            "feature:"
        )
        print(table_scale["mean_geometric"])
        putly.print_terminal_partition(level=4)
        table_ratio = table_scale[names_columns_ratio]
        print(
            "Table of ratios between each original value and the geometric " +
            "mean corresponding to each feature."
        )
        print(table_ratio)
        putly.print_terminal_partition(level=4)
        print(
            "Series of median-ratio scale factors for each observation:"
        )
        print(median_ratios)
        putly.print_terminal_partition(level=4)
        table_scale_preliminary = table_scale[names_columns_scale]
        print(
            "Table of scaled values before organizing the novel product table:"
        )
        print(table_scale_preliminary)
        putly.print_terminal_partition(level=4)
        print(
            "Novel product table of scaled values for observations across " +
            "columns and for features across rows:"
        )
        print(table_scale_clean)
    # Return information.
    return table_scale_clean


##########
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
    # Structure information in NumPy matrix.
    # numpy.expand_dims()
    values_matrix = numpy.reshape(values_array, (-1, 1), order="C")
    # Calculate Rank-Base Inverse Normal.
    values_rank_array = numpy.squeeze(sklearn.preprocessing.quantile_transform(
        values_matrix,
        axis=0,
        n_quantiles=1e+6, # Use one quantile for each sample.
        output_distribution="normal",
        ignore_implicit_zeros=False,
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
        table (object): Pandas data frame table of variables across columns and
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
        table[column_transform] = (
            transform_values_distribution_scale_logarithm(
                values_array=table[column].to_numpy(),
                shift_minimum=1.0, # 0.001 or 1.0
                base=math.e, # e, 2, 10, etc
        ))

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
        putly.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "apply_transformations_to_variable_distribution_scale()"
        )
        print(name_function)
        putly.print_terminal_partition(level=3)
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
                    suffix_rank_inverse="_rank",
                    table=record["table"],
                    report=False,
            ))
            pass
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "drive_transformations_on_multiple_variables_in_cohorts()"
        )
        print(name_function)
        putly.print_terminal_partition(level=3)
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
# End
