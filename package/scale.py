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
        ddof=1, # divisor is (n - 1) for sample standard deviation
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

    This function does not modify the names of the columns for the original
    variables.

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
        putly.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "drive_transform_variables_distribution_scale_standard_z_score()"
        )
        print(name_function)
        putly.print_terminal_partition(level=3)
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
        putly.print_terminal_partition(level=5)
        table_deviation = table_report.aggregate(
            lambda series: series.std(),
            axis="index", # Apply function to each column of table.
        )
        print("Standard deviation")
        print(table_deviation.iloc[0:25])
        putly.print_terminal_partition(level=4)
        print("Summary statistics after standardization.")
        table_mean = table_scale_report.aggregate(
            lambda series: series.mean(),
            axis="index", # Apply function to each column of table.
        )
        print("Mean")
        print(table_mean.iloc[0:25])
        putly.print_terminal_partition(level=5)
        table_deviation = table_scale_report.aggregate(
            lambda series: series.std(),
            axis="index", # Apply function to each column of table.
        )
        print("Standard deviation")
        print(table_deviation.iloc[0:25])
    # Return information.
    return table_scale


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


##########
# Scale adjustment by the DESeq method of median-ratio scaling


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

    Review: 13 September 2024

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
    # For this method, it is necessary to exclude any features that have a
    # geometric mean of zero across samples.
    table["mean_geometric"] = table["mean_geometric"].replace(
        to_replace=0.0,
        value=pandas.NA,
    )
    table.dropna(
        axis="index",
        how="all",
        subset=["mean_geometric"],
        inplace=True,
    )
    # Calculate ratio of each observation of each feature to that feature's
    # geometric mean value across all observations.
    for column in columns:
        table[str(column + "_ratio")] = table.apply(
            lambda row:
                (row[column] / row["mean_geometric"])
                if row["mean_geometric"] != 0.0 else pandas.NA,
            axis="columns", # apply function to each row
        )
        pass
    # Calculate the logarithm of ratios.
    for column in columns:
        # math.log() # optimal for scalar values
        # numpy.log() # optimal for array values
        table[str(column + "_ratio_log")] = table.apply(
            lambda row:
                math.log(row[str(column + "_ratio")])
                if row[str(column + "_ratio")] != 0.0 else pandas.NA,
            axis="columns", # apply function to each row
        )
        pass
    # Copy information in table.
    table = table.copy(deep=True)
    # Return information.
    return table


def filter_table_features_least_change_across_observations(
    table=None,
    name_columns=None,
    name_rows=None,
    count_quantile=None,
    report=None,
):
    """
    Filters rows in table to select features that demonstrate least change
    across observations on the basis of the middle quantile.

    This function is a handy companion to the function below.
    partner.scale.scale_feature_values_between_observations_by_deseq()

    Tables' format and orientation

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

    Review: 2 July 2024

    arguments:
        table (object): Pandas data-frame table of values for observations
            across columns and for features across rows
        name_columns (str): name of single-level index across columns
        name_rows (str): name of single-level index across rows
        count_quantile (int): odd count of quantiles, such as three for
            tertiles
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values for observations across
            columns and for features across rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy names of columns in original table.
    columns = copy.deepcopy(
        table.columns.get_level_values(name_columns).to_list()
    )
    # Determine names of derivative columns.
    columns_ratio = list(map(
        lambda name: str(name + "_ratio"),
        columns,
    ))
    columns_ratio_log = list(map(
        lambda name: str(name + "_ratio_log"),
        columns,
    ))
    # Calculate geometric mean of values for each feature across all
    # observations.
    # Calculate ratio of each observation for each feature to that feature's
    # geometric mean value across all observations.
    table_ratio = (
        calculate_ratios_to_geometric_mean_across_feature_observations(
            table=table,
            columns=columns,
    ))
    # Calculate mean of ratio values across all observations for each feature.
    table_ratio["mean_ratio"] = table_ratio.apply(
        lambda row:
            numpy.nanmean(
                row[columns_ratio].to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
                )
            ),
        axis="columns", # apply function to each row
    )
    table_ratio["mean_ratio_log"] = table_ratio.apply(
        lambda row:
            numpy.nanmean(
                row[columns_ratio_log].to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
                )
            ),
        axis="columns", # apply function to each row
    )
    # Determine indices for quantiles.
    indices = list(range(0, count_quantile, 1))
    index_middle = indices[int(len(indices) // 2)]
    # Determine ordinal sets of features.
    table_ratio["quantile_mean_ratio"] = pandas.qcut(
        table_ratio["mean_ratio"],
        q=count_quantile,
        labels=indices,
    )
    table_ratio["quantile_mean_ratio_log"] = pandas.qcut(
        table_ratio["mean_ratio_log"],
        q=count_quantile,
        labels=indices,
    )
    # Filter table to rows of features in the middle quantile.
    table_middle = table_ratio.loc[
        (table_ratio["quantile_mean_ratio"] == index_middle), :
    ]
    # Determine counts and percentages.
    count_rows_source = table.shape[0]
    count_rows_product = table_middle.shape[0]
    percentage = str(round(
        float(100 * (count_rows_product / count_rows_source)), 2
    ))
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("Report:")
        print("partner")
        print("scale")
        print(
            "filter_table_features_least_change_across_observations()"
        )
        putly.print_terminal_partition(level=4)
        print("Count of quantiles: " + str(count_quantile))
        print("Quantile indices:")
        print(indices)
        print("Middle index: " + str(index_middle))
        putly.print_terminal_partition(level=4)
        print(
            "Count of rows in original, source table: " +
            str(count_rows_source)
        )
        print(
            "Count of rows in novel, product table: " +
            str(count_rows_product) + " (" + percentage + " %)"
        )
        print("Table of middle features that demonstrate least change")
        print(table_middle)
    # Return information.
    return table_middle


def scale_feature_values_between_observations_by_deseq(
    table=None,
    name_columns=None,
    name_rows=None,
    report=None,
):
    """
    Scales values of features between observations by their
    observation-specific median ratio to their geometric mean value across all
    observations.

    This scaling method is defined for comparisons of very many features,
    such as thousands to tens of thousands of genes or proteins, between
    comparatively few observations or samples. The goal of this method is to
    remove the influence of technical confounders so that comparisons between
    observations are more relevant to real differences, such as those due to
    biological regulation of expression of genes and proteins. An important
    assumption of this scaling method is that between observations the vast
    majority of features, such as genes or proteins, will have similar values.

    For this method, it is necessary to exclude any features that have a
    geometric mean of zero across samples.

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

    Review: 13 September 2024

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

    # Calculate median value of feature ratios for each observation.
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
                (row[column] / median_ratios[str(column + "_ratio")])
                if median_ratios[str(column + "_ratio")] != 0.0 else pandas.NA,
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
        print("Report:")
        print("partner")
        print("scale")
        print("scale_feature_values_between_observations_by_deseq()")
        putly.print_terminal_partition(level=4)
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


def describe_variance_across_features_with_least_change(
    table=None,
    name_columns=None,
    name_rows=None,
    count_quantile=None,
    report=None,
):
    """
    Describes the variance of values across features with the least change as
    calculated from the ratio of value for each individual observation relative
    to the geometric mean of values across all observations.

    This function is a handy companion to the function below.
    partner.scale.scale_feature_values_between_observations_by_deseq()

    Tables' format and orientation

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

    Review: 2 July 2024

    arguments:
        table (object): Pandas data-frame table of values for observations
            across columns and for features across rows
        name_columns (str): name of single-level index across columns
        name_rows (str): name of single-level index across rows
        count_quantile (int): odd count of quantiles, such as three for
            tertiles
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Calculate geometric mean of values for each feature across all
    # observations.
    # Calculate ratio of each observation for each feature to that feature's
    # geometric mean value across all observations.
    # Calculate mean of ratio values across all observations for each feature.
    # Filter rows in table to select features that demonstrate least change
    # across observations on the basis of the middle quantile.
    table_middle = (
        filter_table_features_least_change_across_observations(
            table=table,
            name_columns=name_columns,
            name_rows=name_rows,
            count_quantile=count_quantile,
            report=report,
    ))
    # Describe the distribution of values for these sets of middle features
    # that change the least.
    values_middle = table_middle["mean_ratio_log"].to_numpy(
        dtype="float64",
        na_value=numpy.nan,
        copy=True,
    )
    mean_middle = float(numpy.nanmean(values_middle))
    error_middle = float(scipy.stats.sem(values_middle))
    #mean_middle = round(float(numpy.nanmean(values_middle)), 7)
    #error_middle = round(float(scipy.stats.sem(values_middle)), 7)
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("Report:")
        print("partner")
        print("scale")
        print(
            "describe_variance_across_features_with_least_change()"
        )
        putly.print_terminal_partition(level=4)
        print("Variance in the changeability across features with the least")
        print("change across observations")
        print("Mean across middle features: " + str(mean_middle))
        print("Error across middle features: " + str(error_middle))
    # Collect information.
    pail = dict()
    pail["mean"] = mean_middle
    pail["error"] = error_middle
    # Return information.
    return pail


def compare_middle_quantile_feature_sets_by_ratio_to_geometric_mean(
    table_first=None,
    table_second=None,
    name_columns=None,
    name_rows=None,
    count_quantile=None,
    report=None,
):
    """
    Compares the set overlap of features from the middle quantile of their
    ratios relative to the geometric mean across all samples. These features
    demonstrated the least change. The two separate tables could have been
    stratified on the basis of control and intervention samples from the same
    experiment. The goal of this comparison is to evaluate the assumption that
    most features do not change between the experimental groups.

    This function is a handy companion to the function below.
    partner.scale.scale_feature_values_between_observations_by_deseq()

    Tables' format and orientation

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

    Review: 2 July 2024

    arguments:
        table_first (object): Pandas data-frame table of values for
            observations across columns and for features across rows
        table_second (object): Pandas data-frame table of values for
            observations across columns and for features across rows
        name_columns (str): name of single-level index across columns
        name_rows (str): name of single-level index across rows
        count_quantile (int): odd count of quantiles, such as three for
            tertiles
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers of features in the middle tertiles of both
            tables

    """

    # Copy information in table.
    table_first = table_first.copy(deep=True)
    table_second = table_second.copy(deep=True)
    # Calculate geometric mean of values for each feature across all
    # observations.
    # Calculate ratio of each observation for each feature to that feature's
    # geometric mean value across all observations.
    # Calculate mean of ratio values across all observations for each feature.
    # Filter rows in table to select features that demonstrate least change
    # across observations on the basis of the middle quantile.
    table_first_middle = (
        filter_table_features_least_change_across_observations(
            table=table_first,
            name_columns=name_columns,
            name_rows=name_rows,
            count_quantile=count_quantile,
            report=False,
    ))
    table_second_middle = (
        filter_table_features_least_change_across_observations(
            table=table_second,
            name_columns=name_columns,
            name_rows=name_rows,
            count_quantile=count_quantile,
            report=False,
    ))
    # Extract identifiers of features from the middle quantile of each table.
    features_first = copy.deepcopy(
        table_first.index.get_level_values(
            name_rows
        ).unique().to_list()
    )
    features_first_middle = copy.deepcopy(
        table_first_middle.index.get_level_values(
            name_rows
        ).unique().to_list()
    )
    features_second = copy.deepcopy(
        table_second.index.get_level_values(
            name_rows
        ).unique().to_list()
    )
    features_second_middle = copy.deepcopy(
        table_second_middle.index.get_level_values(
            name_rows
        ).unique().to_list()
    )
    # Determine overlap or union from complementary sets.
    list_sets_original = [
        set(features_first),
        set(features_second),
    ]
    list_sets_middle = [
        set(features_first_middle),
        set(features_second_middle),
    ]
    #union_original = set()
    #union_original = set(features_first).union(set(features_second))
    union_original = functools.reduce(
        lambda first, second: first.union(second),
        list_sets_original,
    )
    union_middle = functools.reduce(
        lambda first, second: first.union(second),
        list_sets_middle,
    )
    # Determine counts of features from each group.
    count_first = len(features_first)
    count_first_middle = len(features_first_middle)
    count_second = len(features_second)
    count_second_middle = len(features_second_middle)
    count_union_original = len(union_original)
    count_union_middle = len(union_middle)
    # Determine percentages.
    percentage_first = str(round(
        float(100 * (count_first_middle / count_first)), 2
    ))
    percentage_second = str(round(
        float(100 * (count_second_middle / count_second)), 2
    ))
    # Determine proportion of similarity between the middle sets of features.
    count_union_minimum = int(min([count_first_middle, count_second_middle]))
    count_union_maximum = int(count_first_middle + count_second_middle)
    range_union_middle = int(count_union_maximum - count_union_minimum)
    scale_union_middle = int(count_union_middle - count_union_minimum)
    proportion_union_middle = float(
        scale_union_middle / range_union_middle
    )
    percentage_union_middle = str(round(
        float(100 * proportion_union_middle), 2
    ))
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("Report:")
        print("partner")
        print("scale")
        print(
            "compare_middle_tertile_feature_sets_by_" +
            "ratio_to_geometric_mean()"
        )
        putly.print_terminal_partition(level=4)
        print("Counts of unique features in original tables")
        print("First table: " + str(count_first))
        print("Second table: " + str(count_second))
        print("Union: " + str(count_union_original))
        putly.print_terminal_partition(level=4)
        print("Counts of unique features in middle quantile of tables")
        print(
            "First table: " + str(count_first_middle) +
            " (" + percentage_first + " %)"
        )
        print(
            "Second table: " + str(count_second_middle) +
            " (" + percentage_second + " %)"
        )
        print("Union: " + str(count_union_middle))
        putly.print_terminal_partition(level=4)
        print(
            "Expectable range of unique union features from middle quantiles"
        )
        print(
            "Count of minimum possible unique union features: " +
            str(count_union_minimum)
        )
        print("This count would correspond to 100% similarity.")
        print(
            "Count of maximum possible unique union features: " +
            str(count_union_maximum)
        )
        print("This count would correspond to 0% similarity.")
        putly.print_terminal_partition(level=4)
        print(
            "Percentage of similarity in the features from middle quantiles"
        )
        print("Percentage: " + percentage_union_middle + " %")
        putly.print_terminal_partition(level=4)
    # Collect information.
    pail = dict()
    pail["union_original"] = union_original
    pail["union_middle"] = union_middle
    # Return information.
    return pail


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
# Procedure
# Currently, this module is not executable.
