"""
Supply functionality for Singular Value Decomposition and Principal Components.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of Promiscuity
    (https://github.com/tcameronwaller/promiscuity/).

    Promiscuity supports data analysis in multiple other projects.
    Copyright (C) 2021 Thomas Cameron Waller

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
import scipy
import numpy
import statsmodels.api

# Custom

import promiscuity.utility as utility


#dir()
#importlib.reload()

###############################################################################
# Functionality


# Singular Value Decomposition and Principal Components Analysis


def organize_table_matrix_for_singular_value_decomposition(
    threshold_valid_proportion_per_column=None,
    threshold_column_relative_variance=None,
    table=None,
    report=None,
):
    """
    Organizes a table and relevant information for Singular Value Decomposition
    (SVD).

    Notice that "count_samples" represents the count of samples in the matrix
    after thresholding. This is the count of samples for the actual calculation
    of the Singular Value Decomposition.

    arguments:
        threshold_valid_proportion_per_column (float): minimal proportion of
            a column's rows that must have a valid value
        threshold_column_relative_variance (float): minimal relative variance in
            each column's values
        table (object): Pandas data frame with variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information for the singular value decomposition

    """

    def match_column_variance(
        name=None,
        threshold_column_relative_variance=None,
        variances=None,
    ):
        if (str(name) in variances.keys()):
            variance = variances[name]
            if (
                (not pandas.isna(variance)) and
                (threshold_column_relative_variance <= variance)
            ):
                match = True
            else:
                match = False
        else:
            match = False
        return match

    # Copy information.
    table = table.copy(deep=True)
    # Drop any columns with inadequate valid values across rows.
    if (0.1 <= threshold_valid_proportion_per_column):
        rows = table.shape[0]
        threshold = round(rows*threshold_valid_proportion_per_column)
        table.dropna(
            axis="columns",
            thresh=threshold,
            subset=None,
            inplace=True,
        )
    # Drop any columns with minimal relative variance.
    if (0.1 <= threshold_column_relative_variance):
        series_relative_variance = table.aggregate(
            lambda column: utility.calculate_relative_variance(
                array=column.to_numpy()
            ),
            axis="index", # apply function to each column
        )
        variances = series_relative_variance.to_dict()
        columns = copy.deepcopy(table.columns.to_list())
        columns_variance = list(filter(
            lambda column_trial: match_column_variance(
                name=column_trial,
                threshold_column_relative_variance=(
                    threshold_column_relative_variance
                ),
                variances=variances,
            ),
            columns
        ))
        table = table.loc[
            :, table.columns.isin(columns_variance)
        ]
    # Drop any rows with null values in any columns.
    table.dropna(
        axis="index",
        how="any",
        subset=None,
        inplace=True,
    )

    # Principal components analysis assumptions require at least centering the
    # means (mean = 0) of variables (features).
    # Standardizing the scale of variables (features) is equivalent to
    # calculation on correlation matrix instead of covariance matrix.
    # Standardize scale across variables (features) to mean zero (mean = 0) and
    # standard deviation one (standard deviation = 1).
    table_scale = utility.standardize_table_values_by_column(
        table=table,
        report=report,
    )
    # Copy information.
    index = copy.deepcopy(table_scale.index)
    count_index = len(index.to_list())
    # Organize matrix.
    # Matrix format has samples (cases, observations) across rows (dimension 0)
    # and variables (features) across columns (dimension 1).
    #matrix = numpy.transpose(table_scale.to_numpy())
    matrix = numpy.copy(table_scale.to_numpy())
    count_samples = copy.deepcopy(matrix.shape[0])
    count_variables = copy.deepcopy(matrix.shape[1])
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "organize_table_matrix_for_singular_value_decomposition()"
        )
        utility.print_terminal_partition(level=3)
        print("count matrix samples (dimension 0): " + str(count_samples))
        print("count matrix variables (dimension 1): " + str(count_variables))
        print("count index: " + str(count_index))
    # Compile information.
    pail = dict()
    pail["table_threshold"] = table
    pail["table_threshold_scale"] = table_scale
    pail["index"] = index
    pail["count_index"] = count_index
    pail["matrix"] = matrix
    pail["count_samples"] = count_samples
    pail["count_variables"] = count_variables
    # Return.
    return pail


def sort_decomposition_factors_by_decreasing_singular_values(
    u=None,
    s=None,
    vt=None,
    report=None,
):
    """
    Sorts factors from Singular Value Decomposition (SVD) in order of decreasing
    singular values.

    Sort dimensions of left (u) and right (vt) singular vectors that correspond to the
    singular values (s).

    Sort dimension 1 of matrix "u".
    Sort dimension 0 of matrix "vt".

    arguments:
        u (object): Numpy matrix with left singular vectors as columns
            (dimension 1)
        s (object): Numpy matrix of singular values
        vt (object): Numpy matrix with right singular vectors as rows
            (dimension 0)
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of factors from Singular Value Decomposition (SVD)

    """

    # Copy information.
    u = numpy.copy(u)
    s = numpy.copy(s)
    vt = numpy.copy(vt)
    # Calculate sort indices for Singular Values in decreasing order.
    # Reverse an increasing sort order with either "numpy.flip" or "[::-1]".
    indices_sort_increasing = numpy.argsort(
        s,
        axis=-1,
        kind="stable",
    )
    indices_sort = numpy.flip(
        indices_sort_increasing,
        axis=0,
    )
    # Apply sort order indices.
    u_sort = u[:,indices_sort]
    s_sort = s[indices_sort]
    vt_sort = vt[indices_sort,:]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "sort_decomposition_factors_by_decreasing_singular_values()"
        )
        utility.print_terminal_partition(level=3)
        print("singular values before sort: ")
        print(s)
        utility.print_terminal_partition(level=5)
        print("singular values after sort: ")
        print(s_sort)
    # Compile information.
    pail = dict()
    pail["u"] = u_sort
    pail["s"] = s_sort
    pail["vt"] = vt_sort
    # Return.
    return pail


def calculate_sort_singular_value_decomposition_factors(
    threshold_valid_proportion_per_column=None,
    threshold_column_relative_variance=None,
    table=None,
    report=None,
):
    """
    Calculates the initial, raw factors of a Singular Value Decomposition (SVD).

    Sorts decomposition factors by decreasing singular values.

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    Matrix format: NumPy matrix with samples (cases, observations) across rows
    (dimension 0) and variables (features) across columns (dimension 1)

    Factorize covariance or correlation into direction (eigenvectors) and scale
    (eigenvalues).
    The eigenvalues impart scale or weight to each eigenvector.
    Each eigenvector has its own eigenvalue, and their sort orders mush match.

    loadings = eigenvectors [dot] square_root(eigenvalues)
    Loadings include aspects of both direction (eigenvectors) and scale
    (eigenvalues).

    Singular Value Decomposition assigns Eigenvector direction (sign, positive
    or negative) at random.

    Terminology

    a: original data matrix
    - - shape of matrix "a": (m, n)
    a = u [dot] s_diagonal [dot] vh

    m: count of samples (cases)
    n: count of variables (features)
    k = min(m, n)
    u: unitary matrix with left singular vectors as columns (dimension 1)
    - - shape of matrix "u": (m, m) or (m, k)
    s: singular values
    - - shape of matrix "s": (k,)
    vt: unitary matrix with right singular vectors as rows (dimension 0)
    - - shape of matrix "vh": (n, n) or (k, n)

    arguments:
        threshold_valid_proportion_per_column (float): minimal proportion of
            a column's rows that must have a valid value
        threshold_column_relative_variance (float): minimal relative variance in
            each column's values
        table (object): Pandas data frame with variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the singular value
            decomposition

    """

    # Organize matrix.
    pail_organization = (
        organize_table_matrix_for_singular_value_decomposition(
            threshold_valid_proportion_per_column=(
                threshold_valid_proportion_per_column
            ),
            threshold_column_relative_variance=(
                threshold_column_relative_variance
            ),
            table=table,
            report=report,
    ))

    # Calculate Singular Value Decomposition (SVD).
    u, s, vt = scipy.linalg.svd(
        numpy.copy(pail_organization["matrix"]),
        full_matrices=False, # Full matrices do not convey more information.
        compute_uv=True,
        overwrite_a=False,
        check_finite=True,
        lapack_driver="gesdd",
    )
    # Sort SVD factors in order of decreasing singular values.
    pail_sort = sort_decomposition_factors_by_decreasing_singular_values(
        u=u,
        s=s,
        vt=vt,
        report=report,
    )
    # Calculate the original data matrix as a product of the SVD factors.
    s_diagonal = numpy.diag(s)
    matrix_product = numpy.dot(u, numpy.dot(s_diagonal, vt))
    s_sort_diagonal = numpy.diag(pail_sort["s"])
    matrix_product_sort = numpy.dot(
        pail_sort["u"], numpy.dot(s_sort_diagonal, pail_sort["vt"])
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_singular_value_decomposition_factors()"
        )
        utility.print_terminal_partition(level=4)
        # Compare original matrix to matrix calculation from SVD factors.
        print(
            "Shape of original matrix: " +
            str(pail_organization["matrix"].shape)
        )
        print("Shape of product matrix: " + str(matrix_product.shape))
        print(
            "Shape of sorted product matrix: " + str(matrix_product_sort.shape)
        )
        print("rows (dimension 0): samples (cases, observations)")
        print("columns (dimension 1): variables (features)")
        print("Compare original matrix to product of SVD factors: ")
        print(numpy.allclose(
            pail_organization["matrix"], matrix_product,
            rtol=1e-2,
            atol=1e-3,
            equal_nan=False,
        ))
        print("Compare original matrix to product of sorted SVD factors: ")
        print(numpy.allclose(
            pail_organization["matrix"], matrix_product_sort,
            rtol=1e-2,
            atol=1e-3,
            equal_nan=False,
        ))
        # Describe SVD factor matrices.
        utility.print_terminal_partition(level=4)
        print(
            "Shape of sorted matrix U (left singular vectors): " +
            str(pail_sort["u"].shape)
        )
        print(
            "Shape of sorted matrix S (singular values): " +
            str(pail_sort["s"].shape)
        )
        print(
            "Shape of sorted matrix Vt (transpose right singular vectors): " +
            str(pail_sort["vt"].shape)
        )
        pass
    # Compile information.
    pail = dict()
    pail["table_threshold"] = pail_organization["table_threshold"]
    pail["table_threshold_scale"] = pail_organization["table_threshold_scale"]
    pail["index"] = pail_organization["index"]
    pail["matrix"] = pail_organization["matrix"]
    pail["count_samples"] = pail_organization["count_samples"]
    pail["count_variables"] = pail_organization["count_variables"]
    pail["u_left_singular_vectors_columns"] = pail_sort["u"]
    pail["s_singular_values"] = pail_sort["s"]
    pail["vt_right_singular_vectors_rows"] = pail_sort["vt"]
    # Return.
    return pail


def calculate_principal_component_eigenvalues_from_singular_values(
    singular_values=None,
    count_samples=None,
    report=None,
):
    """
    Calculates Principal Components Analysis (PCA) Eigenvalues from Singular
    Values of Singular Value Decomposition (SVD).

    s = singular value from SVD on original data matrix "a"
    n = count of samples in original data matrix "a"
    l = eigenvalue of covariance matrix of original data matrix "a"

    l = ( (s^2) / (n - 1))

    References:
    "https://towardsdatascience.com/singular-value-decomposition-and-its-"
    + "applications-in-principal-component-analysis-5b7a5f08d0bd"

    arguments:
        singular_values (object): NumPy array of Singular Values
        count_samples (float): count of samples in the original Singular Value
            Decomposition
        report (bool): whether to print reports

    raises:

    returns:
        (object): NumPy array of Eigenvalues

    """

    def divide_by_sample_count(value, count_samples):
        return (value / (count_samples - 1))
    array_divide_by_sample_count = numpy.vectorize(divide_by_sample_count)

    def square_divide_by_sample_count(value, count_samples):
        return ((math.pow(value, 2)) / (count_samples - 1))
    array_square_divide_by_sample_count = numpy.vectorize(
        square_divide_by_sample_count
    )

    # Copy information.
    singular_values = numpy.copy(singular_values)
    # Calculate Eigenvalues first method.
    singular_values_square = numpy.square(singular_values)
    eigenvalues_first = array_divide_by_sample_count(
        singular_values_square, count_samples
    )
    # Calculate Eigenvalues second method.
    eigenvalues_second = array_square_divide_by_sample_count(
        singular_values, count_samples
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_principal_component_eigenvalues_from_singular_values()"
        )
        utility.print_terminal_partition(level=2)
        print("Singular values...")
        print(
            "- - - shape: " + str(singular_values.shape)
        )
        print(singular_values)
        print("Eigenvalues first method...")
        print(
            "- - - shape: " + str(eigenvalues_first.shape)
        )
        print(eigenvalues_first)
        print("Eigenvalues second method...")
        print(
            "- - - shape: " + str(eigenvalues_second.shape)
        )
        print(eigenvalues_second)
    # Return.
    return eigenvalues_second


def calculate_principal_component_explanation_variance_proportions(
    eigenvalues=None,
    report=None,
):
    """
    Calculates the proportion of variance explained by Eigenvectors of Principal
    Components Analysis (PCA).

    Sum of proportional variance explained across all Eigenvectors is one.

    arguments:
        eigenvalues (object): NumPy array of Eigenvalues
        report (bool): whether to print reports

    raises:

    returns:
        (object): NumPy array of proportions of variance explained

    """

    def divide_by_total(value, total):
        return (value / (total))
    array_divide_by_total = numpy.vectorize(divide_by_total)

    # Copy information.
    eigenvalues = numpy.copy(eigenvalues)
    # Calculate total variance across all Eigenvectors.
    variance_total = numpy.sum(eigenvalues)
    # Calculate proportional variance across Eigenvectors.
    variance_proportions = array_divide_by_total(
        eigenvalues, variance_total
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_principal_component_explanation_variance_proportions()"
        )
        utility.print_terminal_partition(level=2)
        print("Eigenvalues...")
        print(eigenvalues)
        print("Variance proportions...")
        print(variance_proportions)
    # Return.
    return variance_proportions


def organize_principal_component_variance_proportion_table(
    variance_proportions=None,
    prefix=None,
    index_name=None,
    report=None,
):
    """
    Organizes a table of proportion of variance explained by each Eigenvector
    and Principal Component factor.

    arguments:
        variance_proportions (object): NumPy array of proportions of variance
            explained
        prefix (str): prefix for names of component columns
        index_name (str): name for table's index column
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    # Copy information.
    utility.print_terminal_partition(level=2)
    print(
        "Report from: " +
        "organize_principal_component_variance_proportion_table()"
    )
    print(variance_proportions)
    # .flatten(order="C").tolist()
    variance_proportions = numpy.copy(variance_proportions)
    variance_proportions_list = variance_proportions.tolist()
    print(variance_proportions_list)
    # Organize information.
    count = 1
    records = list()
    for variance_proportion in variance_proportions:
        record = dict()
        record[index_name] = str(prefix + str(count))
        record["variance_proportion"] = variance_proportion
        records.append(record)
        count += 1
    table = pandas.DataFrame(data=records)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "organize_principal_component_variance_proportion_table()"
        )
        utility.print_terminal_partition(level=3)
        print("Table after organization:")
        print(table)
    # Return.
    return table





# TODO: make this function more of a driver
# 1. calculate SVD factors


def organize_principal_components_by_singular_value_decomposition(
    table=None,
    report=None,
):
    """
    Organizes a Principal Components Analysis (PCA) by Singular Value
    Decomposition (SVD).

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    Relevant dimension: Principal Components represent variance across features

    Reference:
    "https://stats.stackexchange.com/questions/134282/
    relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca"

    arguments:
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the Principal Components
            Analysis (PCA)

    """

    # Threshold and organize original matrix.
    # Calculate and sort Singular Value Decomposition (SVD) factors.
    pail_decomposition = (
        calculate_sort_singular_value_decomposition_factors(
            threshold_valid_proportion_per_column=0.5,
            threshold_column_relative_variance=0.5,
            table=table,
            report=report,
        )
    )
    # Derive factors and summaries for Principal Component Analysis (PCA).
    eigenvalues = (
        calculate_principal_component_eigenvalues_from_singular_values(
            singular_values=pail_decomposition["s_singular_values"],
            count_samples=pail_decomposition["count_samples"],
            report=report,
    ))
    variance_proportions = (
        calculate_principal_component_explanation_variance_proportions(
            eigenvalues=eigenvalues,
            report=report,
    ))
    table_component_variance_proportions = (
        organize_principal_component_variance_proportion_table(
            variance_proportions=variance_proportions,
            prefix="component_",
            index_name="eigenvectors_components",
            report=report,
    ))


    pass






#
