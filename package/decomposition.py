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
    if (0.1 <= threshold_column_relative_variance)
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
    pail["table"] = table
    pail["table_scale"] = table_scale
    pail["index"] = index
    pail["count_index"] = count_index
    pail["matrix"] = matrix
    pail["count_samples"] = count_samples
    pail["count_variables"] = count_variables
    # Return.
    return pail


def calculate_singular_value_decomposition_factors(
    threshold_valid_proportion_per_column=None,
    table=None,
    report=None,
):
    """
    Calculates the initial, raw factors of a Singular Value Decomposition (SVD).

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
    u: unitary matrix with left singular vectors as columns
    - - shape of matrix "u": (m, m) or (m, k)
    s: singular values
    - - shape of matrix "s": (k,)
    vh: unitary matrix with right singular vectors as rows
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

    if False:
        # Calculate Singular Value Decomposition (SVD).
        u, s, vh = scipy.linalg.svd(
            pail_organization["matrix"],
            full_matrices=False, # Full matrices do not convey more information.
            compute_uv=True,
            overwrite_a=False,
            check_finite=True,
            lapack_driver="gesdd",
        )
        # Calculate the original data matrix as a product of the SVD factors.
        s_diagonal = numpy.diag(s)
        matrix_product = numpy.dot(u, numpy.dot(s_diagonal, vh))

        # Report.
        if report:
            utility.print_terminal_partition(level=2)
            print(
                "Report from: " +
                "calculate_singular_value_decomposition_factors()"
            )
            utility.print_terminal_partition(level=2)

            # Original matrix has shape (M, N)
            print(
                "Shape of original matrix: " +
                str(pail_organization["matrix"].shape)
            )
            print("Shape of product matrix: " + str(matrix_product.shape))
            print("rows (dimension 0): samples (cases, observations)")
            print("columns (dimension 1): variables (features)")
            # M: count of samples (cases)
            # N: count of variables (features)
            # K: minimum of M or N
            # Matrix "u" has shape (M, K)
            print("Shape of matrix U (left singular vectors): " + str(u.shape))
            # Matrix "s" has shape (K, )
            # Matrix "s" is basically a one-dimensional array.
            print("Shape of matrix S (singular values): " + str(s.shape))
            # Matrix "vt" has shape (K, N)
            print(
                "Shape of matrix VT (transpose right singular vectors): " +
                str(vt.shape)
            )
            # Compare original matrix to matrix calculation from SVD factors.
            print("Compare original matrix to product of SVD factors: ")
            print(numpy.allclose(
                pail_organization["matrix"], matrix_product,
                rtol=1e-2,
                atol=1e-3,
                equal_nan=False,
            ))

            pass
        # Compile information.
        pail = dict()
        pail["table_valid_scale"] = pail_organization["table_valid_scale"]
        pail["index"] = pail_organization["index"]
        pail["matrix"] = pail_organization["matrix"]
        pail["count_samples"] = pail_organization["count_samples"]
        pail["count_variables"] = pail_organization["count_variables"]
        pail["left_singular_vectors_columns"] = u
        pail["singular_values"] = s
        pail["right_singular_vectors_rows"] = vt
        # Return.
        return pail
    pass




#
