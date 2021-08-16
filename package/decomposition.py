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

#dir()
#importlib.reload()

###############################################################################
# Functionality


# Singular Value Decomposition and Principal Components Analysis


def copy_organize_table_matrix_for_singular_value_decomposition(
    threshold_valid_proportion_per_column=None,
    table=None,
    report=None,
):
    """
    Organizes a table and relevant information for Singular Value Decomposition
    (SVD).

    arguments:
        threshold_valid_proportion_per_column (float): proportion of rows that
            must have a valid value for a column in order to keep the column
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information for the singular value decomposition

    """

    # Copy information.
    table = table.copy(deep=True)
    # Drop any columns with inadequate valid values across rows.
    rows = table.shape[0]
    threshold = round(rows*threshold_valid_proportion_per_column)
    table.dropna(
        axis="columns",
        thresh=threshold,
        subset=None,
        inplace=True,
    )
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
    # Organize matrix.
    # Matrix format has samples (cases, observations) across rows (dimension 0)
    # and variables (features) across columns (dimension 1).
    #matrix = numpy.transpose(table_scale.to_numpy())
    matrix = numpy.copy(table_scale.to_numpy())

    # Compile information.
    pail = dict()
    pail["table_valid"] = table
    pail["table_valid_scale"] = table_scale
    pail["index"] = index
    pail["matrix"] = matrix
    pail["count_samples"] = copy.deepcopy(matrix.shape[0])
    pail["count_variables"] = copy.deepcopy(matrix.shape[1])
    # Return.
    return pail

# TODO: (TCW 16 August 2021) I think this function is sort of complete...
def calculate_initial_raw_singular_value_decomposition_factors(
    threshold_valid_proportion_per_column=None,
    table=None,
    report=None,
):
    """
    Calculates the initial, raw factors of a Singular Value Decomposition (SVD).

    arguments:
        threshold_valid_proportion_per_column (float): proportion of rows that
            must have a valid value for a column in order to keep the column
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the singular value
            decomposition

    """

    # Organize information.
    # Matrix format has samples (cases, observations) across rows (dimension 0)
    # and variables (features) across columns (dimension 1).
    pail_organization = (
        copy_organize_table_matrix_for_singular_value_decomposition(
            threshold_valid_proportion_per_column=(
                threshold_valid_proportion_per_column
            ),
            table=table,
            report=False,
    ))
    # Organize original matrix.

    # Calculate Singular Value Decomposition (SVD).
    # u: unitary matrix with left singular vectors as columns
    # s: singular values
    # vt: unitary matrix with right singular vectors as rows
    u, s, vt = scipy.linalg.svd(
        pail_organization["matrix"],
        full_matrices=False, # Full matrices do not convey more information.
        compute_uv=True,
        overwrite_a=False,
        check_finite=True,
        lapack_driver="gesdd",
    )
    # Calculate the original data matrix as a product of the SVD factors.
    s_diagonal = numpy.diag(s)
    matrix_product = numpy.dot(u, numpy.dot(s_diagonal, vt))

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_initial_raw_singular_value_decomposition_factors()"
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
