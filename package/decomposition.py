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
import sklearn.decomposition
import scipy
import numpy
import statsmodels.api

# Custom

import promiscuity.utility as utility # this import path for subpackage


#dir()
#importlib.reload()

###############################################################################
# Functionality


# Singular Value Decomposition and Principal Components Analysis

# Note: TCW, 26 January 2022
# Principal Component scores and proportional explanation of variance match
# between SKLearn and Singular Value Decomposition (SVD) methods.

# Note: TCW, 17 August 2021
# The SVD method delivers Principal Components that match the SKLearn method.
# I need to confirm that both methods deliver identical loadings.


def organize_table_matrix_for_decomposition(
    table=None,
    report=None,
):
    """
    Organizes a table and relevant information for Singular Value Decomposition
    (SVD).

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    Notice that "count_samples" represents the count of samples in the matrix
    after thresholding. This is the count of samples for the actual calculation
    of the Singular Value Decomposition.

    arguments:
        table (object): Pandas data frame with variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information for the singular value decomposition

    """

    # Copy information.
    table = table.copy(deep=True)

    # Remove columns with inadequate non-missing values or inadequate variance
    # (standard deviation) across rows.
    # This preliminary filter on columns avoids losing table rows for missing
    # values in an unnecessary column.
    table = utility.filter_table_columns_by_nonmissing_variance(
        threshold_valid_proportion_per_column=0.05,
        threshold_column_variance=0.001,
        type_variance="standard_deviation",
        table=table,
        report=report,
    )
    # Drop any rows with missing values in any columns.
    table.dropna(
        axis="index", # drop rows
        how="any",
        subset=None,
        inplace=True,
    )
    # Remove columns with inadequate non-missing values or inadequate variance
    # (standard deviation) across rows.
    # Filter columns again since the removal of missing values changes the
    # variance of columns.
    table = utility.filter_table_columns_by_nonmissing_variance(
        threshold_valid_proportion_per_column=0.05,
        threshold_column_variance=0.001,
        type_variance="standard_deviation",
        table=table,
        report=report,
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
    table_scale.dropna(
        axis="index",
        how="any",
        subset=None,
        inplace=True,
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
    matrix_source=None,
    report=None,
):
    """
    Calculates the initial, raw factors of a Singular Value Decomposition (SVD).

    Sorts decomposition factors by decreasing singular values.

    Matrix format: NumPy matrix with samples (cases, observations) across rows
    (dimension 0) and variables (features) across columns (dimension 1)

    Singular Value Decomposition assigns Eigenvector direction (sign, positive
    or negative) at random.

    Terminology

    m: count of samples (cases, observations) in original data matrix
    n: count of variables (features) in original data matrix
    k = min(m, n)

    A: original data matrix
    - - shape of matrix "A": (m, n)
    A = U <dot product> S_diagonal <dot product> Vt

    U: matrix with left singular vectors as columns (dimension 1)
    - - shape of matrix "U": (m, m) or (m, k)
    Ut: matrix with left singular vectors as rows (dimension 0)
    - - shape of matrix "Ut": (m, m) or (k, m)
    S: singular values
    - - shape of matrix "S": (k,)
    Vt: matrix with right singular vectors as rows (dimension 0)
    - - shape of matrix "Vt": (n, n) or (k, n)
    V: matrix with right singular vectors as columns (dimension 1)
    - - shape of matrix "Vt": (n, n) or (n, k)

    arguments:
        matrix_source (object): NumPy matrix with samples (cases, observations)
            across rows (dimension 0) and variables (features) across columns
            (dimension 1)
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the singular value
            decomposition

    """

    # Copy information.
    matrix_source = numpy.copy(matrix_source)

    # Calculate Singular Value Decomposition (SVD).
    u, s, vt = scipy.linalg.svd(
        matrix_source,
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
    ut_sort = numpy.copy(numpy.transpose(pail_sort["u"]))
    v_sort = numpy.copy(numpy.transpose(pail_sort["vt"]))
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
        print("Shape of original matrix: " + str(matrix_source.shape))
        print("Shape of product matrix: " + str(matrix_product.shape))
        print(
            "Shape of sorted product matrix: " + str(matrix_product_sort.shape)
        )
        print("rows (dimension 0): samples (cases, observations)")
        print("columns (dimension 1): variables (features)")
        print("Compare original matrix to product of SVD factors: ")
        print(numpy.allclose(
            matrix_source, matrix_product,
            rtol=1e-2,
            atol=1e-3,
            equal_nan=False,
        ))
        print("Compare original matrix to product of sorted SVD factors: ")
        print(numpy.allclose(
            matrix_source, matrix_product_sort,
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
            "Shape of sorted matrix Ut (transpose left singular vectors): " +
            str(ut_sort.shape)
        )
        print(
            "Shape of sorted matrix S (singular values): " +
            str(pail_sort["s"].shape)
        )
        print(
            "Shape of sorted matrix Vt (transpose right singular vectors): " +
            str(pail_sort["vt"].shape)
        )
        print(
            "Shape of sorted matrix V (right singular vectors): " +
            str(v_sort.shape)
        )
        pass
    # Compile information.
    pail = dict()
    pail["matrix_source"] = matrix_source
    pail["u_left_singular_vectors_columns"] = pail_sort["u"]
    pail["ut_left_singular_vectors_rows"] = ut_sort
    pail["s_singular_values"] = pail_sort["s"]
    pail["s_singular_values_diagonal"] = s_sort_diagonal
    pail["vt_right_singular_vectors_rows"] = pail_sort["vt"]
    pail["v_right_singular_vectors_columns"] = v_sort
    # Return.
    return pail


def calculate_principal_component_eigenvalues_from_singular_values(
    s_singular_values=None,
    count_samples=None,
    report=None,
):
    """
    Calculates Principal Components Analysis (PCA) Eigenvalues from Singular
    Values of Singular Value Decomposition (SVD).

    s = singular value from SVD on original source data matrix
    n = count of samples in original source data matrix
    l = eigenvalue of covariance matrix of original source data matrix

    l = ( (s^2) / (n - 1))

    arguments:
        s_singular_values (object): Numpy matrix of Singular Values
        count_samples (float): count of samples in the original source matrix
            for Singular Value Decomposition
        count_samples (float): count of samples in the original source matrix
            for Singular Value Decomposition
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
    singular_values = numpy.copy(s_singular_values)
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


def calculate_loadings_from_eigenvalues_eigenvectors(
    eigenvectors=None,
    eigenvalues=None,
    report=None,
):
    """
    Calculates Principal Components Analysis (PCA) loadings from Eigenvectors
    and Eigenvalues.

    Statsmodels erroneously returns "loadings" that have identical values and
    dimensions as the Eigenvectors; however, Eigenvectors and Loadings are not
    equivalent.

    loadings = eigenvectors [dot] square_root(eigenvalues)
    Loadings include aspects of both direction (eigenvectors) and scale
    (eigenvalues).

    arguments:
        eigenvectors (object): NumPy matrix of Eigenvectors
        eigenvalues (object): NumPy array of Eigenvalues
        report (bool): whether to print reports

    raises:

    returns:
        (object): Numpy array of loadings

    """

    # Copy information.
    eigenvectors = numpy.copy(eigenvectors)
    eigenvalues = numpy.copy(eigenvalues)
    # Calculate square roots of Eigenvalues.
    # Organize a diagonal matrix of square roots of Eigenvalues.
    eigenvalues_square_root = numpy.sqrt(eigenvalues)
    eigenvalues_root_diagonal = numpy.diag(eigenvalues_square_root)
    # Calculate loadings.
    loadings = numpy.dot(eigenvectors, eigenvalues_root_diagonal)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_loadings_from_eigenvalues_eigenvectors()"
        )
        utility.print_terminal_partition(level=3)
        print("Shape of loadings: " + str(loadings.shape))
        utility.print_terminal_partition(level=4)
        print("Loadings = Eigenvectors [dot] square_root(diagonal Eigenvalues)")
        print(loadings)
    # Return.
    return loadings


def calculate_loadings_from_decomposition_factors(
    s_singular_values=None,
    vt_right_singular_vectors_rows=None,
    count_samples=None,
    report=None,
):
    """
    Calculates Principal Components Analysis (PCA) loadings from direct factors
    of Singular Value Decomposition.

    arguments:
        s_singular_values (object): Numpy matrix of Singular Values
        vt_right_singular_vectors_rows (object): Numpy matrix
        count_samples (float): count of samples in the original source matrix
            for Singular Value Decomposition
        report (bool): whether to print reports

    raises:

    returns:
        (object): Numpy array of loadings

    """

    def divide_by_sample_count(value, count_samples):
        return (value / math.sqrt(count_samples - 1))
    array_divide_by_sample_count = numpy.vectorize(divide_by_sample_count)

    # Copy information.
    s = numpy.copy(s_singular_values)
    vt = numpy.copy(vt_right_singular_vectors_rows)
    # Calculate loadings.
    quotients = array_divide_by_sample_count(
        s, count_samples
    )
    quotients_diagonal = numpy.diag(quotients)
    loadings = numpy.dot(vt, quotients_diagonal)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_loadings_from_decomposition_factors()"
        )
        utility.print_terminal_partition(level=3)
        print("Shape of loadings: " + str(loadings.shape))
        utility.print_terminal_partition(level=4)
        print("Loadings = (V [dot] (S / square_root(samples - 1)))")
        print(loadings)
    # Return.
    return loadings


# TODO: TCW 17 August 2021
# TODO: organize table with labels of rows and columns
# TODO: feature X PC (which dimension is which?)

def calculate_principal_component_loadings(
    eigenvectors=None,
    eigenvalues=None,
    s_singular_values=None,
    vt_right_singular_vectors_rows=None,
    count_samples=None,
    report=None,
):
    """
    Calculates Principal Components Analysis (PCA) loadings.

    arguments:
        eigenvectors (object): NumPy matrix of Eigenvectors
        eigenvalues (object): NumPy array of Eigenvalues
        s_singular_values (object): Numpy matrix of Singular Values
        vt_right_singular_vectors_rows (object): Numpy matrix
        count_samples (float): count of samples in the original source matrix
            for Singular Value Decomposition
        report (bool): whether to print reports

    raises:

    returns:
        (object): Numpy array of loadings

    """

    # Eigenvectors and Eigenvalues.
    loadings_eigen = calculate_loadings_from_eigenvalues_eigenvectors(
        eigenvectors=eigenvectors,
        eigenvalues=eigenvalues,
        report=report,
    )
    # Raw decomposition factors.
    loadings_factor = calculate_loadings_from_decomposition_factors(
        s_singular_values=s_singular_values,
        vt_right_singular_vectors_rows=vt_right_singular_vectors_rows,
        count_samples=count_samples,
        report=report,
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_principal_component_loadings()"
        )
        utility.print_terminal_partition(level=3)
        print("Shape of Eigen loadings: " + str(loadings_eigen.shape))
        print("Shape of factor loadings: " + str(loadings_factor.shape))
        utility.print_terminal_partition(level=4)
        print("Compare loadings from both methods: ")
        print(numpy.allclose(
            loadings_eigen, loadings_factor,
            rtol=1e-2,
            atol=1e-3,
            equal_nan=False,
        ))
    # Return.
    return loadings_eigen


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
    variance_proportions = numpy.copy(variance_proportions).tolist()
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


def calculate_principal_component_scores_from_factors(
    matrix_source=None,
    u_left_singular_vectors_columns=None,
    s_singular_values_diagonal=None,
    v_right_singular_vectors_columns=None,
    report=None,
):
    """
    Calculates Principal Components Analysis (PCA) scores from factors of
    Singular Value Decomposition (SVD).

    A = U <dot> S_diagonal <dot> Vt
    A / Vt = U <dot> S_diagonal
    A <dot> V = U <dot> S_diagonal
    Principal Component Scores = A <dot> V
    Principal Component Scores = U <dot> S_diagonal

    arguments:
        matrix_source (object): NumPy matrix with samples (cases, observations)
            across rows (dimension 0) and variables (features) across columns
            (dimension 1)
        u_left_singular_vectors_columns (object): Numpy matrix
        s_singular_values_diagonal (object): Numpy matrix
        v_right_singular_vectors_columns (object): Numpy matrix
        report (bool): whether to print reports

    raises:

    returns:
        (object): NumPy matrix with samples (cases, observations) across rows
            (dimension 0) and variables (features) across columns (dimension 1)

    """

    # Copy information.
    a = numpy.copy(matrix_source)
    u = numpy.copy(u_left_singular_vectors_columns)
    s_diagonal = numpy.copy(s_singular_values_diagonal)
    v = numpy.copy(v_right_singular_vectors_columns)

    # Calculate the Principal Component Scores.
    matrix_scores_first = numpy.dot(u, s_diagonal)
    matrix_scores_second = numpy.dot(a, v)

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_principal_component_scores_from_factors()"
        )
        utility.print_terminal_partition(level=4)
        # Compare alternative calculations of the scores.
        print("Shape of source matrix: " + str(matrix_source.shape))
        print(
            "Shape of first scores matrix: " +
            str(matrix_scores_first.shape)
        )
        print(
            "Shape of second scores matrix: " +
            str(matrix_scores_second.shape)
        )
        print("Are the product matrices nearly identical?: ")
        print(numpy.allclose(
            matrix_scores_first, matrix_scores_second,
            rtol=1e-2, # relative tolerance, 1%
            atol=1e-3, # absolute tolerance
            equal_nan=False,
        ))
    # Compile information.
    pail = dict()
    pail["matrix_scores"] = matrix_scores_second
    # Return.
    return pail


def organize_principal_component_scores_table(
    matrix=None,
    index=None,
    index_name=None,
    prefix=None,
    separator=None,
    report=None,
):
    """
    Organizes a table with Principal Components Analysis (PCA) Scores.

    arguments:
        matrix (object): NumPy matrix of Principal Component Scores
        index (object): Pandas series of index from original table
        index_name (str): name of table's index column
        prefix (str): prefix for names of new principal component columns in
            table
        separator (str): separator for names of new columns
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of Principal Component Scores

    """

    # Copy information.
    matrix= numpy.copy(matrix)

    # Organize table.
    count_columns = matrix.shape[1]
    count = 1
    columns = list()
    for component in range(0, count_columns, 1):
        column = str(str(prefix) + str(separator) + str(count))
        columns.append(column)
        count += 1
    table = pandas.DataFrame(
        data=matrix,
        index=index,
        columns=columns,
        dtype="float32",
        copy=True,
    )
    table.rename_axis(
        index=index_name,
        axis="index",
        copy=False,
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "organize_principal_component_scores_table()"
        )
        utility.print_terminal_partition(level=4)
        # Compare alternative calculations of the scores.
        print("Shape of source matrix: " + str(matrix.shape))
        utility.print_terminal_partition(level=4)
        print("table...")
        print(table)
    # Return.
    return table


# TODO: TCW, 27 January 2022
# TODO: the threshold on relative variance might be too strict...


def organize_principal_components_by_singular_value_decomposition(
    table=None,
    index_name=None,
    prefix=None,
    separator=None,
    report=None,
):
    """
    Organizes a Principal Components Analysis (PCA) by Singular Value
    Decomposition (SVD).

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    Relevant dimension: Principal Components represent variance across features

    Principal Components factorize covariance or correlation into direction
    (Eigenvectors) and scale (Eigenvalues).
    The Eigenvalues impart scale or weight to each Eigenvector.
    Each Eigenvector has its own Eigenvalue, and their sort orders mush match.

    Terminology

    m: count of samples (cases, observations) in original data matrix
    n: count of variables (features) in original data matrix
    k = min(m, n)

    A: original data matrix
    - - shape of matrix "A": (m, n)
    A = U <dot product> S_diagonal <dot product> Vt

    U: matrix with left singular vectors as columns (dimension 1)
    - - shape of matrix "U": (m, m) or (m, k)
    Ut: matrix with left singular vectors as rows (dimension 0)
    - - shape of matrix "Ut": (m, m) or (k, m)
    S: singular values
    - - shape of matrix "S": (k,)
    Vt: matrix with right singular vectors as rows (dimension 0)
    - - shape of matrix "Vt": (n, n) or (k, n)
    V: matrix with right singular vectors as columns (dimension 1)
    - - shape of matrix "Vt": (n, n) or (n, k)

    Definitions

    Eigenvalue = ( (S^2) / (m - 1))
    Proportional Variance = Eigenvalue / sum( Eigenvalues )
    Eigenvectors = V columns (right singular vectors)
    - - specifically in this context
    - - left singular vectors can be Eigenvectors in another context
    - - some sources refer to V Eigenvectors as "loadings"
    Loadings include aspects of both direction (Eigenvectors) and scale
    (Eigenvalues).
    Loadings are most useful to understand the contribution of each variable
    (feature) to each Principal Component.
    Loadings matrix is generally square.
    Loadings = Eigenvectors <dot> (Eigenvalues)^0.5
    - - matrix multiplication involves rows of first matrix and columns of
    - - second
    - - Eigenvectors are Vt in this case so that singular fectors are across
    - - rows
    Loadings = Vt <dot> (S / ( (m - 1)^0.5 ))
    - - matrix multiplication involves rows of first matrix and columns of
    - - second
    - - use Vt in this case so that singular fectors are across rows

    A = U <dot> S_diagonal <dot> Vt
    A / Vt = U <dot> S_diagonal
    1 / Vt = V
    A <dot> V = U <dot> S_diagonal
    Principal Component Scores = A <dot> V
    Principal Component Scores = U <dot> S_diagonal

    Reference:

    "https://towardsdatascience.com/singular-value-decomposition-and-its-"
    + "applications-in-principal-component-analysis-5b7a5f08d0bd"

    "https://towardsdatascience.com/principal-component-analysis-for-"
    + "dimensionality-reduction-115a3d157bad

    "https://stats.stackexchange.com/questions/134282/relationship-between-svd-"
    + "and-pca-how-to-use-svd-to-perform-pca"

    "https://stats.stackexchange.com/questions/143905/loadings-vs-eigenvectors-"
    + "in-pca-when-to-use-one-or-another"

    "https://stats.stackexchange.com/questions/126885/methods-to-compute-"
    + "factor-scores-and-what-is-the-score-coefficient-matrix-in"

    arguments:
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        index_name (str): name of table's index column
        prefix (str): prefix for names of new principal component columns in
            table
        separator (str): separator for names of new columns
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the Principal Components
            Analysis (PCA)

    """

    # Threshold and organize original matrix.
    pail_organization = organize_table_matrix_for_decomposition(
        table=table,
        report=report,
    )

    # Calculate and sort Singular Value Decomposition (SVD) factors.
    pail_decomposition = (
        calculate_sort_singular_value_decomposition_factors(
            matrix_source=pail_organization["matrix"],
            report=report,
        )
    )

    # Derive factors and summaries for Principal Component Analysis (PCA).

    # Eigenvalues
    eigenvalues = (
        calculate_principal_component_eigenvalues_from_singular_values(
            s_singular_values=pail_decomposition["s_singular_values"],
            count_samples=pail_organization["count_samples"],
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

    # Loadings
    loadings = calculate_principal_component_loadings(
        eigenvectors=pail_decomposition["vt_right_singular_vectors_rows"],
        eigenvalues=eigenvalues,
        s_singular_values=pail_decomposition["s_singular_values"],
        vt_right_singular_vectors_rows=(
            pail_decomposition["vt_right_singular_vectors_rows"]
        ),
        count_samples=pail_organization["count_samples"],
        report=report,
    )

    # Principal Component Scores
    pail_scores = calculate_principal_component_scores_from_factors(
        matrix_source=pail_organization["matrix"],
        u_left_singular_vectors_columns=(
            pail_decomposition["u_left_singular_vectors_columns"]
        ),
        s_singular_values_diagonal=(
            pail_decomposition["s_singular_values_diagonal"]
        ),
        v_right_singular_vectors_columns=(
            pail_decomposition["v_right_singular_vectors_columns"]
        ),
        report=report,
    )
    table_component_scores = organize_principal_component_scores_table(
        matrix=pail_scores["matrix_scores"],
        index=pail_organization["index"],
        index_name=index_name,
        prefix=prefix,
        separator=separator,
        report=report,
    )

    # Compile information.
    pail = dict()
    pail["table_source_threshold"] = (pail_organization["table_threshold"])
    pail["table_source_threshold_scale"] = (
        pail_organization["table_threshold_scale"]
    )
    pail["count_samples"] = pail_organization["count_samples"]
    pail["count_variables"] = pail_organization["count_variables"]
    pail["matrix_component_scores"] = pail_scores["matrix_scores"]
    pail["table_component_scores"] = table_component_scores
    pail["table_component_variance_proportions"] = (
        table_component_variance_proportions
    )
    # Return.
    return pail


# SKLearn


def calculate_principal_components_by_sklearn(
    matrix_source=None,
    report=None,
):
    """
    Organizes a Principal Components Analysis (PCA) by SKLearn.

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    Relevant dimension: Principal Components represent variance across features

    arguments:
        matrix_source (object): NumPy matrix with samples (cases, observations)
            across rows (dimension 0) and variables (features) across columns
            (dimension 1)
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the Principal Components
            Analysis (PCA)

    """

    # Copy information.
    matrix_source = numpy.copy(matrix_source)
    # Execute principle component analysis.
    pca = sklearn.decomposition.PCA(n_components=None)
    #report = pca.fit_transform(matrix)
    pca.fit(matrix_source)
    # Report extent of variance that each principal component explains.
    variance_ratios = pca.explained_variance_ratio_
    component_numbers = range(len(variance_ratios) + 1)
    variance_series = {
        "components": component_numbers[1:],
        "variance": variance_ratios
    }
    table_component_variance_proportions = pandas.DataFrame(
        data=variance_series
    )
    # Transform data by principal components.
    matrix_scores = pca.transform(matrix_source)

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "calculate_principal_components_by_sklearn()"
        )
        utility.print_terminal_partition(level=4)
        # Compare matrices.
        print("Shape of source matrix: " + str(matrix_source.shape))
        print(
            "Shape of score matrix: " + str(matrix_scores.shape)
        )
        utility.print_terminal_partition(level=4)
        print("Proportional variance for each component...")
        print(table_component_variance_proportions)
    # Compile information.
    pail = dict()
    pail["table_component_variance_proportions"] = (
        table_component_variance_proportions
    )
    pail["matrix_scores"] = matrix_scores
    # Return.
    return pail


def organize_principal_components_by_sklearn(
    table=None,
    index_name=None,
    prefix=None,
    separator=None,
    report=None,
):
    """
    Organizes a Principal Components Analysis (PCA) by SKLearn.

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    Relevant dimension: Principal Components represent variance across features

    Principal Components factorize covariance or correlation into direction
    (Eigenvectors) and scale (Eigenvalues).
    The Eigenvalues impart scale or weight to each Eigenvector.
    Each Eigenvector has its own Eigenvalue, and their sort orders mush match.

    arguments:
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        index_name (str): name of table's index column
        prefix (str): prefix for names of new principal component columns in
            table
        separator (str): separator for names of new columns
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the Principal Components
            Analysis (PCA)

    """

    # Threshold and organize original matrix.
    pail_organization = organize_table_matrix_for_decomposition(
        table=table,
        report=report,
    )

    # Calculate Principal Component Scores in SKLearn.
    pail_sklearn = calculate_principal_components_by_sklearn(
        matrix_source=pail_organization["matrix"],
        report=report,
    )

    table_component_scores = organize_principal_component_scores_table(
        matrix=pail_sklearn["matrix_scores"],
        index=pail_organization["index"],
        index_name=index_name,
        prefix=prefix,
        separator=separator,
        report=report,
    )

    # Compile information.
    pail = dict()
    pail["table_source_threshold"] = (pail_organization["table_threshold"])
    pail["table_source_threshold_scale"] = (
        pail_organization["table_threshold_scale"]
    )
    pail["count_samples"] = pail_organization["count_samples"]
    pail["count_variables"] = pail_organization["count_variables"]
    pail["matrix_component_scores"] = pail_sklearn["matrix_scores"]
    pail["table_component_scores"] = table_component_scores
    pail["table_component_variance_proportions"] = (
        pail_sklearn["table_component_variance_proportions"]
    )
    # Return.
    return pail



# Comparison of SKLearn and SVD methods


def compare_principal_components_methods(
    table=None,
    index_name=None,
    prefix=None,
    separator=None,
    report=None,
):
    """
    Organizes a Principal Components Analysis (PCA) by Singular Value
    Decomposition (SVD).

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    Relevant dimension: Principal Components represent variance across features

    Principal Components factorize covariance or correlation into direction
    (Eigenvectors) and scale (Eigenvalues).
    The Eigenvalues impart scale or weight to each Eigenvector.
    Each Eigenvector has its own Eigenvalue, and their sort orders mush match.

    Singular Value Decomposition assigns Eigenvector direction (sign, positive
    or negative) at random. Hence, take the absolute value of Principal
    Components before comparing matrices between different methods.

    arguments:
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        index_name (str): name of table's index column
        prefix (str): prefix for names of new principal component columns in
            table
        separator (str): separator for names of new columns
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about the Principal Components
            Analysis (PCA)

    """

    # Threshold and organize original matrix.
    pail_organization = organize_table_matrix_for_decomposition(
        table=table,
        report=False,
    )
    # Calculate Principal Component Scores in SKLearn.
    pail_sklearn = organize_principal_components_by_sklearn(
        table=table,
        index_name=index_name,
        prefix=prefix,
        separator=separator,
        report=False,
    )
    # Calculate Prinicipal Component Scores by Singular Value Decomposition.
    pail_svd = organize_principal_components_by_singular_value_decomposition(
        table=table,
        index_name=index_name,
        prefix=prefix,
        separator=separator,
        report=False,
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(
            "Report from: " +
            "compare_principal_components_methods()"
        )
        utility.print_terminal_partition(level=4)
        # Compare original matrix to matrix calculation from SVD factors.
        print(
            "Shape of source matrix: " + str(pail_organization["matrix"].shape)
        )
        utility.print_terminal_partition(level=4)
        print("Shapes of matrices for Principal Component Scores")
        print("SKLearn: " + str(pail_sklearn["matrix_component_scores"].shape))
        print(
            "SVD factor product: " +
            str(pail_svd["matrix_component_scores"].shape)
        )
        utility.print_terminal_partition(level=4)
        print("Compare score matrices between methods...")
        print("SKLearn versus SVD match:")
        print(numpy.allclose(
            numpy.absolute(pail_sklearn["matrix_component_scores"]),
            numpy.absolute(pail_svd["matrix_component_scores"]),
            rtol=1e-2, # relative tolerance, 1%
            atol=1e-3,
            equal_nan=False,
        ))
        utility.print_terminal_partition(level=4)
        print("Compare proportions of variance for principal components...")
        print("SKLearn:")
        print(pail_sklearn["table_component_variance_proportions"])
        utility.print_terminal_partition(level=5)
        print("SVD:")
        print(pail_svd["table_component_variance_proportions"])

    pass



#
