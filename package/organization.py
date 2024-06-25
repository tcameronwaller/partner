"""
Supply functionality for organization of information, especially within
two-dimensional tables.

This module 'organization' is part of the 'partner' package.

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


#dir()
#importlib.reload()

###############################################################################
# Functionality


# Merge, combine, transfer.


def transfer_table_rows_attributes_reference(
    table_main=None,
    column_main_key=None,
    table_reference=None,
    column_reference_key=None,
    columns_reference_transfer=None,
    prefix_reference_main=None,
    suffix_reference_main=None,
    report=None,
):
    """
    Transfers attributes from a reference Pandas data-frame table to a main
    Pandas data-frame table on the basis of a specific factor variable of which
    the individual values serve as keys.

    In this context, the values of the factor variable could be identifiers or
    names for individual samples, groups of samples, groups of experimental
    conditions, or other distinguishing variables such as outcome dependent
    variables and predictor independent variables in regression analyses.

    Notice that the reference table must include information for each and every
    unique value of the key factor variable.

    Review: TCW; 27 March 2024

    arguments:
        table_main (object): Pandas data-frame table
        column_main_key (str): name of column in main table for the key
            factor
        table_reference (object): reference table of attributes corresponding
            to unique values of the key factor in the main table
        column_reference_key (str): name of column in reference table for the
            key factor
        columns_reference_transfer (list<str>): name of columns in reference
            table for attributes corresponding to values of the key factor
        prefix_reference_main (str): prefix to append to names of attributes
            when transfered to the main table
        suffix_reference_main (str): suffix to append to names of attributes
            when transfered to the main table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_reference = table_reference.copy(deep=True)
    # Create dictionary for reference.
    table_reference.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_reference.set_index(
        [column_reference_key],
        append=False,
        drop=True,
        inplace=True,
    )
    reference_attributes = table_reference.to_dict("index")

    # Copy information in table.
    table_main_attribute = table_main.copy(deep=True)
    # Organize information in table.
    table_main_attribute.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Transfer attributes to the table of genetic correlations.
    for column_reference in columns_reference_transfer:
        column_novel = str(
            prefix_reference_main +
            column_reference +
            suffix_reference_main
        )
        table_main_attribute[column_novel] = table_main_attribute.apply(
            lambda row: (
                reference_attributes[row[column_main_key]][column_reference]
            ),
            axis="columns", # apply function to each row
        )
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        count_columns_source = (table_main.shape[1])
        count_columns_product = (table_main_attribute.shape[1])
        print("Count of columns in source table: " + str(count_columns_source))
        print(
            "Count of columns in product table: " +
            str(count_columns_product)
        )
        print("Table")
        print(table_main_attribute)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_main_attribute


# Sort.


def sort_table_rows_primary_secondary_reference(
    table_main=None,
    column_main_1=None,
    column_main_2=None,
    table_reference_1=None,
    table_reference_2=None,
    column_reference_1=None,
    column_reference_2=None,
    column_reference_sequence=None,
    report=None,
):
    """
    Sorts rows within a main Pandas data-frame table on the basis of
    separately-specified sort sequences corresponding to primary and secondary
    factor variables of which the individual values serve as keys.

    In this context, the values of the primary and secondary factor variables
    could be identifiers or names of individual samples, groups of samples,
    groups of experimental conditions, or other distinguishing variables such as
    outcome dependent variables and predictor independent variables in
    regression analyses.

    Notice that the reference tables must include information for each and every
    unique value of the primary and secondary key factor variables.

    The implementation strategy of this function relates closely to the function
    below.
    partner.organization.transfer_table_rows_attributes_reference()

    Review: TCW; 27 March 2024

    arguments:
        table_main (object): Pandas data-frame table
        column_main_1 (str): name of column in main table for primary
            key factor
        column_main_2 (str): name of column in main table for secondary
            key factor
        table_reference_1 (object): reference table for sort sequences
            corresponding to unique values of the primary key factor in the main
            table
        table_reference_2 (object): reference table for sort sequence
            corresponding to unique values of the secondary key factor in the
            main table
        column_reference_1 (str): name of column in primary reference
            table for primary key factor
        column_reference_2 (str): name of column in secondary reference
            table for secondary key factor
        column_reference_sequence (str): name of column in reference tables for
            sort sequence
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_reference_1 = table_reference_1.copy(deep=True)
    table_reference_2 = table_reference_2.copy(deep=True)
    # Create dictionaries for sort sequences of studies.
    table_reference_1.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_reference_2.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_reference_1.set_index(
        [column_reference_1],
        append=False,
        drop=True,
        inplace=True,
    )
    table_reference_2.set_index(
        [column_reference_2],
        append=False,
        drop=True,
        inplace=True,
    )
    sort_primary = table_reference_1.to_dict("index")
    sort_secondary = table_reference_2.to_dict("index")

    # Copy information in table.
    table_main_sort = table_main.copy(deep=True)
    # Organize information in table.
    table_main_sort.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Determine sort sequences for primary and secondary studies.
    table_main_sort["sort_primary_temp_73256"] = table_main_sort.apply(
        lambda row: (
            sort_primary[row[column_main_1]][column_reference_sequence]
        ),
        axis="columns", # apply function to each row
    )
    table_main_sort["sort_secondary_temp_73256"] = table_main_sort.apply(
        lambda row: (
            sort_secondary[row[column_main_2]][column_reference_sequence]
        ),
        axis="columns", # apply function to each row
    )
    # Sort rows within table.
    table_main_sort.sort_values(
        by=["sort_primary_temp_73256", "sort_secondary_temp_73256",],
        axis="index",
        ascending=True,
        inplace=True,
    )
    table_main_sort.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Remove unnecessary columns.
    if True:
        table_main_sort.drop(
            labels=["sort_primary_temp_73256", "sort_secondary_temp_73256"],
            axis="columns",
            inplace=True
        )
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        count_rows_table_source = (table_main.shape[0])
        count_rows_table_product = (table_main_sort.shape[0])
        print("Count of rows in source table: " + str(count_rows_table_source))
        print(
            "Count of rows in product table: " +
            str(count_rows_table_product)
        )
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_main_sort


def sort_table_rows_by_list_indices(
    table=None,
    list_sort=None,
    name_column=None,
    report=None,
):
    """
    Merges a list of indices as a new column to a Pandas data-frame table, sorts
    rows in the table by the new indices, and then removes the sort column.

    It would be useful to implement a conceptually similar function that
    concatenates multiple columns from a secondary table to the primary table
    and then uses those multiple columns for an hierarchical sort.

    arguments:
        table (object): Pandas data-frame table
        list_sort (list<str>): list of string integer indices for sort
        name_column (str): name for new column of indices for sort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table after sort on rows

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Convert strings in list to integers.
    list_sort_integer = copy.deepcopy(list(map(
        lambda value: int(value), list_sort
    )))
    # Check that the count of indices matches the count of rows in the table.
    count_list_elements = len(list_sort_integer)
    count_table_rows = (table.shape[0])
    if (count_list_elements != count_table_rows):
        #list_sort_integer = list_sort_integer[0:count_table_rows:1]
        putly.print_terminal_partition(level=4)
        print(
            "ERROR: Count of indices for sort does not match count of rows!"
        )
        putly.print_terminal_partition(level=4)
    # Introduce list of indices to the data-frame table.
    table[name_column] = list_sort_integer
    # Sort rows in table.
    table.sort_values(
        by=[name_column],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Remove the temporary column of indices for sort.
    table.drop(
        labels=[name_column],
        axis=1, # axis 0: rows; axis 1: columns
        inplace=True,
    )
    # Return information.
    return table


# Filter.

# TODO: TCW; 6 June 2024
# Simplify function 'match_table_row_redundant_interchangeable_pairs' to avoid
# logical errors.

def match_table_row_redundant_interchangeable_pairs(
    table=None,
    name_index=None,
    name_primary=None,
    name_secondary=None,
    row_index=None,
    row_primary=None,
    row_secondary=None,
    match_count=None,
    report=None,
):
    """
    Dependency:
    This function is a dependency of the function below.
    partner.extraction.filter_table_rows_ldsc_correlation()

    Determines whether values of two interchangeable identifiers from a single
    row in a table are redundant with those from a previous row in the original
    table.

    Review: TCW; 4 June 2024

    arguments:
        table (object): Pandas data-frame table from which the row-specific
            values originated
        name_index (str): name of column for sequential index
        name_primary (str): name of column for primary identifier
        name_secondary (str): name of column for secondary identifier
        row_index (int): current row's value of sequential index
        row_primary (str): current row's value of primary identifier
        row_secondary (str): current row's value of secondary identifier
        match_count (str): maximal count of matches to allow before indication
        report (bool): whether to print reports

    raises:

    returns:
        (int): logical binary representation of whether the identifiers in the
            current row are redundant

    """

    # Copy information in table.
    table_check = table.copy(deep=True)
    # Organize information in table.
    table_check.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_check.set_index(
        [name_primary, name_secondary,],
        append=False,
        drop=True,
        inplace=True,
    )
    # Find redundant records in table.
    index_match_1_2 = table_check.index.isin([(row_primary, row_secondary)])
    index_match_2_1 = table_check.index.isin([(row_secondary, row_primary)])
    table_match_1_2 = table_check[index_match_1_2]
    table_match_2_1 = table_check[index_match_2_1]
    values_index_match_1_2 = table_match_1_2[name_index].to_list()
    values_index_match_2_1 = table_match_2_1[name_index].to_list()
    values_index_match = (values_index_match_1_2 + values_index_match_2_1)
    # Determine whether the current combination of primary and secondary studies
    # is either irrelevant if both are identical or redundant with a previous
    # combination in the table.
    # Only keep one self pair regardless.
    indicator = 0
    if (
        (row_primary == row_secondary) and
        (
            (table_match_1_2.shape[0] > 0) or
            (table_match_2_1.shape[0] > 0)
        )
    ):
        if int(row_index) > min(values_index_match):
            indicator = 1
            pass
        pass
    elif (table_match_1_2.shape[0] > 0):
        if int(row_index) > min(values_index_match):
            indicator = 1
            pass
        pass
    elif (table_match_2_1.shape[0] > match_count):
        if int(row_index) > min(values_index_match):
            indicator = 1
            pass
        pass
    # Report.
    if report:
        print("Report.")
        print("row_index: " + str(row_index))
        print(table_match_1_2)
        print(table_match_2_1)
        print(values_index_match)
        print(min(values_index_match))
    return indicator


def filter_sort_table_columns(
    table=None,
    columns_sequence=None,
    report=None,
):
    """
    Filters and sorts the columns within a Pandas data-frame table.

    Review: TCW; 28 March 2024

    arguments:
        table (object): Pandas data-frame table
        columns_sequence (list<str>): identifiers or names of columns to keep
            in sort sequence
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_filter_sort = table.copy(deep=True)

    # Remove unnecessary columns.
    #table_filter_sort.drop(
    #    labels=["test", "blah", "q_significance",],
    #    axis="columns",
    #    inplace=True
    #)

    # Filter and sort table's columns.
    #table_filter_sort = table_filter_sort.loc[
    #    :, table_filter_sort.columns.isin(columns_sequence)
    #]
    table_filter_sort = table_filter_sort.filter(
        items=columns_sequence,
        axis="columns",
    )
    table_filter_sort = table_filter_sort[[*columns_sequence]]

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        count_columns_source = (table.shape[1])
        count_columns_product = (table_filter_sort.shape[1])
        print("Count of columns in source table: " + str(count_columns_source))
        print(
            "Count of columns in product table: " +
            str(count_columns_product)
        )
        print("Table")
        print(table_filter_sort)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_filter_sort


def filter_symmetrical_table_half_diagonal_two_value_types(
    table=None,
    name_index_type_value=None,
    types_value=None,
    report=None,
):
    """
    This function filters two different types of values in a symmetrical table
    to the upper half diagonal. The operation preserves the original multi-level
    indices across rows and columns.

    Format of source table in partial long format:

    group       type_value   group_2_a   group_2_b   group_2_c
    group_1_a   signal       -0.15       -0.20       -0.25
    group_1_a   p_value      0.001       0.001       0.001
    group_1_b   signal       0.15        0.20        0.25
    group_1_b   p_value      0.001       0.001       0.001

    Recommendations for names of indices:
    name_index_type_value="type_value"

    Recommendations for categories within index level "type_value":
    types_value=["signal", "p_value",]

    Review: TCW; 5 June 2024

    arguments:
        table (object): Pandas data-frame table with definitions of indices
            across columns and rows
        name_index_type_value (str): name of index across rows that designates
            which of two types of values corresponds to categorical indices
            across columns and rows
        types_value (list<str>): two types of values within the table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information in table.
    table_source = table.copy(deep=True)

    ##########
    # Separate the different types of values within the table.
    #table_one = table_source.loc[
    #    (table_source[name_index_type_value] == types_value[0]), :
    #].copy(deep=True)
    table_one = table_source[
        table_source.index.get_level_values(
            name_index_type_value
        ).isin([types_value[0]])
    ].copy(deep=True)
    #matrix_one = numpy.copy(table_one.to_numpy())
    table_two = table_source[
        table_source.index.get_level_values(
            name_index_type_value
        ).isin([types_value[1]])
    ].copy(deep=True)

    ##########
    # Filter to the diagonal lower half of a table representing a symmetrical
    # matrix.
    # Create a mask matrix for the lower half diagonal triangle of the table.
    # numpy.triu() # Upper half triangle.
    # numpy.tril() # Lower half triangle.
    matrix_mask_half = numpy.triu(
        numpy.ones(table_one.shape)
    ).astype(numpy.bool_)
    table_one_half = table_one.where(matrix_mask_half)
    table_two_half = table_two.where(matrix_mask_half)

    ##########
    # Combine the half-diagonal tables for the different types of values.
    table_product = pandas.concat(
        [table_one_half, table_two_half,],
        axis="index",
        join="outer",
        ignore_index=False,
        copy=True,
    )

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("porg.filter_symmetrical_table_half_diagonal_two_value_types()")
        count_columns_source = (table_source.shape[1])
        count_rows_source = (table_source.shape[0])
        count_columns_product = (table_product.shape[1])
        count_rows_product = (table_product.shape[0])
        print("Count of columns in source table: " + str(count_columns_source))
        print("Count of rows in source table: " + str(count_rows_source))
        print(
            "Count of columns in product table: " +
            str(count_columns_product)
        )
        print(
            "Count of rows in product table: " +
            str(count_rows_product))
        print("Table")
        print(table_product)
        putly.print_terminal_partition(level=4)

    ##########
    # Return information.
    return table_product



# Shape.


def transform_table_triple_quadruple_index_long_to_wide_partial(
    table_long=None,
    column_index_pivot=None,
    columns_index_stay=None,
    column_value_long=None,
    report=None,
):
    """
    Dependency:
    This function is a dependency of the function below.
    partner.organization.
    transform_table_triple_index_long_to_wide_partial_full()

    Transforms a Pandas data-frame table from full long format with a
    triple-level index across rows to wide formats with a double-level index
    across a first dimension and a single-level index across the second
    dimension.

    This function preserves the original sort sequence of values within each of
    the three or four levels or tiers of indices. This sort capability is a
    major contribution of this function.

    Description: 'long'
    Example format of source table in long format:
    group_1     group_2     type      value
    group_1_a   group_2_a   signal    -0.15
    group_1_a   group_2_a   p_value   0.001
    group_1_a   group_2_a   q_value   0.01
    group_1_a   group_2_b   signal    0.25
    group_1_a   group_2_b   p_value   0.001
    group_1_a   group_2_b   q_value   0.01
    group_1_a   group_2_c   signal    0.35
    group_1_a   group_2_c   p_value   0.001
    group_1_a   group_2_c   q_value   0.01
    group_1_b   group_2_a   signal    0.75
    group_1_b   group_2_a   p_value   0.001
    group_1_b   group_2_a   q_value   0.01
    group_1_b   group_2_b   signal    -0.35
    group_1_b   group_2_b   p_value   0.001
    group_1_b   group_2_b   q_value   0.01
    group_1_b   group_2_c   signal    0.45
    group_1_b   group_2_c   p_value   0.001
    group_1_b   group_2_c   q_value   0.01

    Description: 'wide_partial'
    Example format of table in wide format with single-level index across
    columns and multi-level index across rows:
    group_1     type_value   group_2_a   group_2_b   group_2_c
    group_1_a   signal       -0.15       0.25        0.35
    group_1_a   p_value      0.001       0.001       0.001
    group_1_a   q_value      0.01        0.01        0.01
    group_1_b   signal       0.75        -0.35       0.45
    group_1_b   p_value      0.001       0.001       0.001
    group_1_b   q_value      0.01        0.01        0.01

    Recommendations for names of indices:
    column_index_stay_1="group_primary" (experimental group or category)
    column_index_pivot="group_secondary" (experimental group or category)
    column_index_stay_2="type_value" (type of value, either signal, p-value, or
                                 q-value)

    Limitations:
    1. Notice that the current implementation only supports transformations that
    have unique combinations of the two or three dimensions of the index that
    remain together.
    2. Notice that there can only be two or three dimensions of the index that
    remain together (columns_index_stay).

    Review: TCW; 28 March 2024

    arguments:
        table_long (str): Pandas data-frame table in long format with a
            triple-level index in first three columns for two types of groups
            and one type of values
        column_index_pivot (str): name of column for second-tier index that will
            change to an index across columns in partial wide format
        columns_index_stay (list<str>): names of columns for levels of index
            that will remain together as an index across either rows or columns
        column_value_long (str): name of column for the value in the table in
            full long format
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table in partial wide format

    """

    # Copy information in table.
    table_long_copy = table_long.copy(deep=True)
    # Extract original sequence of values for indices.
    sequence_index_pivot = copy.deepcopy(
        table_long_copy.index.get_level_values(
            column_index_pivot
        ).unique().to_list()
    )
    sequence_index_stay_1 = copy.deepcopy(
        table_long_copy.index.get_level_values(
            columns_index_stay[0]
        ).unique().to_list()
    )
    sequence_index_stay_2 = copy.deepcopy(
        table_long_copy.index.get_level_values(
            columns_index_stay[1]
        ).unique().to_list()
    )
    if (len(columns_index_stay) == 3):
        sequence_index_stay_3 = copy.deepcopy(
            table_long_copy.index.get_level_values(
                columns_index_stay[2]
            ).unique().to_list()
        )
        pass

    # wide_partial

    # Transform table from full long format to partial wide format
    # (wide_partial).
    # Pandas dataframe methods "unstack" and "pivot" both can be useful in this
    # context.
    # Organize information in table.
    table_long_copy.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    if False:
        table_long_copy.set_index(
            [column_index_pivot, columns_index_stay,], # need to use list append
            append=False,
            drop=True,
            inplace=True,
        )
        table_wide_partial = table_long_copy.unstack(
            level=[column_index_pivot,],
        )
    table_wide_partial = table_long_copy.pivot(
        columns=[column_index_pivot,],
        index=columns_index_stay,
        values=column_value_long,
    )
    # Sort table's columns.
    table_wide_partial = table_wide_partial[[*sequence_index_pivot]]
    # Organize information in table.
    table_wide_partial.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Sort table's rows.
    if (len(columns_index_stay) == 2):
        table_wide_partial["sort_row_key"] = table_wide_partial.apply(
            lambda row: str(
                row[columns_index_stay[0]] +
                row[columns_index_stay[1]]
            ),
            axis="columns", # apply function to each row
        )
        sequence = dict()
        counter = 1
        for value_1 in sequence_index_stay_1:
            for value_2 in sequence_index_stay_2:
                sort_key = str(value_1 + value_2)
                if sort_key not in sequence.keys():
                    sequence[sort_key] = counter
                    counter +=1
                    pass
                pass
            pass
    elif (len(columns_index_stay) == 3):
        table_wide_partial["sort_row_key"] = table_wide_partial.apply(
            lambda row: str(
                row[columns_index_stay[0]] +
                row[columns_index_stay[1]] +
                row[columns_index_stay[2]]
            ),
            axis="columns", # apply function to each row
        )
        sequence = dict()
        counter = 1
        for value_1 in sequence_index_stay_1:
            for value_2 in sequence_index_stay_2:
                for value_3 in sequence_index_stay_3:
                    sort_key = str(value_1 + value_2 + value_3)
                    if sort_key not in sequence.keys():
                        sequence[sort_key] = counter
                        counter +=1
                        pass
                    pass
                pass
            pass
        pass
    table_wide_partial.sort_values(
        by=["sort_row_key"],
        axis="index",
        ascending=True,
        inplace=True,
        key=lambda x: x.map(sequence),
    )
    # Remove unnecessary columns.
    table_wide_partial.drop(
        labels=["sort_row_key",],
        axis="columns",
        inplace=True
    )
    # Organize information in table.
    table_wide_partial.set_index(
        columns_index_stay,
        append=False,
        drop=True,
        inplace=True,
    )

    # Return information.
    return table_wide_partial


def transform_table_triple_quadruple_index_long_to_wide(
    table_long=None,
    column_index_pivot=None,
    columns_index_stay=None,
    column_value_long=None,
    report=None,
):
    """
    Transforms a Pandas data-frame table from full long format with a
    triple-level index across rows to wide formats with a double-level index
    across a first dimension and a single-level index across the second
    dimension.

    The first wide format has a single-level index across the dimension of
    columns and a multi-level index across the dimension of rows, consisting of
    columns 0 and 1.

    The second wide format has a multi-level index across the dimension of
    columns, consisting of rows 0 and 1 (compound header), and a single-level
    index across the dimension of rows.

    Multi-level indices across one or both dimensions (rows, columns) of a table
    are useful for filters. It is most often common and convenient to use a
    multi-level index across rows, with the values of two or more columns
    determining a combination of categorical groups that together define the
    values in the other columns of that same row. This format often has the name
    'long format'.

    While the multi-level index is common and intuitive in the dimension across
    rows, in some applications it is more convenient to organize the table with
    the multi-level index in the dimension across columns. This format often has
    the name 'wide format'.

    This function preserves the original sort sequence of values within each of
    the three or four levels or tiers of indices. This sort capability is a
    major contribution of this function.

    Description: 'long'
    Example format of source table in long format:
    group_1     group_2     type      value
    group_1_a   group_2_a   signal    -0.15
    group_1_a   group_2_a   p_value   0.001
    group_1_a   group_2_a   q_value   0.01
    group_1_a   group_2_b   signal    0.25
    group_1_a   group_2_b   p_value   0.001
    group_1_a   group_2_b   q_value   0.01
    group_1_a   group_2_c   signal    0.35
    group_1_a   group_2_c   p_value   0.001
    group_1_a   group_2_c   q_value   0.01
    group_1_b   group_2_a   signal    0.75
    group_1_b   group_2_a   p_value   0.001
    group_1_b   group_2_a   q_value   0.01
    group_1_b   group_2_b   signal    -0.35
    group_1_b   group_2_b   p_value   0.001
    group_1_b   group_2_b   q_value   0.01
    group_1_b   group_2_c   signal    0.45
    group_1_b   group_2_c   p_value   0.001
    group_1_b   group_2_c   q_value   0.01

    Description: 'wide_partial'
    Example format of table in wide format with single-level index across
    columns and multi-level index across rows:
    group_1     type_value   group_2_a   group_2_b   group_2_c
    group_1_a   signal       -0.15       0.25        0.35
    group_1_a   p_value      0.001       0.001       0.001
    group_1_a   q_value      0.01        0.01        0.01
    group_1_b   signal       0.75        -0.35       0.45
    group_1_b   p_value      0.001       0.001       0.001
    group_1_b   q_value      0.01        0.01        0.01

    Description: 'wide_full'
    Example format of table in wide format with multi-level index across
    columns and single-level index across rows:
    Note: notice the use of "index", "index" as a place-holder, which is
    necessary when reading a table in this format from tab-delimited text into
    a Pandas data-frame table.
    index       group_2_a group_2_a group_2_a group_2_b group_2_b group_2_b ...
    index       signal    p_value   q_value   signal    p_value   q_value   ...
    group_1_a   -0.15     0.001     0.01       0.25     0.001     0.01      ...
    group_1_b    0.75     0.001     0.01      -0.35     0.001     0.01      ...

    Recommendations for names of indices:
    column_index_pivot="group_primary" (experimental group or category)
    column_index_stay_1="group_secondary" (experimental group or category)
    column_index_stay_2="type_value" (type of value, either signal, p-value, or
                                 q-value)

    Limitations:
    1. Notice that the current implementation only supports transformations that
    have unique combinations of the two or three dimensions of the index that
    remain together.
    2. Notice that there can only be two or three dimensions of the index that
    remain together (columns_index_stay).

    Review: TCW; 28 March 2024

    arguments:
        table_long (str): Pandas data-frame table in long format with a
            triple-level index in first three columns for two types of groups
            and one type of values
        column_index_pivot (str): name of column for level of index that will
            change to an index across columns in both partial and full wide
            formats
        columns_index_stay (list<str>): names of columns for levels of index
            that will remain together as an index across either rows or columns
        column_value_long (str): name of column for the value in the table in
            full long format
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): Pandas data-frame tables in partial and full wide
            formats

    """

    ##########
    # Transform table from full long to partial and full wide formats.

    # Copy information in table.
    table_long_copy = table_long.copy(deep=True)

    # 1. wide_partial
    table_wide_partial = (
        transform_table_triple_quadruple_index_long_to_wide_partial(
            table_long=table_long_copy,
            column_index_pivot=column_index_pivot,
            columns_index_stay=columns_index_stay,
            column_value_long=column_value_long,
            report=False,
    ))

    # 2. wide_full
    # Transpose table.
    table_wide_full = table_wide_partial.transpose(copy=True)

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("1. Original source table in full long format:")
        print(table_long_copy)
        #print("Column labels:")
        #labels_columns = table_long_copy.columns.to_list()
        #print(labels_columns)
        #print("Row labels:")
        #labels_rows = table_long_copy.index.to_list()
        #print(labels_rows)
        putly.print_terminal_partition(level=4)
        print(
            "The novel product tables below in partial long format and in " +
            "wide format should preserve the original sequences of labels."
        )
        putly.print_terminal_partition(level=4)
        print("2. Novel product table in partial wide format:")
        print(table_wide_partial)
        #print("Column labels:")
        #labels_columns = table_wide_partial.columns.to_list()
        #print(labels_columns)
        #print("Row labels:")
        #labels_rows = table_wide_partial.index.to_list()
        #print(labels_rows)
        putly.print_terminal_partition(level=4)
        print("3. Novel product table in full wide format:")
        print(table_wide_full)
        #print("Column labels:")
        #labels_columns = table_wide_full.columns.to_list()
        #print(labels_columns)
        #print("Row labels:")
        #labels_rows = table_wide_full.index.to_list()
        #print(labels_rows)
        putly.print_terminal_partition(level=4)
    # Collect information.
    pail = dict()
    pail["table_wide_partial"] = table_wide_partial
    pail["table_wide_full"] = table_wide_full
    # Return information.
    return pail

# TODO: TCW; 6 June 2024
# Pay attention to the sort order of the columns from the 'pivot' operation in
# function 'transform_table_quadruple_index_long_to_wide_square'.
# It might be prudent to sort the columns similarly to in the function
# 'transform_table_triple_quadruple_index_long_to_wide_partial'.


def transform_table_quadruple_index_long_to_wide_square(
    table_long=None,
    columns_index_pivot=None,
    columns_index_stay=None,
    column_value_long=None,
    report=None,
):
    """
    Transforms a Pandas data-frame table from full long format with a
    quadruple-level index across rows to wide format with double-level indices
    across both columns and rows.

    This function preserves the original sort sequence of values within each of
    the three or four levels or tiers of indices. This sort capability is a
    major contribution of this function.

    Description: 'long' (with triple levels or tiers of index across rows)
    Example format of source table in long format:
    group_1     group_2     type      value
    group_1_a   group_2_a   signal    -0.15
    group_1_a   group_2_a   p_value   0.001
    group_1_a   group_2_a   q_value   0.01
    group_1_a   group_2_b   signal    0.25
    group_1_a   group_2_b   p_value   0.001
    group_1_a   group_2_b   q_value   0.01
    group_1_a   group_2_c   signal    0.35
    group_1_a   group_2_c   p_value   0.001
    group_1_a   group_2_c   q_value   0.01
    group_1_b   group_2_a   signal    0.75
    group_1_b   group_2_a   p_value   0.001
    group_1_b   group_2_a   q_value   0.01
    group_1_b   group_2_b   signal    -0.35
    group_1_b   group_2_b   p_value   0.001
    group_1_b   group_2_b   q_value   0.01
    group_1_b   group_2_c   signal    0.45
    group_1_b   group_2_c   p_value   0.001
    group_1_b   group_2_c   q_value   0.01

    Limitations:
    1. Notice that the current implementation only supports transformations that
    have unique combinations of the two or three dimensions of the index that
    remain together.

    Review: TCW; 29 March 2024

    arguments:
        table_long (str): Pandas data-frame table in long format with a
            triple-level index in first three columns for two types of groups
            and one type of values
        columns_index_pivot (list<str>): names of columns for levels of index
            that will change to an index across columns
        columns_index_stay (list<str>): names of columns for levels of index
            that will remain together as an index across either rows or columns
        column_value_long (str): name of column for the value in the table in
            full long format
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): Pandas data-frame tables in partial and full wide
            formats

    """

    # Copy information in table.
    table_long_copy = table_long.copy(deep=True)
    # Extract original sequence of values for indices.
    sequence_index_pivot_1 = copy.deepcopy(
        table_long_copy.index.get_level_values(
            columns_index_pivot[0]
        ).unique().to_list()
    )
    sequence_index_pivot_2 = copy.deepcopy(
        table_long_copy.index.get_level_values(
            columns_index_pivot[1]
        ).unique().to_list()
    )
    sequence_index_stay_1 = copy.deepcopy(
        table_long_copy.index.get_level_values(
            columns_index_stay[0]
        ).unique().to_list()
    )
    sequence_index_stay_2 = copy.deepcopy(
        table_long_copy.index.get_level_values(
            columns_index_stay[1]
        ).unique().to_list()
    )

    # Transform table from full long format to square wide format (index is
    # 2 by 2).
    # Pandas dataframe methods "unstack" and "pivot" both can be useful in this
    # context.
    # Organize information in table.
    table_long_copy.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    if False:
        table_long_copy.set_index(
            [column_index_pivot, columns_index_stay,], # need to use list append
            append=False,
            drop=True,
            inplace=True,
        )
        table_wide_partial = table_long_copy.unstack(
            level=[column_index_pivot,],
        )
    table_wide = table_long_copy.pivot(
        columns=columns_index_pivot,
        index=columns_index_stay,
        values=column_value_long,
    )

    # The Pandas pivot transformation preserves the original sort sequence of
    # the index values that shift across columns but does not preserve the
    # original sort sequence of the index values that remain across rows.

    # Notice that the format for accessing the values across rows of the index
    # columns is a bit different due to the multi-dimensional index across
    # columns.

    # Organize information in table.
    table_wide.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Sort table's rows.
    table_wide["sort_row_key"] = table_wide.apply(
        lambda row: str(
                str(row[columns_index_stay[0]][0]) +
                str(row[columns_index_stay[1]][0])
        ),
        axis="columns", # apply function to each row
    )
    sequence = dict()
    counter = 1
    for value_1 in sequence_index_stay_1:
        for value_2 in sequence_index_stay_2:
            sort_key = str(value_1 + value_2)
            if sort_key not in sequence.keys():
                sequence[sort_key] = counter
                counter +=1
                pass
            pass
        pass
    table_wide.sort_values(
        by=["sort_row_key"],
        axis="index",
        ascending=True,
        inplace=True,
        key=lambda x: x.map(sequence),
    )
    # Remove unnecessary columns.
    table_wide.drop(
        labels=["sort_row_key",],
        axis="columns",
        inplace=True
    )
    # Organize information in table.
    table_wide.set_index(
        columns_index_stay,
        append=False,
        drop=True,
        inplace=True,
    )

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("1. Original source table in full long format:")
        print(table_long_copy)
        #print("Column labels:")
        #labels_columns = table_long_copy.columns.to_list()
        #print(labels_columns)
        #print("Row labels:")
        #labels_rows = table_long_copy.index.to_list()
        #print(labels_rows)
        putly.print_terminal_partition(level=4)
        print(
            "The novel product tables below in partial long format and in " +
            "wide format should preserve the original sequences of labels."
        )
        putly.print_terminal_partition(level=4)
        print("2. Novel product table in full wide format:")
        print(table_wide)
        #print("Column labels:")
        #labels_columns = table_wide.columns.to_list()
        #print(labels_columns)
        #print("Row labels:")
        #labels_rows = table_wide.index.to_list()
        #print(labels_rows)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_wide


# Cluster.


def cluster_data_columns(
    data=None,
):
    """
    Clusters features on columns by their similarities across instances on
    rows.

    arguments:
        data (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of values

    """

    data = data.copy(deep=True)
    # Cluster.
    columns = data.columns.to_numpy()#.tolist()
    rows = data.index.to_numpy()#.tolist()
    index_name = data.index.name
    # Plan to cluster across columns.
    # Organize columns across dimension zero.
    matrix = numpy.transpose(data.to_numpy())
    linkage = scipy.cluster.hierarchy.linkage(
        matrix,
        method="average", # "single", "complete", "average"
        metric="euclidean",
        optimal_ordering=True,
    )
    dendrogram = scipy.cluster.hierarchy.dendrogram(
        linkage,
    )
    # Access seriation from dendrogram leaves.
    leaves = dendrogram["leaves"]
    # Sort matrix row and column labels.
    indices = range(0, len(matrix))
    matrix_cluster = list(map(
        lambda index: matrix[leaves[index]],
        indices
    ))
    columns_sort = list(map(
        lambda index: columns[leaves[index]],
        indices
    ))
    # Organize data.
    data_cluster = pandas.DataFrame(
        data=numpy.transpose(matrix_cluster),
        index=rows,
        columns=columns_sort,
    )
    data_cluster.rename_axis(
        index_name,
        axis="index",
        inplace=True,
    )
    # Return information.
    return data_cluster


def cluster_data_rows(
    data=None,
):
    """
    Clusters instances on rows by their similarities across features on
    columns.

    arguments:
        data (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of values

    """

    data = data.copy(deep=True)
    # Cluster.
    columns = data.columns.to_numpy()#.tolist()
    rows = data.index.to_numpy()#.tolist()
    index_name = data.index.name
    # Plan to cluster across columns.
    # Organize rows across dimension zero.
    matrix = data.to_numpy()
    linkage = scipy.cluster.hierarchy.linkage(
        matrix,
        method="average", # "single", "complete", "average"
        metric="euclidean",
        optimal_ordering=True,
    )
    dendrogram = scipy.cluster.hierarchy.dendrogram(
        linkage,
    )
    # Access seriation from dendrogram leaves.
    leaves = dendrogram["leaves"]
    # Sort matrix row and column labels.
    indices = range(0, len(matrix))
    matrix_cluster = list(map(
        lambda index: matrix[leaves[index]],
        indices
    ))
    rows_sort = list(map(
        lambda index: rows[leaves[index]],
        indices
    ))
    # Organize data.
    data_cluster = pandas.DataFrame(
        data=matrix_cluster,
        index=rows_sort,
        columns=columns,
    )
    data_cluster.rename_axis(
        index_name,
        axis="index",
        inplace=True,
    )
    # Return information.
    return data_cluster


def cluster_data_rows_by_group(
    group=None,
    index=None,
    data=None,
):
    """
    Clusters instances on rows by their similarities across features on
    columns.

    arguments:
        group (str): name of column to use for groups
        index (str): name of column to use for index during cluster
        data (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of values

    """

    data = data.copy(deep=True)
    groups = data.groupby(
        level=[group],
    )
    data_collection = pandas.DataFrame()
    for name, data_group in groups:
        data_group = data_group.copy(deep=True)
        data_group.reset_index(
            level=None,
            inplace=True
        )
        data_group.set_index(
            [index],
            append=False,
            drop=True,
            inplace=True
        )
        data_cluster = cluster_data_rows(
            data=data_group,
        )
        data_collection = pandas.concat(
            [data_collection, data_cluster],
            ignore_index=False,
        )
    # Return information.
    return data_collection


# Read and transform tables with multi-level indices.


def read_table_multiindex_columns_transpose(
    path_file_table=None,
    row_index_place_holder=None,
    name_row_index=None,
    name_column_index_1=None,
    name_column_index_2=None,
    report=None,
):
    """
    Dependency:
    This function is a dependency of the function below.
    partner.organization.
    read_table_multiindex_columns_transform_calculate_q_values()

    Reads from file, organizes, and transforms a table as a Pandas data frame.
    The original source table is in a tab-delimited text file that has a
    multi-level index across the dimension of columns consisting of rows 0 and 1
    (compound header) and a single-level index across the dimension of rows.

    Multi-level indices across one or both dimensions (rows, columns) of a table
    are useful for filters. It is most often common and convenient to use a
    multi-level index across rows, with the values of two or more columns
    determining a combination of categorical groups that together define the
    values in the other columns of that same row. This format often has the name
    "long format".

    While the multi-level index is common and intuitive in the dimension across
    rows, in some applications it is far more convenient to organize the
    original data table with the multi-level index in the dimension across
    columns. This format often has the name "wide format".

    This function makes it convenient to convert from wide format to long
    format. This function also calculates Benjamini-Hochberg False Discovery
    Rate q-values from all p-values in the table.

    This function could also serve as a template for further adaptation to
    accommodate multi-level indices across both rows and columns.

    Format of source table in wide format:
    Note: notice the use of "index", "index" as a place-holder.

    index       group_1_a   group_1_a   group_1_b   group_1_b
    index       signal      p_value     signal      p_value
    group_2_a   -0.15       0.001       0.15        0.001
    group_2_b   -0.20       0.001       0.20        0.001
    group_2_c   -0.25       0.001       0.25        0.001

    Format of product table in partial long format:

    group       type_value   group_2_a   group_2_b   group_2_c
    group_1_a   signal       -0.15       -0.20       -0.25
    group_1_a   p_value      0.001       0.001       0.001
    group_1_b   signal       0.15        0.20        0.25
    group_1_b   p_value      0.001       0.001       0.001

    arguments:
        path_file_table (str): path to file for original source table
        row_index_place_holder (str): place holder for index across rows
        name_row_index (str): new name for index across rows
        name_column_index_1 (str): new name for first index across columns
        name_column_index_2 (str): new name for second index across columns;
            one of the two categorical values must be 'p_value'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame tables before and after transposition

    """

    ##########
    # 1. Read table from file.

    # Read information from file.
    table_wide = pandas.read_csv(
        path_file_table,
        sep="\t",
        header=[0,1],
        na_values=["nan", "na", "NAN", "NA",],
    )

    ##########
    # 2. Organize indices and transform table from wide to partial long format.

    # Organize information in table.
    table_wide.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_wide.set_index(
        [(row_index_place_holder, row_index_place_holder),],
        append=False,
        drop=True,
        inplace=True,
    )
    table_wide.index.rename(
        name_row_index,
        inplace=True,
    ) # single-dimensional index
    table_wide.columns.rename(
        name_column_index_1,
        level=0,
        inplace=True,
    ) # multi-dimensional index
    table_wide.columns.rename(
        name_column_index_2,
        level=1,
        inplace=True,
    ) # multi-dimensional index
    # Transpose table.
    table_long_partial = table_wide.transpose(copy=True)

    ##########
    # 3. Report.
    if report:
        print_terminal_partition(level=4)
        print("1. Original source table in wide format after organization:")
        print(table_wide)
        print("Column labels:")
        labels_columns = table_wide.columns.to_list()
        print(labels_columns)
        print("Row labels:")
        labels_rows = table_wide.index.to_list()
        print(labels_rows)
        print_terminal_partition(level=4)
        print("2. Novel product table in partial long format:")
        print(table_long_partial)
        print_terminal_partition(level=4)

    # Collect information.
    pail = dict()
    pail["table_wide"] = table_wide
    pail["table_long_partial"] = table_long_partial
    # Return information.
    return pail


def read_table_multiindex_columns_transform_calculate_q_values(
    path_file_table=None,
    row_index_place_holder=None,
    name_row_index=None,
    name_column_index_1=None,
    name_column_index_2=None,
    filter_half_triangle=None,
    report=None,
):
    """
    Reads from file, organizes, and transforms a table as a Pandas data frame.
    The original source table is in a tab-delimited text file that has a
    multi-level index across the dimension of columns consisting of rows 0 and 1
    (compound header) and a single-level index across the dimension of rows.

    This function calculates Benjamini-Hochberg False Discovery Rate (FDR)
    q-values from all p-values in the original table, unless the argument is
    active to filter a symmetrical table's matrix of values to a half triangle
    along the diagonal axis.

    Multi-level indices across one or both dimensions (rows, columns) of a table
    are useful for filters. It is most often common and convenient to use a
    multi-level index across rows, with the values of two or more columns
    determining a combination of categorical groups that together define the
    values in the other columns of that same row. This format often has the name
    "long format".

    While the multi-level index is common and intuitive in the dimension across
    rows, in some applications it is far more convenient to organize the
    original data table with the multi-level index in the dimension across
    columns. This format often has the name "wide format".

    This function makes it convenient to convert from wide format to long
    format. This function also calculates Benjamini-Hochberg False Discovery
    Rate q-values from all p-values in the table.

    This function could also serve as a template for further adaptation to
    accommodate multi-level indices across both rows and columns.

    Format of source table in wide format:
    Note: notice the use of "index", "index" as a place-holder.

    index       group_1_a   group_1_a   group_1_b   group_1_b
    index       signal      p_value     signal      p_value
    group_2_a   -0.15       0.001       0.15        0.001
    group_2_b   -0.20       0.001       0.20        0.001
    group_2_c   -0.25       0.001       0.25        0.001

    Format of product table in partial long format:

    group       type_value   group_2_a   group_2_b   group_2_c
    group_1_a   signal       -0.15       -0.20       -0.25
    group_1_a   p_value      0.001       0.001       0.001
    group_1_b   signal       0.15        0.20        0.25
    group_1_b   p_value      0.001       0.001       0.001

    Format of product table in complete long format:

    group_1     group_2     type      value
    group_1_a   group_2_a   signal    -0.15
    group_1_a   group_2_a   p_value   0.001
    group_1_a   group_2_b   signal    -0.20
    group_1_a   group_2_b   p_value   0.001
    group_1_a   group_2_c   signal    -0.25
    group_1_a   group_2_c   p_value   0.001
    group_1_b   group_2_a   signal    0.15
    group_1_b   group_2_a   p_value   0.001
    group_1_b   group_2_b   signal    0.20
    group_1_b   group_2_b   p_value   0.001
    group_1_b   group_2_c   signal    0.25
    group_1_b   group_2_c   p_value   0.001

    Format of product table in partial long format with new q-values:

    group       type      group_2_a   group_2_b   group_2_c
    group_1_a   signal    -0.15       -0.20       -0.25
    group_1_a   p_value   0.001       0.001       0.001
    group_1_a   q_value   0.001       0.001       0.001
    group_1_b   signal    0.15        0.20        0.25
    group_1_b   p_value   0.001       0.001       0.001
    group_1_b   q_value   0.001       0.001       0.001

    Recommendations for names of indices:
    row_index_place_holder="group_secondary",
    name_row_index="group_secondary",
    name_column_index_1="group_primary",
    name_column_index_2="type_value",

    Recommendations for categories within name_column_index_2, "type_value":
    signal
    p_value
    q_value (calculated and introduced by this function)

    arguments:
        path_file_table (str): path to file for original source table
        row_index_place_holder (str): place holder for index across rows
        name_row_index (str): new name for index across rows
        name_column_index_1 (str): new name for first index across columns
        name_column_index_2 (str): new name for second index across columns;
            one of the two categorical values must be 'p_value'
        filter_half_triangle (bool): whether to filter the values of a table
            representing a symmetrical matrix to keep only the lower half along
            the diagonal
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame tables before and after transposition and
            with Benjamini-Hochberg False Discovery Rate q-values

    """

    ##########
    # 1. Read table from file.
    # 2. Organize indices and transform table from wide to partial long format.

    source = read_table_multiindex_columns_transpose(
        path_file_table=path_file_table,
        row_index_place_holder=row_index_place_holder,
        name_row_index=name_row_index,
        name_column_index_1=name_column_index_1,
        name_column_index_2=name_column_index_2,
        report=None,
    )
    # Copy information in table.
    table_wide = source["table_wide"].copy(deep=True)
    table_long_partial = source["table_long_partial"].copy(deep=True)

    ##########
    # 3. Filter to the diagonal lower half of a table representing a symmetrical
    # matrix.
    # It is necessary to filter the "signal" values and the "p-value" values
    # separately and then combine them together again.
    if (filter_half_triangle):
        table_long_partial_half = (
            filter_symmetrical_table_half_diagonal_two_value_types(
                table=table_long_partial,
                name_index_type_value=name_column_index_2,
                types_value=["signal", "p_value",],
                report=report,
            )
        )
        # Copy information in table.
        table_long_partial = table_long_partial_half.copy(deep=True)
        pass

    ##########
    # 4. Transform table from partial long to complete long format.

    # Transform table to complete long format.
    # Pandas dataframe methods "stack", "melt", and "wide_to_long", can all be
    # useful in this context.
    # Method "stack" converts to a multi-index series when the column index only
    # has a single level.
    # Method "wide_to_long" assumes that the information about multiple levels
    # in the column index is stored in delimited strings of compound column
    # names.
    # Copy information in table.
    table_long_partial_copy = table_long_partial.copy(deep=True)
    if False:
        table_long_complete = table_long_partial_copy.stack(
            level=name_row_index,
            #future_stack=True,
        )
    table_long_complete = table_long_partial_copy.melt(
        id_vars=None,
        value_vars=None,
        var_name=name_row_index,
        value_name="value",
        ignore_index=False,
    )
    # Organize information in table.
    table_long_complete.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_long_complete.set_index(
        [name_column_index_1, name_row_index, name_column_index_2,],
        append=False,
        drop=True,
        inplace=True,
    )

    ##########
    # 5. Calculate Benjamini-Hochberg False Discovery Rate q-values.
    table_q = pdesc.calculate_table_long_false_discovery_rate_q_values(
        table=table_long_complete,
        name_index_type_value=name_column_index_2,
        type_p_value="p_value",
        threshold=0.05,
        names_indices_rows=[
            name_column_index_1,
            name_row_index,
            name_column_index_2,
        ],
        report=report,
    )

    ##########
    # 6. Transform table to long partial and wide formats.

    # Copy information in table.
    table_long_complete_q_copy = table_long_complete_q.copy(deep=True)
    # Organize information in table.
    table_long_complete_q_copy.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Transform table to partial long format.
    # Pandas dataframe methods "unstack" and "pivot" can be useful in this
    # context.
    if False:
        table_long_complete_q_copy.set_index(
            [name_column_index_1, name_row_index, name_column_index_2,],
            append=False,
            drop=True,
            inplace=True,
        )
        table_long_partial_q = table_long_complete_q_copy.unstack(
            level=[name_row_index,],
        )
    table_long_partial_q = table_long_complete_q_copy.pivot(
        columns=[name_row_index,],
        index=[name_column_index_1, name_column_index_2,],
        values="value",
    )
    # Sort columns in table.
    sequence_row_index = copy.deepcopy(
        table_wide.index.get_level_values(name_row_index).to_list()
    )
    table_long_partial_q = table_long_partial_q[[*sequence_row_index]]
    # Sort rows in table
    table_long_partial_q.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_long_partial_q["sort_row_key"] = table_long_partial_q.apply(
        lambda row: str(
            row[name_column_index_1] +
            row[name_column_index_2]
        ),
        axis="columns", # apply function to each row
    )
    sequence_column_index_1 = copy.deepcopy(
        table_wide.columns.get_level_values(name_column_index_1).to_list()
    )
    sequence_column_index_2 = ["signal", "p_value", "q_value"]
    #values_tier_2 = copy.deepcopy(table_long_partial.columns.to_list())
    #values_tier_3 = ["signal", "p_value", "q_value"]
    sequence = dict()
    counter = 1
    for value_1 in sequence_column_index_1:
        for value_2 in sequence_column_index_2:
            sort_key = str(value_1 + value_2)
            if sort_key not in sequence.keys():
                sequence[sort_key] = counter
                counter +=1
                pass
            pass
        pass
    table_long_partial_q.sort_values(
        by=["sort_row_key"],
        axis="index",
        ascending=True,
        inplace=True,
        key=lambda x: x.map(sequence),
    )
    # Remove unnecessary columns.
    table_long_partial_q.drop(
        labels=["sort_row_key",],
        axis="columns",
        inplace=True
    )
    # Organize information in table.
    table_long_partial_q.set_index(
        [name_column_index_1, name_column_index_2,],
        append=False,
        drop=True,
        inplace=True,
    )
    # Transpose table.
    table_wide_q = table_long_partial_q.transpose(copy=True)

    ##########
    # 7. Report.
    if report:
        print_terminal_partition(level=4)
        print("1. Original source table in wide format after organization:")
        print(table_wide)
        print("Column labels:")
        labels_columns = table_wide.columns.to_list()
        print(labels_columns)
        print("Row labels:")
        labels_rows = table_wide.index.to_list()
        print(labels_rows)
        print_terminal_partition(level=4)
        print(
            "The novel product tables below in partial long format and in " +
            "wide format should preserve the original sequences of labels."
        )
        print(
            "The novel product tables below will also show any half-diagonal " +
            "filter on a table representing a symmetrical matrix."
        )
        print_terminal_partition(level=4)
        print(
            "2. Novel product table in complete long format with newly " +
            "calculated Benjamini-Hochberg False Discovery Rate q-values:"
        )
        print(table_long_complete_q)
        print_terminal_partition(level=4)
        print(
            "3. Novel product table in partial long format with newly " +
            "calculated Benjamini-Hochberg False Discovery Rate q-values:"
        )
        print(table_long_partial_q)
        print_terminal_partition(level=4)
        print(
            "4. Novel product table in wide format with newly " +
            "calculated Benjamini-Hochberg False Discovery Rate q-values:"
        )
        print(table_wide_q)
        print_terminal_partition(level=4)

    # Collect information.
    pail = dict()
    pail["table_wide_q"] = table_wide_q
    pail["table_long_complete_q"] = table_long_complete_q
    pail["table_long_partial_q"] = table_long_partial_q
    # Return information.
    return pail





###############################################################################
# Procedure
# Currently, this module is not executable.

##########
