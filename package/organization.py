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

################################################################################
# Notes

# TODO: TCW; 6 March 2024
# TODO: If the genetic correlation estimate is missing ("nan") with the explanation
# "rg out of bounds", then access the rg, se, p-value, z-score, etc from the
# tab-delimited summary information at the bottom of the report log.


################################################################################
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
import partner.scale as pscl

#dir()
#importlib.reload()

################################################################################
# Functionality


##########
# Indices.


def translate_names_table_indices_columns_rows(
    table=None,
    index_columns_product=None,
    index_rows_source=None,
    index_rows_product=None,
    report=None,
):
    """
    Translate names of explicitly named indices across columns and rows in a
    table.

    The table must have explicitly named indices across columns and rows.

    Review: TCW; 15 October 2024

    arguments:
        table (object): Pandas data-frame table
        index_columns_product (str): name of single-level index across columns
        index_rows_source (str): name of single-level index across rows
        index_rows_product (str): name of single-level index across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_source = table.copy(deep=True)
    table_product = table.copy(deep=True)
    # Organize information in table.
    table_product.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Translate names of columns.
    translations = dict()
    translations[index_rows_source] = index_rows_product
    table_product.rename(
        columns=translations,
        inplace=True,
    )
    table_product.set_index(
        [index_rows_product],
        append=False,
        drop=True,
        inplace=True,
    )
    table_product.columns.rename(
        index_columns_product,
        inplace=True,
    ) # single-dimensional index
    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("Report:")
        print("partner")
        print("organization")
        print("translate_names_table_indices_columns_rows()")
        putly.print_terminal_partition(level=4)
        print("original, source table before translations:")
        print(table_source)
        putly.print_terminal_partition(level=4)
        print("novel, product table after translations:")
        print(table_product)
        pass
    # Return information.
    return table_product


def translate_identifiers_table_indices_columns_rows(
    table=None,
    index_rows=None,
    translations_columns=None,
    translations_rows=None,
    remove_redundancy=None,
    report=None,
):
    """
    Translate names of individual identifiers in a table's indices across
    columns and rows.

    For versatility and convenience, this table does not have explicitly named
    indices across columns or rows.

    This function only translates names that are in the argument references.

    Sometimes translations introduce redundancy in the identifiers within
    indices across columns or rows. This function offers the option to remove
    this redundancy.

    Review: TCW; 25 October 2024

    arguments:
        table (object): Pandas data-frame table
        index_rows (str): name for index across rows, which corresponds to a
            column in the original source table
        translations_columns (dict<str>): translations for names or identifiers
            of columns
        translations_rows (dict<str>): translations for names or identifiers
            of columns
        remove_redundancy (bool): whether to remove columns or rows with
            redundant identifiers after translation
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Define subordinate functions for internal use.
    def determine_translation_identifier(source=None, translations=None):
        if (source in translations.keys()):
            product = translations[source]
        else:
            product = source
        return product

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    translations_columns = copy.deepcopy(translations_columns)
    translations_rows = copy.deepcopy(translations_rows)

    # Translate names of columns.
    #if (len(list(translations_columns.keys())) > 0):
    if (translations_columns is not None):
        table.rename(
            columns=translations_columns,
            inplace=True,
        )
        pass

    # Translate names of rows.
    if (translations_rows is not None):
        index_rows_obsolete = str(index_rows + "_obsolete")
        table[index_rows_obsolete] = table[index_rows]
        table[index_rows] = table.apply(
            lambda row:
                determine_translation_identifier(
                    source=row[index_rows_obsolete],
                    translations=translations_rows,
                ),
            axis="columns", # apply function to each row
        )
        # Remove unnecessary columns.
        table.drop(
            labels=[index_rows_obsolete,],
            axis="columns",
            inplace=True
        )
        pass

    # Determine whether to remove any redundancy in identifiers across columns
    # and rows.
    if remove_redundancy:
        table.drop_duplicates(
            subset=[index_rows,],
            keep="first",
            inplace=True,
        )
        table = table.loc[
            :, ~table.columns.duplicated(
                keep="first",
            )
        ].copy(deep=True)

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        print("function: translate_identifiers_table_indices_columns_rows()")
        putly.print_terminal_partition(level=5)
    # Return information.
    return table



##########
# Transfer, Fill.

# Transfer information between tables.
# This operation differs from a merge or join in that it does not require a
# one to one match of identifiers.
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

    Unlike a merge operation on tables, this function requires the key factor
    variable to become a unique, unambiguous index in the reference table but
    does not require the key factor variable to be unique in the main table.

    Neither the main table nor the reference table ought to have special,
    named indices, as this function would remove them.

    Review: TCW; 1 August 2024

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
    # Transfer attributes from reference table to the main table.
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
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        print("function: transfer_table_rows_attributes_reference()")
        putly.print_terminal_partition(level=5)
        count_columns_source = (table_main.shape[1])
        count_columns_product = (table_main_attribute.shape[1])
        print("count of columns in source table: " + str(count_columns_source))
        print(
            "count of columns in product table: " +
            str(count_columns_product)
        )
        print("table")
        print(table_main_attribute)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_main_attribute


def determine_fill_table_groups_rows(
    table=None,
    column_group=None,
    index_rows=None,
    groups_rows=None,
    report=None,
):
    """
    Determine and fill names of groups to which each row in a table belongs.

    For versatility and convenience, this table does not have explicitly named
    indices across columns or rows.

    Review; TCW; 15 October 2024

    arguments:
        table (object): Pandas data-frame table
        column_group (str): name of new column in product table for names of
            groups
        index_rows (str): name for index across rows, which corresponds to a
            column in the original source table
        groups_rows (dict<list<str>>): names of groups and identifiers of rows
            that belong to each of these groups
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Define subordinate functions for internal use.
    def determine_row_group(identifier=None, groups_rows=None):
        group = ""
        for key in groups_rows.keys():
            if identifier in groups_rows[key]:
                group = key
                pass
            pass
        return group

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    groups_rows = copy.deepcopy(groups_rows)
    columns_source = copy.deepcopy(table.columns.to_list())

    # Determine to which categorical group each sample belongs.
    table[column_group] = table.apply(
        lambda row: determine_row_group(
            identifier=row[index_rows],
            groups_rows=groups_rows,
        ),
        axis="columns", # apply function to each row
    )

    # Filter and sort columns in table.
    columns_sequence = copy.deepcopy(columns_source)
    columns_sequence.remove(index_rows)
    columns_sequence.insert(0, column_group)
    columns_sequence.insert(0, index_rows)
    table = filter_sort_table_columns(
        table=table,
        columns_sequence=columns_sequence,
        report=report,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        print("function: determine_fill_table_groups_rows()")
        putly.print_terminal_partition(level=5)
    # Return information.
    return table



##########
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


def sort_table_rows_by_single_column_reference(
    table=None,
    index_rows=None,
    column_reference=None,
    column_sort_temporary=None,
    reference_sort=None,
):
    """
    Sort the rows in a table by reference to categorical values in a single
    column.

    Original source table must not have an explicitly defined index across
    rows.

    There are multiple convenient options to create the reference indices
    for the sort on the basis of a list of categorical values.

    Option 1.
    reference_sort = dict()
    index = 0
    for name in names_sequence:
        reference_sort[name] = index
        index += 1
        pass

    Option 2.
    reference_sort = dict(zip(names_sequence, range(len(names_sequence))))

    Option 3.
    reference_sort = {key: value for value, key in enumerate(
        names_sequence
        )
    }

    Review: TCW; 16 November 2024

    arguments:
        table (object): Pandas data-frame table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        column_reference (str): name of column to use for match with indices
            in reference for sort
        column_sort_temporary (str): name of temporary column for sort
        reference_sort (dict<int>): indices for sort that match categorical
            values in column in source table

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    reference_sort = copy.deepcopy(reference_sort)
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Sort rows in table by reference.
    table[column_sort_temporary] = table[column_reference].copy(deep=True)
    table[column_sort_temporary] = table[column_sort_temporary].replace(
        reference_sort
    )
    # Sort rows in table.
    table.sort_values(
        by=[column_sort_temporary,],
        axis="index",
        ascending=True,
        inplace=True,
    )
    # Remove unnecessary columns.
    table.drop(
        labels=[column_sort_temporary,],
        axis="columns",
        inplace=True
    )
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Return information.
    return table



##########
# Filter.


def filter_sort_table_columns(
    table=None,
    columns_sequence=None,
    report=None,
):
    """
    Filters and sorts the columns within a Pandas data-frame table.

    Identifiers of columns must only include those that are actually in the
    table's columns.

    Review: TCW; 15 October 2024

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
        count_columns_source = (table.shape[1])
        count_columns_product = (table_filter_sort.shape[1])
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        print("function: filter_sort_table_columns()")
        putly.print_terminal_partition(level=5)
        print("count of columns in source table: " + str(count_columns_source))
        print(
            "count of columns in product table: " +
            str(count_columns_product)
        )
        putly.print_terminal_partition(level=5)
        print("product table:")
        print(table_filter_sort)
        putly.print_terminal_partition(level=5)
    # Return information.
    return table_filter_sort


def filter_select_table_columns_rows_by_identifiers(
    table=None,
    index_rows=None,
    identifiers_columns=None,
    identifiers_rows=None,
    report=None,
):
    """
    Filter and select specific columns and rows from a Pandas data-frame table.

    Review: TCW; 16 October 2024

    arguments:
        table (object): Pandas data-frame table
        index_rows (str): name for index across rows, which corresponds to a
            column in the original source table
        identifiers_columns (list<str>): identifiers or names of columns to
            keep
        identifiers_rows (list<str>): identifiers or names of rows to keep
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information in table.
    table_source = table.copy(deep=True)
    # Copy other information.
    identifiers_columns = copy.deepcopy(identifiers_columns)
    identifiers_rows = copy.deepcopy(identifiers_rows)

    # Filter and sort columns in table.
    columns_sequence = copy.deepcopy(identifiers_columns)
    columns_sequence.insert(0, index_rows)
    table_product = filter_sort_table_columns(
        table=table_source,
        columns_sequence=columns_sequence,
        report=False,
    )

    # Filter rows in table.
    table_product = table_product.loc[
        table_product[index_rows].isin(identifiers_rows), :
    ].copy(deep=True)

    # Report.
    if report:
        count_columns_source = (table_source.shape[1])
        count_columns_product = (table_product.shape[1])
        count_rows_source = (table_source.shape[0])
        count_rows_product = (table_product.shape[0])
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        print("function: filter_select_table_columns_rows_by_identifiers()")
        putly.print_terminal_partition(level=5)
        print("source table")
        print("count of columns: " + str(count_columns_source))
        print("count of rows: " + str(count_rows_source))
        putly.print_terminal_partition(level=5)
        print("product table")
        print("count of columns: " + str(count_columns_product))
        print("count of rows: " + str(count_rows_product))
        putly.print_terminal_partition(level=5)
    # Return information.
    return table_product


def filter_extract_table_row_identifiers_by_columns_categories(
    table=None,
    column_identifier=None,
    name=None,
    columns_categories=None,
    report=None,
):
    """
    Filters rows in a table to extract identifiers for a set of rows
    corresponding to selection by specific categorical values of columns for
    features.

    This function preserves the original sequence of rows from the source
    table.

    arguments:
        table (object): Pandas data-frame table
        column_identifier (str): name of column in table for extraction of
            identifiers
        name (str): name corresponding to selection criteria
        columns_categories (dict<list<str>>): specific columns for features and
            their categorical values by which to filter rows in table
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers from selection of rows from table

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    columns_categories = copy.deepcopy(columns_categories)
    # Filter rows in table by specific criteria.
    # Iterate on columns for features and their categorical values for
    # selection of rows from table.
    table_selection = table.copy(deep=True)
    if (columns_categories is not None):
        for column in columns_categories.keys():
            table_selection = table_selection.loc[(
                table_selection[column].isin(columns_categories[column])
            ), :].copy(deep=True)
            pass
        pass
    # Extract identifiers from selection of rows from table.
    identifiers = copy.deepcopy(
        table_selection[column_identifier].unique().tolist()
    )
    # Report.
    if report:
        count = len(identifiers)
        putly.print_terminal_partition(level=3)
        print("module: exercise.transcriptomics.organize_signal.py")
        function = (
            "filter_extract_table_row_identifiers_by_columns_categories" +
            "()"
        )
        print(str("function: " + function))
        putly.print_terminal_partition(level=5)
        print(str("name: " + name))
        print(str("count: " + str(count)))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return identifiers


def test_extract_organize_values_from_series():
    """
    Test:
    This function tests the function below.
    partner.organization.extract_organize_values_from_series()


    Review: TCW; 24 July 2024

    arguments:

    raises:

    returns:

    """

    # Dictionary.
    dictionary_integer = {
        "a": 1,
        "b": 2,
        "c": 3,
        "d": float("nan"),
        "e": 5,
        "f": 6,
        "g": 7,
        "h": 8,
        "i": 9,
        "j": 10,
    }
    dictionary_float = {
        "a": 0.001,
        "b": 0.005,
        "c": 0.01,
        "d": 0.05,
        "e": 0.1,
        "f": 0.5,
        "g": float("nan"),
        "h": 0.7,
        "i": 0.8,
        "j": 0.9,
    }
    # Series.
    series_integer = pandas.Series(
        data=dictionary_integer,
        index=[
            "a", "b", "c", "d", "e", "f", "g", "h", "i", "j",
        ],
    )
    series_float = pandas.Series(
        data=dictionary_float,
        index=[
            "a", "b", "c", "d", "e", "f", "g", "h", "i", "j",
        ],
    )
    # Call function to test.
    putly.print_terminal_partition(level=3)
    print("module: partner.organization.py")
    print("function: test_extract_organize_values_from_series()")
    print("series type: integer")
    putly.print_terminal_partition(level=5)
    pail_test = extract_organize_values_from_series(
        series=series_integer,
        threshold_low=1,
        threshold_high=9,
        report=True,
    )
    # Call function to test.
    putly.print_terminal_partition(level=3)
    print("module: partner.organization.py")
    print("function: test_extract_organize_values_from_series()")
    print("series type: float")
    putly.print_terminal_partition(level=5)
    pail_test = extract_organize_values_from_series(
        series=series_float,
        threshold_low=0.01,
        threshold_high=0.8,
        report=True,
    )
    pass


def extract_organize_values_from_series(
    series=None,
    threshold_low=None,
    threshold_high=None,
    report=None,
):
    """
    Dependency:
    This function is a dependency of the function below.
    partner.organization.determine_series_signal_validity_threshold()
    partner.organization.fill_missing_values_series()

    Extracts and organizes values from a Pandas series, such as a single column
    or row from a Pandas data-frame table.

    Review: TCW; 24 July 2024

    arguments:
        series (object): Pandas series of values of signal intensity
        threshold_low (float): threshold below which (value <= threshold) all
            values are considered invalid and missing
        threshold_high (float): threshold above which (value > threshold) all
            values are considered invalid and missing
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in series.
    series = series.copy(deep=True)
    # Extract and organize information from series.
    # Values, raw.
    #values_raw = values_raw.astype(numpy.float32)
    values_raw = series.to_numpy(
        dtype="float64",
        na_value=numpy.nan,
        copy=True,
    )
    # Values, nonmissing.
    values_nonmissing = numpy.copy(values_raw[~numpy.isnan(values_raw)])
    # Values, valid.
    if ((threshold_low is None) and (threshold_high is None)):
        values_valid = numpy.copy(values_nonmissing)
    elif ((threshold_low is not None) and (threshold_high is None)):
        values_temporary = numpy.copy(values_raw)
        values_temporary[values_temporary <= threshold_low] = numpy.nan
        values_valid = values_temporary[~numpy.isnan(values_temporary)]
    elif ((threshold_low is None) and (threshold_high is not None)):
        values_temporary = numpy.copy(values_raw)
        values_temporary[values_temporary > threshold_high] = numpy.nan
        values_valid = values_temporary[~numpy.isnan(values_temporary)]
    elif ((threshold_low is not None) and (threshold_high is not None)):
        values_temporary = numpy.copy(values_raw)
        values_temporary[values_temporary <= threshold_low] = numpy.nan
        values_temporary[values_temporary > threshold_high] = numpy.nan
        values_valid = values_temporary[~numpy.isnan(values_temporary)]
        pass
    # Values, logarithm.
    values_log = numpy.copy(values_valid)
    values_log = numpy.log(values_log)
    # Counts.
    #count_nonmissing = numpy.count_nonzero(~numpy.isnan(values_raw))
    #count_nonmissing = int(values_nonmissing.size)
    #count_positive = numpy.count_nonzero(values_nonmissing > 0)
    # Collect information.
    pail = dict()
    pail["values_raw"] = values_raw
    pail["values_nonmissing"] = values_nonmissing
    pail["values_valid"] = values_valid
    pail["values_log"] = values_log
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        print("function: extract_organize_values_from_series()")
        putly.print_terminal_partition(level=5)
        print("original source series:")
        print(series)
        putly.print_terminal_partition(level=5)
        print("values, raw:")
        print(pail["values_raw"])
        putly.print_terminal_partition(level=5)
        print("values, nonmissing:")
        print(pail["values_nonmissing"])
        putly.print_terminal_partition(level=5)
        print("values, nonmissing and within thresholds:")
        print("threshold, low (<=): " + str(threshold_low))
        print("threshold, high (>): " + str(threshold_high))
        print(pail["values_valid"])
        putly.print_terminal_partition(level=4)
        print("values, logarithm:")
        print(values_log)
    # Return information.
    return pail


def determine_series_signal_validity_threshold(
    series=None,
    keys_signal=None,
    threshold_low=None,
    threshold_high=None,
    proportion=None,
    report=None,
):
    """
    Dependency:
    This function is a dependency of the function below.
    exercise.transcriptomics.organize_signal.filter_table_main_rows_signal()

    Determines whether to keep a Pandas series on the basis of the
    proportion of values of signal intensity that are nonmissing and greater
    than zero. The series can come from the row or column in a table.

    Review: TCW; 24 July 2024

    arguments:
        series (object): Pandas series of values of signal intensity
        keys_signal (list<str>): names or identifiers of elements in the series
            that correspond to values of signal intensity
        threshold_low (float): threshold below which (value <= threshold) all
            values are considered invalid and missing
        threshold_high (float): threshold above which (value > threshold) all
            values are considered invalid and missing
        proportion (float): proportion of values of signal intensity that
            must be nonmissing and within low and high thresholds in order to
            keep the series
        report (bool): whether to print reports

    raises:

    returns:
        (int): logical binary representation of whether to keep current row

    """

    # Copy information in series.
    series = series.copy(deep=True)
    # Copy other information.
    keys_signal = copy.deepcopy(keys_signal)
    # Extract and organize information from series.
    pail_values = extract_organize_values_from_series(
        series=series[keys_signal],
        threshold_low=threshold_low,
        threshold_high=threshold_high,
        report=report,
    )
    count_raw = int(pail_values["values_raw"].size)
    count_valid = int(pail_values["values_valid"].size)
    proportion_actual = float(count_valid / count_raw)
    # Determine whether to keep current row from table.
    if (
        (proportion_actual >= proportion)
    ):
        indicator = 1
    else:
        indicator = 0
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        print("function: determine_series_signal_validity_threshold()")
        putly.print_terminal_partition(level=5)
        print("count of total values: " + str(count_raw))
        print(
            "count valid values, nonmissing and threshold: " + str(count_valid)
        )
        putly.print_terminal_partition(level=5)
        print("proportion valid, actual: " + str(proportion_actual))
        print("proportion valid, threshold: " + str(proportion))
        print("indicator: " + str(indicator))
        putly.print_terminal_partition(level=5)
    # Return information.
    return indicator


def segregate_fold_change_values_by_thresholds(
    table=None,
    column_fold_change=None,
    column_significance=None,
    threshold_fold_change=None,
    threshold_significance=None,
    report=None,
):
    """
    Dependency:
    This function is a dependency of the function below.
    partner.plot.plot_scatter_fold_change_volcano()

    Segregates data by values against thresholds on two dimensions.

    Table's format and orientation

    Table has each type of feature oriented across columns with their
    individual observations oriented across rows.

    features        feature_1 feature_2 feature_3 feature_4 feature_5
    observations
    observation_1   ...       ...       ...       ...       ...
    observation_2   ...       ...       ...       ...       ...
    observation_3   ...       ...       ...       ...       ...
    observation_4   ...       ...       ...       ...       ...
    observation_5   ...       ...       ...       ...       ...

    Review: TCW; 3 October 2024

    arguments:
        table (object): Pandas data-frame table of features across columns and
            values for observations across rows
        column_fold_change (str): name of column in table on which to apply the
            threshold for the fold change
        column_significance (str): name of column in table on which to apply
            the threshold for significance, corresponding to the p-value or
            q-value corresponding to the estimate of fold change
        threshold_fold_change (float): value for threshold on fold change
            (fold change > threshold) that is on the same scale, such as
            base-two logarithm, as the actual values themselves
        threshold_significance (float): value for threshold on p-values or
            q-values (p-value or q-value < threshold) that is not on a scale of
            the negative logarithm
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Filter rows in table for segregation of values by thresholds.
    table_pass_any = table.loc[
        (
            (table[column_significance] < threshold_significance) &
            (abs(table[column_fold_change]) > threshold_fold_change)
        ), :
    ].copy(deep=True)
    table_pass_up = table.loc[
        (
            (table[column_significance] < threshold_significance) &
            (table[column_fold_change] > threshold_fold_change)
        ), :
    ].copy(deep=True)
    table_pass_down = table.loc[
        (
            (table[column_significance] < threshold_significance) &
            (table[column_fold_change] < (-1*threshold_fold_change))
        ), :
    ].copy(deep=True)
    # Fail.
    table_fail = table.loc[
        ~table.index.isin(table_pass_any.index), :
    ].copy(deep=True)

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        print("function: segregate_fold_change_values_by_thresholds()")
        putly.print_terminal_partition(level=5)
    # Collect information.
    pail = dict()
    pail["table_pass_any"] = table_pass_any
    pail["table_pass_up"] = table_pass_up
    pail["table_pass_down"] = table_pass_down
    pail["table_fail"] = table_fail
    # Return information.
    return pail


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
    group_1_a   q_value      0.001       0.001       0.001
    group_1_b   signal       0.15        0.20        0.25
    group_1_b   p_value      0.001       0.001       0.001
    group_1_b   q_value      0.001       0.001       0.001

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
    if (len(types_value) > 2):
        table_three = table_source[
            table_source.index.get_level_values(
                name_index_type_value
            ).isin([types_value[2]])
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
    if (len(types_value) > 2):
        table_three_half = table_three.where(matrix_mask_half)

    ##########
    # Combine the half-diagonal tables for the different types of values.
    if (len(types_value) == 2):
        table_product = pandas.concat(
            [table_one_half, table_two_half,],
            axis="index",
            join="outer",
            ignore_index=False,
            copy=True,
        )
    elif (len(types_value) > 2):
        table_product = pandas.concat(
            [table_one_half, table_two_half, table_three_half,],
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


def filter_table_rows_by_quantile_least_variation_across_columns(
    table=None,
    name_columns=None,
    name_rows=None,
    logarithm=None,
    count_quantile=None,
    report=None,
):
    """
    Dependency:
    This function is a dependency of the function below.
    partner.organization.compare_stable_feature_sets_before_after_scale()

    Filters rows in table to select features that demonstrate least change
    across observations on the basis of the minimal quantile.

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

    Review: TCW; 17 September 2024

    arguments:
        table (object): Pandas data-frame table of values for observations
            across columns and for features across rows
        name_columns (str): name of single-level index across columns
        name_rows (str): name of single-level index across rows
        logarithm (bool): whether to transform to logarithmic scale before
            calculation of coefficient of variation
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
    # Calculate the logarithm of values.
    # math.log() # optimal for scalar values
    # numpy.log() # optimal for array values
    if logarithm:
        for column in columns:
            table[column] = numpy.log(
                table[column].to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
            ))
            pass
        pass

    # Calculate coefficient of variation of values for each feature across all
    # observations.
    table["variation_coefficient"] = table.apply(
        lambda row:
            scipy.stats.variation(
                row[columns].to_numpy(
                    dtype="float64",
                    na_value=numpy.nan,
                    copy=True,
                ),
                ddof=1, # divisor is (n - 1) for sample standard deviation
                nan_policy="omit",
            ),
        axis="columns", # apply function to each row
    )
    # Determine indices for quantiles.
    indices = list(range(0, count_quantile, 1))
    index_middle = indices[int(len(indices) // 2)]
    # Determine ordinal sets of features.
    table["quantile_variation"] = pandas.qcut(
        table["variation_coefficient"],
        q=count_quantile,
        labels=indices,
    )
    # Filter table to rows of features in the minimal quantile.
    table_quantile_minimum = table.loc[
        (table["quantile_variation"] == 0), :
    ]
    # Determine counts and percentages.
    count_rows_source = table.shape[0]
    count_rows_product = table_quantile_minimum.shape[0]
    percentage = str(round(
        float(100 * (count_rows_product / count_rows_source)), 2
    ))
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("module: partner.organization.py")
        print(
            "function: " +
            "filter_table_rows_by_quantile_least_variation_across_columns()"
        )
        putly.print_terminal_partition(level=5)
        print("count of quantiles: " + str(count_quantile))
        print("quantile indices:")
        print(indices)
        print("middle index: " + str(index_middle))
        putly.print_terminal_partition(level=5)
        print(
            "count of rows in original, source table: " +
            str(count_rows_source)
        )
        print(
            "count of rows in novel, product table: " +
            str(count_rows_product) + " (" + percentage + " %)"
        )
        print(
            "table of features with minimal coefficient of variation that " +
            "demonstrate least change across observations")
        print(table_quantile_minimum)
    # Return information.
    return table_quantile_minimum


def compare_stable_feature_sets_before_after_scale(
    table_raw=None,
    table_scale=None,
    name_columns=None,
    name_rows=None,
    count_quantile=None,
    report=None,
):
    """
    Compares the set overlap of stable features that have the least coefficient
    of variance from tables of signal values before and after adjustment of
    scale.

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

    Review: 17 September 2024

    arguments:
        table_raw (object): Pandas data-frame table of values for
            observations across columns and for features across rows
        table_scale (object): Pandas data-frame table of values for
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
    table_raw = table_raw.copy(deep=True)
    table_scale = table_scale.copy(deep=True)
    # Filter rows in table to select features that demonstrate least change
    # across observations on the basis of the minimal quantile on the
    # coefficient of variation.
    table_raw_stable = (
        filter_table_rows_by_quantile_least_variation_across_columns(
            table=table_raw,
            name_columns=name_columns,
            name_rows=name_rows,
            logarithm=True,
            count_quantile=count_quantile,
            report=report,
    ))
    table_scale_stable = (
        filter_table_rows_by_quantile_least_variation_across_columns(
            table=table_scale,
            name_columns=name_columns,
            name_rows=name_rows,
            logarithm=False,
            count_quantile=count_quantile,
            report=report,
    ))
    # Extract identifiers of features from the stable quantile of each table.
    features_raw = copy.deepcopy(
        table_raw.index.get_level_values(
            name_rows
        ).unique().to_list()
    )
    features_raw_stable = copy.deepcopy(
        table_raw_stable.index.get_level_values(
            name_rows
        ).unique().to_list()
    )
    features_scale = copy.deepcopy(
        table_scale.index.get_level_values(
            name_rows
        ).unique().to_list()
    )
    features_scale_stable = copy.deepcopy(
        table_scale_stable.index.get_level_values(
            name_rows
        ).unique().to_list()
    )
    # Determine overlap or union from complementary sets.
    list_sets_original = [
        set(features_raw),
        set(features_scale),
    ]
    list_sets_stable = [
        set(features_raw_stable),
        set(features_scale_stable),
    ]
    #union_original = set()
    #union_original = set(features_raw).union(set(features_scale))
    union_original = functools.reduce(
        lambda raw, scale: raw.union(scale),
        list_sets_original,
    )
    union_stable = functools.reduce(
        lambda raw, scale: raw.union(scale),
        list_sets_stable,
    )
    # Determine counts of features from each group.
    count_raw = len(features_raw)
    count_raw_stable = len(features_raw_stable)
    count_scale = len(features_scale)
    count_scale_stable = len(features_scale_stable)
    count_union_original = len(union_original)
    count_union_stable = len(union_stable)
    # Determine percentages.
    percentage_raw = str(round(
        float(100 * (count_raw_stable / count_raw)), 2
    ))
    percentage_scale = str(round(
        float(100 * (count_scale_stable / count_scale)), 2
    ))
    # Determine proportion of similarity between the stable sets of features.
    count_union_minimum = int(min([count_raw_stable, count_scale_stable]))
    count_union_maximum = int(count_raw_stable + count_scale_stable)
    range_union_stable = int(count_union_maximum - count_union_minimum)
    scale_union_stable = int(count_union_stable - count_union_minimum)
    proportion_union_stable = float(
        scale_union_stable / range_union_stable
    )
    percentage_union_stable = str(round(
        float(100 * proportion_union_stable), 2
    ))
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("module: partner.organization.py")
        print(
            "function: " +
            "compare_stable_feature_sets_before_after_scale()"
        )
        putly.print_terminal_partition(level=4)
        print("counts of unique features in original tables")
        print("raw table: " + str(count_raw))
        print("scale table: " + str(count_scale))
        print("union: " + str(count_union_original))
        putly.print_terminal_partition(level=4)
        print("counts of unique features in stable quantile of tables")
        print(
            "raw table: " + str(count_raw_stable) +
            " (" + percentage_raw + " %)"
        )
        print(
            "scale table: " + str(count_scale_stable) +
            " (" + percentage_scale + " %)"
        )
        print("union: " + str(count_union_stable))
        putly.print_terminal_partition(level=4)
        print(
            "expectable range of unique union features from stable quantiles"
        )
        print(
            "count of minimum possible unique union features: " +
            str(count_union_minimum)
        )
        print("This count would correspond to 100% similarity.")
        print("Minimal union size indicates that the two sets were identical.")
        print(
            "count of maximum possible unique union features: " +
            str(count_union_maximum)
        )
        print("This count would correspond to 0% similarity.")
        print(
            "Maximal union size indicates that the two sets were entirely "
            + "different."
        )
        putly.print_terminal_partition(level=4)
        print(
            "percentage of similarity in the features from stable quantiles"
        )
        print("percentage: " + percentage_union_stable + " %")
        putly.print_terminal_partition(level=4)
    # Collect information.
    pail = dict()
    pail["union_original"] = union_original
    pail["union_stable"] = union_stable
    # Return information.
    return pail


# Older functions to apply filters on tables.


def filter_rows_columns_by_threshold_proportion(
    data=None,
    dimension=None,
    threshold=None,
    proportion=None,
):
    """
    Filters either rows or columns by whether a certain proportion of values
    are greater than or equal to a minimal threshold.

    Persistence of a row or column requires at least a specific proportion of
    values beyond a specific threshold.

    Filter rows by consideration of values across columns in each row.
    Filter columns by consideration of values across rows in each column.

    arguments:
        data (object): Pandas data frame of values
        dimension (str): dimension to filter, either "row" or "column"
        threshold (float): minimal signal
        proportion (float): minimal proportion of rows or columns that must
            pass threshold

    raises:

    returns:
        (object): Pandas data frame of values

    """

    def count_true(slice=None, count=None):
        values = slice.tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Determine count from proportion.
    if dimension == "row":
        # Filter rows by consideration of columns for each row.
        columns = data.shape[1]
        count = round(proportion * columns)
    elif dimension == "column":
        # Filter columns by consideration of rows for each columns.
        rows = data.shape[0]
        count = round(proportion * rows)

    # Determine whether values exceed threshold.
    data_threshold = (data >= threshold)
    # Determine whether count of values exceed threshold.
    if dimension == "row":
        axis = "columns"
    elif dimension == "column":
        axis = "index"
    # This aggregation operation produces a series.
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis=axis,
    )

    # Select rows and columns with appropriate values.
    if dimension == "row":
        data_pass = data.loc[data_count, : ]
    elif dimension == "column":
        data_pass = data.loc[:, data_count]
    return data_pass


def filter_rows_columns_by_threshold_outer_proportion(
    data=None,
    dimension=None,
    threshold_high=None,
    threshold_low=None,
    proportion=None,
):
    """
    Filters either rows or columns by whether a certain proportion of values
    are outside of low and high thresholds.

    Persistence of a row or column requires at least a specific proportion of
    values beyond a specific threshold.

    Filter rows by consideration of values across columns in each row.
    Filter columns by consideration of values across rows in each column.

    arguments:
        data (object): Pandas data frame of values
        dimension (str): dimension to filter, either "row" or "column"
        threshold_high (float): value must be greater than this threshold
        threshold_low (float): value must be less than this threshold
        proportion (float): minimal proportion of rows or columns that must
            pass threshold

    raises:

    returns:
        (object): Pandas data frame of values

    """

    def count_true(slice=None, count=None):
        values = slice.values.tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Copy data.
    data = data.copy(deep=True)

    # Determine count from proportion.
    if dimension == "row":
        # Filter rows by consideration of columns for each row.
        columns = data.shape[1]
        count = round(proportion * columns)
    elif dimension == "column":
        # Filter columns by consideration of rows for each columns.
        rows = data.shape[0]
        count = round(proportion * rows)

    # Determine whether values exceed threshold.
    data_threshold = data.applymap(
        lambda value: ((value <= threshold_low) or (threshold_high <= value))
    )
    # Determine whether count of values exceed threshold.
    if dimension == "row":
        axis = "columns"
    elif dimension == "column":
        axis = "index"
    # This aggregation operation produces a series.
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis=axis,
    )

    # Select rows and columns with appropriate values.
    if dimension == "row":
        data_pass = data.loc[data_count, : ]
    elif dimension == "column":
        data_pass = data.loc[:, data_count]
    return data_pass


##########
# Fill missing values across rows.


def fill_missing_values_series(
    fill_missing=None,
    series=None,
    keys_signal=None,
    method=None,
    report=None,
):
    """
    Dependency:
    This function is a dependency of the function below.
    partner.organization.fill_missing_values_table_by_row()

    Determines whether the current row from a table needs fill of missing
    values.

    Within a normal distribution, 99.7% of values occur within 3 standard
    deviations of the mean, either above or below.
    3 standard deviations: 99.73%
    2 standard deviations: 95.45%
    1 standard deviations: 68.27%

    For the 'triple_standard_deviation_below' method, this function transforms
    raw values via natural logarithm before determining their mean and
    standard deviation so that these statistics are representative of a more
    normal distribution. This function calculates the intermediate fill value
    and then inverts the transformation by exponentiation.

    Reference:
    https://en.wikipedia.org/wiki/68-95-99.7_rule

    Note:
    On 27 June 2024, TCW confirmed that the procedure was actually recognizing
    and filling missing values appropriately. TCW also used a separate
    spreadsheet to check the math that determines the fill value via the
    'triple_standard_deviation_below' method. As expected, the fill value
    differed by whether the method transformed by logarithm before calculation
    of mean and standard deviation.

    Review: TCW; 26 July 2024

    arguments:
        fill_missing (int): logical binary representation of whether to fill,
            that is replace, missing values in series
        series (object): Pandas series of values of signal intensity
        keys_signal (list<str>): names or identifiers of elements in the series
            that correspond to values of signal intensity
        method (str): name of method to use for filling missing values, either
            'triple_standard_deviation_below', 'half_minimum', or 'zero'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas series of values of intensity across samples for
            a single protein

    """

    # Copy information in series.
    series = series.copy(deep=True)
    # Copy other information.
    keys_signal = copy.deepcopy(keys_signal)

    # Determine whether the current row has a missing value in the current
    # column.
    check_fill = 0
    value_fill = float("nan")
    if (fill_missing == 1):
        # Fill missing value.
        check_fill = 1
        # Extract and organize information from series.
        pail_values = extract_organize_values_from_series(
            series=series[keys_signal],
            threshold_low=0, # necessary to make logarithm possible
            threshold_high=None,
            report=report,
        )
        # Determine which method to use for fill.
        if (method == "triple_standard_deviation_below"):
            mean = numpy.nanmean(pail_values["values_log"])
            standard_deviation = numpy.nanstd(
                pail_values["values_log"],
                ddof=1, # divisor is (n - 1) for sample standard deviation
            )
            value_fill_log = float(mean - (3 * standard_deviation))
            value_fill = numpy.exp(value_fill_log)
            series_fill = series.replace(
                to_replace=pandas.NA,
                value=value_fill,
            )
        elif (method == "half_minimum"):
            minimum = numpy.nanmin(pail_values["values_nonmissing"])
            value_fill = float(minimum / 2)
            series_fill = series.replace(
                to_replace=pandas.NA,
                value=value_fill,
            )
        elif (method == "zero"):
            value_fill = 0.0
            series_fill = series.replace(
                to_replace=pandas.NA,
                value=value_fill,
            )
        # Report.
        if report:
            if (value_fill <= 0):
                putly.print_terminal_partition(level=1)
                print("**********")
                print("WARNING!")
                print("**********")
                print(
                    "Method selection to fill missing values returned a " +
                    "value less than or equal to zero!"
                )
                print("**********")
                putly.print_terminal_partition(level=1)
            putly.print_terminal_partition(level=4)
            print("fill missing value check: " + str(check_fill))
            print("fill value: " + str(value_fill))
            putly.print_terminal_partition(level=4)
            print("row with fill values: " + str(row_fill))
            putly.print_terminal_partition(level=4)
            pass
        pass
    else:
        series_fill = series
        pass

    # Report.
    if report:
        if (value_fill <= 0):
            putly.print_terminal_partition(level=1)
            print("**********")
            print("WARNING!")
            print("**********")
            print(
                "Method selection to fill missing values returned a value " +
                "less than or equal to zero!"
            )
            print("**********")
            putly.print_terminal_partition(level=1)
        putly.print_terminal_partition(level=4)
        print("fill missing value check: " + str(check_fill))
        print("fill value: " + str(value_fill))
        putly.print_terminal_partition(level=4)
        print("series with fill values: " + str(series_fill))
        putly.print_terminal_partition(level=4)
    # Return information.
    return series_fill


def fill_missing_values_table_by_row(
    table=None,
    columns=None,
    method=None,
    report=None,
):
    """
    Fills missing values in a table row by row with support for multiple
    alternative methods of imputation.

    Within a normal distribution, 99.7% of values occur within 3 standard
    deviations of the mean, either above or below.
    3 standard deviations: 99.73%
    2 standard deviations: 95.45%
    1 standard deviations: 68.27%

    Reference:
    https://en.wikipedia.org/wiki/68-95-99.7_rule

    Table's format and orientation

    Table has each type of feature oriented across rows with their individual
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

    Review: 26 July 2024

    arguments:
        table (object): Pandas data-frame table of values for observations
            across columns and for features across rows
        columns (list<str>): names of columns for which to fill missing values
        method (str): name of method to use for filling missing values, either
            'triple_standard_deviation_below', 'half_minimum', or 'zero'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values for observations across
            columns and for features across rows

    """

    # Copy information in table.
    table_fill = table.copy(deep=True)
    # Copy other information.
    columns = copy.deepcopy(columns)

    # Determine which rows in table have missing values.
    table_fill["match_missing"] = table_fill.apply(
        lambda row:
            1 if any(pandas.isna(row[columns]).to_list()) else 0,
        axis="columns", # apply function to each row
    )
    table_fill_actual = table_fill.loc[
        (table_fill["match_missing"] == 1), :
    ].copy(deep=True)

    # Fill missing values in all relevant columns and across all rows.
    # Apply the function to each row.
    table_fill = table_fill.apply(
        lambda row:
            fill_missing_values_series(
                fill_missing=row["match_missing"],
                series=row,
                keys_signal=columns,
                method=method,
                report=False,
            ),
        axis="columns", # apply function to each row
    )
    # Remove unnecessary columns.
    table_fill.drop(
        labels=["match_missing",],
        axis="columns",
        inplace=True
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        print("function: fill_missing_values_table_by_row()")
        putly.print_terminal_partition(level=5)
        count_rows_missing = (table_fill_actual.shape[0])
        count_rows_total = (table_fill.shape[0])
        proportion_missing = round(
            float(count_rows_missing / count_rows_total), 2
        )
        print(
            "count of rows requiring missing fill: " + str(count_rows_missing)
        )
        print(
            "proportion of rows requiring missing fill: " +
            str(proportion_missing)
        )
        putly.print_terminal_partition(level=4)
        print("table after fill: ")
        print(table_fill)
        putly.print_terminal_partition(level=4)
    # Return information.
    return table_fill


##########
# Combine tables.

# TODO: TCW; 1 July 2024
# From within these two "merge..." functions, call function "filter_sort_table_columns()"
# to clean up the table after the combination.


def merge_columns_tables_supplements_to_main(
    identifier_main=None,
    identifier_supplement=None,
    table_main=None,
    tables_supplements=None,
    report=None,
):
    """
    Merge columns from multiple source tables to columns of a single main table.

    This function does not discard any existing indices of the main and
    supplement tables. This behavior might produce redundant or unnecessary
    columns.

    arguments:
        identifier_main (str): name of column in main table on which
            to merge
        identifier_supplement (str): name of column in supplement tables on
            which to merge
        table_main (object): Pandas data frame of information
        tables_supplements (list<object>): collection of Pandas data frames of
            information
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information

    """

    # Copy information in table.
    table_main = table_main.copy(deep=True)

    # Reset index in original main table.
    table_main.reset_index(
        level=None,
        inplace=True,
        drop=False, # move index to regular columns
    )

    # Collect columns in original main table.
    # This collection must include the index.
    columns_collection = copy.deepcopy(table_main.columns.to_list())

    # Organize main table's index.
    table_main[identifier_main] = table_main[identifier_main].astype(
        "string",
        copy=True,
        errors="raise",
    )
    table_main.set_index(
        identifier_main,
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True,
    )

    # Copy information in table.
    table_merge = table_main.copy(deep=True)

    # Iterate on tables for polygenic scores.
    for table_supplement in tables_supplements:
        # Reset index in original supplement table.
        table_supplement.reset_index(
            level=None,
            inplace=True,
            drop=False, # remove index; do not move to regular columns
        )
        # Collect columns in original supplement table.
        # This collection must include the index.
        columns_supplement = copy.deepcopy(table_supplement.columns.to_list())
        columns_collection.extend(columns_supplement)
        # Organize supplement table's index.
        table_supplement[identifier_supplement] = (
            table_supplement[identifier_supplement].astype(
                "string",
                copy=True,
                errors="raise",
        ))
        table_supplement.set_index(
            identifier_supplement,
            append=False,
            drop=True, # move regular column to index; remove original column
            inplace=True
        )
        # Merge data tables using database-style join.
        # Alternative is to use DataFrame.join().
        table_merge = pandas.merge(
            table_merge, # left table
            table_supplement, # right table
            left_on=None, # "identifier_main"
            right_on=None, # "identifier_supplement"
            left_index=True, # "identifier_main"
            right_index=True, # "identifier_supplement"
            how="outer", # keep union of keys from both tables
            #suffixes=("_main", "_identifiers"), # deprecated?
        )
        pass
    # Organize table's index.
    table_merge.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    #table_merge.set_index(
    #    "identifier_genotype",
    #    append=False,
    #    drop=True, # move regular column to index; remove original column
    #    inplace=True
    #)
    # Determine unique columns from collection.
    columns_collection = list(set(columns_collection))
    # Remove unnecessary columns from transformations on tables.
    # The goal is to avoid any columns from the merge operations themselves.
    #table_merge.drop(
    #    labels=["index_x", "index_y", "index",],
    #    axis="columns",
    #    inplace=True
    #)
    table_merge = table_merge.loc[
        :, table_merge.columns.isin(columns_collection)
    ]
    # Report.
    if report:
        print_terminal_partition(level=2)
        print("report: ")
        print("merge_columns_tables_supplements_to_main()")
        print_terminal_partition(level=3)
        print("table columns: " + str(int(table_merge.shape[1])))
        print("table rows: " + str(int(table_merge.shape[0])))
        print("columns")
        print(table_merge.columns.to_list())
        print(table_merge)
        pass
    # Return information.
    return table_merge


def merge_columns_two_tables(
    identifier_first=None,
    identifier_second=None,
    table_first=None,
    table_second=None,
    preserve_index=None,
    report=None,
):
    """
    Merge columns from two tables.

    This function does not discard any existing indices of the main and
    supplement tables. This behavior might produce redundant or unnecessary
    columns.

    This function also does not check that columns in both tables have unique
    names.

    arguments:
        identifier_first (str): name of column in first table on which to merge
        identifier_second (str): name of column in second table on which to
            merge
        table_first (object): Pandas data frame of information
        tables_second (object): Pandas data frame of information
        preserve_index (bool): whether to preserve the indices across rows in
            the original tables
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information

    """

    # Copy information in table.
    table_first = table_first.copy(deep=True)
    table_second = table_second.copy(deep=True)

    # Reset indices in original tables.
    drop_index = (not preserve_index)
    table_first.reset_index(
        level=None,
        inplace=True,
        drop=drop_index, # do not remove index; move to regular columns
    )
    table_second.reset_index(
        level=None,
        inplace=True,
        drop=drop_index, # do not remove index; move to regular columns
    )

    # Collect columns in original tables.
    # This collection must include the indices.
    columns_collection = copy.deepcopy(table_first.columns.to_list())
    columns_second = copy.deepcopy(table_second.columns.to_list())
    columns_collection.extend(columns_second)

    # Organize first table's index.
    table_first[identifier_first] = table_first[identifier_first].astype(
        "string",
        copy=True,
        errors="raise",
    )
    table_first.set_index(
        identifier_first,
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True,
    )
    # Organize second table's index.
    table_second[identifier_second] = table_second[identifier_second].astype(
        "string",
        copy=True,
        errors="raise",
    )
    table_second.set_index(
        identifier_second,
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True,
    )

    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table = pandas.merge(
        table_first, # left table
        table_second, # right table
        left_on=None, # "identifier_first"
        right_on=None, # "identifier_second"
        left_index=True, # "identifier_first"
        right_index=True, # "identifier_second"
        how="outer", # keep union of keys from both tables
        #suffixes=("_main", "_identifiers"), # deprecated?
    )

    # Organize table's index.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    #table.set_index(
    #    "identifier_genotype",
    #    append=False,
    #    drop=True, # move regular column to index; remove original column
    #    inplace=True
    #)
    # Determine unique columns from collection.
    columns_collection = list(set(columns_collection))
    # Remove unnecessary columns from transformations on tables.
    # The goal is to avoid any columns from the merge operations themselves.
    #table_merge.drop(
    #    labels=["index_x", "index_y", "index",],
    #    axis="columns",
    #    inplace=True
    #)
    table = table.loc[
        :, table.columns.isin(columns_collection)
    ]

    # Report.
    if report:
        putly.print_terminal_partition(level=2)
        print("report: ")
        print("merge_columns_two_tables()")
        putly.print_terminal_partition(level=3)
        print("table columns: " + str(int(table.shape[1])))
        print("table rows: " + str(int(table.shape[0])))
        print("columns")
        print(table.columns.to_list())
        print(table)
        pass
    # Return information.
    return table


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


# Cluster columns or rows in a table.


def cluster_table_columns(
    table=None,
):
    """
    Cluster the sequence of a table's columns by similarity in their values
    across rows.

    Table must have an explicitly defined index across rows to designate that
    remaining columns are those by which to cluster.

    Reference:
    https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    https://stackabuse.com/hierarchical-clustering-with-python-and-scikit-
        learn/
        - description of linkage methods in hierarchical clustering

    Review: TCW; 17 September 2024

    arguments:
        table (object): Pandas data-frame table of floating-point values on
            continuous interval or ratio scales of measurement

    raises:

    returns:
        (object): Pandas data-frame table of values

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Extract names of columns and rows in table.
    #columns = copy.deepcopy(table.columns.to_list())
    columns = copy.deepcopy(table.columns.to_numpy())
    rows = copy.deepcopy(table.index.to_numpy())
    index_name = table.index.name

    # Cluster table's columns.
    # Organize columns across dimension zero.
    matrix = numpy.transpose(table.to_numpy())
    linkage = scipy.cluster.hierarchy.linkage(
        matrix,
        method="average", # "average", "centroid", "ward",
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
    # Organize table.
    table_cluster = pandas.DataFrame(
        data=numpy.transpose(matrix_cluster),
        index=rows,
        columns=columns_sort,
    )
    table_cluster.rename_axis(
        index_name,
        axis="index",
        inplace=True,
    )
    # Return information.
    return table_cluster


def cluster_table_columns_by_external_group(
    table=None,
    indices_rows=None,
    groups_columns=None,
    names_groups_sequence=None,
    report=None,
):
    """
    Cluster the sequence of a table's columns within separate groups by
    similarity in their values across columns.

    A separate dictionary defines the columns in each group, and the product
    table will only include columns that this dictionary assigns to groups.
    This definition of columns in each group must not include the columns that
    define the index across rows.

    A separate list of names of groups enforces a specific sequence.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 15 November 2024

    arguments:
        table (object): Pandas data-frame table of floating-point values on
            continuous interval or ratio scales of measurement
        indices_rows (list<str>): names of columns in source table which
            define an index corresponding to information across rows
        groups_columns (dict<list<str>>): names of columns in groups
        names_groups_sequence (list<str>): names of groups in sequence
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of values

    """

    # Copy information in table.
    table_source = table.copy(deep=True)
    table_temporary = table.copy(deep=True)
    # Copy other information.
    indices_rows = copy.deepcopy(indices_rows)
    groups_columns = copy.deepcopy(groups_columns)
    names_groups_sequence = copy.deepcopy(names_groups_sequence)
    # Organize indices in table.
    table_temporary.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_temporary.set_index(
        indices_rows,
        append=False,
        drop=True,
        inplace=True
    )
    # Create empty table with identical index for collection of columns.
    table_collection = pandas.DataFrame(
        index=table_temporary.index.copy(),
    )
    # Iterate on groups of columns, apply operations, and collect information
    # from each.
    for group in names_groups_sequence:
        # Copy information in table with selection of columns in group.
        table_group = table_temporary.loc[
            :, table_temporary.columns.isin(groups_columns[group])
        ].copy(deep=True)
        table_group = table_group[[*groups_columns[group]]]
        # Cluster table's columns.
        table_cluster = cluster_table_columns(
            table=table_group,
        )
        # Organize indices in table.
        table_cluster.index = pandas.MultiIndex.from_tuples(
            table_cluster.index,
            names=indices_rows,
        )
        # Collect columns from group table with collection table.
        # Merge data tables using database-style join.
        # Use Pandas "DataFrame.join()" to preserve sequence of original source
        # index across rows.
        # Otherwise it would be necessary to sort rows in product table.
        table_collection = table_collection.join(
            table_cluster,
            on=list(table_collection.index.names),
        )
        if False:
            table_collection = pandas.merge(
                table_collection, # left table
                table_cluster, # right table
                left_on=None, # "identifier_first"
                right_on=None, # "identifier_second"
                left_index=True, # "identifier_first"
                right_index=True, # "identifier_second"
                how="outer", # keep union of keys from both tables
                #suffixes=("_main", "_identifiers"), # deprecated?
            )
        pass
    # Copy information in table.
    table_product = table_collection.copy(deep=True)
    # Sort rows in table.
    # This strategy does not work for an index with multiple levels.
    # Attempt to collapse the indices, sort, and then split the indices.
    #table_product = table_product.reindex(table_temporary)
    # Organize indices in table.
    table_product.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        print("function: cluster_table_columns_by_external_group()")
        putly.print_terminal_partition(level=4)
        # Copy information in table.
        table_source_report = table_source.copy(deep=True)
        table_product_report = table_product.copy(deep=True)
        # Organize indices in table.
        table_source_report.set_index(
            indices_rows,
            append=False,
            drop=True,
            inplace=True
        )
        table_product_report.set_index(
            indices_rows,
            append=False,
            drop=True,
            inplace=True
        )
        # Extract information.
        count_columns_source = (table_source_report.shape[1])
        count_rows_source = (table_source_report.shape[0])
        count_columns_product = (table_product_report.shape[1])
        count_rows_product = (table_product_report.shape[0])
        # Collapse and extract multi-level index across rows for comparison.
        values_index_rows_source = copy.deepcopy(
            table_source_report.index.map("|".join).to_list()
        )
        values_index_rows_product = copy.deepcopy(
            table_product_report.index.map("|".join).to_list()
        )
        # Confirm that sets of indices from both sources are inclusive.
        check_inclusion = putly.compare_lists_by_mutual_inclusion(
            list_primary=values_index_rows_source,
            list_secondary=values_index_rows_product,
        )
        # Confirm that sets of indices from both sources are identical across
        # their respective sequences.
        check_identity = putly.compare_lists_by_elemental_identity(
            list_primary=values_index_rows_source,
            list_secondary=values_index_rows_product,
        )
        # Confirm that sets of indices from both sources are equal.
        check_equality = (
            values_index_rows_source == values_index_rows_product
        )
        putly.print_terminal_partition(level=4)
        print("source table before cluster")
        print("count columns: " + str(count_columns_source))
        print("count rows: " + str(count_rows_source))
        putly.print_terminal_partition(level=5)
        print("product table after cluster")
        print("count columns: " + str(count_columns_product))
        print("count rows: " + str(count_rows_product))
        putly.print_terminal_partition(level=4)
        print(
            "confirm that indices across rows are identical between source " +
            "and product tables"
        )
        print(
            "check coherence of indices across rows before and after " +
            "cluster operation"
        )
        print("check inclusion: " + str(check_inclusion))
        print("check identity: " + str(check_identity))
        print("check equality: " + str(check_equality))
        putly.print_terminal_partition(level=5)
        print("product table")
        print(table_product)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_product


def cluster_table_rows(
    table=None,
):
    """
    Cluster the sequence of a table's rows by similarity in their values
    across columns.

    Table must have an explicitly defined index across rows to designate that
    remaining columns are those by which to cluster.

    Reference:
    https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    https://stackabuse.com/hierarchical-clustering-with-python-and-scikit-
        learn/
        - description of linkage methods in hierarchical clustering

    Review: TCW; 17 September 2024

    arguments:
        table (object): Pandas data-frame table of floating-point values on
            continuous interval or ratio scales of measurement

    raises:

    returns:
        (object): Pandas data-frame table of values

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Extract names of columns and rows in table.
    #columns = copy.deepcopy(table.columns.to_list())
    columns = copy.deepcopy(table.columns.to_numpy())
    rows = copy.deepcopy(table.index.to_numpy())
    index_name = table.index.name

    # Cluster table's rows.
    # Organize rows across dimension zero.
    matrix = table.to_numpy()
    linkage = scipy.cluster.hierarchy.linkage(
        matrix,
        method="average", # "average", "centroid", "ward",
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
    # Organize table.
    table_cluster = pandas.DataFrame(
        data=matrix_cluster,
        index=rows_sort,
        columns=columns,
    )
    table_cluster.rename_axis(
        index_name,
        axis="index",
        inplace=True,
    )
    # Return information.
    return table_cluster


def cluster_table_rows_by_group(
    table=None,
    index_rows=None,
    column_group=None,
):
    """
    Cluster the sequence of a table's rows within separate groups by similarity
    in their values across columns.

    Original source table must not have an explicitly defined index across
    rows.

    The split-apply-combine operation disrupts the original sequence of groups
    in the source table. Refer to a sort operation such as in the functions
    below.
    partner.organization.sort_table_rows_by_single_column_reference()

    Review: TCW; 17 October 2024

    arguments:
        table (object): Pandas data-frame table of floating-point values on
            continuous interval or ratio scales of measurement
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        column_group (str): name of column to use for groups

    raises:

    returns:
        (object): Pandas data-frame table of values

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Copy other information.
    columns_sequence = copy.deepcopy(table.columns.to_list())
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table.set_index(
        [index_rows, column_group],
        append=False,
        drop=True,
        inplace=True
    )
    # Split rows within table by factor columns.
    groups = table.groupby(
        level=column_group,
    )
    # Iterate on groups, apply operations, and collect information from each.
    table_collection = pandas.DataFrame()
    for name, table_group in groups:
        # Copy information in table.
        table_group = table_group.copy(deep=True)
        # Cluster table's rows.
        table_cluster = cluster_table_rows(
            table=table_group,
        )
        # Organize indices in table.
        table_cluster.index = pandas.MultiIndex.from_tuples(
            table_cluster.index,
            names=[index_rows, column_group]
        )
        table_cluster.reset_index(
            level=None,
            inplace=True,
            drop=False, # remove index; do not move to regular columns
        )
        # Collect table for current group.
        table_collection = pandas.concat(
            [table_collection, table_cluster],
            ignore_index=True,
        )
        pass
    # Filter and sort columns within table.
    table_collection = filter_sort_table_columns(
        table=table_collection,
        columns_sequence=columns_sequence,
        report=False,
    )
    # Return information.
    return table_collection


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


################################################################################
# Procedure
# Currently, this module is not directly executable.

################################################################################
# End
