"""
Supply functionality for reading, organizing, and processing information from
polygenic scores.

This module 'polygenic_score' is part of the 'partner' package.

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

# TODO: TCW; 10 April 2023
# TODO: I think that this module is now obsolete.
# TODO: Please compare to the Python script at the path below.
# partner/scripts/python/script_collect_standardize_polygenic_scores.py


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
import partner.utility as utility


#dir()
#importlib.reload()

###############################################################################
# Functionality

##########
# Read table of paths and parameters for collection of polygenic scores


def read_source_collection_polygenic_scores(
    path_table=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    This function reads from file a table of information (including paths and
    parameters) about files of a collection of polygenic scores across relevant
    genotypes.


    The column "method" is most relevant as a specification of the format for
    the output of the polygenic scores. It is important to know this format for
    reading the table.

    arguments:
        path_table (str): full path to file for table with paths and parameters
            for relevant polygenic scores
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of parameters for reading from file and
            organizing separate tables for polygenic scores

    """

    # Read information from file.
    table_parameter_scores = pandas.read_csv(
        path_table,
        sep="\t", # ","; "\t"; "\s+"; "\s+|\t+|\s+\t+|\t+\s+"
        header=0,
        dtype={
            "inclusion": "int32",
            "name_file_path_directory_parent": "string",
            "path_directory": "string",
            "name_file": "string",
            "method": "string",
            "name_column_score_old": "string",
            "name_column_score_new": "string",
        },
    )
    table_parameter_scores.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Return information.
    return table_parameter_scores


##########
# Read and organize tables of polygenic scores


def read_source_table_polygenic_score_ldpred2(
    name_file_directory_parent=None,
    path_directory=None,
    name_file=None,
    name_column_score=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    This function reads a private file path from the file system on which it
    runs. Within the user's "~/paths" directory, there must be a collection of
    files in text format with file paths. This function reads the file with the
    name in the parameter "name_file_directory_parent".

    arguments:
        name_file_directory_parent (str): name of file with private path to
            parent directory
        path_directory (str): path beyond parent directory to the file
        name_file (str): name of file
        name_column_score (str): name of column in table for polygenic scores
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of polygenic scores

    """

    # Read private path to parent directory.
    path_home_paths = os.path.expanduser("~/paths")
    path_file_directory_parent = os.path.join(
        path_home_paths, name_file_directory_parent
    )
    path_directory_parent = utility.read_file_text(
        path_file=path_file_directory_parent
    ).rstrip("\n") # remove new line character from string
    # Specify directories and files.
    path_table = os.path.join(
        path_directory_parent, path_directory, name_file
    )
    # Read information from file.
    table = pandas.read_csv(
        path_table,
        sep="\s+", # "\t"; "\s+"; "\s+|\t+|\s+\t+|\t+\s+"
        header=0,
        dtype={
            "FID": "string",
            "IID": "string",
            name_column_score: "float32",
        },
    )
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Compile and return information.
    return table


def organize_table_polygenic_score_ldpred2(
    table=None,
    name_score_old=None,
    name_score_new=None,
    report=None,
):
    """
    Organizes table of polygenic scores.

    arguments:
        table (object): Pandas data frame of information about polygenic scores
        name_score_old (str): old, original name for polygenic score variable
        name_score_new (str): new, novel name for polygenic score variable
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotypes

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Convert all identifiers to type string.
    table["IID"] = table["IID"].astype(
        "string",
        copy=True,
        errors="raise",
    )
    # Replace any empty identifier strings with missing values.
    table["IID"].replace(
        "",
        numpy.nan,
        inplace=True,
    )
    # Remove any records with missing identifiers.
    table.dropna(
        axis="index", # drop rows with missing values in columns
        how="any",
        subset=["IID",],
        inplace=True,
    )
    # Convert identifiers to type string.
    table["IID"] = table["IID"].astype(
        "string",
        copy=True,
        errors="raise",
    )
    table["identifier_genotype"] = table["IID"].astype(
        "string",
        copy=True,
        errors="raise",
    ).copy(deep=True)
    # Translate names of columns.
    #table = table.add_prefix("import_")
    translations = dict()
    translations[name_score_old] = name_score_new
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Convert scores to type float.
    #table[name_score_new] = table[name_score_new].astype("float32")
    table[name_score_new] = pandas.to_numeric(
        table[name_score_new],
        errors="coerce", # force any invalid values to missing or null
        downcast="float",
    )
    # Remove columns.
    table.drop(
        labels=["FID", "IID"],
        axis="columns",
        inplace=True
    )
    # Return information.
    return table


def drive_read_organize_tables_polygenic_scores(
    table_parameter_scores=None,
    filter_inclusion=None,
    report=None,
):
    """
    Drives functions to read and organize source information from file.

    This function extracts from a table information (including paths and
    parameters) about files of a collection of polygenic scores across relevant
    genotypes.

    arguments:
        table_parameter_scores (object): Pandas data frame of parameters for
           reading from file and organizing separate tables for polygenic scores
        filter_inclusion (bool): whether to filter records in polygenic scores
            parameter table by logical binary "inclusion" variable
        report (bool): whether to print reports

    raises:

    returns:
        (list<object>): collection of Pandas data frames of polygenic scores

    """

    # Filter collection of polygenic scores by "inclusion" flag.
    if filter_inclusion:
        table_parameter_scores = table_parameter_scores.loc[
            (
                (table_parameter_scores["inclusion"] == 1)
            ), :
        ]
        pass
    # Extract records from table.
    records = table_parameter_scores.to_dict(
        orient="records",
    )
    # Collect tables of polygenic scores.
    tables_polygenic_scores = list()
    # Iterate on records.
    for record in records:
        # Extract information from record.
        name_file_directory_parent = record["name_file_path_directory_parent"]
        path_directory = record["path_directory"]
        name_file = record["name_file"]
        method = record["method"]
        name_column_score_old = record["name_column_score_old"]
        name_column_score_new = record["name_column_score_new"]
        # Determine appropriate functions for format of polygenic scores.
        if (method == "LDPred2"):
            # Read file.
            table_raw = read_source_table_polygenic_score_ldpred2(
                name_file_directory_parent=name_file_directory_parent,
                path_directory=path_directory,
                name_file=name_file,
                name_column_score=name_column_score_old,
                report=report,
            )
            # Organize table.
            table = organize_table_polygenic_score_ldpred2(
                table=table_raw,
                name_score_old=name_column_score_old,
                name_score_new=name_column_score_new,
                report=report,
            )
            # Collect table.
            tables_polygenic_scores.append(table)
        elif (method == "PRS-CS"):
            print("still need to implement read and organize for PRS-CS...")
            pass
        # Report.
        if report:
            utility.print_terminal_partition(level=2)
            print("report: ")
            print("drive_read_organize_tables_polygenic_scores()")
            utility.print_terminal_partition(level=3)
            print(table)
            print("columns")
            print(table.columns.to_list())
            pass
    # Return information.
    return tables_polygenic_scores



###############################################################################
# Procedure
# Currently, this module is not executable.



##########
