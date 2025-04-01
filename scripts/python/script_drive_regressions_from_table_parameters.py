"""
Drive multiple regressions from a single table of parameters.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This module file is part of the project package directory 'partner'
    (https://github.com/tcameronwaller/partner/).

    Project 'partner' supports data analysis in multiple other projects.
    Copyright (C) 2025 Thomas Cameron Waller

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
# Author: T. Cameron Waller, Ph.D.
# Date, first execution: 28 March 2025
# Date, last execution: 28 March 2025
# Review: TCW; 28 March 2025
################################################################################
# Note

# The specialty of this Python script is to drive multiple regressions from
# parameters within a single table. This script calls versatile functionality
# from the "regression.py" module within the "partner" Python package.


# Eventually, this script ought to execute regressions with parallel processing...


################################################################################
# Installation and importation

# Standard
import sys
import os
import copy

# Relevant
import pandas
import scipy
import numpy

# Custom
import partner.utility as putly
import partner.parallelization as prall
import partner.regression as preg

#dir()
#importlib.reload()

###############################################################################
# Functionality


# Read source information from file.


def define_column_types_table_parameters():
    """
    Defines the types of variables for columns in table of parameters.

    Review: TCW; 31 March 2025

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify types of variables in columns of table.
    types_columns = dict()
    types_columns["execution"] = "string" # "int32"
    types_columns["sequence"] = "string" # "int32"
    types_columns["group"] = "string"
    types_columns["instance"] = "string"
    types_columns["selection_observations"] = "string"
    types_columns["features_continuity_scale"] = "string"
    types_columns["type_regression"] = "string"
    types_columns["formula_text"] = "string"
    types_columns["feature_response"] = "string"
    types_columns["features_predictor_fixed"] = "string"
    types_columns["features_predictor_random"] = "string"
    types_columns["data_path_directory"] = "string"
    types_columns["data_file"] = "string"
    types_columns["review"] = "string"
    types_columns["note"] = "string"
    # Return information.
    return types_columns


def read_source(
    path_file_table_parameters=None,
    report=None,
):
    """
    Read and organize source information about parameters for regressions.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 31 March 2025

    arguments:
        path_file_table_parameters (str): path to source file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters for selection
            of sets from MSigDB

    """

    # Read information from file.

    # Table of parameters for parallel instances.
    types_columns = define_column_types_table_parameters()
    table = pandas.read_csv(
        path_file_table_parameters,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Organize information.
    table["execution"] = pandas.to_numeric(
        table["execution"],
        downcast="integer",
        errors="coerce",
    )
    table["sequence"] = pandas.to_numeric(
        table["sequence"],
        downcast="integer",
        errors="coerce",
    )

    # Collect information.
    records = list()
    for index, row in table.iterrows():
        if (int(row["execution"]) == 1):
            # Collect information and parameters from current row in table.
            record = dict()
            record["execution"] = row["execution"]
            record["sequence"] = row["sequence"]
            record["group"] = str(row["group"]).strip()
            record["instance"] = str(row["instance"]).strip()
            record["name_instance"] = "_".join([
                str(row["group"]).strip(),
                str(row["sequence"]).strip(),
                str(row["instance"]).strip(),
            ])
            record["selection_observations"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["selection_observations"],
                )
            )["features_values"]
            record["features_continuity_scale"] = putly.parse_text_list_values(
                text=row["features_continuity_scale"],
                delimiter=",",
            )
            record["type_regression"] = str(row["type_regression"]).strip()
            record["formula_text"] = str(row["formula_text"]).strip()
            record["feature_response"] = str(row["feature_response"]).strip()
            record["features_predictor_fixed"] = putly.parse_text_list_values(
                text=row["features_predictor_fixed"],
                delimiter=",",
            )
            record["features_predictor_random"] = putly.parse_text_list_values(
                text=row["features_predictor_random"],
                delimiter=",",
            )
            record["data_path_directory"] = str(
                row["data_path_directory"]
            ).strip()
            record["data_file"] = str(row["data_file"]).strip()
            record["review"] = str(row["review"]).strip()
            record["note"] = str(row["note"]).strip()

            # Collect unique names of columns relevant to instance of
            # parameters from current row in table.
            features_relevant = list()
            dictionaries = [
                "selection_observations",
            ]
            for dictionary in dictionaries:
                if record[dictionary] is not None:
                    features_relevant.extend(list(record[dictionary].keys()))
                    pass
                pass
            features_relevant.extend(record["features_continuity_scale"])
            features_relevant.append(record["feature_response"])
            features_relevant.extend(record["features_predictor_fixed"])
            features_relevant.extend(record["features_predictor_random"])
            features_unique = putly.collect_unique_elements(
                elements=features_relevant,
            )
            record["features_relevant"] = copy.deepcopy(features_unique)

            # Collect information and parameters for current row in table.
            records.append(record)
            pass
        pass

    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["records"] = records
    # Report.
    if report:
        # Organize.
        count_records = len(records)
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("parameter table:")
        print(table)
        putly.print_terminal_partition(level=5)
        print("count of records or instances: " + str(count_records))
        putly.print_terminal_partition(level=5)
        print("instance[0]:")
        print(records[0])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail




################################################################################
# Procedure


def execute_procedure(
    path_file_table_parameters=None,
    path_directory_dock=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_file_table_parameters (str): path to source file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parameters.
    if (
        (report is not None) and
        (str(report) != "") and
        (str(report) == "true")
    ):
        report = True
    else:
        report = False
        pass

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_file_table_parameters=path_file_table_parameters,
        report=report,
    )
    for record in pail_source["records"]:
        print(record)
        pass


    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("path_file_table_parameters: " + str(path_file_table_parameters))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass


    ##########
    # Write product information to file.

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_table_parameters = sys.argv[1]
    path_directory_dock = sys.argv[2]
    report = sys.argv[3]

    # Call function for procedure.
    execute_procedure(
        path_file_table_parameters=path_file_table_parameters,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass



#
