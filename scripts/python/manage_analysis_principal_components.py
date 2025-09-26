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
# Date, first execution: 17 September 2025
# Date, last execution or modification: 17 September 2025
# Review: TCW; 17 September 2025
################################################################################
# Note

# The specialty of this Python script is to manage a general and versatile
# design of principal components analysis.

##########
# Note:

# TODO: TCW; 17 September 2025
# This file is currently mostly an empty place holder.
# Look to "drive_regressions_from_table_parameters.py" as a reference.


##########
# Review:

################################################################################
# Installation and importation

# Standard
import sys
# sys.exit() # End execution at this point.
import os
import copy
import textwrap

# Relevant
import pandas
import scipy
import numpy

# Custom
import partner.utility as putly
import partner.parallelization as prall
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
import partner.regression as preg

#dir()
#importlib.reload()

###############################################################################
# Functionality





################################################################################
# Procedure


##########
# Manage principal components analysis.


# arguments
# table
# selection_observations (ex: only older persons)
# groups_observations (ex: females, males)
# features (categorical and continuous)
#


##########
# Execute main procedure.

# path_file_source_table
# path_file_source_list_features


def execute_procedure(
    groups_raw=None,
    path_file_source_table=None,
    path_file_product_table_regressions=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    1. read source information from file
       - source information at this level consists of a table of parameters for
          regressions
       - it is necessary to parse information with hierarchical structure from
          the flat text
    2. drive multiple, concurrent procedures in parallel branches
    3. within each procedural branch, read from file and organize the data
       table consisting of features across observations
       - filter features and observations in table
       - standardize scale of values within specific feature variables
    4. within each procedural branch, execute regression analysis using type of
       model, features, and other parameters from the source table of
       parameters
    5. organize report summary information from the regression analysis and
       write this information to file

    arguments:
        groups_raw (str): names of groups of instances of parameters to include
            in batch for execution
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with parameters for multiple regressions
        path_file_product_table_regressions (str): path to product file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with results of multiple regressions
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
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
    pail_source = read_source_table_parameters(
        groups_raw=groups_raw,
        path_file_source_table_parameters=path_file_source_table_parameters,
        path_directory_dock=path_directory_dock,
        report=report,
    )
    #for record in pail_source["records"]:
    #    print(record)
    #    pass

    # Expand plural response features.
    pail_expansion = expand_plural_response_features(
        instances=pail_source["records"],
        path_directory_dock=path_directory_dock,
        report=report,
    )
    #for instance in pail_expansion["instances"]:
    #    print(instance["feature_response"])
    #    print(instance)
    #    pass

    ##########
    # Organize information.
    count_instances = len(pail_expansion["instances"])
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("path_file_source_table_parameters: " + str(
            path_file_source_table_parameters
        ))
        print("path_directory_product: " + str(path_directory_product))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        print("count of instances: " + str(count_instances))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Control procedure for parallel instances.
    # Define paths to directories.
    path_directory_parallel = os.path.join(
        path_directory_product, "temporary_parallel_branch_products_regress",
    )
    # Remove any previous directories.
    putly.remove_directory(path=path_directory_parallel)
    # Create directories.
    putly.create_directories(
        path=path_directory_parallel,
    )
    control_parallel_instances(
        instances=pail_expansion["instances"],
        path_file_source_table_parameters=path_file_source_table_parameters,
        path_directory_product=path_directory_parallel,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    ##########
    # Collect information from all regressions in set.
    # Prepare table of information as a summary of all regressions in set.
    # Read source information from file.
    pails_parallel = read_source_parallel_branch_products(
        path_directory_parent=path_directory_parallel,
        name_file_child_prefix="branch_",
        name_file_child_suffix=".pickle",
        name_file_child_not="nothing_to_see_here_blargh_actrumphication_317_",
        report=report,
    )

    ##########
    # Organize information within a table.
    table_regressions = preg.organize_summary_table_regressions(
        pails_regression=pails_parallel,
        report=report,
    )

    ##########
    # Write information to file.
    # Collect information.
    # Collections of files.
    pail_write_tables = dict()
    pail_write_tables[str("table_regressions")] = table_regressions
    # Organize information for write.
    # Extract information from path to directory and file.
    path_directory_parent = os.path.dirname(
        path_file_product_table_regressions
    )
    name_suffix_file = os.path.basename(path_file_product_table_regressions)
    name_file, suffix_file = os.path.splitext(name_suffix_file)
    # Create directories.
    putly.create_directories(
        path=path_directory_product,
    )
    putly.create_directories(
        path=path_directory_parent,
    )
    # Write product information to file.
    putly.write_table_to_file(
        table=table_regressions,
        name_file=name_file,
        path_directory=path_directory_parent,
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=suffix_file,
    )
    if False:
        putly.write_table_to_file(
            table=table_regressions,
            name_file=name_file,
            path_directory=path_directory_parent,
            reset_index_rows=False,
            write_index_rows=False,
            write_index_columns=True,
            type="pickle",
            delimiter=None,
            suffix=".pickle",
        )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    groups_raw = sys.argv[1]
    path_file_source_table_parameters = sys.argv[2]
    path_file_product_table_regressions = sys.argv[3]
    path_directory_source = sys.argv[4]
    path_directory_product = sys.argv[5]
    path_directory_dock = sys.argv[6]
    report = sys.argv[7]

    # Call function for procedure.
    execute_procedure(
        groups_raw=groups_raw,
        path_file_source_table_parameters=path_file_source_table_parameters,
        path_file_product_table_regressions=(
            path_file_product_table_regressions
        ),
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass



#
