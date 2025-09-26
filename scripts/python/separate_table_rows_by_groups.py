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
# Date, first execution: 19 September 2025
# Date, last execution or modification: 19 September 2025
# Review: TCW; 19 September 2025
################################################################################
# Note


##########
# Review: TCW; 19 September 2025

################################################################################
# Installation and importation

# Standard
import sys
import os
import copy
import textwrap
import math

# Relevant
import pandas
import scipy
import numpy

# Custom
import partner.utility as putly
#import partner.parallelization as prall
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
#import partner.regression as preg
import partner.plot as pplot


#dir()
#importlib.reload()

###############################################################################
# Functionality


def parse_text_parameters(
    path_file_source_table=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    column_group=None,
    groups=None,
    prefix_file_product=None,
    suffix_file_product=None,
    report=None,
):
    """
    Parse parameters from text.

    arguments:
        path_file_source_table (str): path to source file for table in text
            format
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        column_group (str): name of column in table for groups of rows
        groups (list<str>): names of groups for rows in table
        prefix_file_product (str):
        suffix_file_product (str):
        report (str): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information
    """

    # Bundle information.
    pail = dict()
    # Parse information.
    pail["path_file_source_table"] = str(path_file_source_table).strip()
    pail["path_directory_source"] = str(path_directory_source).strip()
    pail["path_directory_product"] = str(path_directory_product).strip()
    pail["path_directory_dock"] = str(path_directory_dock).strip()
    pail["column_group"] = str(column_group).strip()
    pail["groups"] = putly.parse_text_list_values(
        text=groups,
        delimiter=",",
    )
    pail["groups"] = putly.collect_unique_items(
        items=pail["groups"],
    )
    pail["prefix_file_product"] = str(prefix_file_product).strip()
    pail["suffix_file_product"] = str(suffix_file_product).strip()
    if (
        (report is not None) and
        (str(report) != "") and
        (str(report) != "none") and
        (str(report) == "true")
    ):
        pail["report"] = True
    else:
        pail["report"] = False
        pass

    # Report.
    if pail["report"]:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: operate_sets.py")
        print("function: parse_text_parameters()")
        putly.print_terminal_partition(level=5)
        print("parameters:")
        print(pail)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def define_column_types_table_source_effects(
    columns_source_text=None,
    columns_source_number=None,
):
    """
    Defines the types of variables in columns of table.

    Review: TCW; 12 September 2025
    Review: TCW; 23 July 2025

    arguments:
        columns_source_text (list<str>): names of relevant columns in table
        columns_source_number (list<str>): names of relevant columns in table

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify types of variables in columns of table.
    types_columns = dict()
    for column in columns_source_text:
        types_columns[column] = "string"
        pass
    for column in columns_source_number:
        types_columns[column] = "float32"
        pass
    # Return information.
    return types_columns


def read_source(
    path_file_source_table=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 19 September 2025

    arguments:
        path_file_source_table (str): path to source file for table in text
            format
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (str): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Define paths to parent directories.
    #path_directory_source = os.path.join()

    # Define types of variables in columns of table.
    #types_columns = define_column_types_table_source_effects(
    #    columns_source_text=[],
    #    columns_source_number=[],
    #)

    # Read information from file.
    table = pandas.read_csv(
        path_file_source_table,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Bundle information.
    pail = dict()
    pail["table"] = table

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: merge_filter_gsea_enrichments.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print(list(pail.keys()))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


################################################################################
# Procedure


##########
# Call main procedure.


def execute_procedure(
    path_file_source_table=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    column_group=None,
    groups=None,
    prefix_file_product=None,
    suffix_file_product=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review: TCW; 19 September 2025

    arguments:
        path_file_source_table (str): path to source file for table in text
            format
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        column_group (str): name of column in table for groups of rows
        groups (list<str>): names of groups for rows in table
        prefix_file_product (str):
        suffix_file_product (str):
        report (str): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parse parameters.
    pail_parameters = parse_text_parameters(
        path_file_source_table=path_file_source_table,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        column_group=column_group,
        groups=groups,
        prefix_file_product=prefix_file_product,
        suffix_file_product=suffix_file_product,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        # Organize.
        path_directory_source = pail_parameters["path_directory_source"]
        path_directory_product = pail_parameters["path_directory_product"]
        path_directory_dock = pail_parameters["path_directory_dock"]
        # Report.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: merge_filter_gsea_enrichments.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("path_directory_source: " + str(path_directory_source))
        print("path_directory_product: " + str(path_directory_product))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_file_source_table=pail_parameters["path_file_source_table"],
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        report=pail_parameters["report"],
    )
    #print(pail_source.keys())
    print(pail_source["table"])

    # Filter groups.
    groups_all = copy.deepcopy(
        pail_source["table"][column_group].unique().tolist()
    )
    groups_available = list(filter(
        lambda item: item in groups_all, pail_parameters["groups"]
    ))

    # Separate rows from table by groups.
    pail_separation = dict()
    for group in groups_available:
        # Copy information.
        table = pail_source["table"].copy(deep=True)
        # Separate rows in table.
        table_product = table.loc[
            (table[column_group] == group), :
        ].copy(deep=True)

        # Name.
        prefix_file_product = pail_parameters["prefix_file_product"]
        suffix_file_product = pail_parameters["suffix_file_product"]
        if (
            (prefix_file_product is not None) and
            (prefix_file_product != "") and
            (len(prefix_file_product) > 0) and
            (prefix_file_product != "none") and
            (suffix_file_product is not None) and
            (suffix_file_product != "") and
            (len(suffix_file_product) > 0) and
            (suffix_file_product != "none")
        ):
            name_table_product = str(
                prefix_file_product + group + suffix_file_product
            )
        elif (
            (prefix_file_product is not None) and
            (prefix_file_product != "") and
            (len(prefix_file_product) > 0) and
            (prefix_file_product != "none")
        ):
            name_table_product = str(
                prefix_file_product + group
            )
        elif (
            (suffix_file_product is not None) and
            (suffix_file_product != "") and
            (len(suffix_file_product) > 0) and
            (suffix_file_product != "none")
        ):
            name_table_product = str(
                group + suffix_file_product
            )
            pass
        # Collect information.
        pail_separation[name_table_product] = table_product
        pass

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Lists.
    # Tables.
    pail_write_tables = pail_separation

    ##########
    # Write product information to file.
    # Lists.
    # Tables.
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=pail_parameters["path_directory_product"],
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )

    pass


# Execute program process in Python.

if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_source_table = sys.argv[1]
    path_directory_source = sys.argv[2]
    path_directory_product = sys.argv[3]
    path_directory_dock = sys.argv[4]
    column_group = sys.argv[5]
    groups = sys.argv[6]
    prefix_file_product = sys.argv[7]
    suffix_file_product = sys.argv[8]
    report = sys.argv[9]

    # Call function for procedure.
    execute_procedure(
        path_file_source_table=path_file_source_table,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        column_group=column_group,
        groups=groups,
        prefix_file_product=prefix_file_product,
        suffix_file_product=suffix_file_product,
        report=report,
    )

    pass



#
