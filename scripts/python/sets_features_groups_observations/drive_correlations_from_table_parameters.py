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
# Date, initialization: 3 December 2025
# Date, review or revision: 3 December 2025
################################################################################
# Note

# The specialty of this Python script is to drive multiple analyses of pairwise
# correlations between sets of features across groups of observations.

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

import utility_special as sutly
import plot_special as splot
import correlate_sets_features_tables_plot_charts as scorr

#dir()
#importlib.reload()

###############################################################################
# Functionality


def parse_text_parameters(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_parameters=None,
    name_batch=None,
    categories_batch=None,
    report=None,
):
    """
    Parse parameters from text.

    arguments:

        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with one row for each instance of
            parameters that define correlations
        name_batch (str): name for a set or group of categories that designate
            instances of parameters in a batch for execution
        categories_batch (list<str>): names of categories that designate sets
            or groups of instances of parameters in a batch for execution
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information
    """

    # Bundle information.
    pail = dict()

    # Parse information.

    # Paths to directories.
    pail["path_directory_source"] = str(path_directory_source).strip()
    pail["path_directory_product"] = str(path_directory_product).strip()
    pail["path_directory_dock"] = str(path_directory_dock).strip()

    # Paths to files.
    pail["path_file_source_table_parameters"] = str(
        path_file_source_table_parameters
    ).strip()

    # Names and categories.
    # It is problematic to pass any white space in parameters from a script in
    # Bash. Designate the hash symbol "#" as a substitute for white space.
    # It is also problematic to pass an empty string in parameters from a
    # script in Bash. Designate the word "none" as a substitute for missing or
    # empty.
    # Iterate on individual names that could be empty or missing.
    categories = {
        "name_batch": name_batch,
    }
    for key_category in categories.keys():
        # Determine whether parameter has a valid value that is not none.
        if (
            (str(categories[key_category]).strip().lower() != "none")
        ):
            # Parse value.
            pail[key_category] = str(
                categories[key_category]
            ).strip().replace("#", " ")
        else:
            pail[key_category] = ""
            pass
        pass

    # Simple list of names and categories.
    pail["categories_batch"] = putly.parse_text_list_values(
        text=categories_batch,
        delimiter=",",
    )

    # Boolean, true or false.
    # Iterate on individual of Boolean designations.
    designations = {
        "report": report,
    }
    for key_designation in designations.keys():
        # Determine whether parameter has a valid value.
        if (
            (designations[key_designation] is not None) and
            (len(str(designations[key_designation])) > 0) and
            (str(designations[key_designation]) != "") and
            (str(designations[key_designation]).strip().lower() != "none") and
            (str(designations[key_designation]) == "true")
        ):
            # Designation is true.
            pail[key_designation] = True
        else:
            # Designation is false.
            pail[key_designation] = False
            pass
        pass

    # Report.
    if pail["report"]:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "drive_correlations_from_table_parameters.py"
        )
        print(str("module: " + module))
        function = str(
            "parse_text_parameters()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        print("parameters:")
        print(pail)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


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
    types_columns["category"] = "string"
    types_columns["name"] = "string"
    types_columns["abbreviation"] = "string"
    types_columns["directories_source_table_features_observations"] = (
        "string"
    )
    types_columns["name_file_source_table_features_observations"] = "string"
    types_columns["directories_source_list_features_first"] = "string"
    types_columns["name_file_source_list_features_first"] = "string"
    types_columns["directories_source_list_features_second"] = "string"
    types_columns["name_file_source_list_features_second"] = "string"
    types_columns["directories_source_table_groups_observations"] = "string"
    types_columns["name_file_source_table_groups_observations"] = "string"
    types_columns["column_identifier_observation"] = "string"
    types_columns["column_name_observation"] = "string"
    types_columns["proportion_nonmissing_observations"] = "string"
    types_columns["type_correlation"] = "string"
    types_columns["match_pairs"] = "string"
    types_columns["prefix_match_first"] = "string"
    types_columns["suffix_match_first"] = "string"
    types_columns["prefix_match_second"] = "string"
    types_columns["suffix_match_second"] = "string"
    types_columns["intersect_features"] = "string"
    types_columns["threshold_filter_first"] = "string"
    types_columns["threshold_filter_second"] = "string"
    types_columns["threshold_type"] = "string"
    types_columns["threshold_value"] = "string"
    types_columns["threshold_proportion"] = "string"
    types_columns["threshold_optimization"] = "string"
    types_columns["threshold_optimization_count"] = "string"
    types_columns["sort_match_pairs_diagonal"] = "string"
    types_columns["cluster_features_first"] = "string"
    types_columns["cluster_features_second"] = "string"
    types_columns["sort_other_features"] = "string"
    types_columns["plot_scale_minimum"] = "string"
    types_columns["plot_scale_center"] = "string"
    types_columns["plot_scale_maximum"] = "string"
    types_columns["report"] = "string"
    types_columns["date_review"] = "string"
    types_columns["note"] = "string"

    # Return information.
    return types_columns


def read_source_table_parameters(
    name_batch=None,
    categories_batch=None,
    path_file_source_table_parameters=None,
    path_directory_dock=None,
    filter_instances_parameters=None,
    report=None,
):
    """
    Read and organize source information about parameters for regressions.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 3 December 2025

    arguments:
        name_batch (str): name for a set or group of categories that designate
            instances of parameters in a batch for execution
        categories_batch (list<str>): names of categories that designate sets
            or groups of instances of parameters in a batch for execution
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with one row for each instance of
            parameters that define correlations
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        filter_instances_parameters (bool): whether to filter instances of
            parameters
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Read information from file.

    # Table of parameters for parallel instances.
    types_columns = define_column_types_table_parameters()
    table = pandas.read_csv(
        path_file_source_table_parameters,
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

    # Filter rows in table by names of categories.
    if filter_instances_parameters:
        table = table.loc[(
            (table["execution"] == 1) &
            (table["category"].isin(categories_batch))
        ), :].copy(deep=True)
        pass

    # Parse and extract information for distinct, individual instances of
    # parameters.
    # Collect information.
    records = list()
    for index, row in table.iterrows():
        # Collect information and parameters from current row in table.
        record = dict()

        # Identity for the instance of parameters.
        record["execution"] = int(row["execution"])
        record["sequence"] = int(row["sequence"])
        record["category"] = str(row["category"]).strip()
        record["name"] = str(row["name"]).strip() # name for instance of parameters
        record["abbreviation"] = str(row["abbreviation"]).strip()
        record["name_combination"] = "_".join([
            str(row["category"]).strip(),
            str(row["sequence"]).strip(),
            str(row["name"]).strip(),
        ])

        # Paths to directories and files.
        # Iterate on aliases.
        aliases_temporary = [
            str("table_features_observations"),
            str("list_features_first"),
            str("list_features_second"),
            str("table_groups_observations"),
        ]
        for alias_temporary in aliases_temporary:
            # Directories.
            alias_temporary_directories = str(
                "directories_source_" + alias_temporary
            )
            record[alias_temporary_directories] = (
                putly.parse_text_list_values(
                    text=str(row[alias_temporary_directories]).strip(),
                    delimiter=",",
            ))
            # Files.
            alias_temporary_file = str(
                "name_file_source_" + alias_temporary
            )
            record[alias_temporary_file] = str(
                row[alias_temporary_file]
            ).strip()
            # Paths to directories and files.
            pail_path = putly.extract_organize_path_directory_file(
                name_file=record[alias_temporary_file],
                directories_path=record[alias_temporary_directories],
                name_parent="dock",
                path_directory_parent=path_directory_dock,
                report=report,
            )
            record[str("path_file_source_" + alias_temporary)] = (
                pail_path["path_file"]
            )
            pass

        # Names and categories.
        # It is problematic to pass any white space in parameters from a script in
        # Bash. Designate the hash symbol "#" as a substitute for white space.
        # It is also problematic to pass an empty string in parameters from a
        # script in Bash. Designate the word "none" as a substitute for missing or
        # empty.
        # Iterate on individual names that could be empty or missing.
        keys_categories = [
            "column_identifier_observation",
            "column_name_observation",
            "type_correlation",
            "prefix_match_first",
            "suffix_match_first",
            "prefix_match_second",
            "suffix_match_second",
            "threshold_type",
            "date_review",
            "note",
        ]
        for key_category in keys_categories:
            # Determine whether parameter has a valid value that is not none.
            if (
                (str(row[key_category]).strip().lower() != "none")
            ):
                # Parse value.
                record[key_category] = str(
                    row[key_category]
                ).strip().replace("#", " ")
            else:
                record[key_category] = ""
                pass
            pass

        # Numbers.
        # Iterate on individual names that could be empty or missing.
        # Float.
        keys_numbers = [
            "proportion_nonmissing_observations",
            "threshold_value",
            "threshold_proportion",
            "plot_scale_minimum",
            "plot_scale_center",
            "plot_scale_maximum",
        ]
        for key_number in keys_numbers:
            # Determine whether parameter has a valid value that is not none.
            if (
                (len(str(row[key_number]).strip()) > 0) and
                (str(row[key_number]).strip().lower() != "none")
            ):
                # Parse value.
                record[key_number] = float(
                    str(row[key_number]).strip()
                )
            else:
                record[key_number] = None
                pass
            pass
        # Integer.
        keys_numbers = [
            "threshold_optimization_count",
        ]
        for key_number in keys_numbers:
            # Determine whether parameter has a valid value that is not none.
            if (
                (len(str(row[key_number]).strip()) > 0) and
                (str(row[key_number]).strip().lower() != "none")
            ):
                # Parse value.
                record[key_number] = int(
                    str(row[key_number]).strip()
                )
            else:
                record[key_number] = None
                pass
            pass

        # Boolean, true or false.
        # Iterate on individual of Boolean designations.
        designations = [
            "match_pairs",
            "intersect_features",
            "threshold_filter_first",
            "threshold_filter_second",
            "threshold_optimization",
            "sort_match_pairs_diagonal",
            "cluster_features_first",
            "cluster_features_second",
            "sort_other_features",
            "report",
        ]
        for key_designation in designations:
            # Determine whether parameter has a valid value.
            if (
                (row[key_designation] is not None) and
                (len(str(row[key_designation])) > 0) and
                (str(row[key_designation]) != "") and
                (str(row[key_designation]).strip().lower() != "none") and
                (str(row[key_designation]) == "true")
            ):
                # Designation is true.
                record[key_designation] = True
            else:
                # Designation is false.
                record[key_designation] = False
                pass
            pass

        # Collect information and parameters for current row in table.
        records.append(record)
        pass

    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["records"] = records

    # Report.
    if report:
        # Organize information.
        count_records = len(records)
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_correlations_from_table_parameters.py")
        print("function: read_source_table_parameters()")
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


##########
# Manage parallelization.


def control_parallel_instance(
    instance=None,
    parameters=None,
):
    """
    Control procedure for a single instance in parallel with others.

    Review: TCW; 3 December 2025

    arguments:
        instance (dict): parameters specific to current instance
            ...
            review (str): notes about review of instance in table of parameters
            note (str): notes about instance in table of parameters

        parameters (dict): parameters common to all instances

            name_batch (str): name for a set or group of categories that
                designate instances of parameters in a batch for execution
            categories_batch (list<str>): names of categories that designate
                sets or groups of instances of parameters in a batch for
                execution
            path_file_source_table_parameters (str): path to source file in
                text format as a table with tab delimiters between columns and
                newline delimiters between rows, with one row for each instance
                of parameters that define correlations
            path_directory_source (str): path to directory for procedure's
                source directories and files
            path_directory_product (str): path to directory for procedure's
                product directories and files
            path_directory_dock (str): path to dock directory for procedure's
                source and product directories and files
            report (bool): whether to print reports


    raises:

    returns:

    """

    ##########
    # Copy information.
    instance = copy.deepcopy(instance)

    ##########
    # Extract and organize parameters.

    # Parameters specific to each instance.
    # instance[]
    # ...

    # Parameters common across all instances.
    # parameters[]
    # ...

    ##########
    # Control procedure with split branch for parallelization.
    # Determine whether current instance of parameters belongs to the batch for
    # execution.
    if (
        (int(instance["execution"]) == 1) and
        (instance["category"] in parameters["categories_batch"])
    ):

        # The parent directory passed to this procedure corresponds to the name
        # for a batch of instances. Each instance of parameters is independent
        # and needs to have its own product directory.
        # Define paths to directories.
        path_directory_product_child = os.path.join(
            path_directory_product, instance["name_combination"],
        )
        # Create directories.
        putly.create_directories(
            path=path_directory_product_child,
        )

        # Control procedure for current instance of parameters.
        scorr.control_procedure(
            path_directory_source=parameters["path_directory_source"],
            path_directory_product=path_directory_product_child,
            path_directory_dock=parameters["path_directory_dock"],

            path_file_source_table_features_observations=(
                instance["path_file_source_table_features_observations"]
            ),
            path_file_source_list_features_first=(
                instance["path_file_source_list_features_first"]
            ),
            path_file_source_list_features_second=(
                instance["path_file_source_list_features_second"]
            ),
            path_file_source_table_groups_observations=(
                instance["path_file_source_table_groups_observations"]
            ),
            column_identifier_observation=(
                instance["column_identifier_observation"]
            ),
            column_name_observation=instance["column_name_observation"],
            proportion_nonmissing_observations=(
                instance["proportion_nonmissing_observations"]
            ),
            type_correlation=instance["type_correlation"],
            match_pairs=instance["match_pairs"],
            prefix_match_first=instance["prefix_match_first"],
            suffix_match_first=instance["suffix_match_first"],
            prefix_match_second=instance["prefix_match_second"],
            suffix_match_second=instance["suffix_match_second"],
            intersect_features=instance["intersect_features"],
            threshold_filter_first=instance["threshold_filter_first"],
            threshold_filter_second=instance["threshold_filter_second"],
            threshold_type=instance["threshold_type"],
            threshold_value=instance["threshold_value"],
            threshold_proportion=instance["threshold_proportion"],
            threshold_optimization=instance["threshold_optimization"],
            threshold_optimization_count=(
                instance["threshold_optimization_count"]
            ),
            sort_match_pairs_diagonal=instance["sort_match_pairs_diagonal"],
            cluster_features_first=instance["cluster_features_first"],
            cluster_features_second=instance["cluster_features_second"],
            sort_other_features=instance["sort_other_features"],
            plot_scale_minimum=instance["plot_scale_minimum"],
            plot_scale_center=instance["plot_scale_center"],
            plot_scale_maximum=instance["plot_scale_maximum"],
            report=instance["report"],
        )

    pass


def control_parallel_instances(
    name_batch=None,
    categories_batch=None,
    instances=None,
    path_file_source_table_parameters=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    report=None,
):
    """
    Control procedure for parallel instances.

    Review: TCW; 3 December 2025

    arguments:

        name_batch (str): name for a set or group of categories that designate
            instances of parameters in a batch for execution
        categories_batch (list<str>): names of categories that designate sets
            or groups of instances of parameters in a batch for execution
        instances (list<dict>): parameters to control individual instances in
            parallel
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with one row for each instance of
            parameters that define correlations
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

    # Collect parameters common across all instances.
    parameters = dict()
    parameters["name_batch"] = name_batch
    parameters["categories_batch"] = categories_batch
    parameters["path_file_source_table_parameters"] = (
        path_file_source_table_parameters
    )
    parameters["path_directory_source"] = path_directory_source
    parameters["path_directory_product"] = path_directory_product
    parameters["path_directory_dock"] = path_directory_dock
    parameters["report"] = report

    # Execute procedure iteratively with parallelization across instances.
    if True:
        prall.drive_procedure_parallel(
            function_control=(
                control_parallel_instance
            ),
            instances=instances,
            parameters=parameters,
            cores=7,
            report=True,
        )
    else:
        # Execute procedure directly for testing.
        control_parallel_instance(
            instance=instances[0],
            parameters=parameters,
        )
    pass


##########
# Call main procedure.


def execute_procedure(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_parameters=None,
    name_batch=None,
    categories_batch=None,
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
    4. within each procedural branch, execute correlation analysis using
       instances of parameters from the source table of parameters
    5. organize and write to file report tables, plot charts, and other summary
       information from the correlation analysis

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with one row for each instance of
            parameters that define correlations
        name_batch (str): name for a set or group of categories that designate
            instances of parameters in a batch for execution
        categories_batch (list<str>): names of categories that designate sets
            or groups of instances of parameters in a batch for execution
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parse parameters.
    pail_parameters = parse_text_parameters(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_parameters=path_file_source_table_parameters,
        name_batch=name_batch,
        categories_batch=categories_batch,
        report=report,
    )

    ##########
    # Read source information from file.
    pail_source = read_source_table_parameters(
        name_batch=pail_parameters["name_batch"],
        categories_batch=pail_parameters["categories_batch"],
        path_file_source_table_parameters=(
            pail_parameters["path_file_source_table_parameters"]
        ),
        path_directory_dock=pail_parameters["path_directory_dock"],
        filter_instances_parameters=True,
        report=pail_parameters["report"],
    )
    #for record in pail_source["records"]:
    #    print(record)
    #    pass

    ##########
    # Organize information.
    count_instances = len(pail_source["records"])
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_correlations_from_table_parameters.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("path_file_source_table_parameters: " + str(
            pail_parameters["path_file_source_table_parameters"]
        ))
        print("path_directory_product: " + str(
            pail_parameters["path_directory_product"]
        ))
        print("path_directory_dock: " + str(
            pail_parameters["path_directory_dock"]
        ))
        putly.print_terminal_partition(level=5)
        print("count of instances: " + str(count_instances))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Control procedure for parallel instances.
    control_parallel_instances(
        name_batch=pail_parameters["name_batch"],
        categories_batch=pail_parameters["categories_batch"],
        instances=pail_source["records"],
        path_file_source_table_parameters=(
            pail_parameters["path_file_source_table_parameters"]
        ),
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        report=pail_parameters["report"],
    )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_source = sys.argv[1]
    path_directory_product = sys.argv[2]
    path_directory_dock = sys.argv[3]
    path_file_source_table_parameters = sys.argv[4]
    name_batch = sys.argv[5]
    categories_batch = sys.argv[6]
    report = sys.argv[7]

    # Call function for procedure.
    execute_procedure(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_parameters=path_file_source_table_parameters,
        name_batch=name_batch,
        categories_batch=categories_batch,
        report=report,
    )

    pass



#
