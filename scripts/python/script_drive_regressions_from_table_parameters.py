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
            data_path_directory = putly.parse_text_list_values(
                text=row["data_path_directory"],
                delimiter=",",
            )
            record["data_path_directory"] = os.path.join(
                *data_path_directory, # 'splat' operator unpacks list items
            )
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



##########
# Control procedure within branch for parallelization.


def control_procedure_part_branch(
    execution=None,
    sequence=None,
    group=None,
    instance=None,
    name_instance=None,
    selection_observations=None,
    features_continuity_scale=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    features_relevant=None,
    data_path_directory=None,
    data_file=None,
    review=None,
    note=None,
    path_file_table_parameters=None,
    path_directory_dock=None,
    report=None,
):
    """
    Execute module's main behavior for a single parallel branch.

    arguments:
        execution (int): logical binary indicator of whether to execute and
            handle the parameters for the current instance
        sequence (int): sequential index for instance's name and sort order
        group (str): categorical group of instances
        instance (str): name or designator of instance
        name_instance (str): compound name for instance of parameters
        selection_observations (dict<list<str>>): names of columns in data
            table for feature variables and their categorical values by
            which to filter rows for observations in data table
        features_continuity_scale (list<str>): names of columns in data table
            for feature variables with values on quantitative, continuous scale
            of measurement, interval or ratio, for which to standardize the
            scale by z score
        type_regression (str): name of type of regression model to use
        formula_text (str): human readable formula for regression model,
            treated as a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        features_predictor_random (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with random effects
        features_relevant (list<str>): names of columns in data table for
            feature variables that are relevant to the current instance of
            parameters
        data_path_directory (str): path to parent directory in which to find
            the file for the data table of information about features and
            observations
        data_file (str): path to file for the data table of information about
            features and observations
        review (str): notes about review of instance in table of parameters
        note (str): notes about instance in table of parameters
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
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: control_procedure_part_branch()")
        putly.print_terminal_partition(level=5)
        print("instance: " + instance)
        putly.print_terminal_partition(level=5)
        pass

    pass


##########
# Manage parallelization.


def control_parallel_instance(
    instance=None,
    parameters=None,
):
    """
    Control procedure for a single instance in parallel with others.

    Review: TCW; 31 March 2025

    arguments:
        instance (dict): parameters specific to current instance
            execution (int): logical binary indicator of whether to execute
                and handle the parameters for the current instance
            sequence (int): sequential index for instance's name and sort
                order
            group (str): categorical group of instances
            instance (str): name or designator of instance
            name_instance (str): compound name for instance of parameters
            selection_observations (dict<list<str>>): names of columns in data
                table for feature variables and their categorical values by
                which to filter rows for observations in data table
            features_continuity_scale (list<str>): names of columns in data
                table for feature variables with values on quantitative,
                continuous scale of measurement, interval or ratio, for which
                to standardize the scale by z score
            type_regression (str): name of type of regression model to use
            formula_text (str): human readable formula for regression model,
                treated as a note for clarification
            feature_response (str): name of column in data table for feature
                variable to include in regression model as response dependent
                variable
            features_predictor_fixed (list<str>): names of columns in data
                table for feature variables to include in regression model as
                predictor independent variables with fixed effects
            features_predictor_random (list<str>): names of columns in data
                table for feature variables to include in regression model as
                predictor independent variables with random effects
            features_relevant (list<str>): names of columns in data table for
                feature variables that are relevant to the current instance of
                parameters
            data_path_directory (str): path to parent directory in which to
                find the file for the data table of information about features
                and observations
            data_file (str): path to file for the data table of information
                about features and observations
            review (str): notes about review of instance in table of parameters
            note (str): notes about instance in table of parameters

        parameters (dict): parameters common to all instances
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
    # Copy information.
    instance_record = copy.deepcopy(instance)

    ##########
    # Extract parameters.

    # Extract parameters specific to each instance.
    execution = instance_record["execution"]
    sequence = instance_record["sequence"]
    group = instance_record["group"]
    instance = instance_record["instance"]
    name_instance = instance_record["name_instance"]
    selection_observations = instance_record["selection_observations"]
    features_continuity_scale = instance_record["features_continuity_scale"]
    type_regression = instance_record["type_regression"]
    formula_text = instance_record["formula_text"]
    feature_response = instance_record["feature_response"]
    features_predictor_fixed = instance_record["features_predictor_fixed"]
    features_predictor_random = instance_record["features_predictor_random"]
    features_relevant = instance_record["features_relevant"]
    data_path_directory = instance_record["data_path_directory"]
    data_file = instance_record["data_file"]
    review = instance_record["review"]
    note = instance_record["note"]

    # Extract parameters common across all instances.
    path_file_table_parameters = parameters["path_file_table_parameters"]
    path_directory_dock = parameters["path_directory_dock"]
    report = parameters["report"]

    ##########
    # Control procedure with split branch for parallelization.
    control_procedure_part_branch(
        execution=execution,
        sequence=sequence,
        group=group,
        instance=instance,
        name_instance=name_instance,
        selection_observations=selection_observations,
        features_continuity_scale=features_continuity_scale,
        type_regression=type_regression,
        formula_text=formula_text,
        feature_response=feature_response,
        features_predictor_fixed=features_predictor_fixed,
        features_predictor_random=features_predictor_random,
        features_relevant=features_relevant,
        data_path_directory=data_path_directory,
        data_file=data_file,
        review=review,
        note=note,
        path_file_table_parameters=path_file_table_parameters,
        path_directory_dock=path_directory_dock,
        report=report,
    )
    pass


def control_parallel_instances(
    instances=None,
    path_file_table_parameters=None,
    path_directory_dock=None,
    report=None,
):
    """
    Control procedure for parallel instances.

    Review: TCW; 31 March 2025

    arguments:
        instances (list<dict>): parameters to control individual instances in
            parallel
        path_file_table_parameters (str): path to source file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports


    raises:

    returns:

    """

    # Collect parameters common across all instances.
    parameters = dict()
    parameters["path_file_table_parameters"] = path_file_table_parameters
    parameters["path_directory_dock"] = path_directory_dock
    parameters["report"] = report

    # Execute procedure iteratively with parallelization across instances.
    if False:
        prall.drive_procedure_parallel(
            function_control=(
                control_parallel_instance
            ),
            instances=instances,
            parameters=parameters,
            cores=5,
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
    # Control procedure for parallel instances.
    control_parallel_instances(
        instances=pail_source["records"],
        path_file_table_parameters=path_file_table_parameters,
        path_directory_dock=path_directory_dock,
        report=report,
    )

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
