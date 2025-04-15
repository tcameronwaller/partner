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
# Date, last execution or modification: 10 April 2025
# Review: TCW; 10 April 2025
################################################################################
# Note

# The specialty of this Python script is to drive multiple regressions from
# parameters within a single table. This script calls versatile functionality
# from the "regression.py" module within the "partner" Python package.

# Useful functionality for preparing the table of data for regression.
# pandas.get_dummies(groups).values

##########
# Review: TCW; 3 April 2025
# - On 3 April 2025, TCW confirmed that extracted values from linear OLS
#   and discrete generalized Logit regression models, intercept, and
#   predictors in the summary table match those reported directly in the
#   summary from the implementations of the respective regression models on
#   demonstration data.

################################################################################
# Installation and importation

# Standard
import sys
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


##########
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
    types_columns["name"] = "string"
    #types_columns["name_combination"] = "string"
    types_columns["selection_observations"] = "string"
    types_columns["type_regression"] = "string"
    types_columns["formula_text"] = "string"
    types_columns["feature_response"] = "string"
    types_columns["features_predictor_fixed"] = "string"
    types_columns["features_predictor_random"] = "string"
    types_columns["groups_random"] = "string"
    types_columns["features_continuity_scale"] = "string"
    types_columns["identifier_observations"] = "string"
    types_columns["method_scale"] = "string"
    types_columns["data_path_directory"] = "string"
    types_columns["data_file"] = "string"
    types_columns["review"] = "string"
    types_columns["note"] = "string"
    # Return information.
    return types_columns


# TODO: TCW; 7 April 2025
# I don't think that "identifier_observations" should be in the "features_regression" list...
# but it does need to be in the "features_relevant" list

def read_source_table_parameters(
    path_file_table_parameters=None,
    path_directory_dock=None,
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
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

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
        # Collect information and parameters from current row in table.
        record = dict()
        record["execution"] = int(row["execution"])
        record["sequence"] = row["sequence"]
        record["group"] = str(row["group"]).strip()
        record["name"] = str(row["name"]).strip() # name for instance of parameters
        record["name_combination"] = "_".join([
            str(row["group"]).strip(),
            str(row["sequence"]).strip(),
            str(row["name"]).strip(),
        ])
        record["selection_observations"] = (
            putly.parse_extract_text_keys_values_semicolon_colon_comma(
                text=row["selection_observations"],
            )
        )["features_values"]
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
        if (
            (row["groups_random"] is not None) and
            (len(str(row["groups_random"]).strip()) > 0) and
            (str(row["groups_random"]).strip().lower() != "none")
        ):
            record["groups_random"] = str(row["groups_random"]).strip()
        else:
            record["groups_random"] = None
            pass
        if (
            (row["features_continuity_scale"] is not None) and
            (len(str(row["features_continuity_scale"])) > 0) and
            (str(row["features_continuity_scale"]).strip().lower() != "none")
        ):
            record["features_continuity_scale"] = putly.parse_text_list_values(
                text=row["features_continuity_scale"],
                delimiter=",",
            )
        else:
            record["features_continuity_scale"] = list()
            pass
        if (
            (row["identifier_observations"] is not None) and
            (len(str(row["identifier_observations"])) > 0) and
            (str(row["identifier_observations"]).strip().lower() != "none")
        ):
            record["identifier_observations"] = str(
                row["identifier_observations"]
            ).strip()
        else:
            record["identifier_observations"] = None
            pass
        record["method_scale"] = str(row["method_scale"]).strip()
        record["data_path_directory"] = putly.parse_text_list_values(
            text=row["data_path_directory"],
            delimiter=",",
        )
        record["data_file"] = str(row["data_file"]).strip()
        record["review"] = str(row["review"]).strip()
        record["note"] = str(row["note"]).strip()

        # Collect unique names of columns relevant to instance of
        # parameters from current row in table.
        features_regression = list()
        features_regression.append(record["feature_response"])
        features_regression.extend(record["features_predictor_fixed"])
        features_regression.extend(record["features_predictor_random"])
        #if (record["identifier_observations"] is not None):
        #    features_regression.insert(
        #        0,
        #        record["identifier_observations"],
        #    )
        #    pass
        features_regression = putly.collect_unique_elements(
            elements=features_regression,
        )
        record["features_regression"] = copy.deepcopy(features_regression)
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
        features_relevant.extend(features_regression)
        if (record["groups_random"] is not None):
            features_relevant.insert(
                0,
                record["groups_random"],
            )
            pass
        if (record["identifier_observations"] is not None):
            features_relevant.insert(
                0,
                record["identifier_observations"],
            )
            pass
        features_relevant = putly.collect_unique_elements(
            elements=features_relevant,
        )
        record["features_relevant"] = copy.deepcopy(features_relevant)
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
        print("module: script_drive_regressions_from_table_parameters.py")
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


def read_source_table_data(
    name_file_table_data=None,
    directories_path_data=None,
    path_directory_dock=None,
    report=None,
):
    """
    Read and organize source information about data for regression.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 1 April 2025

    arguments:
        name_file_table_data (str): name of file for the table of data with
            features and observations for regression
        directories_path_data (list<str>): names of directories in path at
            which to find the file for the table of data with features and
            observations for regression
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters for selection
            of sets from MSigDB

    """

    # Define paths to directories and files.
    if ("dock" in directories_path_data):
        directories_path_data = list(filter(
            lambda directory: (directory != "dock"),
            directories_path_data
        ))
        pass
    path_directory_data = os.path.join(
        path_directory_dock,
        *directories_path_data, # 'splat' operator unpacks list items
    )
    path_file_table = os.path.join(
        path_directory_data, name_file_table_data,
    )

    # Read information from file.

    # Table of data with values for observations of features.
    table = pandas.read_csv(
        path_file_table,
        sep="\t",
        header=0,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )

    # Report.
    if report:
        # Organize.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: read_source_table_data()")
        putly.print_terminal_partition(level=5)
        print("data table:")
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


def read_source_parallel_branch_products(
    path_directory_parent=None,
    name_file_child_prefix=None,
    name_file_child_suffix=None,
    name_file_child_not=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_directory_parent (str): path to parent directory in which to find
            child files
        name_file_child_prefix (str): prefix in name by which to recognize
            relevant child files within parent directory
        name_file_child_suffix (str): suffix in name by which to recognize
            relevant child files within parent directory
        name_file_child_not (str): character string in names of files to exclude
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict<str>>): collections of information from files of parallel
            branch procedures

    """

    # Extract and filter complete paths to child files within parent directory.
    paths = putly.extract_filter_child_file_names_paths(
        path_directory=path_directory_parent,
        name_file_prefix=name_file_child_prefix,
        name_file_suffix=name_file_child_suffix,
        name_file_not=name_file_child_not,
        report=report,
    )
    # Iterate on paths to files.
    # Collect information from files.
    pails = list()
    for path in paths:
        # Read information from file.
        pail = putly.read_object_from_file_pickle(
            path_file=path,
        )
        # Collect information.
        pails.append(pail)
        pass
    # Return information.
    return pails


##########
# Organize information from regressions in a table.


def organize_summary_table_regressions(
    pails_regression=None,
    report=None,
):
    """
    Organize information from multiple regressions within a table for summary.

    arguments:
        pails_regression (list<dict<str>>): collections of information from
            files of parallel branch procedures
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of information that summarizes
            multiple regression analyses

    """

    # pail (dict): collection of information about a regression analysis
    #     entries_parameter_instance (list<str>): names and sequence of entries
    #         with information about the instance of parameters for the
    #         regression
    #     entries_model (list<str>): names and sequence of entries with
    #         information about the regression model as a whole
    #     entries_intercept_predictors (list<str>): names and sequence of
    #         entries with with information about model's intercept and
    #         predictors
    #     record (dict<str>): record for collection and assembly as rows within
    #         a Pandas data-frame table, using the lists of entries to sort the
    #         columns in the table

    # When creating a Pandas data-frame table from a list of dictionary
    # records, if one or more records have key entries that do not match those
    # from other records, then the operation automatically introduces missing
    # values.

    # Copy information.
    pails_regression = copy.deepcopy(pails_regression)

    # Collect information about all entries that the records represent
    # collectively.
    entries = dict()
    entries["entries_parameter_instance"] = list() # expect to be same for all
    #                                                 records
    entries["entries_model"] = list() # composition depends on types of
    #                                    regression
    entries["entries_intercept_predictors"] = list() # composition depends on
    #                                                   count of predictors in
    #                                                   regression model
    # Example.
    # For linear ordinary least squares regression, "entries_model" includes
    # "r_square", "r_square_adjust", and "condition", all of which are not in
    # "entries_model" for logistic regression. In contrast, for logistic
    # regression, "entries_model" includes "r_square_pseudo", which is not in
    # "entries_model" for linear ordinary least squares regression.
    # Here is the preferrable sequence of columns.
    # "r_square", "r_square_adjust", "log_likelihood", "akaike", "bayes",
    # "condition",
    # A simple solution is to insert missing values for the respective types of
    # regression.
    # Collect records.
    records = list()
    for pail in pails_regression:
        for entry_type in entries.keys():
            if (pail[entry_type] is not None):
                # Filter to keep novel entries only.
                entries_novel = list(filter(
                    lambda entry: (entry not in entries[entry_type]),
                    pail[entry_type]
                ))
                entries[entry_type].extend(entries_novel)
            pass
        # Collect records.
        records.append(copy.deepcopy(pail["record"]))
        pass

    # Create table.
    table = pandas.DataFrame(data=records)
    # Organize information in table.
    entries_total = list()
    entries_total.extend(entries["entries_parameter_instance"])
    entries_total.extend(entries["entries_model"])
    entries_total.extend(entries["entries_intercept_predictors"])
    # Filter and sort columns in table.
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=entries_total,
        report=report,
    )
    # Sort rows in table.
    table.sort_values(
        by=["sequence"],
        axis="index",
        ascending=True,
        inplace=True,
    )
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    # Return information.
    return table


################################################################################
# Procedure

# TODO: TCW; 3 April 2025
# Now integrate the functionality for logistic regression.



##########
# Control procedure within branch for parallelization.


# TODO: TCW; 10 April 2025
# Switch to use the new function in "regression.py" that organizes information
# about regression analysis in a text summary

def control_procedure_part_branch(
    execution=None,
    sequence=None,
    group=None,
    name=None,
    name_combination=None,
    selection_observations=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    groups_random=None,
    features_regression=None,
    features_continuity_scale=None,
    features_relevant=None,
    identifier_observations=None,
    method_scale=None,
    data_path_directory=None,
    data_file=None,
    review=None,
    note=None,
    path_file_table_parameters=None,
    path_directory_product=None,
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
        name (str): name or designator for instance of parameters
        name_combination (str): compound name for instance of parameters
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
        groups_random (str): name of column in data table for identifiers or
            designations of groups of observations for which to allow random
            effects in the regression model
        features_regression (list<str>): names of columns in data table for
            feature variables that are relevant to the actual regression model
        features_relevant (list<str>): names of columns in data table for
            feature variables that are relevant to the current instance of
            parameters
        identifier_observations (str): name of column in data table for unique
            identifiers of observations across rows
        method_scale (str): name of method to use to adjust the scale of values
            for features across observations, either 'z_score' or 'unit_range'
        data_path_directory (list<str>): names of directories in path at which
            to find the file for the table of data with features and
            observations for regression
        data_file (str): name of file for the table of data with features and
            observations for regression
        review (str): notes about review of instance in table of parameters
        note (str): notes about instance in table of parameters
        path_file_table_parameters (str): path to source file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Read source information from file.
    table = read_source_table_data(
        name_file_table_data=data_file,
        directories_path_data=data_path_directory,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    ##########
    # Evaluate information in table.
    # This function is useful to evaluate the variance of values for features
    # across observations after adjusting their scale to the common unit range.
    # This evaluation enables a manual, supervised approach to a method of
    # feature selection that is called 'variance thresholding'.
    if False:
        table = preg.evaluate_table_data(
            table=table,
            selection_observations=selection_observations,
            features_relevant=features_relevant,
            features_regression=features_regression,
            features_continuity_scale=features_continuity_scale,
            index_columns_source="features",
            index_columns_product="features",
            index_rows_source=identifier_observations,
            index_rows_product="observations",
            adjust_scale=True,
            method_scale=method_scale, # 'z_score' or 'unit_range'
            explicate_indices=True,
            report=report,
        )
        pass


    ##########
    # Organize information in table.
    if True:
        table = porg.prepare_table_features_observations_for_analysis(
            table=table,
            selection_observations=selection_observations,
            features_relevant=features_relevant,
            features_essential=features_regression,
            features_continuity_scale=features_continuity_scale,
            index_columns_source="features",
            index_columns_product="features",
            index_rows_source=identifier_observations,
            index_rows_product="observations",
            remove_missing=True,
            remove_redundancy=True,
            adjust_scale=True,
            method_scale=method_scale, # 'z_score' or 'unit_range'
            explicate_indices=True,
            report=False,
        )
        pass

    ##########
    # Check parameters and table of data for performing regression analysis.
    # It only makes sense to compare the variance of features if there has not
    # been scale standardization by z-score.
    if True:
        pail_check = preg.check_parameters_table_data_regression(
            table=table,
            name_combination=name_combination,
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            features_predictor_random=features_predictor_random,
            groups_random=groups_random,
            threshold_features_variance=0.01,
            threshold_observations_count=5,
            measure_variance="variance",
            report=report,
        )
    else:
        pail_check = dict()
        pail_check["check_overall"] = True
        pass

    ##########
    # Perform regression analysis.
    record_extra = dict()
    pail_regression = preg.collect_organize_record_regression_analysis(
        table=table,
        index_columns="features",
        index_rows="observations",
        sequence=sequence,
        group=group,
        name=name,
        name_combination=name_combination,
        check_overall=pail_check["check_overall"],
        type_regression=type_regression,
        formula_text=formula_text,
        feature_response=feature_response,
        features_predictor_fixed=features_predictor_fixed,
        features_predictor_random=features_predictor_random,
        groups_random=groups_random,
        record_extra=record_extra,
        delimiter_list_items=",", # delimiter to flatten list items in strings
        report=report,
    )

    #print(str(pail_regression["table_summary"]))
    summary_text = str(
        str(pail_regression["record"]["name_combination"]) +
        textwrap.dedent("""\

            ----------
        """) +
        str(pail_regression["record"]["formula_text"]) +
        textwrap.dedent("""\

            ----------
        """) +
        str(pail_regression["record"]["type_regression"]) +
        textwrap.dedent("""\

            ----------
            ----------
            ----------
        """) +
        str(pail_regression["table_summary"]) +
        textwrap.dedent("""\

            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------
        """)
    )
    print(summary_text)

    ##########
    # Write information to file.
    # Write for each parallel instance of regression.
    # A subsequent procedure will read the information from file and collect it
    # within a summary for all instances of regression.

    # Collect information.
    # Collections of files.
    pail_write_text = dict()
    pail_write_text[str("branch_" + name_combination)] = summary_text
    pail_write_objects = dict()
    pail_write_objects[str("branch_" + name_combination)] = pail_regression
    # Write product information to file.
    putly.write_character_strings_to_file_text(
        pail_write=pail_write_text,
        path_directory=path_directory_product,
    )
    putly.write_objects_to_file_pickle(
        pail_write=pail_write_objects,
        path_directory=path_directory_product,
    )

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: control_procedure_part_branch()")
        putly.print_terminal_partition(level=5)
        print("name for instance of parameters: " + name)
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
            name (str): name or designator for instance of parameters
            name_combination (str): compound name for instance of parameters
            selection_observations (dict<list<str>>): names of columns in data
                table for feature variables and their categorical values by
                which to filter rows for observations in data table
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
            groups_random (str): name of column in data table for identifiers
                or designations of groups of observations for which to allow
                random effects in the regression model
            features_regression (list<str>): names of columns in data table for
                feature variables that are relevant to the actual regression
                model
            features_continuity_scale (list<str>): names of columns in data
                table for feature variables with values on quantitative,
                continuous scale of measurement, interval or ratio, for which
                to standardize the scale by z score
            features_relevant (list<str>): names of columns in data table for
                feature variables that are relevant to the current instance of
                parameters
            identifier_observations (str): name of column in data table for
                unique identifiers of observations across rows
            method_scale (str): name of method to use to adjust the scale of
                values for features across observations, either 'z_score' or
                'unit_range'
            data_path_directory (list<str>): names of directories in path at
                which to find the file for the table of data with features and
                observations for regression
            data_file (str): name of file for the data table of information
                about features and observations
            review (str): notes about review of instance in table of parameters
            note (str): notes about instance in table of parameters

        parameters (dict): parameters common to all instances
            path_file_table_parameters (str): path to source file in text format as
                a table with tab delimiters between columns and newline delimiters
                between rows, with parameters for multiple regressions
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
    instance_record = copy.deepcopy(instance)

    ##########
    # Extract parameters.

    # Extract parameters specific to each instance.
    execution = instance_record["execution"]
    sequence = instance_record["sequence"]
    group = instance_record["group"]
    name = instance_record["name"]
    name_combination = instance_record["name_combination"]
    selection_observations = instance_record["selection_observations"]
    type_regression = instance_record["type_regression"]
    formula_text = instance_record["formula_text"]
    feature_response = instance_record["feature_response"]
    features_predictor_fixed = instance_record["features_predictor_fixed"]
    features_predictor_random = instance_record["features_predictor_random"]
    groups_random = instance_record["groups_random"]
    features_regression = instance_record["features_regression"]
    features_continuity_scale = instance_record["features_continuity_scale"]
    features_relevant = instance_record["features_relevant"]
    identifier_observations = instance_record["identifier_observations"]
    method_scale = instance_record["method_scale"]
    data_path_directory = instance_record["data_path_directory"]
    data_file = instance_record["data_file"]
    review = instance_record["review"]
    note = instance_record["note"]

    # Extract parameters common across all instances.
    path_file_table_parameters = parameters["path_file_table_parameters"]
    path_directory_product = parameters["path_directory_product"]
    path_directory_dock = parameters["path_directory_dock"]
    report = parameters["report"]

    ##########
    # Control procedure with split branch for parallelization.
    if (int(execution) == 1):
        control_procedure_part_branch(
            execution=execution,
            sequence=sequence,
            group=group,
            name=name,
            name_combination=name_combination,
            selection_observations=selection_observations,
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            features_predictor_random=features_predictor_random,
            groups_random=groups_random,
            features_regression=features_regression,
            features_continuity_scale=features_continuity_scale,
            features_relevant=features_relevant,
            identifier_observations=identifier_observations,
            method_scale=method_scale,
            data_path_directory=data_path_directory,
            data_file=data_file,
            review=review,
            note=note,
            path_file_table_parameters=path_file_table_parameters,
            path_directory_product=path_directory_product,
            path_directory_dock=path_directory_dock,
            report=report,
        )
    pass


def control_parallel_instances(
    instances=None,
    path_file_table_parameters=None,
    path_directory_product=None,
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
    parameters["path_file_table_parameters"] = path_file_table_parameters
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
            cores=5,
            report=True,
        )
    else:
        # Execute procedure directly for testing.
        control_parallel_instance(
            instance=instances[1],
            parameters=parameters,
        )
    pass


##########
# Call main procedure.


def execute_procedure(
    path_file_table_parameters=None,
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
        path_file_table_parameters (str): path to source file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
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
        path_file_table_parameters=path_file_table_parameters,
        path_directory_dock=path_directory_dock,
        report=report,
    )
    #for record in pail_source["records"]:
    #    print(record)
    #    pass


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
        print("path_directory_product: " + str(path_directory_product))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass


    ##########
    # Control procedure for parallel instances.
    # Define paths to directories.
    path_directory_parallel = os.path.join(
        path_directory_product, "temporary_parallel_branch_products",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_parallel,
    )
    control_parallel_instances(
        instances=pail_source["records"],
        path_file_table_parameters=path_file_table_parameters,
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
    table_regressions = organize_summary_table_regressions(
        pails_regression=pails_parallel,
        report=report,
    )

    ##########
    # Write information to file.
    # Collect information.
    # Collections of files.
    pail_write_tables = dict()
    pail_write_tables[str("table_regressions")] = table_regressions

    ##########
    # 5. Write product information to file.
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=path_directory_product,
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=path_directory_product,
        reset_index_rows=None,
        write_index_rows=None,
        write_index_columns=None,
        type="pickle",
        delimiter=None,
        suffix=".pickle",
    )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_table_parameters = sys.argv[1]
    path_directory_product = sys.argv[2]
    path_directory_dock = sys.argv[3]
    report = sys.argv[4]

    # Call function for procedure.
    execute_procedure(
        path_file_table_parameters=path_file_table_parameters,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass



#
