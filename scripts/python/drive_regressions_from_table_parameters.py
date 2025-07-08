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
# Date, last execution or modification: 14 May 2025
# Review: TCW; 14 May 2025
################################################################################
# Note

# The specialty of this Python script is to drive multiple regressions from
# parameters within a single table. This script calls versatile functionality
# from the "regression.py" module within the "partner" Python package.

# Useful functionality for preparing the table of data for regression.
# pandas.get_dummies(groups).values



# TODO: TCW; 30 May 2025
# It will be a common occurrence to need to run regressions with and without
# transformation (logarithm, exponent, etc) and scale standardization (to
# happen after any transformation; mean center, z score, etc).
# Implement a new column in the parameter table for names of continuous
# features for which to apply log transformation. This transformation will
# replace the values in the original column without changing the name of the
# column.

# TODO: TCW; 1 July 2025
# confirm that the paths are valid and the primary and secondary data files exist
# read the primary and secondary data tables



##########
# Note: TCW; 3 July 2025
# Previously, for linear regression with mixed effects, I had not been
# introducing a constant intercept as a predictor with fixed effect. It was
# necessary to add this intercept explicitly to the data table and to the
# regression model. After I included this constant intercept with fixed effect,
# the results of regression changed, substantially in some cases. I think that
# before this change, the regression model was not able adequately to permit
# each group to have its own, random intercept. Rather, I think that the model
# only described the variance between groups in the "Group Var" row of the
# summary table.

##########
# Note: TCW; 2 July 2025
# With 6 parallel threads on local machine halyard, driving ordinary least
# squares regression with a simple regression model across 15,000 gene features
# requires about 1 hour and 30 minutes of time to run.

##########
# Note: TCW; 1 July 2025
# The primary data table orients features across columns and observations
# across rows.
# The secondary data table orients features across rows and observations
# across columns.

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
    types_columns["method_scale"] = "string"
    types_columns["identifier_observations_primary"] = "string"
    types_columns["identifier_features_secondary"] = "string"
    types_columns["identifier_merge"] = "string"
    types_columns["path_directory_response"] = "string"
    types_columns["name_file_list_response"] = "string"
    types_columns["path_directory_table_primary"] = "string"
    types_columns["name_file_table_primary"] = "string"
    types_columns["path_directory_table_secondary"] = "string"
    types_columns["name_file_table_secondary"] = "string"
    types_columns["review"] = "string"
    types_columns["note"] = "string"
    # Return information.
    return types_columns


def extract_organize_path_directory_file(
    name_file=None,
    directories_path=None,
    name_parent=None,
    path_directory_parent=None,
    report=None,
):
    """
    Organize from raw parameters the path to a directory and file.

    To specify the parent directory as the directory in which to find the file,
    include only the name of the parent directory.

    Review: TCW; 1 July 2025

    arguments:
        name_file (str): name of file
        directories_path (list<str>): names of directories in a path
        name_parent (str): name of parent directory as origin of path
        path_directory_parent (str): path to parent directory
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Organize information.
    name_file = str(name_file).strip()
    # Determine whether parameters are valid.
    if (
        (len(name_file) > 0) and
        (name_file != "none") and
        (len(directories_path) > 0)
    ):
        validity = True
        # Define paths to directories and files.
        if (name_parent in directories_path):
            directories_path = list(filter(
                lambda directory: (directory != name_parent),
                directories_path
            ))
            pass
        path_directory = os.path.join(
            path_directory_parent,
            *directories_path, # 'splat' operator unpacks list items
        )
        path_file = os.path.join(
            path_directory, name_file,
        )
        # Determine whether the path points to a file that exists.
        existence_directory = os.path.exists(path_directory)
        existence_file = os.path.exists(path_file)
    else:
        validity = False
        path_directory = None
        path_file = None
        existence_directory = False
        existence_file = False
        pass

    # Collect information.
    pail = dict()
    pail["validity"] = validity
    pail["path_directory"] = path_directory
    pail["path_file"] = path_file
    pail["existence_directory"] = existence_directory
    pail["existence_file"] = existence_file

    # Report.
    if report:
        # Organize.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_regressions_from_table_parameters.py")
        print("function: extract_organize_path_directory_file()")
        putly.print_terminal_partition(level=5)
        print("path to file:")
        print(path_file)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source_table_parameters(
    groups_raw=None,
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
        groups_raw (str): names of groups of instances of parameters to include
            in batch for execution
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

    # Extract names of groups.
    groups = putly.parse_text_list_values(
        text=groups_raw,
        delimiter=",",
    )

    # Filter rows in table by names of groups.
    table = table.loc[(
        table["group"].isin(groups)
    ), :].copy(deep=True)


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
        record["method_scale"] = str(row["method_scale"]).strip()

        # Identifiers across rows in tables.
        identifiers = [
            "identifier_observations_primary",
            "identifier_features_secondary",
            "identifier_merge",
        ]
        for identifier in identifiers:
            if (
                (row[identifier] is not None) and
                (len(str(row[identifier])) > 0) and
                (str(row[identifier]).strip().lower() != "none")
            ):
                record[identifier] = str(row[identifier]).strip()
            else:
                record[identifier] = None
                pass
            pass
        # Directories and files.
        record["directories_path_response"] = putly.parse_text_list_values(
            text=row["path_directory_response"],
            delimiter=",",
        )
        record["name_file_list_response"] = str(
            row["name_file_list_response"]
        ).strip()
        record["directories_path_table_primary"] = putly.parse_text_list_values(
            text=row["path_directory_table_primary"],
            delimiter=",",
        )
        record["name_file_table_primary"] = str(
            row["name_file_table_primary"]
        ).strip()
        record["directories_path_table_secondary"] = putly.parse_text_list_values(
            text=row["path_directory_table_secondary"],
            delimiter=",",
        )
        record["name_file_table_secondary"] = str(
            row["name_file_table_secondary"]
        ).strip()
        record["review"] = str(row["review"]).strip()
        record["note"] = str(row["note"]).strip()

        # Collect unique names of columns relevant to instance of
        # parameters from current row in table.
        features_regression = list()
        if (record["feature_response"] != "file_list_response"):
            features_regression.append(record["feature_response"])
            pass
        features_regression.extend(record["features_predictor_fixed"])
        features_regression.extend(record["features_predictor_random"])
        #if (record["identifier_observations_primary"] is not None):
        #    features_regression.insert(
        #        0,
        #        record["identifier_observations_primary"],
        #    )
        #    pass
        features_regression = putly.collect_unique_items(
            items=features_regression,
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
        identifiers = [
            "identifier_observations_primary",
            "identifier_features_secondary",
            "identifier_merge",
        ]
        for identifier in identifiers:
            if (record[identifier] is not None):
                features_relevant.insert(0, record[identifier])
                pass
            pass
        features_relevant = putly.collect_unique_items(
            items=features_relevant,
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
        (dict): collection of source information about parameters

    """

    # Define paths to directories and files.
    pail_path = extract_organize_path_directory_file(
        name_file=name_file_table_data,
        directories_path=directories_path_data,
        name_parent="dock",
        path_directory_parent=path_directory_dock,
        report=report,
    )
    # Determine whether parameters point path to a file exists.
    if (pail_path["existence_file"]):
        # Read information from file.
        # Table of data with values for observations of features.
        table = pandas.read_csv(
            pail_path["path_file"],
            sep="\t",
            header=0,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        table = None
        pass

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


def read_organize_source_data(
    features_relevant=None,
    identifier_observations_primary=None,
    identifier_features_secondary=None,
    identifier_merge=None,
    directories_path_table_primary=None,
    name_file_table_primary=None,
    directories_path_table_secondary=None,
    name_file_table_secondary=None,
    path_directory_dock=None,
    report=None,
):
    """
    Read and organize source information in tables of data.

    The primary data table orients features across columns and observations
    across rows.
    The secondary data table orients features across rows and observations
    across columns.

    Review: TCW; 2 July 2025

    arguments:
        features_relevant (list<str>): names of columns in data table for
            feature variables that are relevant to the current instance of
            parameters
        identifier_observations_primary (str): name of column in primary data
            table for unique identifiers of observations across rows
        identifier_features_secondary (str): name of column in secondary data
            table after transposition for unique identifiers of features across
            rows
        identifier_merge (str): name of column in primary data table and in
            secondary data table after transposition by which to merge
        directories_path_table_primary (list<str>): names of directories in
            file system path at which to find the file for the table of data
            with observations of features for regression
        name_file_table_primary (str): name of file for the table of data with
            observations of features for regression
        directories_path_table_secondary (list<str>): names of directories in
            file system path at which to find the file for the table of data
            with observations of features for regression
        name_file_table_secondary (str): name of file for the table of data
            with observations of features for regression
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Copy information.
    features_relevant = copy.deepcopy(features_relevant)

    ##########
    # Read source information from file.
    table_primary = read_source_table_data(
        name_file_table_data=name_file_table_primary,
        directories_path_data=directories_path_table_primary,
        path_directory_dock=path_directory_dock,
        report=report,
    )
    table_secondary = read_source_table_data(
        name_file_table_data=name_file_table_secondary,
        directories_path_data=directories_path_table_secondary,
        path_directory_dock=path_directory_dock,
        report=report,
    )
    if (
        (table_primary is not None) and
        (table_secondary is not None)
    ):
        # Filter features according to current instance of parameters.
        if (len(features_relevant) > 0):
            table_secondary = table_secondary.loc[
                table_secondary[identifier_features_secondary].isin(
                    features_relevant
                ), :
            ].copy(deep=True)
            pass
        # Remove unnecessary columns.
        if ("index" in table_secondary.columns.to_list()):
            table_secondary.drop(
                labels=["index",],
                axis="columns",
                inplace=True
            )
        # Organize indices of table.
        table_secondary = (
            porg.explicate_table_indices_columns_rows_single_level(
                table=table_secondary,
                index_columns=identifier_observations_primary,
                index_rows=identifier_features_secondary,
                explicate_indices=True,
                report=report,
        ))
        # Transpose table.
        table_secondary = table_secondary.transpose(copy=True)
        # Organize indices of table.
        table_secondary.reset_index(
            level=None,
            inplace=True,
            drop=False, # remove index; do not move to regular columns
        )
        table_secondary.columns.rename(
            None,
            inplace=True,
        ) # single-dimensional index
        # Merge features and observations from primary and secondary tables.
        table_merge = porg.merge_columns_two_tables(
            identifier_first=identifier_merge,
            identifier_second=identifier_merge,
            table_first=table_primary,
            table_second=table_secondary,
            preserve_index=False,
            report=report,
        )
        # Remove from list of relevant features the name of identifiers for
        # features that the table now orients across columns after
        # transposition.
        if (identifier_features_secondary in features_relevant):
            features_relevant = list(filter(
                lambda feature: (feature != identifier_features_secondary),
                features_relevant
            ))
            pass
        pass
    else:
        table_merge = table_primary
        pass

    # Bundle product information.
    pail = dict()
    pail["table_primary"] = table_primary
    pail["table_secondary"] = table_secondary
    pail["table_merge"] = table_merge
    pail["features_relevant"] = features_relevant

    # Report.
    if report:
        # Organize.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: read_organize_source_data()")
        putly.print_terminal_partition(level=5)
        print("data table:")
        print(table_merge)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# TODO: TCW; 30 June 2025
# consider moving this more general function to the regression.py module for
# easier accessibility.
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
# Expand plural response features.


def expand_plural_response_features(
    instances=None,
    path_directory_dock=None,
    report=None,
):
    """
    Consider whether instances of parameters include any definitions of plural
    response features. If so, then create new instances of parameters to
    represent each of the plural response features.

    Review: TCW; 01 July 2025

    arguments:
        instances (list<dict>): parameters to control individual instances in
            parallel
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports


    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Copy information.
    instances_source = copy.deepcopy(instances)
    # Determine count of source instances.
    count_instances_source = int(len(instances_source))

    # Evaluate instances of parameters to consider whether any define plural
    # response features that require expansion.
    # If any plural response features are found, replace the source instance
    # with product instances for each response feature.
    # Collect novel product instances and combine with the original source
    # instances.
    instances_original = list()
    instances_novel = list()
    counter = int(count_instances_source + 1)
    for instance in instances_source:
        # Copy information.
        instance = copy.deepcopy(instance)
        # Organize information.
        feature_response = str(instance["feature_response"]).strip()
        # Determine whether there are plural response features.
        if (
            (feature_response == "file_list_response")
        ):
            # Define paths to directories and files.
            pail_path = extract_organize_path_directory_file(
                name_file=instance["name_file_list_response"],
                directories_path=instance["directories_path_response"],
                name_parent="dock",
                path_directory_parent=path_directory_dock,
                report=report,
            )
            # Determine whether parameters point path to a file exists.
            if (pail_path["existence_file"]):
                # Read information from file.
                features_response = putly.read_file_text_list(
                    path_file=pail_path["path_file"],
                    delimiter="\n",
                    unique=True,
                )
            else:
                features_response = None

            # Determine whether to iterate on plural response features.
            if (
                (features_response is not None) and
                (len(features_response) > 0)
            ):
                for feature in features_response:
                    # Copy information.
                    instance_plurality = copy.deepcopy(instance)
                    # Specify sequence.
                    instance_plurality["sequence"] = counter
                    # Specify name.
                    instance_plurality["name_combination"] = "_".join([
                        str(instance_plurality["group"]).strip(),
                        str(instance_plurality["sequence"]).strip(),
                        str(instance_plurality["name"]).strip(),
                    ])
                    # Specify response feature.
                    instance_plurality["feature_response"] = feature
                    # Include response feature in selections of features that
                    # are relevant and part of the regression model.
                    instance_plurality["features_regression"].insert(
                        0, feature,
                    )
                    instance_plurality["features_relevant"].insert(
                        0, feature,
                    )
                    # Collect novel instance of parameters.
                    instances_novel.append(copy.deepcopy(instance_plurality))
                    # Increment counter for sequence.
                    counter += 1
                    pass
                pass

        else:
            # Collect and preserve original instance of parameters.
            instances_original.append(copy.deepcopy(instance))
            pass

    # Copy information.
    instances_product = copy.deepcopy(instances_original)
    # Combine original instances with novel instances.
    instances_product.extend(instances_novel)

    # Collect information.
    pail = dict()
    pail["instances"] = instances_product

    # Report.
    if report:
        # Organize information.
        count_instances_novel = len(instances_novel)
        count_instances_product = len(instances_product)
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_regressions_from_table_parameters.py")
        print("function: expand_plural_response_features()")
        putly.print_terminal_partition(level=5)
        print(
            str(
                "count of novel instances from plural expansion for " +
                "response features: "
            ) +
            str(count_instances_novel)
        )
        putly.print_terminal_partition(level=5)
        print("example of a novel instance: ")
        print(instances_novel[0])
        putly.print_terminal_partition(level=5)
        print(
            str(
                "count of product instances: "
            ) +
            str(count_instances_product)
        )
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
    method_scale=None,
    identifier_observations_primary=None,
    identifier_features_secondary=None,
    identifier_merge=None,
    directories_path_response=None,
    name_file_list_response=None,
    directories_path_table_primary=None,
    name_file_table_primary=None,
    directories_path_table_secondary=None,
    name_file_table_secondary=None,
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
        method_scale (str): name of method to use to adjust the scale of values
            for features across observations, either 'z_score' or 'unit_range'
        identifier_observations_primary (str): name of column in primary data
            table for unique identifiers of observations across rows
        identifier_features_secondary (str): name of column in secondary data
            table after transposition for unique identifiers of features across
            rows
        identifier_merge (str): name of column in primary data table and in
            secondary data table after transposition by which to merge
        directories_path_response (list<str>): names of directories in file
            system path at which to find the file for the table of data with
            observations of features for regression
        name_file_list_response (str): name of file for the table of data with
            observations of features for regression
        directories_path_table_primary (list<str>): names of directories in
            file system path at which to find the file for the table of data
            with observations of features for regression
        name_file_table_primary (str): name of file for the table of data with
            observations of features for regression
        directories_path_table_secondary (list<str>): names of directories in
            file system path at which to find the file for the table of data
            with observations of features for regression
        name_file_table_secondary (str): name of file for the table of data
            with observations of features for regression
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
    # Read and organize source information from file.
    pail_source = read_organize_source_data(
        features_relevant=features_relevant,
        identifier_observations_primary=identifier_observations_primary,
        identifier_features_secondary=identifier_features_secondary,
        identifier_merge=identifier_merge,
        directories_path_table_primary=directories_path_table_primary,
        name_file_table_primary=name_file_table_primary,
        directories_path_table_secondary=directories_path_table_secondary,
        name_file_table_secondary=name_file_table_secondary,
        path_directory_dock=path_directory_dock,
        report=report,
    )
    # pail_source["features_relevant"]

    # Copy information.
    features_predictor_fixed_no_intercept = copy.deepcopy(
        features_predictor_fixed
    )
    features_predictor_random_no_intercept = copy.deepcopy(
        features_predictor_random
    )
    features_regression_no_intercept = copy.deepcopy(features_regression)
    features_relevant_no_intercept = copy.deepcopy(
        pail_source["features_relevant"]
    )
    # Organize predictor features without intercept.
    features_predictor_fixed_no_intercept = list(filter(
        lambda feature: (feature not in ["intercept",]),
        features_predictor_fixed_no_intercept
    ))
    features_predictor_random_no_intercept = list(filter(
        lambda feature: (feature not in ["intercept",]),
        features_predictor_random_no_intercept
    ))
    features_regression_no_intercept = list(filter(
        lambda feature: (feature not in ["intercept",]),
        features_regression_no_intercept
    ))
    features_relevant_no_intercept = list(filter(
        lambda feature: (feature not in ["intercept",]),
        features_relevant_no_intercept
    ))

    ##########
    # Evaluate information in table.
    # This function is useful to evaluate the variance of values for features
    # across observations after adjusting their scale to the common unit range.
    # This evaluation enables a manual, supervised approach to a method of
    # feature selection that is called 'variance thresholding'.
    if False:
        table = preg.evaluate_table_data(
            table=pail_source["table_merge"],
            selection_observations=selection_observations,
            features_relevant=features_relevant_no_intercept,
            features_regression=features_regression_no_intercept,
            features_continuity_scale=features_continuity_scale,
            index_columns_source="features",
            index_columns_product="features",
            index_rows_source=identifier_observations_primary,
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
            table=pail_source["table_merge"],
            selection_observations=selection_observations,
            features_relevant=features_relevant_no_intercept,
            features_essential=features_regression_no_intercept,
            features_continuity_scale=features_continuity_scale,
            index_columns_source="features",
            index_columns_product="features",
            index_rows_source=identifier_observations_primary,
            index_rows_product="observations",
            remove_missing=True,
            remove_redundancy=True,
            adjust_scale=True,
            method_scale=method_scale, # 'z_score' or 'unit_range'
            explicate_indices=True,
            report=True,
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
            features_predictor_fixed=features_predictor_fixed_no_intercept,
            features_predictor_random=features_predictor_random_no_intercept,
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
        method_fit=None,
        record_extra=record_extra,
        delimiter_list_items=",", # delimiter to flatten list items in strings
        report=report,
    )

    ##########
    # Prepare text summary of regression analysis.
    description_analysis = str(
        "Linear Regression with Mixed Effects: Fixed Slopes and Random " +
        " Slopes and Intercepts; Standard Model"
    )
    formula_text = str(
        formula_text
    )
    description_response = str(
        feature_response
    )
    description_groups_random = str(
        groups_random
    )
    description_predictor = textwrap.dedent("""\
        fixed effects:
           {fixed_effects}
        random effects, intercepts:
           {random_effects}
        random effects, slope coefficients:
    """).format(
        fixed_effects=features_predictor_fixed,
        random_effects=features_predictor_random,
    )
    summary_text = preg.prepare_text_summary_regression_anova(
        title="Regressions from table of parameters",
        description_analysis=description_analysis,
        formula_text=formula_text,
        description_response=description_response,
        description_groups_random=description_groups_random,
        description_predictor=description_predictor,
        summary_1=str(pail_regression["table_summary"]),
        summary_2="",
        summary_3="",
        report=report,
    )

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
            method_scale (str): name of method to use to adjust the scale of
                values for features across observations, either 'z_score' or
                'unit_range'
            identifier_observations_primary (str): name of column in primary
                data table for unique identifiers of observations across rows
            identifier_features_secondary (str): name of column in secondary
                data table after transposition for unique identifiers of
                features across rows
            identifier_merge (str): name of column in primary data table and in
                secondary data table after transposition by which to merge
            directories_path_response (list<str>): names of directories in file
                system path at which to find the file for the table of data
                with observations of features for regression
            name_file_list_response (str): name of file for the table of data
                with observations of features for regression
            directories_path_table_primary (list<str>): names of directories in
                file system path at which to find the file for the table of
                data with observations of features for regression
            name_file_table_primary (str): name of file for the table of data
                with observations of features for regression
            directories_path_table_secondary (list<str>): names of directories
                in file system path at which to find the file for the table of
                data with observations of features for regression
            name_file_table_secondary (str): name of file for the table of data
                with observations of features for regression
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
    method_scale = instance_record["method_scale"]
    identifier_observations_primary = (
        instance_record["identifier_observations_primary"]
    )
    identifier_features_secondary = (
        instance_record["identifier_features_secondary"]
    )
    identifier_merge = (
        instance_record["identifier_merge"]
    )
    directories_path_response = (
        instance_record["directories_path_response"]
    )
    name_file_list_response = instance_record["name_file_list_response"]
    directories_path_table_primary = (
        instance_record["directories_path_table_primary"]
    )
    name_file_table_primary = instance_record["name_file_table_primary"]
    directories_path_table_secondary = (
        instance_record["directories_path_table_secondary"]
    )
    name_file_table_secondary = instance_record["name_file_table_secondary"]
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
            method_scale=method_scale,
            identifier_observations_primary=identifier_observations_primary,
            identifier_features_secondary=identifier_features_secondary,
            identifier_merge=identifier_merge,
            directories_path_response=directories_path_response,
            name_file_list_response=name_file_list_response,
            directories_path_table_primary=directories_path_table_primary,
            name_file_table_primary=name_file_table_primary,
            directories_path_table_secondary=directories_path_table_secondary,
            name_file_table_secondary=name_file_table_secondary,
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
            instance=instances[0],
            parameters=parameters,
        )
    pass


##########
# Call main procedure.


def execute_procedure(
    groups_raw=None,
    path_file_table_parameters=None,
    path_file_table_results=None,
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
        path_file_table_parameters (str): path to source file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
        path_file_table_results (str): path to product file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
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
        path_file_table_parameters=path_file_table_parameters,
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
        print("path_file_table_parameters: " + str(path_file_table_parameters))
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
    path_directory_parent = os.path.dirname(path_file_table_results)
    name_suffix_file = os.path.basename(path_file_table_results)
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
    path_file_table_parameters = sys.argv[2]
    path_file_table_results = sys.argv[3]
    path_directory_source = sys.argv[4]
    path_directory_product = sys.argv[5]
    path_directory_dock = sys.argv[6]
    report = sys.argv[7]

    # Call function for procedure.
    execute_procedure(
        groups_raw=groups_raw,
        path_file_table_parameters=path_file_table_parameters,
        path_file_table_results=path_file_table_results,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass



#
