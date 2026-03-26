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
# Date, revision or review: 20 February 2026
# Date, revision or review: 14 May 2025
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
import math
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

import utility_special as sutly

#dir()
#importlib.reload()

###############################################################################
# Functionality



def parse_text_parameters(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
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

    # Paths to directories and files.
    pail["path_directory_dock"] = str(path_directory_dock).strip()
    pail["path_directory_dock_pail"] = str(path_directory_dock_pail).strip()
    pail["path_directory_source"] = str(path_directory_source).strip()
    pail["path_directory_product"] = str(path_directory_product).strip()
    pail["path_file_source_table_parameters"] = str(
        path_file_source_table_parameters
    ).strip()

    # Names and categories.
    # It is problematic to pass any white space in parameters from a script in
    # Bash. Designate the hash symbol "#" as a substitute for white space.
    # It is also problematic to pass an empty string in parameters from a
    # script in Bash. Designate the word "none" as a substitute for missing or
    # empty.
    # Iterate on individual parameters for names and categories.
    names = {
        "name_batch": name_batch,
    }
    for key_name in names.keys():
        # Determine whether parameter has a valid value that is not none.
        if (
            (str(names[key_name]).strip().lower() != "none")
        ):
            # Parse value.
            pail[key_name] = str(
                names[key_name]
            ).strip().replace("#", " ")
        else:
            pail[key_name] = ""
            pass
        pass

    # Lists.
    # Iterate on individual parameters for names and categories.
    lists = {
        "categories_batch": categories_batch,
    }
    for key_list in lists.keys():
        # Parse information.
        pail[key_list] = putly.parse_text_list_values(
            text=lists[key_list],
            delimiter=",",
        )
        pail[key_list] = putly.collect_unique_items(
            items=pail[key_list],
        )
        pass

    # Booleans, true or false.
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
            "drive_regressions_from_table_parameters.py"
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
    types_columns["category"] = "string"
    types_columns["name"] = "string"
    types_columns["abbreviation"] = "string"
    #types_columns["name_combination"] = "string"
    types_columns["selection_observations"] = "string"
    types_columns["type_regression"] = "string"
    types_columns["formula_text"] = "string"
    types_columns["feature_response"] = "string"
    types_columns["features_predictor_fixed"] = "string"
    types_columns["features_predictor_random"] = "string"
    types_columns["feature_groups_random"] = "string"
    types_columns["features_continuity_scale"] = "string"
    types_columns["method_scale"] = "string"
    types_columns["identifier_observations_primary"] = "string"
    types_columns["identifier_features_secondary"] = "string"
    types_columns["identifier_merge"] = "string"
    types_columns["directories_source_series_response_predictor"] = "string"
    types_columns["name_file_source_series_response_predictor"] = "string"
    types_columns["directories_source_table_primary"] = "string"
    types_columns["name_file_source_table_primary"] = "string"
    types_columns["directories_source__table_secondary"] = "string"
    types_columns["name_file_source_table_secondary"] = "string"
    types_columns["date_review"] = "string"
    types_columns["note"] = "string"
    # Return information.
    return types_columns


def read_source_table_parameters(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_parameters=None,
    report=None,
):
    """
    Read source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review or revision: TCW; 31 December 2025
    Review or revision: TCW; 31 March 2025

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
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Bundle information.
    pail = dict()

    # Read information from file.

    # Table of parameters for parallel instances.
    types_columns = define_column_types_table_parameters()
    pail["table"] = pandas.read_csv(
        path_file_source_table_parameters,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_correlations_from_table_parameters.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("parameter table:")
        print(pail["table"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source_table_features_plural(
    path_directory_dock=None,
    string_directories=None,
    string_file=None,
    report=None,
):
    """
    Read from file information about series of features if a valid file exists.

    Date, revision or review: 25 February 2026

    arguments:
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        string_directories (str): list of directories in text format with comma
            (',') delimiters
        string_file (str): name of file
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Determine whether there exists a file with valid information about a
    # series of features.
    # Determine whether there are plural response features.
    if (
        (string_directories != "none") and
        (string_file != "none")
    ):
        # Read from file the series of features.
        # Directories.
        directories = (
            putly.parse_text_list_values(
                text=str(string_directories).strip(),
                delimiter=",",
        ))
        # Files.
        name_file = str(string_file).strip()
        # Paths to directories and files.
        pail_path = putly.extract_organize_path_directory_file(
            name_file=name_file,
            directories_path=directories,
            name_parent="dock",
            path_directory_parent=path_directory_dock,
            report=False,
        )
        if (pail_path["existence_file"]):
            path_file_source_features = pail_path["path_file"]
        else:
            path_file_source_features = None
        pass
        # Determine whether a valid file exists.
        if (path_file_source_features is not None):
            # Read information from file.
            table_features_plural = pandas.read_csv(
                path_file_source_features,
                sep="\t",
                header=0,
                dtype={
                    "response": "string",
                    "predictor": "string",
                },
                na_values=[
                    "nan", "na", "NAN", "NA",
                    "<nan>", "<na>", "<NAN>", "<NA>",
                ],
                encoding="utf-8",
            )
        else:
            # Instance does not include enough information to expland
            # plural series of features.
            table_features_plural = None
            pass
        pass
    else:
        table_features_plural = None
        pass

    # Report.
    if report:
        # Organize information.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_regressions_from_table_parameters.py")
        print("function: read_source_table_features_plural()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_features_plural


def control_string_replace_expand_plural_features(
    path_directory_dock=None,
    table=None,
    name_series_response=None,
    name_series_predictor=None,
    report=None,
):
    """
    Consider whether instances of parameters include any definitions of plural
    response features. If so, then create new instances of parameters to
    represent each of the plural response features.

    Date, revision or review: 25 February 2026

    arguments:
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        table (object): Pandas data-frame table of parameters
        name_series_response (str): name to designate a series of features for
            use in plural expansion of the response variable
        name_series_predictor (str): name to designate a series of features for
            use in plural expansion of the predictor variable
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information.
    table = table.copy(deep=True)
    # Parse and extract information for distinct, individual instances of
    # parameters.
    # Collect information.
    instances_original = list()
    instances_novel = list()
    counter = 1
    for index, row_parameter in table.iterrows():
        # Copy information.
        row_parameter = copy.deepcopy(row_parameter)
        # Determine whether there exists a file with valid information about a
        # series of features.
        alias_directories = str("directories_source_series_response_predictor")
        alias_file = str("name_file_source_series_response_predictor")
        table_features_plural = read_source_table_features_plural(
            path_directory_dock=path_directory_dock,
            string_directories=row_parameter[alias_directories],
            string_file=row_parameter[alias_file],
            report=False,
        )
        # Determine whether there are plural response features.
        if (table_features_plural is not None):
            # Collect information.
            instances_plurality = list()
            # Iterate on series of features.
            for index, row_feature in table_features_plural.iterrows():
                # Copy information.
                row_feature = copy.deepcopy(row_feature)
                instance_plurality = copy.deepcopy(row_parameter)
                # Replace key phrases with variable names of features.
                # Iterate on individual parameters for names and categories.
                strings_features = [
                    "feature_response",
                    "feature_groups_random",
                    "features_predictor_fixed",
                    "features_predictor_random",
                    "features_continuity_scale",
                ]
                for string_features in strings_features:
                    if (
                        ("response" in row_feature.keys()) and
                        (pandas.notna(row_feature["response"])) and
                        (row_feature["response"] is not None) and
                        (len(str(row_feature["response"])) > 0) and
                        (str(row_feature["response"]) != "") and
                        (str(
                            row_feature["response"]
                        ).strip().lower() != "none")
                    ):
                        instance_plurality[string_features] = str(
                            instance_plurality[string_features]
                        ).replace(
                            name_series_response,
                            row_feature["response"],
                        )
                    if (
                        ("predictor" in row_feature.keys()) and
                        (pandas.notna(row_feature["predictor"])) and
                        (row_feature["predictor"] is not None) and
                        (len(str(row_feature["predictor"])) > 0) and
                        (str(row_feature["predictor"]) != "") and
                        (str(
                            row_feature["predictor"]
                        ).strip().lower() != "none")
                    ):
                        instance_plurality[string_features] = str(
                            instance_plurality[string_features]
                        ).replace(
                            name_series_predictor,
                            row_feature["predictor"],
                        )
                        pass
                    pass
                # Specify sequence.
                instance_plurality["sequence"] = counter
                # Specify name.
                instance_plurality["name_combination"] = "_".join([
                    str(instance_plurality["category"]).strip(),
                    str(instance_plurality["sequence"]).strip(),
                    str(instance_plurality["name"]).strip(),
                ])
                # Collect instance of parameters.
                instances_plurality.append(instance_plurality)
                # Increment counter for sequence.
                counter += 1
                pass
            # Collect information.
            instances_novel.extend(instances_plurality)
        else:
            # Copy information.
            instance_original = copy.deepcopy(row_parameter)
            # Specify sequence.
            instance_original["sequence"] = counter
            # Specify name.
            instance_original["name_combination"] = "_".join([
                str(instance_original["category"]).strip(),
                str(instance_original["sequence"]).strip(),
                str(instance_original["name"]).strip(),
            ])
            # Collect instance of parameters.
            instances_original.append(instance_original)
            # Increment counter for sequence.
            counter += 1
            pass
        pass
    # Copy information.
    instances_product = copy.deepcopy(instances_original)
    # Combine original instances with novel instances.
    instances_product.extend(instances_novel)
    # Organize table.
    table_instances = pandas.DataFrame(data=instances_product)
    # Sort rows within table.
    table_instances.sort_values(
        by=["sequence",],
        axis="index",
        ascending=True,
        inplace=True,
    )
    table_instances.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    # Report.
    if report:
        # Organize information.
        count_instances_original = len(instances_original)
        count_instances_novel = len(instances_novel)
        count_instances_product = len(instances_product)
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_regressions_from_table_parameters.py")
        print("function: control_string_replace_expand_plural_features()")
        putly.print_terminal_partition(level=5)
        print("count original: " + str(count_instances_original))
        print("count novel: " + str(count_instances_novel))
        print("count product: " + str(count_instances_product))
        putly.print_terminal_partition(level=5)
        #if (len(instances_novel) > 0):
        #    print("example of a novel instance: ")
        #    print(instances_novel[0])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_instances


def extract_parameter_instances(
    path_directory_dock=None,
    table=None,
    name_batch=None,
    categories_batch=None,
    name_series_response=None,
    name_series_predictor=None,
    report=None,
):
    """
    Extract and organize instances of parameters from table.


    Date, revision or review: 20 February 2026
    Date, revision or review: 31 December 2025
    Date, revision or review: 31 March 2025

    arguments:
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        table (object): Pandas data-frame table of parameters
        name_batch (str): name for a set or group of categories that designate
            instances of parameters in a batch for execution
        categories_batch (list<str>): names of categories that designate sets
            or groups of instances of parameters in a batch for execution
        name_series_response (str): name to designate a series of features for
            use in plural expansion of the response variable
        name_series_predictor (str): name to designate a series of features for
            use in plural expansion of the predictor variable
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Copy information.
    table = table.copy(deep=True)

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

        # Names and categories.
        # Designate the hash symbol "#" as a substitute for white space.
        # Designate the word "none" as a substitute for missing or empty.
        # Iterate on individual parameters for names and categories.
        names = [
            "type_regression",
            "formula_text",
            "feature_response",
            "feature_groups_random",
            "method_scale",
            "identifier_observations_primary",
            "identifier_features_secondary",
            "identifier_merge",
            "date_review",
            "note",
        ]
        for name in names:
            # Determine whether parameter has a valid value that is not none.
            if (
                (row[name] is not None) and
                (len(str(row[name]).strip()) > 0) and
                (str(row[name]).strip().lower() != "none")
            ):
                # Parse value.
                record[name] = str(row[name]).strip().replace("#", " ")
            else:
                # record[name] = ""
                record[name] = None
                pass
            pass

        # Lists.
        # Iterate on individual parameters for names and categories.
        lists = [
            "features_predictor_fixed",
            "features_predictor_random",
            "features_continuity_scale",
        ]
        for list_key in lists:
            # Parse information.
            list_extraction = putly.parse_text_list_values(
                text=row[list_key],
                delimiter=",",
            )
            record[list_key] = putly.collect_unique_items(
                items=list_extraction,
            )
            pass

        # Dictionaries.
        record["selection_observations"] = (
            putly.parse_extract_text_keys_values_semicolon_colon_comma(
                text=row["selection_observations"],
            )
        )["features_values"]

        # Paths to directories and files.
        # Iterate on aliases.
        aliases_temporary = [
            str("series_response_predictor"),
            str("table_primary"),
            str("table_secondary"),
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
            if (pail_path["existence_file"]):
                record[str("path_file_source_" + alias_temporary)] = (
                    pail_path["path_file"]
                )
            else:
                record[str("path_file_source_" + alias_temporary)] = None
            pass

        ##########
        # Collect unique names of columns relevant to instance of
        # parameters from current row in table.
        # Features for variables in the regression models.
        features_regression = list()
        features_regression.append(record["feature_response"])
        features_regression.extend(record["features_predictor_fixed"])
        features_regression.extend(record["features_predictor_random"])
        #if (record["identifier_observations_primary"] is not None):
        #    features_regression.insert(
        #        0,
        #        record["identifier_observations_primary"],
        #    )
        #    pass
        # Exclude key flag words for series or lists of features for response
        # or predictor variables.
        features_regression = list(filter(
            lambda item: (
                item not in [name_series_response, name_series_predictor,]
            ),
            features_regression
        ))
        features_regression = putly.collect_unique_items(
            items=features_regression,
        )
        record["features_regression"] = copy.deepcopy(features_regression)
        # Features for variables otherwise relevant to the analysis.
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
        if (record["feature_groups_random"] is not None):
            features_relevant.insert(
                0,
                record["feature_groups_random"],
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
        # Exclude key flag words for series or lists of features for response
        # or predictor variables.
        features_relevant = list(filter(
            lambda item: (
                item not in [name_series_response, name_series_predictor,]
            ),
            features_relevant
        ))
        features_relevant = putly.collect_unique_items(
            items=features_relevant,
        )
        record["features_relevant"] = copy.deepcopy(features_relevant)
        # Collect information and parameters for current row in table.
        records.append(record)
        pass

    # Bundle information.
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
        print("module: drive_regressions_from_table_parameters.py")
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


def expand_assemble_plural_features(
    instance=None,
    name_series_response=None,
    name_series_predictor=None,
    counter=None,
    report=None,
):
    """
    Create new instances of parameters to represent each of the plural features
    for response or predictor variables.

    Date, revision or review: 20 February 2026

    arguments:
        instance (list<dict>): parameters to control individual instances in
            parallel
        name_series_response (str): name to designate a series of features for
            use in plural expansion of the response variable
        name_series_predictor (str): name to designate a series of features for
            use in plural expansion of the predictor variable
        counter (int): current value of count index for instance
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Copy information.
    instance_source = copy.deepcopy(instance)

    # Expand and assemble instances with plural series of features for response
    # or predictor variables.
    # Determine whether parameters point path to a file exists.
    alias = str("path_file_source_series_response_predictor")
    if (instance[alias] is not None):
        # Read information from file.
        #features_response = putly.read_file_text_list(
        #    path_file=instance["path_file_source_series_response_predictor"],
        #    delimiter="\n",
        #    unique=True,
        #)
        table_features_plural = pandas.read_csv(
            instance[alias],
            sep="\t",
            header=0,
            dtype={
                "response": "string",
                "predictor": "string",
            },
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        # Instance does not include enough information to expland plural series
        # of features.
        table_features_plural = None
        pass
    # Determine whether to iterate on plural features for response or predictor
    # variables.
    if (
        (table_features_plural is not None) and
        (len(table_features_plural) > 0)
    ):
        # Collect information
        instances_novel = list()
        # Iterate on features.
        for index, row in table_features_plural.iterrows():
            # Collect information.
            # Copy information.
            record = copy.deepcopy(instance_source)
            # Specify sequence.
            record["sequence"] = counter
            # Specify name.
            record["name_combination"] = "_".join([
                str(record["category"]).strip(),
                str(record["sequence"]).strip(),
                str(record["name"]).strip(),
            ])
            # Specify response feature.
            if (record["feature_response"] == name_series_response):
                # Collect information.
                record["feature_response"] = row["response"]
                # Include response feature in selections of features that
                # are relevant and part of the regression model.
                record["features_regression"].insert(
                    0, row["response"],
                )
                record["features_relevant"].insert(
                    0, row["response"],
                )
                record["features_regression"] = putly.collect_unique_items(
                    items=record["features_regression"],
                )
                record["features_relevant"] = putly.collect_unique_items(
                    items=record["features_relevant"],
                )
                pass
            if (name_series_response in record["features_continuity_scale"]):
                # Define alias for current list.
                alias = "features_continuity_scale"
                # Determine position at which the flag appears in the list.
                position = record[alias].index(
                    name_series_response
                )
                # Insert the plural feature in place of the flag.
                record[alias].insert(
                    position,
                    row["response"]
                )
                # Exclude key flag words for series or lists of features for
                # response or predictor variables.
                record[alias] = list(filter(
                    lambda item: (
                        item not in [
                            name_series_response,
                        ]
                    ),
                    record[alias]
                ))
                record[alias] = putly.collect_unique_items(
                    items=record[alias],
                )
                # Include response feature in selections of features that
                # are relevant and part of the regression model.
                record["features_relevant"].insert(
                    0, row["response"],
                )
                record["features_relevant"] = putly.collect_unique_items(
                    items=record["features_relevant"],
                )
                pass
            # Specify predictor feature.
            if (name_series_predictor in record["features_predictor_fixed"]):
                # Define alias for current list.
                alias = "features_predictor_fixed"
                # Determine position at which the flag appears in the list.
                position = record[alias].index(
                    name_series_predictor
                )
                # Insert the plural feature in place of the flag.
                record[alias].insert(
                    position,
                    row["predictor"]
                )
                # Exclude key flag words for series or lists of features for
                # response or predictor variables.
                record[alias] = list(filter(
                    lambda item: (
                        item not in [
                            name_series_response, name_series_predictor,
                        ]
                    ),
                    record[alias]
                ))
                record[alias] = putly.collect_unique_items(
                    items=record[alias],
                )
                # Include response feature in selections of features that
                # are relevant and part of the regression model.
                record["features_regression"].insert(
                    0, row["predictor"],
                )
                record["features_relevant"].insert(
                    0, row["predictor"],
                )
                record["features_regression"] = putly.collect_unique_items(
                    items=record["features_regression"],
                )
                record["features_relevant"] = putly.collect_unique_items(
                    items=record["features_relevant"],
                )
                pass
            if (name_series_predictor in record["features_predictor_random"]):
                # Define alias for current list.
                alias = "features_predictor_random"
                # Determine position at which the flag appears in the list.
                position = record[alias].index(
                    name_series_predictor
                )
                # Insert the plural feature in place of the flag.
                record[alias].insert(
                    position,
                    row["predictor"]
                )
                # Exclude key flag words for series or lists of features for
                # response or predictor variables.
                record[alias] = list(filter(
                    lambda item: (
                        item not in [
                            name_series_response, name_series_predictor,
                        ]
                    ),
                    record[alias]
                ))
                record[alias] = putly.collect_unique_items(
                    items=record[alias],
                )
                # Include response feature in selections of features that
                # are relevant and part of the regression model.
                record["features_regression"].insert(
                    0, row["predictor"],
                )
                record["features_relevant"].insert(
                    0, row["predictor"],
                )
                record["features_regression"] = putly.collect_unique_items(
                    items=record["features_regression"],
                )
                record["features_relevant"] = putly.collect_unique_items(
                    items=record["features_relevant"],
                )
                pass
            if (name_series_predictor in record["features_continuity_scale"]):
                # Define alias for current list.
                alias = "features_continuity_scale"
                # Determine position at which the flag appears in the list.
                position = record[alias].index(
                    name_series_predictor
                )
                # Insert the plural feature in place of the flag.
                record[alias].insert(
                    position,
                    row["predictor"]
                )
                # Exclude key flag words for series or lists of features for
                # response or predictor variables.
                record[alias] = list(filter(
                    lambda item: (
                        item not in [
                            name_series_response,
                        ]
                    ),
                    record[alias]
                ))
                record[alias] = putly.collect_unique_items(
                    items=record[alias],
                )
                # Include response feature in selections of features that
                # are relevant and part of the regression model.
                record["features_relevant"].insert(
                    0, row["predictor"],
                )
                record["features_relevant"] = putly.collect_unique_items(
                    items=record["features_relevant"],
                )
                pass

            # Exclude key flag words for series or lists of features for response
            # or predictor variables.
            record["features_regression"] = list(filter(
                lambda item: (
                    item not in [name_series_response, name_series_predictor,]
                ),
                record["features_regression"]
            ))
            record["features_relevant"] = list(filter(
                lambda item: (
                    item not in [name_series_response, name_series_predictor,]
                ),
                record["features_relevant"]
            ))
            # Collect novel instance of parameters.
            instances_novel.append(copy.deepcopy(record))
            # Increment counter for sequence.
            counter += 1
            pass
        pass
    else:
        instances_novel = None
        pass

    # Bundle information.
    pail = dict()
    pail["instances"] = instances_novel
    pail["counter"] = counter

    # Report.
    if report:
        # Organize information.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_regressions_from_table_parameters.py")
        print("function: expand_assemble_plural_features()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def control_expand_assemble_plural_features(
    instances=None,
    name_series_response=None,
    name_series_predictor=None,
    report=None,
):
    """
    Consider whether instances of parameters include any definitions of plural
    response features. If so, then create new instances of parameters to
    represent each of the plural response features.

    Date, revision or review: 20 February 2026
    Date, revision or review: 02 January 2026
    Date, revision or review: 01 July 2025

    arguments:
        instances (list<dict>): parameters to control individual instances in
            parallel
        name_series_response (str): name to designate a series of features for
            use in plural expansion of the response variable
        name_series_predictor (str): name to designate a series of features for
            use in plural expansion of the predictor variable
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
    counter = 1
    for instance in instances_source:
        # Copy information.
        instance = copy.deepcopy(instance)
        # Organize information.
        feature_response = str(instance["feature_response"]).strip()
        features_regression = instance["features_regression"]
        # Determine whether there are plural response features.
        if (
            (feature_response == name_series_response) or
            (name_series_response in features_regression) or
            (name_series_predictor in features_regression)
        ):
            # Control the expansion and assembly of instances with plural
            # series of features for response or predictor variables.
            pail_plurality = expand_assemble_plural_features(
                instance=instance,
                name_series_response=name_series_response,
                name_series_predictor=name_series_predictor,
                counter=counter,
                report=False,
            )
            # Determine whether to include novel instances from expansion.
            if (pail_plurality["instances"] is not None):
                # Collect novel instances of parameters.
                instances_novel.extend(copy.deepcopy(
                    pail_plurality["instances"]
                ))
                # Increment counter for sequence.
                #counter += (int(len(pail_plurality["instances"])))
                counter = (int(pail_plurality["counter"]))
                pass
            pass
        else:
            # Specify sequence.
            instance["sequence"] = counter
            # Specify name.
            instance["name_combination"] = "_".join([
                str(instance["category"]).strip(),
                str(instance["sequence"]).strip(),
                str(instance["name"]).strip(),
            ])
            # Collect and preserve original instance of parameters.
            instances_original.append(copy.deepcopy(instance))
            # Increment counter for sequence.
            counter += 1
            pass
        pass
    # Copy information.
    instances_product = copy.deepcopy(instances_original)
    # Combine original instances with novel instances.
    instances_product.extend(instances_novel)

    # Organize table.
    table_instances = pandas.DataFrame(data=instances_product)
    # Sort rows within table.
    table_instances.sort_values(
        by=["sequence",],
        axis="index",
        ascending=True,
        inplace=True,
    )
    table_instances.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    # Collect information.
    pail = dict()
    pail["table"] = table_instances
    pail["instances"] = instances_product

    # Report.
    if report:
        # Organize information.
        count_instances_original = len(instances_original)
        count_instances_novel = len(instances_novel)
        count_instances_product = len(instances_product)
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_regressions_from_table_parameters.py")
        print("function: expand_plural_response_features()")
        putly.print_terminal_partition(level=5)
        print("count original: " + str(count_instances_original))
        print("count novel: " + str(count_instances_novel))
        print("count product: " + str(count_instances_product))
        putly.print_terminal_partition(level=5)
        if (len(instances_novel) > 0):
            print("example of a novel instance: ")
            print(instances_novel[0])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source_table_data(
    path_file_source_table=None,
    report=None,
):
    """
    Read and organize source information about data for regression.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 1 April 2025

    arguments:
        path_file_source_table (str): path to source file in text format as a
            table with tab delimiters between columns and newline delimiters
            between rows, with observations of features for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Determine whether parameters point path to a file exists.
    if (path_file_source_table is not None):
        # Read information from file.
        # Table of data with values for observations of features.
        table = pandas.read_csv(
            path_file_source_table,
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
        print("module: drive_regressions_from_table_parameters.py")
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
    path_file_source_table_primary=None,
    path_file_source_table_secondary=None,
    report=None,
):
    """
    Read and organize source information in tables of data.

    The primary data table orients features across columns and observations
    across rows.
    The secondary data table orients features across rows and observations
    across columns.

    Review or revision: TCW; 02 January 2026
    Review or revision: TCW; 02 July 2025

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
        path_file_source_table_primary (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with observations of features for
            regression
        path_file_source_table_secondary (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with observations of features for
            regression
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
        path_file_source_table=path_file_source_table_primary,
        report=report,
    )
    table_secondary = read_source_table_data(
        path_file_source_table=path_file_source_table_secondary,
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
        print("module: drive_regressions_from_table_parameters.py")
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


################################################################################
# Procedure


##########
# Control procedure within branch for parallelization.

# TODO: TCW; 5 March 2026
# TODO: TCW; 22 December 2025
# It doesn't make much sense to check the variance of features after transforming
# to z-score standard scale. Instead, transform to "unit_range" in a special,
# temporary table for the sake of the check.


def control_procedure_part_branch(
    execution=None,
    sequence=None,
    category=None,
    name=None,
    name_combination=None,
    selection_observations=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    feature_groups_random=None,
    features_regression=None,
    features_continuity_scale=None,
    features_relevant=None,
    method_scale=None,
    identifier_observations_primary=None,
    identifier_features_secondary=None,
    identifier_merge=None,
    path_file_source_series_response_predictor=None,
    directories_source_series_response_predictor=None,
    name_file_source_series_response_predictor=None,
    path_file_source_table_primary=None,
    directories_source_table_primary=None,
    name_file_source_table_primary=None,
    path_file_source_table_secondary=None,
    directories_source_table_secondary=None,
    name_file_source_table_secondary=None,
    date_review=None,
    note=None,
    path_file_source_table_parameters=None,
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
        category (str): categorical group of instances
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
        feature_groups_random (str): name of column in data table for
            identifiers or designations of groups of observations for which to
            allow random effects in the regression model
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
        path_file_source_series_response_predictor (str): path to source file
            in text format as a table with newline and tab delimiters, with
            names of features for response and predictor variables
        directories_source_series_response_predictor (list<str>): names of
            directories in file system path at which to find the file for the
            table of features for use as response and predictor variables
        name_file_source_series_response_predictor (str): name of file for the
            table of features for use as response and predictor variables
        path_file_source_table_primary (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with observations of features for
            regression
        directories_source_table_primary (list<str>): names of directories in
            file system path at which to find the file for the table of data
            with observations of features for regression
        name_file_source_table_primary (str): name of file for the table of
            data with observations of features for regression
        path_file_source_table_secondary (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with observations of features for
            regression
        directories_source_table_secondary (list<str>): names of directories in
            file system path at which to find the file for the table of data
            with observations of features for regression
        name_file_source_table_secondary (str): name of file for the table of
            data with observations of features for regression
        date_review (str): date of latest review of parameters
        note (str): notes about instance in table of parameters
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with parameters for multiple regressions
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
        path_file_source_table_primary=path_file_source_table_primary,
        path_file_source_table_secondary=path_file_source_table_secondary,
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
        table = preg.prepare_table_features_observations_for_analysis(
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
        table_check_raw = preg.prepare_table_features_observations_for_analysis(
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
            adjust_scale=False,
            method_scale="none", # 'z_score' or 'unit_range'
            explicate_indices=True,
            report=True,
        )
        pail_check_1 = preg.check_parameters_table_data_regression(
            table=table_check_raw,
            name_combination=name_combination,
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed_no_intercept,
            features_predictor_random=features_predictor_random_no_intercept,
            feature_groups_random=feature_groups_random,
            threshold_features_variance=0.01,
            threshold_observations_count=5,
            measure_variance="coefficient_variation",
            report=report,
        )
        pail_check_2 = preg.check_parameters_table_data_regression(
            table=table,
            name_combination=name_combination,
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed_no_intercept,
            features_predictor_random=features_predictor_random_no_intercept,
            feature_groups_random=feature_groups_random,
            threshold_features_variance=0.01,
            threshold_observations_count=5,
            measure_variance="standard_deviation", # values might be on z-score standard scale
            report=report,
        )
    else:
        pail_check_1 = dict()
        pail_check_2 = dict()
        pail_check_1["check_overall"] = True
        pail_check_2["check_overall"] = True
        pass
    # Determine whether both tables pass checks.
    if (pail_check_1["check_overall"] and pail_check_2["check_overall"]):
        check_overall = True
    else:
        check_overall = False
        pass
    #check_overall = True

    ##########
    # Perform regression analysis.
    record_extra = dict()
    pail_regression = preg.collect_organize_record_regression_analysis(
        table=table,
        index_columns="features",
        index_rows="observations",
        sequence=sequence,
        category=category,
        name=name,
        name_combination=name_combination,
        check_overall=check_overall,
        type_regression=type_regression,
        formula_text=formula_text,
        feature_response=feature_response,
        features_predictor_fixed=features_predictor_fixed,
        features_predictor_random=features_predictor_random,
        #features_predictor_random=list(),
        feature_groups_random=feature_groups_random,
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
    description_feature_groups_random = str(
        feature_groups_random
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
        description_feature_groups_random=description_feature_groups_random,
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
    pail_write_text[str("instance_" + name_combination)] = summary_text
    pail_write_objects = dict()
    pail_write_objects[str("instance_" + name_combination)] = pail_regression
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
        print("module: drive_regressions_from_table_parameters.py")
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

    Review or revision: TCW; 02 January 2026
    Review or revision: TCW; 31 March 2025

    arguments:
        instance (dict): parameters specific to current instance
            execution (int): logical binary indicator of whether to execute
                and handle the parameters for the current instance
            sequence (int): sequential index for instance's name and sort
                order
            category (str): categorical group of instances
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
            feature_groups_random (str): name of column in data table for
                identifiers or designations of groups of observations for which
                to allow random effects in the regression model
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
            path_file_source_series_response_predictor (str): path to source
                file in text format as a table with newline and tab delimiters,
                with names of features for response and predictor variables
            directories_source_series_response_predictor (list<str>): names of
                directories in file system path at which to find the file for
                the table of features for use as response and predictor
                variables
            name_file_source_series_response_predictor (str): name of file for
                the table of features for use as response and predictor
                variables
            path_file_source_table_primary (str): path to source file in
                text format as a table with tab delimiters between columns and
                newline delimiters between rows, with observations of features
                for regression
            directories_source_table_primary (list<str>): names of directories
                in file system path at which to find the file for the table of
                data with observations of features for regression
            name_file_source_table_primary (str): name of file for the table of
                data with observations of features for regression
            path_file_source_table_secondary (str): path to source file in
                text format as a table with tab delimiters between columns and
                newline delimiters between rows, with observations of features
                for regression
            directories_source_table_secondary (list<str>): names of directories
                in file system path at which to find the file for the table of
                data with observations of features for regression
            name_file_source_table_secondary (str): name of file for the table
                of data with observations of features for regression
            date_review (str): date of latest review of parameters
            note (str): notes about instance in table of parameters

        parameters (dict): parameters common to all instances
            path_file_source_table_parameters (str): path to source file in
                text format as a table with tab delimiters between columns and
                newline delimiters between rows, with parameters for multiple
                regressions
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
    category = instance_record["category"]
    name = instance_record["name"]
    name_combination = instance_record["name_combination"]
    selection_observations = instance_record["selection_observations"]
    type_regression = instance_record["type_regression"]
    formula_text = instance_record["formula_text"]
    feature_response = instance_record["feature_response"]
    features_predictor_fixed = instance_record["features_predictor_fixed"]
    features_predictor_random = instance_record["features_predictor_random"]
    feature_groups_random = instance_record["feature_groups_random"]
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
    path_file_source_series_response_predictor = (
        instance_record["path_file_source_series_response_predictor"]
    )
    directories_source_series_response_predictor = (
        instance_record["directories_source_series_response_predictor"]
    )
    name_file_source_series_response_predictor = (
        instance_record["name_file_source_series_response_predictor"]
    )
    path_file_source_table_primary = (
        instance_record["path_file_source_table_primary"]
    )
    directories_source_table_primary = (
        instance_record["directories_source_table_primary"]
    )
    name_file_source_table_primary = (
        instance_record["name_file_source_table_primary"]
    )
    path_file_source_table_secondary = (
        instance_record["path_file_source_table_secondary"]
    )
    directories_source_table_secondary = (
        instance_record["directories_source_table_secondary"]
    )
    name_file_source_table_secondary = (
        instance_record["name_file_source_table_secondary"]
    )
    date_review = instance_record["date_review"]
    note = instance_record["note"]

    # Extract parameters common across all instances.
    path_file_source_table_parameters = (
        parameters["path_file_source_table_parameters"]
    )
    path_directory_product = parameters["path_directory_product"]
    path_directory_dock = parameters["path_directory_dock"]
    report = parameters["report"]

    ##########
    # Control procedure with split branch for parallelization.
    if (int(execution) == 1):
        control_procedure_part_branch(
            execution=execution,
            sequence=sequence,
            category=category,
            name=name,
            name_combination=name_combination,
            selection_observations=selection_observations,
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            features_predictor_random=features_predictor_random,
            feature_groups_random=feature_groups_random,
            features_regression=features_regression,
            features_continuity_scale=features_continuity_scale,
            features_relevant=features_relevant,
            method_scale=method_scale,
            identifier_observations_primary=identifier_observations_primary,
            identifier_features_secondary=identifier_features_secondary,
            identifier_merge=identifier_merge,
            path_file_source_series_response_predictor=(
                path_file_source_series_response_predictor
            ),
            directories_source_series_response_predictor=(
                directories_source_series_response_predictor
            ),
            name_file_source_series_response_predictor=(
                name_file_source_series_response_predictor
            ),
            path_file_source_table_primary=path_file_source_table_primary,
            directories_source_table_primary=directories_source_table_primary,
            name_file_source_table_primary=name_file_source_table_primary,
            path_file_source_table_secondary=path_file_source_table_secondary,
            directories_source_table_secondary=(
                directories_source_table_secondary
            ),
            name_file_source_table_secondary=name_file_source_table_secondary,
            date_review=date_review,
            note=note,
            path_file_source_table_parameters=(
                path_file_source_table_parameters
            ),
            path_directory_product=path_directory_product,
            path_directory_dock=path_directory_dock,
            report=report,
        )
    pass


def control_parallel_instances(
    instances=None,
    path_file_source_table_parameters=None,
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
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with parameters for multiple regressions
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
    parameters["path_file_source_table_parameters"] = (
        path_file_source_table_parameters
    )
    parameters["path_directory_product"] = path_directory_product
    parameters["path_directory_dock"] = path_directory_dock
    parameters["report"] = report

    # Execute procedure iteratively with parallelization across instances.
    if True:
        prall.drive_procedure_parallel(
            function_control=(
                control_parallel_instance
            ),
            instances=instances, # instances[5:25]
            parameters=parameters,
            cores=8,
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
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
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
    4. within each procedural branch, execute regression analysis using type of
       model, features, and other parameters from the source table of
       parameters
    5. organize report summary information from the regression analysis and
       write this information to file

    arguments:
        name_batch (str): name for current batch of parameters
        categories_batch (list<str>): names of categories that designate sets
            or groups of instances of parameters in a batch for execution
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with parameters for multiple regressions
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
    # Parse parameters.
    pail_parameters = parse_text_parameters(
        path_directory_dock=path_directory_dock,
        path_directory_dock_pail=path_directory_dock_pail,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_file_source_table_parameters=path_file_source_table_parameters,
        name_batch=name_batch,
        categories_batch=categories_batch,
        report=report,
    )

    ##########

    # Read source information from file.
    pail_source = read_source_table_parameters(
        path_directory_dock=pail_parameters["path_directory_dock"],
        path_directory_dock_pail=pail_parameters["path_directory_dock_pail"],
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_file_source_table_parameters=(
             pail_parameters["path_file_source_table_parameters"]
        ),
        report=pail_parameters["report"],
    )
    # Organize and filter table of parameters by categories in current batch.
    pail_batch = sutly.organize_filter_table_parameters(
        table=pail_source["table"],
        name_batch=pail_parameters["name_batch"],
        categories_batch=pail_parameters["categories_batch"],
        filter_parameters=True,
        report=pail_parameters["report"],
    )

    # Determine whether table of parameters includes any references to files
    # with series of features for response or predictors.
    # Filter rows in table of parameters.
    alias_directories = str("directories_source_series_response_predictor")
    alias_file = str("name_file_source_series_response_predictor")
    table_parameters_series = pail_batch["table"].loc[
        (
            (pail_batch["table"][alias_directories] != "none") &
            (pail_batch["table"][alias_file] != "none")
        ), :
    ].copy(deep=True)
    # Replace variable flags in character strings. This strategy makes it
    # possible to define more complex models, such as those with interaction
    # terms.
    if (len(table_parameters_series) > 0):
        table_parameters = control_string_replace_expand_plural_features(
            path_directory_dock=pail_parameters["path_directory_dock"],
            table=pail_batch["table"],
            name_series_response="series_response",
            name_series_predictor="series_predictor",
            report=pail_parameters["report"],
        )
    else:
        table_parameters = pail_source["table"]
        pass

    # Extract from table information about instances of parameters.
    pail_extraction = extract_parameter_instances(
        path_directory_dock=pail_parameters["path_directory_dock"],
        table=table_parameters,
        name_batch=pail_parameters["name_batch"],
        categories_batch=pail_parameters["categories_batch"],
        report=pail_parameters["report"],
    )
    #for record in pail_extraction["records"]:
    #    print(record)
    #    pass

    # Expand plural response features.
    pail_expansion = control_expand_assemble_plural_features(
        instances=pail_extraction["records"],
        name_series_response="series_response",
        name_series_predictor="series_predictor",
        report=pail_parameters["report"],
    )
    #for instance in pail_expansion["instances"]:
    #    print(instance["feature_response"])
    #    print(instance)
    #    pass

    ##########
    # Organize information.
    count_instances = len(pail_expansion["instances"])
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_regressions_from_table_parameters.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        putly.print_terminal_partition(level=5)
        print(str("name_batch: " + name_batch))
        print(str(
            "path_file_source_table_parameters: " +
            str(pail_parameters["path_file_source_table_parameters"])
        ))
        print(str(
            "path_directory_product: " +
            str(pail_parameters["path_directory_product"])
        ))
        print(str(
            "path_directory_dock: " +
            str(pail_parameters["path_directory_dock"])
        ))
        putly.print_terminal_partition(level=5)
        print("count of instances: " + str(count_instances))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Text.
    # Lists.
    # Tables.
    pail_write_tables = dict()
    pail_write_tables[str("table_parameters_expansion")] = (
        pail_expansion["table"]
    )
    ##########
    # Write product information to file.
    # Extract information from path to directory and file.
    #path_directory = os.path.dirname(path_file)
    #name_suffix_file = os.path.basename(path_file)
    #name_file, suffix_file = os.path.splitext(name_suffix_file)
    # Define paths to directories.
    path_directory_tables = os.path.join(
        path_directory_product, "tables",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_tables,
    )
    # Text.
    # Lists.
    # Tables.
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=path_directory_tables,
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )

    ##########
    # Control procedure for parallel instances.
    # Define paths to directories.
    name_directory = str(
        "regression_parallel_instances_" +
        pail_parameters["name_batch"]
    )
    path_directory_parallel = os.path.join(
        path_directory_product, name_directory,
    )
    if True:
        # Remove any previous directories.
        putly.remove_directory(path=path_directory_parallel)
        # Create directories.
        putly.create_directories(
            path=path_directory_parallel,
        )
        control_parallel_instances(
            instances=pail_expansion["instances"],
            path_file_source_table_parameters=(
                pail_parameters["path_file_source_table_parameters"]
            ),
            path_directory_product=path_directory_parallel,
            path_directory_dock=pail_parameters["path_directory_dock"],
            report=pail_parameters["report"],
        )
    pass

    ##########
    # Collect information from all regressions in set.
    # Prepare table of information as a summary of all regressions in set.
    # Read source information from file.
    pails_parallel = read_source_parallel_branch_products(
        path_directory_parent=path_directory_parallel,
        name_file_child_prefix="instance_",
        name_file_child_suffix=".pickle",
        name_file_child_not="nothing_to_see_here_blargh_actrumphication_317_",
        report=pail_parameters["report"],
    )

    ##########
    # Organize information within a table.
    table_regressions = preg.organize_summary_table_regressions(
        pails_regression=pails_parallel,
        report=pail_parameters["report"],
    )

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Text.
    # Lists.
    # Tables.
    pail_write_tables = dict()
    pail_write_tables[str("table_regressions")] = table_regressions
    ##########
    # Write product information to file.
    # Extract information from path to directory and file.
    #path_directory = os.path.dirname(path_file)
    #name_suffix_file = os.path.basename(path_file)
    #name_file, suffix_file = os.path.splitext(name_suffix_file)
    # Define paths to directories.
    path_directory_tables = os.path.join(
        path_directory_product, "tables",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_tables,
    )
    # Text.
    # Lists.
    # Tables.
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=path_directory_tables,
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_dock = sys.argv[1]
    path_directory_dock_pail = sys.argv[2]
    path_directory_source = sys.argv[3]
    path_directory_product = sys.argv[4]
    path_file_source_table_parameters = sys.argv[5]
    name_batch = sys.argv[6]
    categories_batch = sys.argv[7]
    report = sys.argv[8]

    # Call function for procedure.
    execute_procedure(
        path_directory_dock=path_directory_dock,
        path_directory_dock_pail=path_directory_dock_pail,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_file_source_table_parameters=path_file_source_table_parameters,
        name_batch=name_batch,
        categories_batch=categories_batch,
        report=report,
    )

    pass



#
