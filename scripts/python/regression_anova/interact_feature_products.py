"""
This module script is part of a series of studies of age, exercise, and dietary
intervention with omega-3 fatty acids in skeletal muscle and subcutaneous
adipose of healthy adult humans.

This module is part of the 'phenotypes' subpackage within
the 'age_exercise' package.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This file is part of the project package directory 'age_exercise'
    (https://github.com/tcameronwaller/age_exercise/).

    Project 'age_exercise' supports data analysis with team in endocrinology.
    Copyright (C) 2025 Thomas Cameron Waller

    The code within project 'age_exercise' is free software: you can
    redistribute it and/or modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation, either version 3 of
    the GNU General Public License, or (at your option) any later version.

    The code within project 'age_exercise' is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'age_exercise'. If not, see <http://www.gnu.org/licenses/>.
"""

################################################################################
# Author: T. Cameron Waller, Ph.D.
# Date, initialization: 24 February 2026
# Date, revision or review: 3 March 2026
################################################################################



###############################################################################
# Notes

# Note: TCW; 2 March 2026
# If the subsequent regression model will use quantitative features on the
# standard z-score scale, it is important to adjust the scale before the
# calculation of interaction terms and using the same groups of observations.
# The calculation of the standard z-score depends on the group of observations.
# If filtering or stratification of observations happens after calculation of
# the standard z-score, then the scale of the quantitative feature will no
# longer have mean of zero and standard deviation of one.


###############################################################################
# Installation and importation

# Standard
import sys
#print(sys.path)
import os
import math
import statistics
import pickle
import copy
import random
import itertools
import time
from datetime import datetime
#import dateutil # requires explicit installation

# Relevant
import numpy
import scipy.stats
import pandas
pandas.options.mode.chained_assignment = None # default = "warn"
pandas.set_option('future.no_silent_downcasting', True) # set option to suppress warnings

# Custom
import partner.utility as putly
import partner.extraction as pextr
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
#import partner.regression as preg
import partner.decomposition as pdecomp
import partner.plot as pplot
import partner.parallelization as prall
import utility_special as sutly

###############################################################################
# Functionality


##########
# Parse parameters.


def parse_text_parameters(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_features_observations=None,
    path_file_source_table_groups_observations=None,
    path_file_source_table_parameters=None,
    column_identifier_observation=None,
    column_name_observation=None,
    name_batch=None,
    categories_batch=None,
    scale_before_interaction=None,
    method_scale=None,
    report=None,
):
    """
    Parse parameters from text.

    arguments:
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_directory_dock_pail (str): path to pail directory for procedure's
            source and product directories and files
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_file_source_table_features_observations (object): path to source
            file in text format as a table
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with parameters
        name_batch (str): name for current batch of parameters
        categories_batch (list<str>): names of categories that designate sets
            or groups of instances of parameters in a batch for execution
        report (str): whether to print reports

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
    pail["path_file_source_table_features_observations"] = str(
        path_file_source_table_features_observations
    ).strip()
    pail["path_file_source_table_groups_observations"] = str(
        path_file_source_table_groups_observations
    ).strip()
    pail["path_file_source_table_parameters"] = str(
        path_file_source_table_parameters
    ).strip()

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

    # Names and categories.
    # It is problematic to pass any white space in parameters from a script in
    # Bash. Designate the hash symbol "#" as a substitute for white space.
    # It is also problematic to pass an empty string in parameters from a
    # script in Bash. Designate the word "none" as a substitute for missing or
    # empty.
    # Iterate on individual parameters for names and categories.
    names = {
        "name_batch": name_batch,
        "column_identifier_observation": column_identifier_observation,
        "column_name_observation": column_name_observation,
        "method_scale": method_scale,
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

    # Booleans, true or false.
    # Iterate on individual of Boolean designations.
    designations = {
        "scale_before_interaction": scale_before_interaction,
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
        print("module: interact_feature_products.py")
        print("function: parse_text_parameters()")
        putly.print_terminal_partition(level=5)
        print("parameters:")
        print(pail)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Read source information.


def read_source(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_features_observations=None,
    path_file_source_table_groups_observations=None,
    path_file_source_table_parameters=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Date, review or revision: 18 February 2026

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_file_source_table_features_observations (object): path to source
            file in text format as a table
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with parameters
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Determine whether paths point to a directory or file that exist.
    #existence_directory = os.path.exists(path_directory)

    # Bundle information.
    pail = dict()

    # Read information from file.
    pail["table_features_observations"] = pandas.read_csv(
        path_file_source_table_features_observations,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_groups_observations"] = pandas.read_csv(
        path_file_source_table_groups_observations,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_parameters"] = pandas.read_csv(
        path_file_source_table_parameters,
        sep="\t",
        header=0,
        #dtype=types_columns,
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
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("interact_feature_products.py")
        print("module: " + name_module)
        name_function = str("read_source()")
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Transform scale for features with values on a quantitative, continuous scale
# of measurement.


def transform_scale_features_by_parameters(
    table_features_observations=None,
    table_parameters=None,
    column_identifier_observation=None,
    column_name_observation=None,
    scale_before_interaction=None,
    method_scale=None,
    report=None,
):
    """
    Transform the scale for features with values on a quantitative, continuous
    scale of measurement.

    Date, revision or review: 3 March 2026

    arguments:
        table_features_observations (object): Pandas data-frame table of
            subjects, samples, and their attribute features across observations
        table_parameters (object): Pandas data-frame table
        ...
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information.
    table_main = table_features_observations.copy(deep=True)
    table_parameters = table_parameters.copy(deep=True)

    # Collect identifiers of features for which to transform scale.
    # Collect information.
    features_scale = list()
    # Iterate on rows in table of parameters.
    for index, row_parameter in table_parameters.iterrows():
        # Parse information.
        features_scale_row = putly.parse_text_list_values(
            text=row_parameter["features_continuity_scale"],
            delimiter=",",
        )
        # Collect information.
        features_scale.extend(features_scale_row)
        pass
    # Extract identifiers of features.
    features_observations = copy.deepcopy(
        table_main.columns.unique().tolist()
    )
    # Filter to features with available signals.
    features_selection = list(filter(
        lambda feature: (feature in features_observations),
        copy.deepcopy(features_scale)
    ))
    # Collect unique features.
    features_selection = putly.collect_unique_items(
        items=features_selection,
    )

    # Transform scale of values for the features before calculating difference
    # between instances.
    # Standardize scale of values for observations of features.
    table_main = pscl.manage_transform_scale_feature_by_table_columns(
        table=table_main,
        features_continuity_scale=features_selection,
        adjust_scale=scale_before_interaction,
        method_scale=method_scale,
        report=report,
    )

    # Bundle information.
    pail = dict()
    pail["table_features_observations"] = table_main
    pail["features_scale"] = features_selection

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("interact_feature_products.py")
        print("module: " + name_module)
        name_function = str(
            "transform_scale_features_by_parameters()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Calculate interaction products.


def calculate_interaction_products_between_feature_pairs(
    table_features_observations=None,
    table_parameters=None,
    report=None,
):
    """
    Calculate multiplication products of values between pairs of features as a
    representation of their interactions for regression analysis.

    Date, revision or review: 25 February 2026

    arguments:
        table_features_observations (object): Pandas data-frame table of
            subjects, samples, and their attribute features across observations
        table_parameters (object): Pandas data-frame table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Define subordinate functions for internal use.
    def calculate_interaction_product(
        series=None,
        factor_one=None,
        factor_two=None,
    ):
        # Determine whether there is adequate information for division.
        if (
            (factor_one is not None) and
            (factor_two is not None) and
            (factor_one in series.keys()) and
            (factor_two in series.keys()) and
            (pandas.notna(series[factor_one])) and
            (pandas.notna(series[factor_two]))
            #(not math.isnan(series[factor_one])) and
            #(not math.isnan(series[factor_two]))
        ):
            # Manage types of values.
            factor_one = float(series[factor_one])
            factor_two = float(series[factor_two])
            # There is adequate information for multiplication.
            product = (factor_one * factor_two)
        else:
            # There is inadequate information for division.
            product = float("nan")
            pass
        # Return information.
        return product

    # Copy information.
    table_main = table_features_observations.copy(deep=True)
    table_parameters = table_parameters.copy(deep=True)

    # Available features.
    # Extract identifiers of features.
    features_observations = copy.deepcopy(
        table_main.columns.unique().tolist()
    )
    # Collect unique features.
    features_observations = putly.collect_unique_items(
        items=features_observations,
    )
    # Filter rows in table for parameters about unique pairs of features.
    table_parameters.drop_duplicates(
        subset=["feature_first", "feature_second",],
        keep="first",
        inplace=True,
    )
    # Filter rows in table for parameters about features with available
    # observations.
    table_parameters = table_parameters.loc[
        (
            (table_parameters["feature_first"].isin(features_observations)) &
            (table_parameters["feature_second"].isin(features_observations))
        ), :
    ].copy(deep=True)
    # Extract names of new features representing interaction products.
    # Copy information.
    features_interaction = copy.deepcopy(
        table_parameters["name_interaction"].tolist()
    )
    # Collect unique features.
    features_interaction = putly.collect_unique_items(
        items=features_interaction,
    )
    # Iterate on rows in table for parameters.
    for index, row_parameter in table_parameters.iterrows():
        # Calculate features for interaction products.
        table_main[row_parameter["name_interaction"]] = table_main.apply(
            lambda row: calculate_interaction_product(
                series=row,
                factor_one=row_parameter["feature_first"],
                factor_two=row_parameter["feature_second"],
            ),
            axis="columns", # apply function to each row
        )
        pass

    # Bundle information.
    pail = dict()
    pail["table_features_observations"] = table_main
    pail["features_interaction"] = features_interaction

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("7_interact_feature_products.py")
        print("module: " + name_module)
        name_function = str(
            "calculate_interaction_products_between_feature_pairs()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


###############################################################################
# Procedure


def execute_procedure(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_features_observations=None,
    path_file_source_table_groups_observations=None,
    path_file_source_table_parameters=None,
    column_identifier_observation=None,
    column_name_observation=None,
    name_batch=None,
    categories_batch=None,
    scale_before_interaction=None,
    method_scale=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_directory_dock_pail (str): path to pail directory for procedure's
            source and product directories and files
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_file_source_table_features_observations (object): path to source
            file in text format as a table
        path_file_source_table_parameters (str): path to source file in text
            format as a table with tab delimiters between columns and newline
            delimiters between rows, with parameters
        name_batch (str): name for current batch of parameters
        categories_batch (list<str>): names of categories that designate sets
            or groups of instances of parameters in a batch for execution
        report (str): whether to print reports

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
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        path_file_source_table_parameters=(
            path_file_source_table_parameters
        ),
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        name_batch=name_batch,
        categories_batch=categories_batch,
        scale_before_interaction=scale_before_interaction,
        method_scale=method_scale,
        report=report,
    )

    ##########
    # Organize information.
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: 7_interact_feature_products.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print(
            "path_directory_dock: " +
            str(pail_parameters["path_directory_dock"])
        )
        print(
            "path_directory_dock_pail: " +
            str(pail_parameters["path_directory_dock_pail"])
        )
        print(
            "path_directory_source: " +
            str(pail_parameters["path_directory_source"])
        )
        print(
            "path_directory_product: " +
            str(pail_parameters["path_directory_product"])
        )
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_directory_dock=pail_parameters["path_directory_dock"],
        path_directory_dock_pail=pail_parameters["path_directory_dock_pail"],
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_file_source_table_features_observations=(
            pail_parameters["path_file_source_table_features_observations"]
        ),
        path_file_source_table_groups_observations=(
            pail_parameters["path_file_source_table_groups_observations"]
        ),
        path_file_source_table_parameters=(
            pail_parameters["path_file_source_table_parameters"]
        ),
        report=pail_parameters["report"],
    )
    #pail_source["table_features_observations"]
    #pail_source["table_groups_observations"]
    #pail_source["table_parameters"]

    ##########
    # Organize parameters.

    # Organize parameters for interaction products between pairs of features.
    pail_instances = sutly.organize_filter_table_parameters(
        table=pail_source["table_parameters"],
        name_batch=pail_parameters["name_batch"],
        categories_batch=pail_parameters["categories_batch"],
        filter_parameters=True,
        report=False,
    )

    # Organize parameters about groups of observations.
    pail_groups = sutly.organize_parameters_groups_observations(
        table=pail_source["table_groups_observations"],
        column_name="abbreviation",
        report=pail_parameters["report"],
    )
    #pail_groups["table"]
    #pail_groups["names_groups_observations_sequence"]
    #pail_groups["categories_groups"]
    #pail_groups["records"]
    pail_observations = sutly.organize_parameters_further_groups_observations(
        table_observations=pail_source["table_features_observations"],
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_name_observation=pail_parameters["column_name_observation"],
        column_identifier_groups_observations=pail_parameters["column_identifier_observation"],
        instances_groups_observations=pail_groups["records"],
        key_name="abbreviation",
        names_groups_observations_sequence=(
            pail_groups["names_groups_observations_sequence"]
        ),
        report=pail_parameters["report"],
    )
    #pail_observations["table_observations_selection"]
    #pail_observations["observations_selection"]
    #pail_observations["translations_observations"]
    #pail_observations["names_groups_observations_sequence"]
    #pail_observations["groups_observations"]

    ##########
    # Filter observations.
    # It is important for any filters and stratifications of observations to
    # happen before calculations for transformation of scale and for
    # interaction products between pairs of features.

    # Extract identifiers of features.
    features_selection = copy.deepcopy(
        pail_source["table_features_observations"].columns.unique().tolist()
    )
    # Filter columns for features and rows for observations in main table.
    table_main = sutly.filter_table_columns_features_rows_observations(
        table=pail_source["table_features_observations"],
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_name_observation=pail_parameters["column_name_observation"],
        column_identifier_groups_observations=(
            pail_parameters["column_identifier_observation"]
        ),
        columns_categories=pail_groups["categories_groups"],
        features_selection=features_selection,
        observations_selection=pail_observations["observations_selection"],
        filter_table_main=True,
        report=pail_parameters["report"],
    )

    ##########
    # Transform features on quantitative, continuous scales of measurement to
    # standard z-score scale.
    pail_scale = transform_scale_features_by_parameters(
        table_features_observations=table_main,
        table_parameters=pail_instances["table"],
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_name_observation=pail_parameters["column_name_observation"],
        scale_before_interaction=pail_parameters["scale_before_interaction"],
        method_scale=pail_parameters["method_scale"],
        report=pail_parameters["report"],
    )
    #pail_scale["table_features_observations"]

    ##########
    # Calculate interaction products between pairs of features.

    # Calculate interaction products between pairs of features on quantitative,
    # continuous scales of measurement, either ratio or interval scales.
    pail_interaction = calculate_interaction_products_between_feature_pairs(
        table_features_observations=pail_scale["table_features_observations"],
        table_parameters=pail_instances["table"],
        report=pail_parameters["report"],
    )


    ##########
    # Bundle information.
    # Bundles of information for files.
    # Texts.
    #pail_write_texts = dict()
    # Objects.
    #pail_write_objects = dict()
    # Lists.
    pail_write_lists = dict()
    pail_write_lists[str("features_interaction")] = (
        pail_interaction["features_interaction"]
    )
    # Tables.
    pail_write_tables = dict()
    pail_write_tables[str("table_interaction")] = (
        pail_interaction["table_features_observations"]
    )

    ##########
    # Write product information to file.

    # Define paths to directories.
    path_directory_lists = os.path.join(
        path_directory_product, "lists",
    )
    path_directory_tables = os.path.join(
        path_directory_product, "tables",
    )
    path_directory_charts = os.path.join(
        path_directory_product, "charts",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_lists,
    )
    putly.create_directories(
        path=path_directory_tables,
    )
    putly.create_directories(
        path=path_directory_charts,
    )
    # Text.
    # Objects.
    # Lists.
    putly.write_lists_to_file_text(
        pail_write=pail_write_lists,
        path_directory=path_directory_lists,
        delimiter="\n",
    )
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
    # Plot charts and write to file.


    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_dock = sys.argv[1]
    path_directory_dock_pail = sys.argv[2]
    path_directory_source = sys.argv[3]
    path_directory_product = sys.argv[4]
    path_file_source_table_features_observations = sys.argv[5]
    path_file_source_table_groups_observations = sys.argv[6]
    path_file_source_table_parameters = sys.argv[7]
    column_identifier_observation = sys.argv[8]
    column_name_observation = sys.argv[9]
    name_batch = sys.argv[10]
    categories_batch = sys.argv[11]
    scale_before_interaction = sys.argv[12]
    method_scale = sys.argv[13]
    report = sys.argv[14]

    # Call function for procedure.
    execute_procedure(
        path_directory_dock=path_directory_dock,
        path_directory_dock_pail=path_directory_dock_pail,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        path_file_source_table_parameters=(
            path_file_source_table_parameters
        ),
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        name_batch=name_batch,
        categories_batch=categories_batch,
        scale_before_interaction=scale_before_interaction,
        method_scale=method_scale,
        report=report,
    )

    pass




###############################################################################
# End
