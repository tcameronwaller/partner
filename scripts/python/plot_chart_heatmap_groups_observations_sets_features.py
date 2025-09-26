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
# Date, first execution: 24 September 2025
# Date, last execution or modification: 24 September 2025
# Review: TCW; 24 September 2025
################################################################################
# Note


# The functions below might be helpful when organizing the signals between
# groups of observations.

#    partner.description.describe_features_from_table_columns_by_groups_rows()

#    partner.organization.
#    extract_array_values_from_column_by_separate_tables_rows()

#    partner.organization.
#    extract_array_values_from_table_column_by_groups_rows()




# Note: TCW; 25 September 2025
# It is possible to include features that do not explicitly belong to any of
# the allocation sets.

# Note: TCW; 25 September 2025
# The selection of features that the heatmap represents ultimately depends on
# the intersection of the three types of parameters below.
# 1. features with available in the table of signals (table_signals)
# 2. selection of features (features_selection)
# 3. special features for constraint on cluster (features_cluster_one,
# features_cluster_two)

# Note: TCW; 25 September 2025
# It is optional to designate one or two special sets of features within which
# to constrain the cluster operation on the features dimension. Since the
# cluster operation cannot accommodate any overlap between these special sets
# of features, it is necessary for the cluster sets to be mutually exclussive.
# If the user designates one or two of these special sets for cluster, then the
# total inclusive selection of features becomes the intersection of the custer
# set with the distinct selection set.


# Recent example of usage:
# /.../pails_process/omega3/2025-09-22_heterogeneity_candidate_adipose_fibrosis

##########
# Review: TCW; 26 September 2025

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
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams.update({'figure.max_open_warning': 0})
import matplotlib.pyplot
import matplotlib.lines

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


# Organize raw parameters.

def parse_text_parameters(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_observations=None,
    path_file_source_table_features=None,
    path_file_source_table_signals=None,
    path_file_source_table_groups_observations=None,
    path_file_source_table_sets_features=None,
    path_file_source_list_features_selection=None,
    path_file_source_list_features_cluster_one=None,
    path_file_source_list_features_cluster_two=None,
    column_identifier_observation=None,
    column_identifier_feature=None,
    column_identifier_signal=None,
    column_name_feature=None,
    transpose_table_signals=None,
    allow_replicate_observations=None,
    report=None,
):
    """
    Parse parameters from text.

    arguments:

    TODO: update documentation

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
    pail["path_file_source_table_observations"] = str(
        path_file_source_table_observations
    ).strip()
    pail["path_file_source_table_features"] = str(
        path_file_source_table_features
    ).strip()
    pail["path_file_source_table_signals"] = str(
        path_file_source_table_signals
    ).strip()
    pail["path_file_source_table_groups_observations"] = str(
        path_file_source_table_groups_observations
    ).strip()
    pail["path_file_source_table_sets_features"] = str(
        path_file_source_table_sets_features
    ).strip()
    pail["path_file_source_list_features_selection"] = str(
        path_file_source_list_features_selection
    ).strip()
    pail["path_file_source_list_features_cluster_one"] = str(
        path_file_source_list_features_cluster_one
    ).strip()
    pail["path_file_source_list_features_cluster_two"] = str(
        path_file_source_list_features_cluster_two
    ).strip()
    # Names of columns.
    pail["column_identifier_observation"] = str(
        column_identifier_observation
    ).strip()
    pail["column_identifier_feature"] = str(
        column_identifier_feature
    ).strip()
    pail["column_identifier_signal"] = str(
        column_identifier_signal
    ).strip()
    pail["column_name_feature"] = str(
        column_name_feature
    ).strip()

    # Boolean, true or false.
    if (
        (transpose_table_signals is not None) and
        (str(transpose_table_signals) != "") and
        (str(transpose_table_signals) != "none") and
        (str(transpose_table_signals) == "true")
    ):
        pail["transpose_table_signals"] = True
    else:
        pail["transpose_table_signals"] = False
        pass
    if (
        (allow_replicate_observations is not None) and
        (str(allow_replicate_observations) != "") and
        (str(allow_replicate_observations) != "none") and
        (str(allow_replicate_observations) == "true")
    ):
        pail["allow_replicate_observations"] = True
    else:
        pail["allow_replicate_observations"] = False
        pass
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
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        print("function: parse_text_parameters()")
        putly.print_terminal_partition(level=5)
        print("parameters:")
        print(pail)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Read source information.


def define_column_types_table_source(
    columns_source_text=None,
    columns_source_number=None,
):
    """
    Defines the types of variables in columns of table.

    Review: TCW; 5 May 2025

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
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_observations=None,
    path_file_source_table_features=None,
    path_file_source_table_signals=None,
    path_file_source_table_groups_observations=None,
    path_file_source_table_sets_features=None,
    path_file_source_list_features_selection=None,
    path_file_source_list_features_cluster_one=None,
    path_file_source_list_features_cluster_two=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 24 September 2025

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files

        ...


        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Bundle information.
    pail = dict()

    # Read information from file.
    pail["table_observations"] = pandas.read_csv(
        path_file_source_table_observations,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_features"] = pandas.read_csv(
        path_file_source_table_features,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_signals"] = pandas.read_csv(
        path_file_source_table_signals,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_groups"] = pandas.read_csv(
        path_file_source_table_groups_observations,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_sets"] = pandas.read_csv(
        path_file_source_table_sets_features,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["features_selection"] = putly.read_file_text_list(
        path_file=path_file_source_list_features_selection,
        delimiter="\n",
        unique=True,
    )
    pail["features_cluster_one"] = putly.read_file_text_list(
        path_file=path_file_source_list_features_cluster_one,
        delimiter="\n",
        unique=True,
    )
    pail["features_cluster_two"] = putly.read_file_text_list(
        path_file=path_file_source_list_features_cluster_two,
        delimiter="\n",
        unique=True,
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Organize parameters and information about observations and features.


def organize_parameters_groups_observations(
    table=None,
    column_name=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 24 September 2025

    arguments:
        table (object): Pandas data-frame table
        column_name (str): name of column to use for names of groups
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): bundle of information

    """

    # Copy information.
    table = table.copy(deep=True)

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

    # Filter rows in table.
    table_execution = table.loc[
        (table["execution"] == 1), :
    ].copy(deep=True)
    # Sort rows in table.
    table_execution.sort_values(
        by=["sequence",],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Organize indices in table.
    table_execution.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Extract names of groups.
    groups_sequence = copy.deepcopy(
        table_execution[column_name].unique().tolist()
    )

    # Collect information.
    categories_groups = list()
    records = list()
    for index, row in table_execution.iterrows():
        # Collect information and parameters from current row in table.
        record = dict()
        record["execution"] = int(row["execution"])
        record["sequence"] = row["sequence"]
        record["category"] = str(row["category"]).strip()
        record["name"] = str(row["name"]).strip() # name for group
        record["name_combination"] = "_".join([
            str(row["sequence"]).strip(),
            str(row["category"]).strip(),
            str(row["name"]).strip(),
        ])
        record["abbreviation"] = str(row["abbreviation"]).strip() # name for group
        record["selection_observations"] = (
            putly.parse_extract_text_keys_values_semicolon_colon_comma(
                text=row["selection_observations"],
            )
        )["features_values"]
        record["review"] = str(row["review"]).strip()
        record["note"] = str(row["note"]).strip()

        # Collect unique names of columns for features relevant to instance of
        # parameters from current row in table.
        categories_groups_instance = list()
        dictionaries = [
            "selection_observations",
        ]
        for dictionary in dictionaries:
            if record[dictionary] is not None:
                categories_groups_instance.extend(copy.deepcopy(list(
                    record[dictionary].keys()
                )))
                pass
            pass
        #categories_groups.extend(record["any_others"])
        categories_groups_instance = putly.collect_unique_items(
            items=categories_groups_instance,
        )
        record["categories_groups_instance"] = copy.deepcopy(
            categories_groups_instance
        )
        categories_groups.extend(copy.deepcopy(
            categories_groups_instance
        ))
        # Collect information and parameters for current row in table.
        records.append(record)
        pass

    # Names of columns for relevant features.
    categories_groups = putly.collect_unique_items(
        items=categories_groups,
    )

    # Bundle information.
    pail = dict()
    pail["table"] = table_execution
    pail["names_groups_observations_sequence"] = groups_sequence
    pail["categories_groups"] = categories_groups
    pail["records"] = records

    # Report.
    if report:
        # Organize.
        count_records = len(records)
        # Report.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        print("function: organize_parameters_groups_observations()")
        putly.print_terminal_partition(level=5)
        print("count of records: " + str(count_records))
        print("sequence of groups:")
        print(groups_sequence)
        pass
    # Return information.
    return pail


def combine_filter_sets_features(
    features_available=None,
    features_selection=None,
    features_cluster_one=None,
    features_cluster_two=None,
    report=None,
):
    """
    Combine and filter custom sets of features.

    Review: TCW; 25 September 2025

    arguments:
        features_available (list<str>): identifiers of features with available
            signals
        features_selection (list<str>): identifiers of features in custom
            selection for inclusion and representation
        features_cluster_one (list<str>): identifiers of features for
            constraint on cluster operation
        features_cluster_two (list<str>): identifiers of features for
            constraint on cluster operation
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): bundle of information

    """

    # Copy information.
    features_available = copy.deepcopy(features_available)
    features_selection = copy.deepcopy(features_selection)
    features_cluster_one = copy.deepcopy(features_cluster_one)
    features_cluster_two = copy.deepcopy(features_cluster_two)

    # Filter to features with available signals.
    features_selection_available = list(filter(
        lambda feature: (feature in features_available),
        features_selection
    ))
    features_cluster_one_available = list(filter(
        lambda feature: (feature in features_available),
        features_cluster_one
    ))
    features_cluster_two_available = list(filter(
        lambda feature: (feature in features_available),
        features_cluster_two
    ))

    # Take union of sets for constraint of cluster operation.
    features_cluster_union = putly.combine_sets_items_union_unique(
        sets_items=[
            features_cluster_one_available,
            features_cluster_two_available,
        ],
        report=None,
    )

    # Take intersection of features in custom selection and in union of
    # constraints for cluster operation.
    if (
        (features_cluster_union is not None) and
        (len(features_cluster_union) > 0)
    ):
        features_selection_intersection = (
            putly.combine_sets_items_intersection_unique(
                items_first=features_selection_available,
                items_second=features_cluster_union,
                report=None,
        ))
    else:
        features_selection_intersection = features_selection_available

    # Filter to features in selection.
    features_cluster_one_selection = list(filter(
        lambda feature: (feature in features_selection_intersection),
        features_cluster_one_available
    ))
    features_cluster_two_selection = list(filter(
        lambda feature: (feature in features_selection_intersection),
        features_cluster_two_available
    ))

    # Bundle information.
    pail = dict()
    # Filter to unique features.
    pail["features_available"] = putly.collect_unique_items(
        items=features_available,
    )
    pail["features_selection"] = putly.collect_unique_items(
        items=features_selection_intersection,
    )
    pail["features_cluster_one"] = putly.collect_unique_items(
        items=features_cluster_one_selection,
    )
    pail["features_cluster_two"] = putly.collect_unique_items(
        items=features_cluster_two_selection,
    )

    # Report.
    if report:
        # Organize.
        count_available = len(features_available)
        count_selection_source = len(features_selection)
        count_cluster_one_source = len(features_cluster_one)
        count_cluster_two_source = len(features_cluster_two)
        count_cluster_union = len(features_cluster_union)
        count_intersection = len(features_selection_intersection)
        count_selection_product = len(pail["features_selection"])
        count_cluster_one_product = len(pail["features_cluster_one"])
        count_cluster_two_product = len(pail["features_cluster_two"])

        # Report.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        print("function: combine_filter_sets_features()")
        putly.print_terminal_partition(level=5)
        print("counts of features in sets")
        putly.print_terminal_partition(level=5)
        print("features_available: " + str(count_available))
        print("--- source ---")
        print("features_selection: " + str(count_selection_source))
        print("features_cluster_one: " + str(count_cluster_one_source))
        print("features_cluster_two: " + str(count_cluster_two_source))
        print("--- intermediate ---")
        print("features_cluster_union: " + str(count_cluster_union))
        print("features_intersection: " + str(count_intersection))
        print("--- product ---")
        print("features_selection: " + str(count_selection_product))
        print("features_cluster_one: " + str(count_cluster_one_product))
        print("features_cluster_two: " + str(count_cluster_two_product))
        pass
    # Return information.
    return pail


def read_organize_parameters_sets_features(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    table=None,
    column_name=None,
    features_available=None,
    features_selection=None,
    features_cluster_one=None,
    features_cluster_two=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 24 September 2025

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        table (object): Pandas data-frame table
        column_name (str): name of column to use for names of sets
        features_available (list<str>): identifiers of features with available
            signals
        features_selection (list<str>): identifiers of features in custom
            selection for inclusion and representation
        features_cluster_one (list<str>): identifiers of features for
            constraint on cluster operation
        features_cluster_two (list<str>): identifiers of features for
            constraint on cluster operation
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): bundle of information

    """

    # Copy information.
    table = table.copy(deep=True)
    features_available = copy.deepcopy(features_available)
    features_selection = copy.deepcopy(features_selection)
    features_cluster_one = copy.deepcopy(features_cluster_one)
    features_cluster_two = copy.deepcopy(features_cluster_two)

    # Filter sets of features for custom selection and constraint of cluster
    # operation.
    pail_sets = combine_filter_sets_features(
        features_available=features_available,
        features_selection=features_selection,
        features_cluster_one=features_cluster_one,
        features_cluster_two=features_cluster_two,
        report=report,
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

    # Filter rows in table.
    table_execution = table.loc[
        (table["execution"] == 1), :
    ].copy(deep=True)
    # Sort rows in table.
    table_execution.sort_values(
        by=["sequence",],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Organize indices in table.
    table_execution.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Extract names of groups.
    sets_sequence = copy.deepcopy(
        table_execution[column_name].unique().tolist()
    )

    # Collect information.
    sets_features = dict()
    features_sets_union = list()
    records = list()
    for index, row in table_execution.iterrows():
        # Collect information and parameters from current row in table.
        record = dict()
        record["execution"] = int(row["execution"])
        record["sequence"] = row["sequence"]
        record["category"] = str(row["category"]).strip()
        record["name"] = str(row["name"]).strip() # name for group
        record["name_combination"] = "_".join([
            str(row["sequence"]).strip(),
            str(row["category"]).strip(),
            str(row["name"]).strip(),
        ])
        record["abbreviation"] = str(row["abbreviation"]).strip() # name for group
        record["directories_path_sets"] = (
            putly.parse_text_list_values(
                text=row["path_directory_sets_features"],
                delimiter=",",
        ))
        record["names_files_sets"] = (
            putly.parse_text_list_values(
                text=row["names_files_sets_features"],
                delimiter=",",
        ))
        record["review"] = str(row["review"]).strip()
        record["note"] = str(row["note"]).strip()

        # Read, organize, and collect sets of features.
        # A single set can be the union of identifiers from multiple separate
        # files.
        features_set = list()
        for name_file in record["names_files_sets"]:
            # Define paths to directories and files.
            pail_path = putly.extract_organize_path_directory_file(
                name_file=name_file,
                directories_path=record["directories_path_sets"],
                name_parent="dock",
                path_directory_parent=path_directory_dock,
                report=report,
            )
            # Determine whether parameters point path to a file exists.
            if (pail_path["existence_file"]):
                # Read information from file.
                features = putly.read_file_text_list(
                    path_file=pail_path["path_file"],
                    delimiter="\n",
                    unique=True,
                )
                # Collect features in set.
                features_set.extend(features)
                pass
            pass
        # Filter to features with available signals.
        features_set_available = list(filter(
            lambda feature: (feature in features_available),
            features_set
        ))
        # Collect unique features in set.
        record["features_set"] = putly.collect_unique_items(
            items=features_set_available,
        )
        sets_features[record[column_name]] = copy.deepcopy(
            record["features_set"]
        )
        features_sets_union.extend(record["features_set"])
        # Collect information and parameters for current row in table.
        records.append(record)
        pass

    # Selection of total inclusive features from all sets.
    features_sets_union = putly.collect_unique_items(
        items=features_sets_union,
    )

    # Bundle information.
    pail = dict()
    pail.update(pail_sets)
    pail["table"] = table_execution
    pail["names_sets_features_sequence"] = sets_sequence
    pail["sets_features"] = sets_features
    pail["features_sets_union"] = features_sets_union
    pail["records"] = records

    # Report.
    if report:
        # Organize.
        count_records = len(records)
        count_features_sets_union = len(features_sets_union)
        # Report.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        print("function: read_organize_parameters_sets_features()")
        putly.print_terminal_partition(level=5)
        print("count of records: " + str(count_records))
        print("sequence of sets:")
        print(sets_sequence)
        print(
            "count of total selection of features: " +
            str(count_features_sets_union)
        )
        #print("sets of features:")
        #print(sets_features)
        pass
    # Return information.
    return pail


def organize_parameters_further_groups_observations(
    table_observations=None,
    column_identifier_observation=None,
    column_identifier_signal=None,
    instances_groups_observations=None,
    key_name=None,
    names_groups_observations_sequence=None,
    report=None,
):
    """
    Organize parameters and information about observations.

    Review: TCW; 25 September 2025

    arguments:
        table_observations (object): Pandas data-frame table
        column_identifier_observation (str): name of column in source table
        column_identifier_signal (str): name of column in source table
        instances_groups_observations (list<dict>): multiple instances, each
            with parameters for the name and selection of a group of
            observations
            execution (int): logical binary indicator of whether to execute and
                handle the parameters for the current instance
            sequence (int): sequential index for instance's name and sort order
            category (str): categorical group of instances
            name (str): name or designator for instance of parameters
            name_combination (str): compound name for instance of parameters
            abbreviation (str): name abbreviation
            selection_observations (dict<list<str>>): names of columns in table
                for feature variables and their categorical values by which to
                filter rows for observations in table
            review (str):
            note (str):
        key_name (str): name of entry to use for names of groups
        names_groups_observations_sequence (list<str>): names of groups for
            observations in specific sequence
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Define subordinate functions for internal use.
    def alias_filter_extract(
        table=None,
        column_identifier=None,
        name=None,
        columns_categories=None,
        report=None,
    ):
        identifiers = (
            porg.filter_extract_table_row_identifiers_by_columns_categories(
                table=table,
                column_identifier=column_identifier,
                name=name,
                columns_categories=columns_categories,
                report=report,
        ))
        return identifiers

    # Copy information.
    table_observations = table_observations.copy(deep=True)
    instances_groups_observations = copy.deepcopy(
        instances_groups_observations
    )
    names_groups_observations_sequence = copy.deepcopy(
        names_groups_observations_sequence
    )

    # Extract and organize information about samples in groups.
    # Collect information.
    groups_observations = dict()
    observations_selection = list()
    # Iterate across instances of parameters.
    instances_execution = list(filter(
        lambda instance: (int(instance["execution"]) == 1),
        instances_groups_observations
    ))
    for instance in instances_execution:
        if (
            (names_groups_observations_sequence is not None) and
            (len(names_groups_observations_sequence) > 0) and
            (instance[key_name] in names_groups_observations_sequence)
        ):
            # Filter and extract identifiers of cohort sample observations
            # corresponding to selection criteria for current instance.
            observations = (
                alias_filter_extract(
                    table=table_observations,
                    column_identifier=column_identifier_signal,
                    name=instance[key_name],
                    columns_categories=instance["selection_observations"],
                    report=report,
            ))
            # Collect information.
            groups_observations[instance[key_name]] = observations
            observations_selection.extend(observations)
            pass
        pass
    # Collect unique names of sample observations.
    observations_selection = putly.collect_unique_items(
        items=observations_selection,
    )
    # Filter rows in table for selection of features.
    if (len(observations_selection) > 0):
        table_observations_selection = table_observations.loc[
            table_observations[column_identifier_signal].isin(
                observations_selection
            ), :
        ].copy(deep=True)
        pass

    # Bundle information.
    pail = dict()
    pail["table_observations_selection"] = table_observations_selection
    pail["observations_selection"] = observations_selection
    pail["translations_observations"] = None
    pail["names_groups_observations_sequence"] = (
        names_groups_observations_sequence
    )
    pail["groups_observations"] = groups_observations

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        print("function: organize_parameters_further_groups_observations()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def organize_parameters_further_sets_features(
    table_features=None,
    column_identifier_feature=None,
    column_name_feature=None,
    features_selection=None,
    features_cluster_one=None,
    features_cluster_two=None,
    features_sets_union=None,
    sets_features=None,
    names_sets_features_sequence=None,
    prefix_name_feature=None,
    report=None,
):
    """
    Organize parameters and information about features.

    Review: TCW; 25 September 2025

    arguments:
        table_features (object): Pandas data-frame table
        column_identifier_feature (str): name of column in source table
        column_name_feature (str): name of column in source table
        features_selection (list<str>): identifiers of features for which to
            include signals across observations
        sets_features (dict<list<str>>): names of sets and identifiers
            of features that belong to each of these sets
        names_sets_features_sequence (list<str>): names of sets for features in
            specific sequence
        prefix_name_feature (str): prefix for names of features
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): bundle of information

    """

    # Copy information.
    table_features = table_features.copy(deep=True)
    features_selection = copy.deepcopy(features_selection)
    features_selection_translation = copy.deepcopy(features_selection)
    features_cluster_one = copy.deepcopy(features_cluster_one)
    features_cluster_two = copy.deepcopy(features_cluster_two)
    features_sets_union = copy.deepcopy(features_sets_union)
    sets_features = copy.deepcopy(sets_features)
    names_sets_features_sequence = copy.deepcopy(names_sets_features_sequence)

    # Filter rows in table for selection of features.
    if (len(features_selection) > 0):
        table_features_selection = table_features.loc[
            table_features[column_identifier_feature].isin(
                features_selection
            ), :
        ].copy(deep=True)
        pass
    table_features_selection = table_features_selection.loc[
        (table_features_selection[column_name_feature].str.len() > 0), :
    ].copy(deep=True)
    table_features_selection = table_features_selection.loc[
        (table_features_selection[column_identifier_feature].str.len() > 0), :
    ].copy(deep=True)
    # Extract information for translation of names of columns.
    table_translations = table_features_selection.filter(
        items=[column_identifier_feature, column_name_feature,],
        axis="columns",
    )
    series_translations = pandas.Series(
        table_translations[column_name_feature].to_list(),
        index=table_translations[column_identifier_feature],
    )
    translations_features = series_translations.to_dict()

    # Translate names of genes in set for selection.
    if (translations_features is not None):
        features_selection_translation = list(map(
            lambda feature: (translations_features[feature]),
            features_selection_translation
        ))
        pass

    # Append prefix to names of columns in table for genes.
    features_selection_prefix = list()
    for name in features_selection_translation:
        name_prefix = str(prefix_name_feature + name)
        features_selection_prefix.append(name_prefix)
        pass
    translations_features_prefix = copy.deepcopy(translations_features)
    for key in translations_features_prefix.keys():
        name = str(translations_features_prefix[key])
        name_prefix = str(prefix_name_feature + name)
        translations_features_prefix[key] = name_prefix
        pass

    # Determine sets of features for constraint on cluster operation.
    if (
        (len(features_cluster_one) > 0) and
        (len(features_cluster_two) > 0)
    ):
        sets_features_cluster = dict()
        sets_features_cluster["one"] = copy.deepcopy(
            features_cluster_one
        )
        sets_features_cluster["two"] = copy.deepcopy(
            features_cluster_two
        )
        names_sets_features_sequence_cluster = ["one", "two",]
    elif (
        (len(features_cluster_one) > 0)
    ):
        sets_features_cluster = dict()
        sets_features_cluster["one"] = copy.deepcopy(
            features_cluster_one
        )
        names_sets_features_sequence_cluster = ["one",]
    elif (
        (len(features_sets_union) > 0)
    ):
        sets_features_cluster = dict()
        sets_features_cluster["union"] = copy.deepcopy(
            features_sets_union
        )
        names_sets_features_sequence_cluster = ["union",]
        pass

    # Bundle information.
    pail = dict()
    pail["table_features_selection"] = table_features_selection
    pail["features_selection"] = features_selection
    pail["features_selection_translation"] = features_selection_translation
    pail["features_selection_prefix"] = features_selection_prefix
    pail["translations_features"] = translations_features
    pail["translations_features_prefix"] = translations_features_prefix
    pail["names_sets_features_sequence"] = names_sets_features_sequence
    pail["names_sets_features_sequence_cluster"] = (
        names_sets_features_sequence_cluster
    )
    pail["sets_features"] = sets_features
    pail["sets_features_cluster"] = sets_features_cluster

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        print("function: organize_parameters_further_sets_features()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Prepare derivative, deliverable, product tables for individual signals and
# their means corresponding to features in sets across observations in groups.


def organize_preliminary_information_to_prepare_tables_signal(
    index_features=None,
    index_observations=None,
    features_all=None,
    observations_all=None,
    features_selection=None,
    observations_selection=None,
    sets_features=None,
    sets_features_cluster=None,
    groups_observations=None,
    names_sets_features_sequence=None,
    names_sets_features_sequence_cluster=None,
    names_groups_observations_sequence=None,
    translations_features=None,
    translations_observations=None,
    report=None,
):
    """
    Organize preliminary information for subsequent use in preparation of
    derivative tables.

    Features could correspond to columns in an original source table, or they
    could correspond to rows. Likewise, observations could correspond either
    to rows or columns. This function handles features and observations
    equivalently, so the difference is semantic.

    Review: 6 December 2024

    arguments:
        index_features (str): name for index corresponding to features across
            columns in the original source table
        index_observations (str): name for index corresponding to observations
            across rows in the original source table
        features_all (list<str>): identifiers of all features in the original
            source table
        observations_all (list<str>): identifiers of all observations in the
            original source table
        features_selection (list<str>): identifiers of features for which to
            include and describe values of signal intensity across observations
        observations_selection (list<str>): identifiers of observations for
            which to include and describe values of signal intensity across
            features
        sets_features (dict<list<str>>): names of sets and identifiers
            of features that belong to each of these sets; clustering of
            features will have constraint within these sets
        sets_features_cluster (dict<list<str>>): names of sets and identifiers
            of features that belong to each of these sets; alternative option
            since cluster operation cannot accommodate overlapping sets
        groups_observations (dict<list<str>>): names of groups and identifiers
            of observations that belong to each of these groups; clustering of
            observations will have constraint within these groups
        names_sets_features_sequence (list<str>): names of sets for
            features in specific sequence
        names_sets_features_sequence_cluster (list<str>): names of sets for
            features in specific sequence; alternative option since cluster
            operation cannot accommodate overlapping sets
        names_groups_observations_sequence (list<str>): names of groups for
            observations in specific sequence
        translations_features (dict<str>): translations for names or
            identifiers of features
        translations_observations (dict<str>): translations for names or
            identifiers of observations
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Organize and collect information.
    pail = dict()

    # Copy information.
    pail["index_features"] = index_features
    pail["index_observations"] = index_observations
    pail["features_all"] = copy.deepcopy(features_all)
    pail["observations_all"] = copy.deepcopy(observations_all)
    pail["features_selection"] = copy.deepcopy(features_selection)
    pail["observations_selection"] = copy.deepcopy(observations_selection)
    pail["sets_features"] = copy.deepcopy(sets_features)
    pail["sets_features_cluster"] = copy.deepcopy(sets_features_cluster)
    pail["groups_observations"] = copy.deepcopy(groups_observations)
    pail["names_sets_features_sequence"] = copy.deepcopy(
        names_sets_features_sequence
    )
    pail["names_sets_features_sequence_cluster"] = copy.deepcopy(
        names_sets_features_sequence_cluster
    )
    pail["names_groups_observations_sequence"] = copy.deepcopy(
        names_groups_observations_sequence
    )
    pail["translations_features"] = copy.deepcopy(translations_features)
    pail["translations_observations"] = copy.deepcopy(
        translations_observations
    )

    ##########
    # Prepare source information.
    # Prepare reference to sort rows or columns in a subsequent table.
    # Option 1.
    if True:
        # Features.
        sequence_sets_features = dict()
        index = 0
        for name in pail["names_sets_features_sequence"]:
            sequence_sets_features[name] = index
            index += 1
            pass
        # Observations.
        sequence_groups_observations = dict()
        index = 0
        for name in pail["names_groups_observations_sequence"]:
            sequence_groups_observations[name] = index
            index += 1
            pass
    # Option 2.
    if False:
        sequence_groups_observations = dict(zip(
            pail["names_groups_observations_sequence"],
            range(len(pail["names_groups_observations_sequence"]))
        ))
    # Option 3.
    if False:
        sequence_groups_observations = {
            key: value for value, key in enumerate(
                pail["names_groups_observations_sequence"]
            )
        }
    pail["sequence_sets_features"] = sequence_sets_features
    pail["sequence_groups_observations"] = sequence_groups_observations
    # Ensure that features are unique.
    features_selection_unique = putly.collect_unique_items(
        items=pail["features_selection"],
    )
    pail["features_selection_unique"] = features_selection_unique
    # Ensure that all features are in the source table.
    features_available = list(filter(
        lambda feature: (feature in pail["features_all"]),
        pail["features_selection_unique"]
    ))
    pail["features_available"] = features_available
    # Translate identifiers of features.
    features_available_translation = copy.deepcopy(
        features_available
    )
    if (pail["translations_features"] is not None):
        features_available_translation = list(map(
            lambda feature: (pail["translations_features"][feature]),
            features_available_translation
        ))
        pass
    # Ensure that features are unique after translation.
    features_available_translation = putly.collect_unique_items(
        items=features_available_translation,
    )
    pail["features_available_translation"] = features_available_translation

    # Ensure that observations are unique.
    observations_selection_unique = putly.collect_unique_items(
        items=pail["observations_selection"],
    )
    pail["observations_selection_unique"] = observations_selection_unique
    # Ensure that all observations are in the source table.
    observations_available = list(filter(
        lambda observation: (observation in pail["observations_all"]),
        observations_selection_unique
    ))
    pail["observations_available"] = observations_available
    # Translate identifiers of features.
    observations_available_translation = copy.deepcopy(
        observations_available
    )
    if (pail["translations_observations"] is not None):
        observations_available_translation = list(map(
            lambda observation: (
                pail["translations_observations"][observation]
            ),
            observations_available_translation
        ))
        pass
    # Ensure that features are unique after translation.
    observations_available_translation = putly.collect_unique_items(
        items=observations_available_translation,
    )
    pail["observations_available_translation"] = (
        observations_available_translation
    )

    # Return information.
    return pail


def prepare_tables_signals_groups_observations_sets_features(
    table=None,
    index_observations=None,
    index_features=None,
    observations_selection=None,
    features_selection=None,
    groups_observations=None,
    sets_features=None,
    sets_features_cluster=None,
    names_groups_observations_sequence=None,
    names_sets_features_sequence=None,
    names_sets_features_sequence_cluster=None,
    translations_observations=None,
    translations_features=None,
    transpose_source_table=None,
    allow_replicate_observations=None,
    report=None,
):
    """
    Prepare derivative tables of information about measurement signal
    intensities corresponding to individual features in sets across individual
    observations in groups.

    Any translations of identifiers or names of features and observations occur
    after the selection of features and observations. Hence identifiers or
    names for selection of features and observations must match those in the
    original source table.

    Each observation must only belong to a single group. That is, the groups of
    observations must be exclusive.

    Each feature can belong to multiple sets. That is, the sets of features are
    not exclusive.

    By intentional design, this function does not apply any transformation of
    scale or normalization of distribution to the values of signal intensity.

    ----------
    Format of source table (name: "table_source")
    ----------
    Format of source table is in wide format with floating-point values of
    signal intensities for measurements corresponding to individual features
    across columns and distinct individual observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. For versatility, this table does not have
    explicitly defined indices across rows or columns.
    ----------
    features        feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observation
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------
    Notice that it is also possible for the original source table to have a
    format corresponding to the simple transposition of the format in the
    description above. In this case, it is optional to transpose the table
    before further modification.

    ----------
    Format of product table 1 (name: "table_product_1")
    ----------
    Format of product table 1 is in wide format with floating-point values of
    signal intensities for measurements corresponding to individual features
    across columns and distinct individual observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. For versatility, this table does not have
    explicitly defined indices across rows or columns. This novel product table
    shares its format with the original source table. The difference between
    them is that this novel product table only includes only a specific
    selection of features and observations from the original source table.
    ----------
    features        feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observation
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    ----------
    Format of product table 2 (name: "table_product_2")
    ----------
    Format of product table 2 is in wide format with floating-point values of
    signal intensities for measurements corresponding to individual features
    across columns and distinct individual observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns. A special column gives identifiers or names corresponding to each
    observation across rows, and another special column provides names of
    categorical groups for these observations. For versatility, this table does
    not have explicitly defined indices across rows or columns.
    ----------
    features        group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_2   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_3   0.001     0.001     0.001     0.001     0.001
    ----------

    ----------
    Format of product tables 3, 4, and 5 (name: "table_product_3",
    "table_product_4", "table_product_5")
    ----------
    Format of product tables 3, 4, and 5 is in wide format with floating-point
    values of signal intensities for measurements corresponding to individual
    features across columns and distinct individual observations across rows. A
    special header row gives identifiers or names corresponding to each feature
    across columns. A special column gives identifiers or names corresponding
    to each observation across rows, and another special column provides names
    of categorical groups for these observations. For versatility, these tables
    do not have explicitly defined indices across rows or columns.

    Product tables 3, 4, and 5 are derivations from product table 2. This
    derivation includes standardization of values of signal intensities and
    clustering of columns for features and rows for observations by their
    respective values of signal intensities. The standardization transforms by
    z score the values of signal intensity for each feature such that these
    values have a mean of zero and standard deviation of one across all
    observations. This standardization simplifies the scale and distribution of
    values for subsequent visual representation on charts, especially heatmaps.

    The difference between product tables 3, 4, and 5 is whether the clustering
    of columns for features occurs with constrait by sets and whether the
    cluster of rows for observations occurs with constraint by groups.
    Clustering for table 3 has constraint by sets of features and by groups of
    observations. Clustering for table 4 has constraint by groups of
    observations only. Clustering for table 5 does not have any constraints.
    ----------
    features        group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_2   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_3   0.001     0.001     0.001     0.001     0.001
    ----------


    ----------
    Format of product table 6 (name: "table_product_6")
    ----------
    Format of product table 6 is in partial long format with floating-point
    values of statistics corresponding to type of descriptive statistic across
    columns and features and groups across rows.
    Product table 6 is a derivation from product table 2.
    ----------
    detail    group   mean standard_error standard_deviation median interqua...
    feature
    feature_1 group_1 0.01 0.001          0.001              0.015  0.5
    feature_1 group_2 0.01 0.001          0.001              0.015  0.5
    feature_1 group_3 0.01 0.001          0.001              0.015  0.5
    feature_1 group_4 0.01 0.001          0.001              0.015  0.5
    feature_2 group_1 0.01 0.001          0.001              0.015  0.5
    feature_2 group_2 0.01 0.001          0.001              0.015  0.5
    feature_2 group_3 0.01 0.001          0.001              0.015  0.5
    feature_2 group_4 0.01 0.001          0.001              0.015  0.5
    ----------

    ----------
    Format of product table 7 (name: "table_product_7")
    ----------
    Format of product table 7 is in partial long format with floating-point
    values of statistics corresponding to type of descriptive statistic across
    columns and features and groups across rows.
    Product table 7 is a derivation from product table 3.
    Product table 7 shares its format with product table 6. The difference
    between them is that product table 7 is a derivation after z score
    standardization of signals for each feature across all observations.
    ----------
    detail    group   mean standard_error standard_deviation median interqua...
    feature
    feature_1 group_1 0.01 0.001          0.001              0.015  0.5
    feature_1 group_2 0.01 0.001          0.001              0.015  0.5
    feature_1 group_3 0.01 0.001          0.001              0.015  0.5
    feature_1 group_4 0.01 0.001          0.001              0.015  0.5
    feature_2 group_1 0.01 0.001          0.001              0.015  0.5
    feature_2 group_2 0.01 0.001          0.001              0.015  0.5
    feature_2 group_3 0.01 0.001          0.001              0.015  0.5
    feature_2 group_4 0.01 0.001          0.001              0.015  0.5
    ----------

    ----------
    Format of product table 8 (name: "table_product_8")
    ----------
    Format of product table 8 is in wide format with floating-point values of
    a single, specific type of descriptive statistics (usually either mean or
    median) corresponding to features across rows and groups of observations
    across columns.
    Product table 8 is a derivation from product table 7.
    ----------
    group     group_1 group_2 group_3 group_4
    feature
    feature_1 0.01    0.001   0.001   0.015
    feature_2 0.01    0.001   0.001   0.015
    feature_3 0.01    0.001   0.001   0.015
    feature_4 0.01    0.001   0.001   0.015
    feature_5 0.01    0.001   0.001   0.015
    ----------

    Review: 25 September 2025
    Review: 24 December 2024

    arguments:
        table (object): Pandas data-frame table for values of signal intensity
            corresponding to features across columns and observations across
            rows (whether before or after any necessary transposition)
        index_observations (str): name for index corresponding to observations
            across rows in the original source table
        index_features (str): name for index corresponding to features across
            columns in the original source table
        observations_selection (list<str>): identifiers of observations for
            which to include and describe values of signal intensity across
            features
        features_selection (list<str>): identifiers of features for which to
            include and describe values of signal intensity across observations
        groups_observations (dict<list<str>>): names of groups and identifiers
            of observations that belong to each of these groups; clustering of
            observations will have constraint within these groups
        sets_features (dict<list<str>>): names of sets and identifiers
            of features that belong to each of these sets
        sets_features_cluster (dict<list<str>>): names of sets and identifiers
            of features that belong to each of these sets; alternative option
            since cluster operation cannot accommodate overlapping sets
        names_groups_observations_sequence (list<str>): names of groups for
            observations in specific sequence
        names_sets_features_sequence (list<str>): names of sets for features in
            specific sequence
        names_sets_features_sequence_cluster (list<str>): names of sets for
            features in specific sequence; alternative option since cluster
            operation cannot accommodate overlapping sets
        translations_observations (dict<str>): translations for names or
            identifiers of observations
        translations_features (dict<str>): translations for names or
            identifiers of features
        transpose_source_table (bool): whether to transpose the original source
            table before further transformation
        allow_replicate_observations (bool): whether to allow replicate
            observations or to require groups to be mutually exclusive, such
            that any individual observation can only belong to one group
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Copy information.
    table_source = table.copy(deep=True)

    ##########
    # Optional preliminary transposition.
    # Copy information.
    table_source_format = table_source.copy(deep=True)
    # Determine whether to apply optional transposition.
    if (transpose_source_table):
        # Organize indices in table.
        table_source_format.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table_source_format.set_index(
            index_features,
            append=False,
            drop=True,
            inplace=True
        )
        table_source_format.columns.rename(
            index_observations,
            inplace=True,
        ) # single-dimensional index
        # Transpose table.
        table_source_format = table_source_format.transpose(copy=True)
        # Organize indices in table.
        table_source_format.reset_index(
            level=None,
            inplace=True,
            drop=False, # remove index; do not move to regular columns
        )
        table_source_format.columns.rename(
            None,
            inplace=True,
        ) # single-dimensional index
        pass

    ##########
    # Organize preliminary information.
    features_all = copy.deepcopy(table_source_format.columns.to_list())
    features_all.remove(index_observations)
    observations_all = copy.deepcopy(
        table_source_format[index_observations].unique().tolist()
    )
    pail = organize_preliminary_information_to_prepare_tables_signal(
        index_features=index_features,
        index_observations=index_observations, # assigned in new tables
        features_all=features_all,
        observations_all=observations_all,
        features_selection=features_selection,
        observations_selection=observations_selection,
        sets_features=sets_features,
        sets_features_cluster=sets_features_cluster,
        groups_observations=groups_observations,
        names_sets_features_sequence=names_sets_features_sequence,
        names_sets_features_sequence_cluster=(
            names_sets_features_sequence_cluster
        ),
        names_groups_observations_sequence=(
            names_groups_observations_sequence
        ),
        translations_features=translations_features,
        translations_observations=translations_observations,
        report=report,
    )

    ##########
    # Tables 1, 2, 3, 4, and 5 for individual observations.

    ##########
    # Prepare product table 1.
    # Filter specific features and observations from table.
    table_selection = porg.filter_select_table_columns_rows_by_identifiers(
        table=table_source_format,
        index_rows=pail["index_observations"],
        identifiers_columns=pail["features_available"],
        identifiers_rows=pail["observations_available"],
        report=False,
    )
    # Copy information.
    table_product_1 = table_selection.copy(deep=True)
    # Translate names of features and observations.
    table_product_1_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product_1,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=pail["translations_observations"],
            remove_redundancy=True,
            report=False,
    ))

    ##########
    # Prepare product table 2.
    # Determine and fill groups of observations.
    # These next two functions create a new column named "group" and assign
    # categorical names corresponding to specific groups of observations.
    if (allow_replicate_observations):
        # Each observation has potential to belong to multiple groups.
        table_group = porg.determine_fill_table_groups_rows_with_replicates(
            table=table_selection,
            index_rows=pail["index_observations"],
            column_group="group",
            groups_rows=pail["groups_observations"],
            report=False,
        )
        # Create a new index for rows in the table.
        table_group["index_observations_old"] = (
            table_group[pail["index_observations"]]
        )
        table_group[pail["index_observations"]] = table_group.apply(
            lambda row: str(
                row["group"] + "_" + row["index_observations_old"]
            ),
            axis="columns", # apply function to each row
        )
        table_group.drop(
            labels=["index_observations_old",],
            axis="columns",
            inplace=True
        )
    else:
        # Each observation can only belong to a single group.
        table_group = porg.determine_fill_table_groups_rows(
            table=table_selection,
            index_rows=pail["index_observations"],
            column_group="group",
            groups_rows=pail["groups_observations"],
            report=False,
        )
        pass
    # Sort rows in table by groups.
    table_group = porg.sort_table_rows_by_single_column_reference(
        table=table_group,
        index_rows=pail["index_observations"],
        column_reference="group",
        column_sort_temporary="sort_temporary",
        reference_sort=pail["sequence_groups_observations"],
    )
    # Copy information.
    table_product_2 = table_group.copy(deep=True)
    # Translate names of features and observations.
    table_product_2_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product_2,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=pail["translations_observations"],
            remove_redundancy=True,
            report=False,
    ))

    ##########
    # Prepare product table 3.
    # Calculate z scores of values for each feature across observations to
    # standardize their scales and distributions.
    table_scale = pscl.transform_standard_z_score_by_table_columns(
        table=table_group,
        columns=pail["features_available"],
        report=False,
    )
    # Cluster rows in table by groups of observations.
    table_product_3 = porg.cluster_table_rows_by_group(
        table=table_scale,
        index_rows=pail["index_observations"],
        column_group="group",
    )
    # Sort rows in table by groups.
    table_product_3 = porg.sort_table_rows_by_single_column_reference(
        table=table_product_3,
        index_rows=pail["index_observations"],
        column_reference="group",
        column_sort_temporary="sort_temporary",
        reference_sort=pail["sequence_groups_observations"],
    )

    # Cluster columns in table by sets of features.
    # There is an error if there is any overlap between sets of features.
    # Support the option of a separate set of features that is simpler without
    # any overlap between sets.
    table_product_3 = porg.cluster_table_columns_by_external_group(
        table=table_product_3,
        indices_rows=[pail["index_observations"], "group",],
        groups_columns=pail["sets_features_cluster"],
        names_groups_sequence=pail["names_sets_features_sequence_cluster"],
        report=False,
    )

    # Translate names of features and observations.
    table_product_3_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product_3,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=pail["translations_observations"],
            remove_redundancy=True,
            report=False,
    ))

    ##########
    # Prepare product table 4.
    # Copy information.
    table_product_4 = table_scale.copy(deep=True)
    # Cluster rows in table by groups of observations.
    table_product_4 = porg.cluster_table_rows_by_group(
        table=table_product_4,
        index_rows=pail["index_observations"],
        column_group="group",
    )
    # Sort rows in table by groups.
    table_product_4 = porg.sort_table_rows_by_single_column_reference(
        table=table_product_4,
        index_rows=pail["index_observations"],
        column_reference="group",
        column_sort_temporary="sort_temporary",
        reference_sort=pail["sequence_groups_observations"],
    )
    # Organize indices in table.
    table_product_4.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_product_4.set_index(
        [pail["index_observations"], "group"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Cluster columns in table.
    table_product_4 = porg.cluster_table_columns(
        table=table_product_4,
    )
    table_product_4.index = pandas.MultiIndex.from_tuples(
        table_product_4.index,
        names=[pail["index_observations"], "group"]
    )
    # Organize indices in table.
    table_product_4.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Translate names of features and observations.
    table_product_4_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product_4,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=pail["translations_observations"],
            remove_redundancy=True,
            report=False,
    ))

    ##########
    # Prepare product table 5.
    # Copy information.
    table_product_5 = table_scale.copy(deep=True)
    # Organize indices in table.
    table_product_5.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_product_5.set_index(
        [pail["index_observations"], "group"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Cluster rows in table.
    table_product_5 = porg.cluster_table_rows(
        table=table_product_5,
    )
    table_product_5.index = pandas.MultiIndex.from_tuples(
        table_product_5.index,
        names=[pail["index_observations"], "group"]
    )
    # Cluster columns in table.
    table_product_5 = porg.cluster_table_columns(
        table=table_product_5,
    )
    table_product_5.index = pandas.MultiIndex.from_tuples(
        table_product_5.index,
        names=[pail["index_observations"], "group"]
    )
    # Organize indices in table.
    table_product_5.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Translate names of features and observations.
    table_product_5_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_product_5,
            index_rows=pail["index_observations"],
            translations_columns=pail["translations_features"],
            translations_rows=pail["translations_observations"],
            remove_redundancy=True,
            report=False,
    ))


    # Extract available groups in table.
    #groups_observations_sequence_available = (
    #    pail["names_groups_observations_sequence"]
    #)
    groups_observations_sequence_available = copy.deepcopy(
            table_product_2_translation["group"].unique().tolist()
    )

    ##########
    # Tables 6, 7, and 8 for summaries using descriptive statistics on groups
    # of observations.

    ##########
    # Prepare product table 6.
    # Calculate descriptive statistics for each feature across observations.
    #table_product_6 = pdesc.describe_table_features_by_groups(
    #    table=table_product_2_translation,
    #    column_group="group",
    #    columns_features=pail["features_available_translation"],
    #    index_feature=pail["index_features"],
    #    translations_feature=None,
    #    threshold_observations=5,
    #    digits_round=3,
    #    report=False,
    #)
    table_product_6 = (
        pdesc.describe_features_from_table_columns_by_groups_rows(
            table_group=table_product_2_translation,
            index_columns=pail["index_features"],
            index_rows=pail["index_observations"],
            index_features=pail["index_features"],
            column_group="group",
            groups_sequence=groups_observations_sequence_available,
            columns_features=pail["features_available_translation"],
            translations_feature=None,
            key_group="group",
            threshold_observations=5,
            digits_round=3,
            ttest_one=None,
            ttest_two=None,
            ttest_three=None,
            report=False,
    ))

    ##########
    # Prepare product table 7.
    #table_product_7 = pdesc.describe_table_features_by_groups(
    #    table=table_product_3_translation,
    #    column_group="group",
    #    columns_features=pail["features_available_translation"],
    #    index_feature=pail["index_features"],
    #    translations_feature=None,
    #    threshold_observations=5,
    #    digits_round=3,
    #    report=False,
    #)
    table_product_7 = (
        pdesc.describe_features_from_table_columns_by_groups_rows(
            table_group=table_product_3_translation,
            index_columns=pail["index_features"],
            index_rows=pail["index_observations"],
            index_features=pail["index_features"],
            column_group="group",
            groups_sequence=groups_observations_sequence_available,
            columns_features=pail["features_available_translation"],
            translations_feature=None,
            key_group="group",
            threshold_observations=5,
            digits_round=3,
            ttest_one=None,
            ttest_two=None,
            ttest_three=None,
            report=False,
    ))

    ##########
    # Prepare product table 8.
    pail_table_8 = (
        pdesc.prepare_table_signal_summaries_features_observation_groups(
            table=table_scale,
            index_columns=pail["index_features"],
            column_group="group",
            columns_features=pail["features_available"],
            index_rows=pail["index_observations"],
            names_groups_observations_sequence=(
                groups_observations_sequence_available
            ),
            summary="mean", # "mean", "median",
            translations_features=pail["translations_features"],
            report=report,
    ))

    ##########
    # Collect information.
    pail_return = dict()
    pail_return["table_1"] = table_product_1
    pail_return["table_1_translation"] = table_product_1_translation
    pail_return["table_2"] = table_product_2
    pail_return["table_2_translation"] = table_product_2_translation
    pail_return["table_3"] = table_product_3
    pail_return["table_3_translation"] = table_product_3_translation
    pail_return["table_4"] = table_product_4
    pail_return["table_4_translation"] = table_product_4_translation
    pail_return["table_5"] = table_product_5
    pail_return["table_5_translation"] = table_product_5_translation
    pail_return["table_6"] = table_product_6
    pail_return["table_7"] = table_product_7
    pail_return["table_8"] = pail_table_8["table_summary"]
    pail_return["table_8_translation"] = (
       pail_table_8["table_summary_translation"]
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        function = str(
            "prepare_tables_signals_individual_features_sets_observations_" +
            "groups()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail_return


# Prepare table for allocation of features to sets.


def assemble_table_features_sets_allocation(
    features_selection=None,
    sets_features=None,
    names_sets_features_sequence=None,
    group_other=None,
    index_features=None,
    report=None,
):
    """
    Assemble table with logical binary designations of allocation for a
    selection of features to sets.

    ----------
    Format of product table
    ----------
    sets        set_1 set_2 set_3
    features
    feature_1   1     0     0
    feature_2   0     1     0
    feature_3   0     0     1
    feature_4   0     1     0
    feature_5   1     0     0
    ----------

    Review: 12 December 2024

    arguments:
        features_selection (list<str>): identifiers of features to allocate to
            sets
        sets_features (dict<list<str>>): names of sets and identifiers
            of features that belong to each of these sets
        names_sets_features_sequence (list<str>): names of sets for
            features in specific sequence
        group_other (str): name for group to which to assign all features, even
            without allocation to any other groups
        index_features (str): name for index corresponding to features across
            rows in the allocation table
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy information.
    features_selection = copy.deepcopy(features_selection)
    sets_features = copy.deepcopy(sets_features)
    names_sets_features_sequence = copy.deepcopy(
        names_sets_features_sequence
    )

    # Determine whether to include column for other sets.
    if (
        (group_other is not None) and
        (len(group_other) > 0)
    ):
        names_sets_features_sequence.append(group_other)
        pass

    ##########
    # Determine set allocations for each relevant gene.
    # Iterate on relevant genes.
    # Collect information about set allocations for each gene.
    records = list()
    for feature in features_selection:
        # Collect information.
        record = dict()
        record[index_features] = feature
        if (
            (group_other is not None) and
            (len(group_other) > 0)
        ):
            record[group_other] = 1
        # Iterate on sets for allocation.
        for name_set in sets_features.keys():
            if (feature in sets_features[name_set]):
                record[name_set] = 1
                if (
                    (group_other is not None) and
                    (len(group_other) > 0)
                ):
                    record[group_other] = 0
            else:
                record[name_set] = 0
            pass
        # Collect information.
        records.append(record)
        pass
    # Create pandas data-frame table.
    table = pandas.DataFrame(data=records)

    # Filter and sort columns within table.
    columns_sequence = copy.deepcopy(names_sets_features_sequence)
    columns_sequence.insert(0, index_features)
    if (
        (group_other is not None) and
        (len(group_other) > 0)
    ):
        columns_sequence.append("other")
        pass
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=columns_sequence,
        report=False,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: age_exercise.proteomics")
        print("module: organize_subject.py")
        function = str(
            "assemble_table_features_sets_allocation()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        print("table with allocations of features to sets")
        print(table)
        # Determine whether to describe column for other sets.
        if (
            (group_other is not None) and
            (len(group_other) > 0)
        ):
            # Filter rows within table.
            table_nonother = table.loc[
                (table[group_other] == 0), :
            ].copy(deep=True)
            putly.print_terminal_partition(level=5)
            print(table_nonother)
            pass
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


def sort_check_table_rows_sequence_custom(
    table=None,
    index_rows=None,
    values_sort_sequence=None,
    report=None,
):
    """
    Sort rows in a table by a custom sequence of values corresponding to the
    table's index across rows. Check that the table's rows after sort match the
    custom sequence.

    arguments:
        table (object): Pandas data-frame table
        index_rows (str): name for index across rows by which to sort rows
        values_sort_sequence (list<str>): values corresponding to index across
            rows in custom sequence by which to sort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information.
    table = table.copy(deep=True)
    values_sort_sequence = copy.deepcopy(values_sort_sequence)

    # Prepare indices for sort.
    reference_sort = dict(zip(
        values_sort_sequence, range(len(values_sort_sequence))
    ))

    # Sort rows in table reference indices.
    table_sort = porg.sort_table_rows_by_single_column_reference(
        table=table,
        index_rows=index_rows,
        column_reference=index_rows,
        column_sort_temporary="sort_temporary_982416",
        reference_sort=reference_sort,
    )

    # Report.
    if report:
        # Extract index across rows for comparison.
        values_index_rows = copy.deepcopy(
            table_sort[index_rows].to_list()
        )
        # Count identifiers of features.
        count_sequence = len(values_sort_sequence)
        count_sort = len(values_index_rows)
        # Confirm that sets of indices from both sources are inclusive.
        check_inclusion = putly.compare_lists_by_mutual_inclusion(
            list_primary=values_sort_sequence,
            list_secondary=values_index_rows,
        )
        # Confirm that sets of indices from both sources are identical across
        # their respective sequences.
        check_identity = putly.compare_lists_by_elemental_identity(
            list_primary=values_sort_sequence,
            list_secondary=values_index_rows,
        )
        # Confirm that sets of indices from both sources are equal.
        check_equality = (
            values_sort_sequence == values_index_rows
        )

        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        function = str(
            "sort_check_table_rows_sequence_custom()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=4)
        print(
            "count of values in sort sequence: " +
            str(count_sequence)
        )
        print(
            "count of values of index in table: " +
            str(count_sort)
        )
        putly.print_terminal_partition(level=5)
        print(
            "confirm that values of index across rows are identical to " +
            "reference"
        )
        print("check inclusion: " + str(check_inclusion))
        print("check identity: " + str(check_identity))
        print("check equality: " + str(check_equality))
        putly.print_terminal_partition(level=5)
        print("product table after sort")
        print(table_sort)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_sort


def prepare_table_sets_features_allocation_match_table_signal(
    table_signal=None,
    index_features=None,
    indices_observations=None,
    sets_features=None,
    names_sets_features_sequence=None,
    translations_features=None,
    report=None,
):
    """
    Prepare tables for allocation of features to sets that match with a table
    of signals.

    Optional translation of the identifiers or names of features occurs before
    assembly of the allocation table and before comparison and sort to match
    the signal table.

    ----------
    Format of source table for signals (name: "table_signal")
    ----------
    features        group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_2   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_3   0.001     0.001     0.001     0.001     0.001
    ----------

    Review: 30 December 2024

    arguments:
        table_signal (object): Pandas data-frame table of values of signal
            intensity corresponding to features across columns and observations
            across rows
        index_features (str): name for index corresponding to features across
            columns in the signal table and across rows in the allocation table
        indices_observations (str): names for indices corresponding to
            observations across rows in the signal table
        features_selection (list<str>): identifiers of features for which to
            include and describe values of signal intensity across observations
        sets_features (dict<list<str>>): names of sets and identifiers
            of features that belong to each of these sets
        names_sets_features_sequence (list<str>): names of sets for
            features in specific sequence
        translations_features (dict<str>): translations for names or
            identifiers of features
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Copy information.
    table_signal = table_signal.copy(deep=True)
    sets_features = copy.deepcopy(sets_features)
    names_sets_features_sequence = copy.deepcopy(
        names_sets_features_sequence
    )

    # Extract identifiers of features from the columns in the table of signals.
    features_signal = copy.deepcopy(table_signal.columns.to_list())
    for index in indices_observations:
        features_signal.remove(index)
        pass
    features_signal_unique = putly.collect_unique_items(
        items=features_signal,
    )

    # Determine whether to translate identifiers of features in sets to match
    # the identifiers of features from the table of signals.
    if (translations_features is not None):
        sets_features_translation = dict()
        for group_features in sets_features.keys():
            features_group = sets_features[group_features]
            features_group_translation = list(map(
                lambda feature: (translations_features[feature]),
                features_group
            ))
            sets_features_translation[group_features] = (
                features_group_translation
            )
            pass
    else:
        sets_features_translation = copy.deepcopy(sets_features)
        pass

    # Assemble table with logical binary designations of allocation for a
    # selection of features to sets.
    table_allocation = assemble_table_features_sets_allocation(
        features_selection=features_signal_unique,
        sets_features=sets_features_translation,
        names_sets_features_sequence=names_sets_features_sequence,
        group_other=None, # None, "", or "other"
        index_features=index_features,
        report=report,
    )
    # Sort the sequence of features across rows in the set allocation table to
    # match the sequence of features across columns in the signal table.
    # Extract identifiers of features in original sequence from the columns in
    # the table of signals.
    features_sequence = copy.deepcopy(table_signal.columns.to_list())
    for index in indices_observations:
        features_sequence.remove(index)
        pass
    table_allocation_sort = sort_check_table_rows_sequence_custom(
        table=table_allocation,
        index_rows=index_features,
        values_sort_sequence=features_sequence,
        report=report,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        function = str(
            "prepare_table_sets_features_allocation_match_table_signal()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_allocation_sort


def manage_prepare_tables(
    table_observations=None,
    table_features=None,
    table_signals=None,
    column_identifier_observation=None,
    column_identifier_feature=None,
    column_identifier_signal=None,
    column_name_feature=None,
    instances_groups_observations=None,
    key_name=None,
    names_groups_observations_sequence=None,
    categories_groups=None,
    features_selection=None,
    features_cluster_one=None,
    features_cluster_two=None,
    features_sets_union=None,
    sets_features=None,
    names_sets_features_sequence=None,
    transpose_table_signals=None,
    allow_replicate_observations=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 25 September 2025

    arguments:

    TODO: update documentation

    key_name (str): name of entry to use for names of groups

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table_observations = table_observations.copy(deep=True)
    table_features = table_features.copy(deep=True)
    table_signals = table_signals.copy(deep=True)
    instances_groups_observations = copy.deepcopy(
        instances_groups_observations
    )
    names_groups_observations_sequence = copy.deepcopy(
        names_groups_observations_sequence
    )
    categories_groups = copy.deepcopy(categories_groups)
    features_selection = copy.deepcopy(features_selection)
    features_cluster_one = copy.deepcopy(features_cluster_one)
    features_cluster_two = copy.deepcopy(features_cluster_two)
    features_sets_union = copy.deepcopy(features_sets_union)
    sets_features = copy.deepcopy(sets_features)
    names_sets_features_sequence = copy.deepcopy(names_sets_features_sequence)

    # Organize parameters further for observations and features.

    pail_observations = organize_parameters_further_groups_observations(
        table_observations=table_observations,
        column_identifier_observation=column_identifier_observation,
        column_identifier_signal=column_identifier_signal,
        instances_groups_observations=instances_groups_observations,
        key_name=key_name,
        names_groups_observations_sequence=names_groups_observations_sequence,
        report=report,
    )
    #pail_observations["table_observations_selection"]
    #pail_observations["observations_selection"]
    #pail_observations["translations_observations"]
    #pail_observations["names_groups_observations_sequence"]
    #pail_observations["groups_observations"]

    pail_features = organize_parameters_further_sets_features(
        table_features=table_features,
        column_identifier_feature=column_identifier_feature,
        column_name_feature=column_name_feature,
        features_selection=features_selection,
        features_cluster_one=features_cluster_one,
        features_cluster_two=features_cluster_two,
        features_sets_union=features_sets_union,
        sets_features=sets_features,
        names_sets_features_sequence=names_sets_features_sequence,
        prefix_name_feature="",
        report=report,
    )
    #pail_features["table_features_selection"]
    #pail_features["features_selection"]
    #pail_features["features_selection_translation"]
    #pail_features["features_selection_prefix"]
    #pail_features["translations_features"]
    #pail_features["translations_features_prefix"]
    #pail_features["names_sets_features_sequence"]
    #pail_features["sets_features"]

    # Note: TCW; 25 September 2025
    # The identifiers extracted for observations in the two functions below
    # must match.
    # organize_parameters_further_groups_observations()
    # prepare_tables_signals_groups_observations_sets_features()

    ##########
    # Prepare basic tables.
    # Prepare tables for signals.
    pail_tables_signal = (
        prepare_tables_signals_groups_observations_sets_features(
            table=table_signals,
            index_observations=column_identifier_signal,
            index_features=column_identifier_feature,
            observations_selection=pail_observations["observations_selection"],
            features_selection=pail_features["features_selection"],
            groups_observations=pail_observations["groups_observations"],
            sets_features=pail_features["sets_features"],
            sets_features_cluster=pail_features["sets_features_cluster"],
            names_groups_observations_sequence=(
                pail_observations["names_groups_observations_sequence"]
            ),
            names_sets_features_sequence=(
                pail_features["names_sets_features_sequence"]
            ),
            names_sets_features_sequence_cluster=(
                pail_features["names_sets_features_sequence_cluster"]
            ),
            translations_observations=(
                pail_observations["translations_observations"]
            ),
            translations_features=pail_features["translations_features"],
            transpose_source_table=transpose_table_signals,
            allow_replicate_observations=allow_replicate_observations,
            report=False,
        )
    )

    # Extract available groups in table.
    #groups_observations_sequence_available = (
    #    pail_samples["names_groups_samples_sequence"]
    #)
    #groups_observations_sequence_available = copy.deepcopy(
    #        pail_tables["table_2_translation"]["group"].unique().tolist()
    #)

    # Prepare tables for allocation of features (genes) to sets with sequence
    # of features that matches those in tables of signals.
    table_allocation_3 = (
        prepare_table_sets_features_allocation_match_table_signal(
            table_signal=pail_tables_signal["table_3"],
            index_features=column_identifier_feature,
            indices_observations=[column_identifier_signal, "group",],
            sets_features=pail_features["sets_features"],
            names_sets_features_sequence=(
                pail_features["names_sets_features_sequence"]
            ),
            translations_features=None,
            report=report,
        )
    )
    table_allocation_4 = (
        prepare_table_sets_features_allocation_match_table_signal(
            table_signal=pail_tables_signal["table_4"],
            index_features=column_identifier_feature,
            indices_observations=[column_identifier_signal, "group",],
            sets_features=pail_features["sets_features"],
            names_sets_features_sequence=(
                pail_features["names_sets_features_sequence"]
            ),
            translations_features=None,
            report=report,
        )
    )
    table_allocation_5 = (
        prepare_table_sets_features_allocation_match_table_signal(
            table_signal=pail_tables_signal["table_5"],
            index_features=column_identifier_feature,
            indices_observations=[column_identifier_signal, "group",],
            sets_features=pail_features["sets_features"],
            names_sets_features_sequence=(
                pail_features["names_sets_features_sequence"]
            ),
            translations_features=None,
            report=report,
        )
    )
    # Prepare allocation to feature sets in sort sequence that matches table of
    # mean signals.
    features_sequence = copy.deepcopy(
        pail_tables_signal["table_8"][column_identifier_feature].tolist()
    )
    table_allocation_8 = sort_check_table_rows_sequence_custom(
        table=table_allocation_3,
        index_rows=column_identifier_feature,
        values_sort_sequence=features_sequence,
        report=report,
    )

    # Bundle information.
    pail = dict()
    pail.update(pail_observations)
    pail.update(pail_features)
    pail.update(pail_tables_signal)
    pail["table_allocation_3"] = table_allocation_3 # cluster with constraint by sets of features and by groups of observations
    pail["table_allocation_4"] = table_allocation_4 # cluster with constraint by groups of observations only
    pail["table_allocation_5"] = table_allocation_5 # cluster without constraint
    pail["table_allocation_8"] = table_allocation_8

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        print("function: manage_prepare_tables()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Create plot charts.
# 1.
# function: plot_heatmap_signal_features_observations_labels()
# signal: individual signals or summaries (such as mean or median) for small
# counts of features across observations or groups of observations
# design: explicit labels for individual features on vertical ordinate axis and
# for individual observations or groups on horizontal abscissa axis
# 2.
# function: plot_heatmap_signal_features_sets_observations_labels()
# signal: individual signals or summaries (such as mean or median) for small or
# large counts of features which belong to sets across observations or groups
# of observations
# design: graphical representation of allocation of each feature to overlapping
# sets and explicit labels for individual observations or groups
# 3.
# function: plot_heatmap_signal_features_sets_observations_groups()
# signal: individual signals for large counts of features which belong to sets
# across observations which belong to groups
# design: graphical representation of allocation of each feature to overlapping
# sets and graphical representation of observations in mutually exclusive
# groups


def extract_prepare_table_signals_categories_for_heatmap(
    table=None,
    format_table=None,
    index_columns=None,
    index_rows=None,
    column_group=None,
    transpose_table=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_maximum=None,
    labels_ordinate_categories=None,
    labels_abscissa_categories=None,
    report=None,
):
    """
    Extract from table and prepare signals and categorical labels for heatmap.

    Format of source table

    Option 1.
    Format of source table is in wide format with floating-point values of
    signal intensities or a single, specific type of descriptive statistics
    (usually either mean or median) corresponding to features across rows and
    observations or groups of observations across columns.
    ----------
    observation_or_group group_1 group_2 group_3 group_4
    feature
    feature_1            0.01    0.001   0.001   0.015
    feature_2            0.01    0.001   0.001   0.015
    feature_3            -0.01   0.001   0.001   0.015
    feature_4            -0.01   0.001   0.001   0.015
    feature_5            -0.01   0.001   0.001   0.015
    ----------

    Option 2.
    Format of source table is in wide format with floating-point values of
    signal intensities corresponding to features across columns and distinct
    observations across rows. A special column gives identifiers corresponding
    to each observation across rows. Another special column provides names
    of categorical groups of observations.
    ----------
    observation     group   feature_1 feature_2 feature_3 feature_4 feature_5
    observation_1   group_1 0.001     0.001     0.001     0.001     0.001
    observation_2   group_1 0.001     0.001     0.001     0.001     0.001
    observation_3   group_2 0.001     0.001     0.001     0.001     0.001
    observation_4   group_2 0.001     0.001     0.001     0.001     0.001
    observation_5   group_3 0.001     0.001     0.001     0.001     0.001
    ----------

    For versatility, the source table does not have explicitly defined indices
    across columns or rows.

    This function preserves the original sequence of features. This function
    also preserves the original sequence of groups and observations within
    groups.

    Review: 27 December 2024

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            for features across sample observations or groups of sample
            observations
        format_table (int): value 1 for features across rows and observations
            or groups of observations across columns, value 2 for features
            across columns and observations across rows with potential special
            column for groups of observations
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        column_group (str): name of column in table to use for groups
        transpose_table (bool): whether to transpose matrix from table
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_missing_fill (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
        labels_ordinate_categories (list<str>): optional, explicit labels for
            ordinate or vertical axis
        labels_abscissa_categories (list<str>): optional, explicit labels for
            abscissa or horizontal axis
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information for figure

    """

    # Copy information in table.
    table_signal = table.copy(deep=True)
    # Organize indices in table.
    table_signal.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_signal.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    if (format_table == 1):
        table_signal.set_index(
            index_rows,
            append=False,
            drop=True,
            inplace=True
        )
    elif ((format_table == 2) and (column_group is None)):
        table_signal.set_index(
            index_rows,
            append=False,
            drop=True,
            inplace=True
        )
    elif ((format_table == 2) and (column_group is not None)):
        table_signal.set_index(
            [index_rows, column_group,],
            append=False,
            drop=True,
            inplace=True
        )
        pass
    # Extract categorical names from original indices in table.
    #labels_columns = copy.deepcopy(table_signal.columns.to_list())
    #labels_rows = table.index.to_list()
    #labels_rows = table_signal["group_primary"].to_list()
    labels_index_columns = copy.deepcopy(
        table_signal.columns.get_level_values(
            index_columns
        ).to_list()
    )
    labels_index_rows = copy.deepcopy(
        table_signal.index.get_level_values(
            index_rows
        ).to_list()
    )
    if ((format_table == 2) and (column_group is not None)):
        labels_group_rows = copy.deepcopy(
            table_signal.index.get_level_values(
                column_group
            ).to_list()
        )
        labels_group_unique = putly.collect_unique_elements(
            elements=labels_group_rows,
        )
    else:
        labels_group_rows = None
        labels_group_unique = None
        pass
    # Extract counts of individual observations in each group.
    if ((format_table == 2) and (column_group is not None)):
        groups_counts = table_signal.groupby(
            level=column_group,
        ).size().to_dict()
        pass
    else:
        groups_counts = None
        pass

    # Transpose table.
    if (transpose_table):
        table_signal = table_signal.transpose(copy=True)
        pass
    # Extract values.
    #matrix_signal = numpy.transpose(numpy.copy(table_signal.to_numpy()))
    matrix_signal = numpy.copy(table_signal.to_numpy())
    # Extract minimal and maximal values.
    if (
        (not constrain_signal_values) and
        (value_minimum is None) and
        (value_maximum is None)
    ):
        round_offset = (numpy.nanmin(matrix_signal) * 0.10)
        value_minimum = round((numpy.nanmin(matrix_signal) - round_offset), 3)
        value_maximum = round((numpy.nanmax(matrix_signal) + round_offset), 3)
    # Fill missing values.
    if fill_missing:
        matrix_signal = numpy.nan_to_num(
            matrix_signal,
            copy=True,
            nan=value_missing_fill,
            posinf=value_maximum, # or + 1.0 for correlations
            neginf=value_minimum, # or - 1.0 for correlations
        )
        pass
    # Constrain values.
    if constrain_signal_values:
        matrix_signal[matrix_signal < value_minimum] = value_minimum
        matrix_signal[matrix_signal > value_maximum] = value_maximum
        pass
    # Determine labels for axes.
    if (transpose_table):
        if (
            (format_table == 2) and
            (column_group is not None) and
            (
                (labels_abscissa_categories is None) or
                (len(labels_abscissa_categories) < 2)
            )
        ):
            labels_abscissa_categories = labels_index_rows # vertical_axis
            #labels_abscissa_categories = labels_group_rows # vertical_axis
            labels_abscissa_groups = labels_group_unique # vertical_axis

        elif (
            (
                (labels_abscissa_categories is None) or
                (len(labels_abscissa_categories) < 2)
            )
        ):
            labels_abscissa_categories = labels_index_rows # vertical axis
            labels_abscissa_groups = None
        if (
            (labels_ordinate_categories is None) or
            (len(labels_ordinate_categories) < 2)
        ):
            labels_ordinate_categories = labels_index_columns # vertical axis
            labels_ordinate_groups = None
    else:
        if (
            (labels_abscissa_categories is None) or
            (len(labels_abscissa_categories) < 2)
        ):
            labels_abscissa_categories = labels_index_columns # horizontal axis
            labels_abscissa_groups = None,
        if (
            (format_table == 2) and
            (column_group is not None) and
            (
                (labels_ordinate_categories is None) or
                (len(labels_ordinate_categories) < 2)
            )
        ):
            labels_ordinate_categories = labels_index_rows # vertical_axis
            #labels_ordinate_categories = labels_group_rows # vertical_axis
            labels_ordinate_groups = labels_group_unique # vertical_axis
        elif (
            (
                (labels_ordinate_categories is None) or
                (len(labels_ordinate_categories) < 2)
            )
        ):
            labels_ordinate_categories = labels_index_rows # vertical axis
            labels_ordinate_groups = None
        pass

    # Define discrete numerical representation of categorical groups.
    if ((format_table == 2) and (column_group is not None)):
        groups_indices = dict()
        indices_groups = dict()
        index = 0
        for name in labels_group_unique:
            groups_indices[name] = index
            indices_groups[index] = name
            index += 1
            pass
        integers_group_rows = list(map(
            lambda name: groups_indices[name], labels_group_rows
        ))
        #groups_representation = dict()
        #groups_representation["names"] = labels_group_rows
        #groups_representation["integers"] = integers_group_rows
        #table_groups_representation = pandas.DataFrame(
        #    data=groups_representation,
        #)
        # Organize the integer representations of discrete categorical groups
        # as a matrix.
        matrix_group_integers = numpy.array(integers_group_rows).reshape(
            1, len(integers_group_rows)
        )
    else:
        groups_indices = None
        indices_groups = None
        integers_group_rows = None
        matrix_group_integers = None
        pass

    ##########
    # 9. Collect information.
    pail = dict()
    pail["matrix_signal"] = matrix_signal
    pail["matrix_group_integers"] = matrix_group_integers
    pail["value_minimum"] = value_minimum
    pail["value_maximum"] = value_maximum
    pail["labels_index_columns"] = labels_index_columns
    pail["labels_index_rows"] = labels_index_rows
    pail["labels_group_rows"] = labels_group_rows
    pail["integers_group_rows"] = integers_group_rows
    pail["labels_group_unique"] = labels_group_unique
    pail["groups_indices"] = groups_indices
    pail["indices_groups"] = indices_groups
    pail["groups_counts"] = groups_counts
    pail["labels_ordinate_categories"] = labels_ordinate_categories
    pail["labels_ordinate_groups"] = labels_ordinate_groups
    pail["labels_abscissa_categories"] = labels_abscissa_categories
    pail["labels_abscissa_groups"] = labels_abscissa_groups

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.plot.py")
        function = str(
            "extract_prepare_table_signals_categories_for_heatmap" +
            "()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=4)
        print("source table of signals:")
        print("rows in table: " + str(table_signal.shape[0]))
        print("columns in table: " + str(table_signal.shape[1]))
        putly.print_terminal_partition(level=4)
        print("matrix of signals:")
        count_rows = copy.deepcopy(matrix_signal.shape[0])
        count_columns = copy.deepcopy(matrix_signal.shape[1])
        print("matrix rows (dimension 0): " + str(count_rows))
        print("matrix columns (dimension 1): " + str(count_columns))
        putly.print_terminal_partition(level=4)
        print("abscissa, horizontal axis labels:")
        #print(labels_abscissa_categories)
        print("ordinate, vertical axis labels:")
        #print(labels_ordinate_categories)
        putly.print_terminal_partition(level=4)
        print("labels_group_rows:")
        print(labels_group_rows)
        print("integers_group_rows:")
        print(integers_group_rows)
        putly.print_terminal_partition(level=4)

    # Return information.
    return pail


def determine_size_axis_labels_categories(
    count_labels=None,
):
    """
    Determine appropriate size of font for text labels on an axis to represent
    explicit names of features or categories.

    The maximal count of labels is 200.

    These sizes of font correspond to definitions in the function below.
    package: partner
    module: plot.py
    function: define_font_properties()

    Review: 24 October 2024

    arguments:
        count_labels (int): count of labels for axis

    raises:

    returns:
        (str): text designation of size for font in labels

    """

    # Determine appropriate size of font for text labels of explicit names of
    # features or categories along axis.
    if (150 <= count_labels and count_labels < 201):
        size_label = "nineteen"
    elif (120 <= count_labels and count_labels < 150):
        size_label = "eighteen"
    elif (100 <= count_labels and count_labels < 120):
        size_label = "seventeen"
    elif (80 <= count_labels and count_labels < 100):
        size_label = "sixteen"
    elif (60 <= count_labels and count_labels < 80):
        size_label = "fifteen"
    elif (40 <= count_labels and count_labels < 60):
        size_label = "fourteen"
    elif (30 <= count_labels and count_labels < 40):
        size_label = "thirteen"
    elif (20 <= count_labels and count_labels < 30):
        size_label = "twelve"
    elif (10 <= count_labels and count_labels < 20):
        size_label = "eleven"
    elif (1 <= count_labels and count_labels < 10):
        size_label = "ten"
    else:
        size_label = None
        pass
    # Return information.
    return size_label


def plot_heatmap_signal_features_observations_labels(
    table=None,
    format_table=None,
    index_columns=None,
    index_rows=None,
    transpose_table=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_maximum=None,
    title_ordinate=None,
    title_abscissa=None,
    title_bar=None,
    labels_ordinate_categories=None,
    labels_abscissa_categories=None,
    size_title_ordinate=None,
    size_title_abscissa=None,
    size_title_bar=None,
    size_label_ordinate=None,
    size_label_abscissa=None,
    size_label_bar=None,
    show_labels_ordinate=None,
    show_labels_abscissa=None,
    show_scale_bar=None,
    aspect=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Heatmap.

    features of this chart design...
    labels of categorical groups on both axes: True
    labels of significance on individual cells: False
    clustering: False

    Format of source table

    Format of source table is in wide format with floating-point values of
    signal intensities or a single, specific type of descriptive statistics
    (usually either mean or median) corresponding to features across rows and
    observations or groups of observations across columns.
    ----------
    observation_or_group group_1 group_2 group_3 group_4
    feature
    feature_1            0.01    0.001   0.001   0.015
    feature_2            0.01    0.001   0.001   0.015
    feature_3            -0.01   0.001   0.001   0.015
    feature_4            -0.01   0.001   0.001   0.015
    feature_5            -0.01   0.001   0.001   0.015
    ----------

    For versatility, the source table does not have explicitly defined indices
    across columns or rows.

    This function preserves the original sequence of features. This function
    also preserves the original sequence of groups and observations within
    groups.

    MatPlotLib color maps.
    https://matplotlib.org/stable/tutorials/colors/colormaps.html

    Review: 30 December 2024

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            for features in rows across sample observations or groups of
            sample observations in columns
        format_table (int): value 1 for features across rows and observations
            or groups of observations across columns, value 2 for features
            across columns and observations across rows with potential special
            column for groups of observations
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        transpose_table (bool): whether to transpose matrix from table
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_missing_fill (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
        title_ordinate (str): title for ordinate vertical axis
        title_abscissa (str): title for abscissa horizontal axis
        title_bar (str): title for scale bar
        labels_ordinate_categories (list<str>): optional, explicit labels for
            ordinate or vertical axis
        labels_abscissa_categories (list<str>): optional, explicit labels for
            abscissa or horizontal axis
        size_title_ordinate (str): font size
        size_title_abscissa (str): font size
        size_title_bar (str): font size
        size_label_ordinate (str): font size
        size_label_abscissa (str): font size
        size_label_bar (str): font size
        show_labels_ordinate (bool): whether to show on vertical ordinate axis
            of plot chart explicit text labels for individual categories
        show_labels_abscissa (bool): whether to show on horizontal abscissa
            axis of plot chart explicit text labels for individual categories
        show_scale_bar (bool): whether to create scale bar
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): MatPlotLib figure object

    """

    ##########
    # Prepare information for figure.
    pail = extract_prepare_table_signals_categories_for_heatmap(
        table=table,
        format_table=format_table, # 1: features in rows, observations in columns
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=None,
        transpose_table=transpose_table,
        fill_missing=fill_missing,
        value_missing_fill=value_missing_fill,
        constrain_signal_values=constrain_signal_values,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        labels_ordinate_categories=labels_ordinate_categories,
        labels_abscissa_categories=labels_abscissa_categories,
        report=report,
    )

    ##########
    # Create figure.
    figure = pplot.initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Create axes.
    #axes_main = matplotlib.pyplot.axes()
    axes_main = figure.add_subplot(111)
    # Keep axes, ticks, and labels, but remove border.
    for position in ['right', 'top', 'bottom', 'left']:
        matplotlib.pyplot.gca().spines[position].set_visible(False)
    # Adjust margins.
    figure.subplots_adjust(
        left=0.02,
        right=0.99,
        top=0.98,
        bottom=0.15,
    )

    # Plot values as a grid of color on continuous scale.
    # This function represents values acros matrix dimension 0 as vertical
    # rows.
    # This function represents values across matrix dimension 1 as horizontal
    # columns.
    # Diverging color maps: "PRGn", "PRGn_r", "PiYG", "PiYG_r",
    # Diverging color maps: "PuOr", "PuOr_r",
    # Diverging color maps: "PuOr", "PuOr_r", "RdBu", "RdBu_r", "BrBG",
    # Sequential color maps: "Reds", "Reds_r", "Oranges", "Oranges_r",
    # site: https://montoliu.naukas.com/2021/11/18/color-blindness-purple-and-
    #     orange-are-the-solution/
    image_main = axes_main.imshow(
        pail["matrix_signal"],
        cmap=matplotlib.colormaps["PuOr"], # binary, Reds, RdBu_r, PuOr, PuOr_r
        vmin=pail["value_minimum"],
        vmax=pail["value_maximum"],
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set titles for axes.
    if (len(title_ordinate) > 0):
        axes_main.set_ylabel(
            ylabel=title_ordinate,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_ordinate]
        )
        pass
    if (len(title_abscissa) > 0):
        axes_main.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
        pass
    # Set tick parameters for axes.
    axes_main.tick_params(
        axis="both", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=5.0, # 5.0
        width=3.5, # 3.5
        pad=10, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=False,
        bottom=False,
        left=False,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=False,
    )
    # Manage labels for categories along the axes.
    # Determine feasibility and appropriate font size for representing labels
    # of individual categorical features or observations along the axes.
    count_labels_ordinate = len(pail["labels_ordinate_categories"])
    count_labels_abscissa = len(pail["labels_abscissa_categories"])
    if (size_label_ordinate is None):
        size_label_ordinate = determine_size_axis_labels_categories(
            count_labels=count_labels_ordinate,
        )
        pass
    if (size_label_abscissa is None):
        size_label_abscissa = determine_size_axis_labels_categories(
            count_labels=count_labels_abscissa,
        )
        pass
    # Determine whether to show labels for features or categories along the
    # vertical ordinate axis.
    if (
        (show_labels_ordinate) and
        (size_label_ordinate is not None) and
        (pail["labels_ordinate_categories"] is not None) and
        (count_labels_ordinate > 1)
    ):
        # Set tick positions and labels on vertical ordinate axis.
        axes_main.set_yticks(
            numpy.arange(pail["matrix_signal"].shape[0]),
        )
        axes_main.set_yticklabels(
            pail["labels_ordinate_categories"],
            #minor=False,
            ha="right", # horizontal alignment
            va="center", # vertical alignment
            alpha=1.0,
            rotation=0.0,
            rotation_mode="anchor",
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_label_ordinate]
        )
        # Set tick parameters for vertical ordinate axis.
        axes_main.tick_params(
            axis="y", # "y", "x", or "both"
            which="both", # "major", "minor", or "both"
            length=5.0, # 5.0
            width=3.5, # 3.5
            pad=10, # 7.5
            direction="out",
            color=colors["black"],
            labelcolor=colors["black"],
            left=True,
            right=False,
            labelleft=True,
            labelright=False,
        )
        pass
    # Determine whether to show labels for features or categories along the
    # horizontal abscissa axis.
    if (
        (show_labels_abscissa) and
        (size_label_abscissa is not None) and
        (pail["labels_abscissa_categories"] is not None) and
        (count_labels_abscissa > 1)
    ):
        # Set tick positions and labels on horizontal abscissa axis.
        axes_main.set_xticks(
            numpy.arange(pail["matrix_signal"].shape[1]),
        )
        axes_main.set_xticklabels(
            pail["labels_abscissa_categories"],
            #minor=False,
            ha="left", # horizontal alignment
            va="top", # vertical alignment
            alpha=1.0,
            rotation=-60,
            rotation_mode="anchor",
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_label_abscissa]
        )
        # Set tick parameters for horizontal abscissa axis.
        axes_main.tick_params(
            axis="x", # "y", "x", or "both"
            which="both", # "major", "minor", or "both"
            length=5.0, # 5.0
            width=3.5, # 3.5
            pad=17, # 7.5 - 17.5
            direction="out",
            color=colors["black"],
            labelcolor=colors["black"],
            top=False,
            bottom=True,
            labeltop=False,
            labelbottom=True,
        )
        pass

    # Create legend for scale of color grid.
    if show_scale_bar:
        bar = axes_main.figure.colorbar(
            image_main,
            orientation="vertical",
            ax=axes_main,
            location="right",
            shrink=0.5, # 0.7; factor for dimensions of the Scale Bar.
        )
        if (len(title_bar) > 0):
            bar.ax.set_ylabel(
                title_bar,
                rotation=-90,
                va="bottom",
                labelpad=5, # 5
                alpha=1.0,
                backgroundcolor=colors["white"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_title_bar],
            )
        bar.ax.tick_params(
            axis="both",
            which="both", # major, minor, or both
            direction="out",
            length=7.5, # 5.0, 7.5
            width=3, # 2.5, 5.0
            color=colors["black"],
            pad=5, # 5, 7
            labelsize=fonts["values"][size_label_bar]["size"],
            labelcolor=colors["black"],
        )

    # Return figure.
    return figure


def plot_heatmap_signal_features_sets_observations_labels(
    table_signal=None,
    table_feature_sets=None,
    format_table_signal=None,
    index_columns=None,
    index_rows=None,
    transpose_table=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_maximum=None,
    title_ordinate=None,
    title_abscissa=None,
    title_bar=None,
    labels_abscissa_categories=None,
    size_title_ordinate=None,
    size_title_abscissa=None,
    size_title_bar=None,
    size_label_feature_set=None,
    size_label_abscissa=None,
    size_label_bar=None,
    show_labels_abscissa=None,
    show_scale_bar=None,
    aspect=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Heat map.

    Notice that this chart design does not show explicit categorical labels for
    features on the vertical ordinate axis or for observations or groups of
    observations on the horizontal abscissa axis.

    Format of source table of signals
    Format of source table is in wide format with floating-point values of
    signal intensities or a single, specific type of descriptive statistics
    (usually either mean or median) corresponding to features across rows and
    observations or groups of observations across columns.
    ----------
    observation_or_group group_1 group_2 group_3 group_4
    feature
    feature_1            0.01    0.001   0.001   0.015
    feature_2            0.01    0.001   0.001   0.015
    feature_3            -0.01   0.001   0.001   0.015
    feature_4            -0.01   0.001   0.001   0.015
    feature_5            -0.01   0.001   0.001   0.015
    ----------

    Format of source table of feature allocations to sets
    ----------
    feature     set_1 set_2 set_3 set_4 set_5
    feature_1   1     0     0     0     1
    feature_2   1     1     0     0     0
    feature_3   0     1     1     0     0
    feature_4   0     0     1     1     0
    feature_5   0     0     0     1     1
    ----------

    For versatility, the source tables do not have explicitly defined indices
    across columns or rows.

    This function preserves the original sequence of features. This function
    also preserves the original sequence of groups and observations within
    groups.

    This function assumes that the table has many more observations than
    features, and for this reason, the design orients features across the
    vertical axis and observations across the horizontal axis. With a landscape
    aspect ratio, this design allows more space for the horizontal axis.

    ----------

    Reference:

    Review: 30 December 2024

    arguments:
        table_signal (object): Pandas data-frame table of values of signal
            intensity for features in columns across sample observations or
            groups of sample observations in rows
        table_feature_sets (object): Pandas data-frame table of indications of
            allocation of features to sets in a sort sequence that matches the
            sequence of features across columns in table of signals
        format_table_signal (int): value 1 for features across rows and
            observations or groups of observations across columns, value 2 for
            features across columns and observations across rows with potential
            special column for groups of observations
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        transpose_table (bool): whether to transpose matrix from table
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_missing_fill (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
        title_ordinate (str): title for ordinate vertical axis
        title_abscissa (str): title for abscissa horizontal axis
        title_bar (str): title for scale bar
        labels_abscissa_categories (list<str>): optional, explicit labels for
            abscissa or horizontal axis
        size_title_ordinate (str): font size
        size_title_abscissa (str): font size
        size_title_bar (str): font size
        size_label_feature_set (str): font size
        size_label_abscissa (str): font size
        size_label_bar (str): font size
        show_labels_abscissa (bool): whether to show on horizontal abscissa
            axis of plot chart explicit text labels for individual categories
        show_scale_bar (bool): whether to create scale bar
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): MatPlotLib figure object

    """

    ##########
    # Prepare information for figure.
    pail = extract_prepare_table_signals_categories_for_heatmap(
        table=table_signal,
        format_table=format_table_signal, # 1: features in rows, observations in columns
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=None,
        transpose_table=transpose_table,
        fill_missing=fill_missing,
        value_missing_fill=value_missing_fill,
        constrain_signal_values=constrain_signal_values,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        labels_ordinate_categories=None,
        labels_abscissa_categories=labels_abscissa_categories,
        report=report,
    )
    # Copy information in table.
    table_feature_sets = table_feature_sets.copy(deep=True)
    # Organize indices in table.
    table_feature_sets.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_feature_sets.set_index(
        index_rows,
        append=False,
        drop=True,
        inplace=True
    )
    table_feature_sets.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    # Extract categorical names from original indices in table.
    labels_feature_sets_index_columns = copy.deepcopy(
        table_feature_sets.columns.get_level_values(
            index_columns
        ).to_list()
    )
    labels_feature_sets_index_rows = copy.deepcopy(
        table_feature_sets.index.get_level_values(
            index_rows
        ).to_list()
    )
    # Extract values.
    #matrix_signal = numpy.transpose(numpy.copy(table_signal.to_numpy()))
    matrix_feature_sets = numpy.copy(table_feature_sets.to_numpy())

    ##########
    # Create figure.

    ##########
    # Initialize figure.
    figure = pplot.initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Initialize grid within figure.
    #matplotlib.gridspec.GridSpec()
    #figure.add_gridspec()
    #sharey=True, # sharey="row"
    #sharex=True, # sharex="col",
    grid = matplotlib.gridspec.GridSpec(
        nrows=1,
        ncols=4,
        wspace=0.005, # horizontal width space between grid blocks for subplots
        hspace=0.005, # vertical height space between grid blocks for subplots
        width_ratios=(40,55,1,4),
        height_ratios=(100,),
    )
    grid.update(
        wspace=0.005, # horizontal width space between grid blocks for subplots
        hspace=0.005, # vertical height space between grid blocks for subplots
    )
    # Initialize axes within grid within figure.
    axes_set = figure.add_subplot(grid[0,0]) # first row, first column
    axes_main = figure.add_subplot(grid[0,1]) # first row, second column
    axes_space = figure.add_subplot(grid[0,2]) # first row, third column
    axes_bar = figure.add_subplot(grid[0,3]) # first row, fourt column
    axes_set.clear()
    axes_main.clear()
    axes_space.clear()
    axes_bar.clear()
    # Set axes to empty as a space holder.
    axes_space.axis("off")
    # Keep axes, ticks, and labels, but remove border.
    for position in ['right', 'top', 'bottom', 'left']:
        matplotlib.pyplot.gca().spines[position].set_visible(False)
    # Adjust margins.
    # Method "tight_layout()" does not function properly with object
    # "GridSpec()".
    #grid.tight_layout(
    #    figure,
    #    #pad=1.0,
    #    #h_pad=1.0,
    #    #w_pad=1.0,
    #    rect=[0,0.05,1.0,1.0], # left, bottom, right, top
    #)
    #grid.update(
    #    wspace=0.005, # horizontal width space between grid blocks for subplots
    #    hspace=0.005, # vertical height space between grid blocks for subplots
    #)
    figure.subplots_adjust(
        left=0.01,
        right=0.90,
        top=0.99,
        bottom=0.20,
    )

    ##########
    # axes: main
    # Plot values as a grid of color on continuous scale.
    # This function represents values acros matrix dimension 0 as vertical
    # rows.
    # This function represents values across matrix dimension 1 as horizontal
    # columns.
    # Diverging color maps: "PRGn", "PRGn_r", "PiYG", "PiYG_r",
    # Diverging color maps: "PuOr", "PuOr_r",
    # Diverging color maps: "PuOr", "PuOr_r", "RdBu", "RdBu_r", "BrBG",
    # Sequential color maps: "Reds", "Reds_r", "Oranges", "Oranges_r",
    # site: https://montoliu.naukas.com/2021/11/18/color-blindness-purple-and-
    #     orange-are-the-solution/
    image_main = axes_main.imshow(
        pail["matrix_signal"],
        cmap=matplotlib.colormaps["PuOr"], # binary, Reds, RdBu_r, PuOr, PuOr_r
        vmin=pail["value_minimum"],
        vmax=pail["value_maximum"],
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set titles for axes.
    if (len(title_ordinate) > 0):
        axes_main.set_ylabel(
            ylabel=title_ordinate,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_ordinate]
        )
        pass
    if (len(title_abscissa) > 0):
        axes_main.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
        pass
    # Set tick parameters for axes.
    axes_main.tick_params(
        axis="both", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=5.0, # 5.0
        width=3.5, # 3.5
        pad=10, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=False,
        bottom=False,
        left=False,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=False,
    )
    # Manage labels for categories along the axes.
    # Determine feasibility and appropriate font size for representing labels
    # of individual categorical features or observations along the axes.
    count_labels_abscissa = len(pail["labels_abscissa_categories"])
    if (size_label_abscissa is None):
        size_label_abscissa = determine_size_axis_labels_categories(
            count_labels=count_labels_abscissa,
        )
        pass
    # Determine whether to show labels for features or categories along the
    # horizontal abscissa axis.
    if (
        (show_labels_abscissa) and
        (size_label_abscissa is not None) and
        (pail["labels_abscissa_categories"] is not None) and
        (count_labels_abscissa > 1)
    ):
        # Set tick positions and labels on horizontal abscissa axis.
        axes_main.set_xticks(
            numpy.arange(pail["matrix_signal"].shape[1]),
        )
        axes_main.set_xticklabels(
            pail["labels_abscissa_categories"],
            #minor=False,
            ha="left", # horizontal alignment
            va="top", # vertical alignment
            alpha=1.0,
            rotation=-60,
            rotation_mode="anchor",
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_label_abscissa]
        )
        # Set tick parameters for horizontal abscissa axis.
        axes_main.tick_params(
            axis="x", # "y", "x", or "both"
            which="both", # "major", "minor", or "both"
            length=5.0, # 5.0
            width=3.5, # 3.5
            pad=10, # 7.5
            direction="out",
            color=colors["black"],
            labelcolor=colors["black"],
            top=False,
            bottom=True,
            labeltop=False,
            labelbottom=True,
        )
        pass

    ##########
    # axes: set
    # Define color map for discrete, binary integer representation of
    # allocation to sets.
    # https://matplotlib.org/3.1.1/gallery/color/named_colors.html
    # "white", "black"
    # "white", "dimgray"
    color_map_set = matplotlib.colors.ListedColormap([
        "white", "dimgray"
    ])
    # Plot values as a grid of color on discrete scale.
    image = axes_set.imshow(
        matrix_feature_sets,
        cmap=color_map_set,
        vmin=0,
        vmax=1,
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set tick positions and labels on axes.
    axes_set.set_xticks(
        numpy.arange(matrix_feature_sets.shape[1]),
    )
    axes_set.set_xticklabels(
        labels_feature_sets_index_columns,
        #minor=False,
        ha="center", # horizontal alignment
        va="top", # vertical alignment
        alpha=1.0,
        rotation=90, # negative: clockwise; positive: count-clockwise
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"][size_label_feature_set],
    )
    # Set tick parameters for axes.
    axes_set.tick_params(
        axis="y", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        left=False,
        right=False,
        labelleft=False,
        labelright=False,
    )
    axes_set.tick_params(
        axis="x", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=5.0, # 5.0
        width=3.5, # 3.5
        pad=10, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=False,
        bottom=True,
        labeltop=False,
        labelbottom=True,
    )

    ##########
    # axes: bar
    # Create legend for scale of color grid.
    # Notice that use of the "cax" argument causes to ignore the "shrink"
    # argument.
    if show_scale_bar:
        # Create scale bar.
        bar_main = axes_bar.figure.colorbar(
            image_main,
            orientation="vertical",
            cax=axes_bar,
            location="right",
            #shrink=0.5, # 0.7; factor for dimensions of the Scale Bar.
        )
        if (len(title_bar) > 0):
            bar_main.ax.set_ylabel(
                title_bar,
                rotation=-90,
                loc="center",
                va="center",
                labelpad=15, # 5
                alpha=1.0,
                backgroundcolor=colors["white"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_title_bar],
            )
        bar_main.ax.tick_params(
            axis="both",
            which="both", # major, minor, or both
            direction="out",
            length=7.5, # 5.0, 7.5
            width=3, # 2.5, 5.0
            color=colors["black"],
            pad=3, # 5, 7
            labelsize=fonts["values"][size_label_bar]["size"],
            labelcolor=colors["black"],
        )
        pass

    ##########
    # Return figure.
    return figure


def plot_heatmap_signal_features_sets_observations_groups(
    table_signal=None,
    table_feature_sets=None,
    format_table_signal=None,
    index_columns=None,
    index_rows=None,
    column_group=None,
    transpose_table=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_maximum=None,
    title_ordinate=None,
    title_abscissa=None,
    title_bar=None,
    size_title_ordinate=None,
    size_title_abscissa=None,
    size_title_bar=None,
    size_label_feature_set=None,
    size_label_legend_observation_group=None,
    size_label_bar=None,
    show_scale_bar=None,
    aspect=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Heat map.

    Notice that this chart design does not show explicit categorical labels for
    features on the vertical ordinate axis or for observations or groups of
    observations on the horizontal abscissa axis.

    Format of source table of signals
    Format of source table is in wide format with floating-point values of
    signal intensities corresponding to features across columns and distinct
    observations across rows. A special column gives identifiers corresponding
    to each observation across rows. Another special column provides names
    of categorical groups of observations.
    ----------
    observation     group   feature_1 feature_2 feature_3 feature_4 feature_5
    observation_1   group_1 0.001     0.001     0.001     0.001     0.001
    observation_2   group_1 0.001     0.001     0.001     0.001     0.001
    observation_3   group_2 0.001     0.001     0.001     0.001     0.001
    observation_4   group_2 0.001     0.001     0.001     0.001     0.001
    observation_5   group_3 0.001     0.001     0.001     0.001     0.001
    ----------

    Format of source table of feature allocations to sets
    ----------
    feature     set_1 set_2 set_3 set_4 set_5
    feature_1   1     0     0     0     1
    feature_2   1     1     0     0     0
    feature_3   0     1     1     0     0
    feature_4   0     0     1     1     0
    feature_5   0     0     0     1     1
    ----------

    For versatility, the source tables do not have explicitly defined indices
    across columns or rows.

    This function preserves the original sequence of features. This function
    also preserves the original sequence of groups and observations within
    groups.

    This function assumes that the table has many more observations than
    features, and for this reason, the design orients features across the
    vertical axis and observations across the horizontal axis. With a landscape
    aspect ratio, this design allows more space for the horizontal axis.

    ----------

    Reference:
    https://www.kaggle.com/code/sgalella/correlation-heatmaps-with-
        hierarchical-clustering
    https://www.kaggle.com/code/chinoysen/beginners-guide-to-heatmaps-cluster-
        heatmap
    https://seaborn.pydata.org/generated/seaborn.clustermap.html
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.heatmap.html
    https://gist.github.com/peterk87/5505691
    https://matplotlib.org/stable/tutorials/colors/colormaps.html

    Maybe useful for problem-solving:
    https://matplotlib.org/stable/gallery/subplots_axes_and_figures/
        subplots_demo.html
        - subplots and gridspec
    http://www.acgeospatial.co.uk/colour-bar-for-discrete-rasters-with-matplotlib/
        - color maps for discrete values
    https://matplotlib.org/stable/users/explain/colors/colorbar_only.html
        #colorbar-with-arbitrary-colors
        - color bar with custom dimensions of discrete, categorical intervals

    https://matplotlib.org/stable/gallery/color/colormap_reference.html
        - color maps in MatPlotLib

    https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
        - categorical color bars
        - helpful reference for figuring out how to set the thresholds between
          the discrete levels of the color bar to correspond properly to the
          columns of the main heatmap...
    https://stackoverflow.com/questions/9707676/defining-a-discrete-colormap-for-imshow
    https://stackoverflow.com/questions/7229971/2d-grid-data-visualization-in-python
    https://matplotlib.org/stable/users/explain/colors/colormapnorms.html
        - normalization of information for color maps???

    Review: 30 December 2024

    arguments:
        table_signal (object): Pandas data-frame table of values of signal
            intensity for features in columns across sample observations or
            groups of sample observations in rows
        table_feature_sets (object): Pandas data-frame table of indications of
            allocation of features to sets in a sort sequence that matches the
            sequence of features across columns in table of signals
        format_table_signal (int): value 1 for features across rows and
            observations or groups of observations across columns, value 2 for
            features across columns and observations across rows with potential
            special column for groups of observations
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        column_group (str): name of column in table to use for groups
        transpose_table (bool): whether to transpose matrix from table
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_missing_fill (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
        title_ordinate (str): title for ordinate vertical axis
        title_abscissa (str): title for abscissa horizontal axis
        title_bar (str): title for scale bar
        size_title_ordinate (str): font size
        size_title_abscissa (str): font size
        size_title_bar (str): font size
        size_label_feature_set (str): font size
        size_label_legend_observation_group (str): font size
        size_label_bar (str): font size
        show_scale_bar (bool): whether to create scale bar
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): MatPlotLib figure object

    """

    ##########
    # Prepare information for figure.
    pail = extract_prepare_table_signals_categories_for_heatmap(
        table=table_signal,
        format_table=format_table_signal, # 1: features in rows, observations in columns
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=column_group,
        transpose_table=transpose_table,
        fill_missing=fill_missing,
        value_missing_fill=value_missing_fill,
        constrain_signal_values=constrain_signal_values,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        report=report,
    )
    # Copy information in table.
    table_feature_sets = table_feature_sets.copy(deep=True)
    # Organize indices in table.
    table_feature_sets.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_feature_sets.set_index(
        index_columns,
        append=False,
        drop=True,
        inplace=True
    )
    table_feature_sets.columns.rename(
        "feature_sets",
        inplace=True,
    ) # single-dimensional index
    # Extract categorical names from original indices in table.
    labels_feature_sets_index_columns = copy.deepcopy(
        table_feature_sets.columns.get_level_values(
            "feature_sets"
        ).to_list()
    )
    labels_feature_sets_index_rows = copy.deepcopy(
        table_feature_sets.index.get_level_values(
            index_columns
        ).to_list()
    )
    # Extract values.
    #matrix_signal = numpy.transpose(numpy.copy(table_signal.to_numpy()))
    matrix_feature_sets = numpy.copy(table_feature_sets.to_numpy())

    ##########
    # Create figure.

    ##########
    # Initialize figure.
    figure = pplot.initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Initialize grid within figure.
    #matplotlib.gridspec.GridSpec()
    #figure.add_gridspec()
    #sharey=True, # sharey="row"
    #sharex=True, # sharex="col",
    grid = matplotlib.gridspec.GridSpec(
        nrows=4,
        ncols=2,
        wspace=0.005, # horizontal width space between grid blocks for subplots
        hspace=0.005, # vertical height space between grid blocks for subplots
        width_ratios=(15,85), # first column 1/10th width of second column
        height_ratios=(90,3,5,2), # first row 30 times the height of second row
    )
    grid.update(
        wspace=0.005, # horizontal width space between grid blocks for subplots
        hspace=0.005, # vertical height space between grid blocks for subplots
    )
    # Initialize axes within grid within figure.
    axes_set = figure.add_subplot(grid[0,0]) # first row, first column
    axes_main = figure.add_subplot(grid[0,1]) # first row, second column
    axes_space_first = figure.add_subplot(grid[1,0]) # second row, first column
    axes_group = figure.add_subplot(grid[1,1]) # second row, second column
    axes_space_second = figure.add_subplot(grid[2,1]) # third row, second column
    axes_bar = figure.add_subplot(grid[3,1]) # fourth row, second column
    axes_set.clear()
    axes_main.clear()
    axes_space_first.clear()
    axes_group.clear()
    axes_space_second.clear()
    axes_bar.clear()
    # Set axes to empty as a space holder.
    axes_space_first.axis("off")
    axes_space_second.axis("off")
    #axes_set.axis("off")
    #axes_group.axis("off")
    #axes_bar.axis("off")
    #axes_main.axis("off")
    # Keep axes, ticks, and labels, but remove border.
    for position in ['right', 'top', 'bottom', 'left']:
        matplotlib.pyplot.gca().spines[position].set_visible(False)
    # Adjust margins.
    # Method "tight_layout()" does not function properly with object
    # "GridSpec()".
    #grid.tight_layout(
    #    figure,
    #    #pad=1.0,
    #    #h_pad=1.0,
    #    #w_pad=1.0,
    #    rect=[0,0.05,1.0,1.0], # left, bottom, right, top
    #)
    #grid.update(
    #    wspace=0.005, # horizontal width space between grid blocks for subplots
    #    hspace=0.005, # vertical height space between grid blocks for subplots
    #)
    figure.subplots_adjust(
        left=0.02,
        right=0.98,
        top=0.98,
        bottom=0.05,
    )

    ##########
    # axes: main
    # Plot values as a grid of color on continuous scale.
    # This function represents values acros matrix dimension 0 as vertical
    # rows.
    # This function represents values across matrix dimension 1 as horizontal
    # columns.
    # Diverging color maps: "PRGn", "PRGn_r", "PiYG", "PiYG_r",
    # Diverging color maps: "PuOr", "PuOr_r",
    # Diverging color maps: "PuOr", "PuOr_r", "RdBu", "RdBu_r", "BrBG",
    # Sequential color maps: "Reds", "Reds_r", "Oranges", "Oranges_r",
    # site: https://montoliu.naukas.com/2021/11/18/color-blindness-purple-and-
    #     orange-are-the-solution/
    image_main = axes_main.imshow(
        pail["matrix_signal"],
        cmap=matplotlib.colormaps["PuOr"], # binary, Reds, RdBu_r, PuOr, PuOr_r
        vmin=pail["value_minimum"],
        vmax=pail["value_maximum"],
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set titles for axes.
    if (len(title_ordinate) > 0):
        axes_main.set_ylabel(
            ylabel=title_ordinate,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_ordinate]
        )
        pass
    if (len(title_abscissa) > 0):
        axes_main.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
        pass
    # Manage tick gradations, labels, and titles on axes.
    axes_main.tick_params(
        axis="both", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        left=False,
        right=False,
        labelleft=False,
        labelright=False,
    )

    ##########
    # axes: set
    # Define color map for discrete, binary integer representation of
    # allocation to sets.
    # https://matplotlib.org/3.1.1/gallery/color/named_colors.html
    # "white", "black"
    # "white", "dimgray"
    color_map_set = matplotlib.colors.ListedColormap([
        "white", "dimgray"
    ])
    # Plot values as a grid of color on discrete scale.
    image_set = axes_set.imshow(
        matrix_feature_sets,
        cmap=color_map_set,
        vmin=0,
        vmax=1,
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set tick positions and labels on axes.
    axes_set.set_xticks(
        numpy.arange(matrix_feature_sets.shape[1]),
    )
    axes_set.set_xticklabels(
        labels_feature_sets_index_columns,
        #minor=False,
        ha="center", # horizontal alignment
        va="top", # vertical alignment
        alpha=1.0,
        rotation=90, # negative: clockwise; positive: count-clockwise
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"][size_label_feature_set],
    )
    # Set tick parameters for axes.
    axes_set.tick_params(
        axis="y", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        left=False,
        right=False,
        labelleft=False,
        labelright=False,
    )
    axes_set.tick_params(
        axis="x", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=5.0, # 5.0
        width=3.5, # 3.5
        pad=10, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=False,
        bottom=True,
        labeltop=False,
        labelbottom=True,
    )

    ##########
    # axes: group
    # Define color map for discrete, integer representation of categorical
    # groups.
    #discrete_minimum = 0
    #discrete_maximum = len(pail["labels_group_unique"])
    discrete_minimum = numpy.nanmin(pail["matrix_group_integers"])
    discrete_maximum = numpy.nanmax(pail["matrix_group_integers"])
    color_map_group = matplotlib.pyplot.get_cmap(
        "tab10", # "Set1", "Set2", "Dark2", "tab10",
        ((discrete_maximum - discrete_minimum) + 1)
    )
    # Plot values as a grid of color on discrete scale.
    image_group = axes_group.imshow(
        pail["matrix_group_integers"],
        cmap=color_map_group,
        vmin=discrete_minimum,
        vmax=discrete_maximum,
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(
        #    0.0,
        #    (pail["matrix_group_integers"].shape[1]),
        #    0.0,
        #    (3*(pail["matrix_group_integers"].shape[0]))
        #),
    )
    axes_group.tick_params(
        top=False,
        bottom=False,
        left=False,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=False,
    )
    # Create legend.
    # See nested function "create_legend_elements()" within main function
    # "plot_scatter_factor_groups()".
    handles = [matplotlib.patches.Patch(
        color=color_map_group(i),
        label=pail["indices_groups"][i]
    ) for i in pail["indices_groups"].keys()]
    axes_group.legend(
        handles=handles,
        prop=fonts["properties"][size_label_legend_observation_group],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=4,
    )

    ##########
    # axes: bar
    # Create legend for scale of color grid.
    # Notice that use of the "cax" argument causes to ignore the "shrink"
    # argument.
    if show_scale_bar:
        # Set definitive scale to avoid conflicts with global objects.
        #color_map_main = matplotlib.pyplot.cm.ScalarMappable(
        #    cmap="PuOr",
        #    norm=matplotlib.pyplot.Normalize(
        #        vmin=pail["value_minimum"],
        #        vmax=pail["value_maximum"],
        #    ),
        #)
        #color_map_main.set_array(pail["matrix_signal"])
        # Create scale bar.
        bar_main = axes_bar.figure.colorbar(
            image_main, # image_main, color_map_main
            orientation="horizontal",
            cax=axes_bar,
            location="bottom",
            #shrink=0.5, # 0.7; factor for dimensions of the Scale Bar.
        )
        if (len(title_bar) > 0):
            bar_main.ax.set_xlabel(
                title_bar,
                rotation=0,
                loc="center",
                va="bottom",
                labelpad=20, # 5
                alpha=1.0,
                backgroundcolor=colors["white"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_title_bar],
            )
        bar_main.ax.tick_params(
            axis="both",
            which="both", # major, minor, or both
            direction="out",
            length=7.5, # 5.0, 7.5
            width=3, # 2.5, 5.0
            color=colors["black"],
            pad=3, # 5, 7
            labelsize=fonts["values"][size_label_bar]["size"],
            labelcolor=colors["black"],
        )
        pass

    ##########
    # Return figure.
    return figure


# Manage and write plot charts.

# TODO: TCW; 26 September 2025
# Create a heatmap that has labeled features while having individual (not mean)
# values for observations.


def plot_heatmap_features_sets_observations_groups(
    table_signal=None,
    table_feature=None,
    index_columns=None,
    index_rows=None,
    column_group=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 27 December 2024

    arguments:
        table_signal (object): Pandas data-frame table of floating-point values
            of a signal corresponding features in columns across observations
            in rows
        table_feature (object): Pandas data-frame table of indications of
            allocation of genes to sets in a sort sequence that matches the
            sequence of genes across columns in table of signals
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        column_group (str): name of column in table to use for groups
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table_signal = table_signal.copy(deep=True)
    table_signal_extract = table_signal.copy(deep=True)
    table_feature = table_feature.copy(deep=True)
    # Organize indices in table.
    table_signal_extract.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_signal_extract.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    table_signal_extract.set_index(
        [index_rows, column_group,],
        append=False,
        drop=True,
        inplace=True
    )
    # Extract minimal and maximal values of signal intensity.
    matrix = numpy.copy(table_signal_extract.to_numpy())
    #value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
    #value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)
    round_offset = (numpy.nanmin(matrix) * 0.10)
    value_minimum = round((numpy.nanmin(matrix) - round_offset), 3)
    value_maximum = round((numpy.nanmax(matrix) + round_offset), 3)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = plot_heatmap_signal_features_sets_observations_groups(
        table_signal=table_signal,
        table_feature_sets=table_feature,
        format_table_signal=2, # 2: features in columns; observations, groups in rows
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=column_group,
        transpose_table=True,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        title_ordinate="",
        title_abscissa="",
        title_bar="individual signal (z-score)",
        size_title_ordinate="ten",
        size_title_abscissa="ten",
        size_title_bar="thirteen",
        size_label_feature_set="fifteen",
        size_label_legend_observation_group="fourteen",
        size_label_bar="fifteen",
        show_scale_bar=True,
        aspect="portrait",
        fonts=fonts,
        colors=colors,
        report=report,
    )
    # Return information.
    return figure


def plot_heatmap_features_sets_observations_labels(
    table_signal=None,
    table_feature=None,
    index_columns=None,
    index_rows=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 26 December 2024

    arguments:
        table_signal (object): Pandas data-frame table of floating-point values
            of a signal corresponding features in columns across observations
            in rows
        table_feature (object): Pandas data-frame table of indications of
            allocation of genes to sets in a sort sequence that matches the
            sequence of genes across columns in table of signals
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table_signal = table_signal.copy(deep=True)
    table_signal_extract = table_signal.copy(deep=True)
    table_feature = table_feature.copy(deep=True)
    # Organize indices in table.
    table_signal_extract.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_signal_extract.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    table_signal_extract.set_index(
        [index_rows,],
        append=False,
        drop=True,
        inplace=True
    )
    # Extract minimal and maximal values of signal intensity.
    matrix = numpy.copy(table_signal_extract.to_numpy())
    #value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
    #value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)
    round_offset = (numpy.nanmin(matrix) * 0.10)
    value_minimum = round((numpy.nanmin(matrix) - round_offset), 3)
    value_maximum = round((numpy.nanmax(matrix) + round_offset), 3)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = plot_heatmap_signal_features_sets_observations_labels(
        table_signal=table_signal,
        table_feature_sets=table_feature,
        format_table_signal=1, # 1: features in rows, observations or groups in columns
        index_columns=index_columns,
        index_rows=index_rows,
        transpose_table=False,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        title_ordinate="",
        title_abscissa="",
        title_bar="Mean Signal (z-score)",
        labels_abscissa_categories=None,
        size_title_ordinate="eight",
        size_title_abscissa="eight",
        size_title_bar="nine",
        size_label_feature_set="eleven",
        size_label_abscissa="nine",
        size_label_bar="eleven",
        show_labels_abscissa=True,
        show_scale_bar=True, # whether to show scale bar on individual figures
        aspect="square", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )
    # Return information.
    return figure


def plot_heatmap_features_observations_labels(
    table=None,
    index_columns=None,
    index_rows=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review: TCW; 27 December 2024

    arguments:
        table (object): Pandas data-frame table of floating-point values of a
            single, specific type of descriptive statistics (usually either
            mean or median) corresponding to groups of observations across
            columns and features across rows
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table = table.copy(deep=True)
    table_extract = table.copy(deep=True)
    # Organize indices in table.
    table_extract.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_extract.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    table_extract.set_index(
        [index_rows],
        append=False,
        drop=True,
        inplace=True
    )
    # Extract minimal and maximal values of signal intensity.
    matrix = numpy.copy(table_extract.to_numpy())
    #value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
    #value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)
    round_offset = (numpy.nanmin(matrix) * 0.10)
    value_minimum = round((numpy.nanmin(matrix) - round_offset), 3)
    value_maximum = round((numpy.nanmax(matrix) + round_offset), 3)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = plot_heatmap_signal_features_observations_labels(
        table=table,
        format_table=1, # 1: features in rows, observations or groups in columns
        index_columns=index_columns,
        index_rows=index_rows,
        transpose_table=False,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        title_ordinate="",
        title_abscissa="",
        title_bar="Mean Signal (z-score)",
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        size_title_ordinate="eight",
        size_title_abscissa="eight",
        size_title_bar="ten",
        size_label_ordinate="eighteen", # determine automatically if "None"; "fifteen"
        size_label_abscissa="eleven", # determine automatically if "None"
        size_label_bar="twelve",
        show_labels_ordinate=True,
        show_labels_abscissa=True,
        show_scale_bar=True,
        aspect="portrait", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )
    # Return information.
    return figure


def manage_create_write_plot_charts(
    table_heatmap_individual_1=None,
    table_heatmap_individual_2=None,
    table_heatmap_individual_3=None,
    table_heatmap_mean_set=None,
    table_heatmap_mean_label=None,
    table_allocation_1=None,
    table_allocation_2=None,
    table_allocation_3=None,
    table_allocation_4=None,
    index_features=None,
    index_observations=None,
    heatmap_individual=None,
    heatmap_mean=None,
    path_directory_parent=None,
    report=None,
):
    """
    Plot chart representations of values of signal intensity for features
    across sample observations or groups of sample observations.

    arguments:
        table_heatmap_individual_1 (object): Pandas data-frame table of values
            of signal intensity for features across columns and sample
            observations in groups across rows
        table_heatmap_individual_2 (object): Pandas data-frame table of values
            of signal intensity for features across columns and sample
            observations in groups across rows
        table_heatmap_individual_3 (object): Pandas data-frame table of values
            of signal intensity for features across columns and sample
            observations in groups across rows
        table_heatmap_mean_set (object): Pandas data-frame table of
            descriptive statistics for values of signal intensity for features
            across rows and groups of sample observations across columns
        table_heatmap_mean_label (object): Pandas data-frame table of
            descriptive statistics for values of signal intensity for features
            across rows and groups of sample observations across columns
        table_allocation_1 (object): Pandas data-frame table of indications of
            allocation of genes to sets in a sort sequence that matches the
            sequence of genes across columns in table
            'table_heatmap_individual_1' and the sequence of genes across rows
            in table 'table_heatmap_mean'
        table_allocation_2 (object): Pandas data-frame table of indications of
            allocation of genes to sets in a sort sequence that matches the
            sequence of genes across columns in table
            'table_heatmap_individual_2'
        table_allocation_3 (object): Pandas data-frame table of indications of
            allocation of genes to sets in a sort sequence that matches the
            sequence of genes across columns in table
            'table_heatmap_individual_3'
        table_allocation_4 (object): Pandas data-frame table of indications of
            allocation of genes to sets in a sort sequence that matches the
            sequence of genes across columns in table
            'table_heatmap_mean_set'
        index_features (str): name for index corresponding to features
        index_observations (str): name for index corresponding to observations
        heatmap_individual (bool): whether to create heatmap chart for
            individual values of signal intensity
        heatmap_mean (bool): whether to create heatmap chart for means of
            signal intensity
        path_directory_parent (str): path to directory for procedure's product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:


    """

    # Copy information.
    table_heatmap_individual_1 = table_heatmap_individual_1.copy(deep=True)
    table_heatmap_individual_2 = table_heatmap_individual_2.copy(deep=True)
    table_heatmap_individual_3 = table_heatmap_individual_3.copy(deep=True)
    table_heatmap_mean_set = table_heatmap_mean_set.copy(deep=True)
    table_heatmap_mean_label = table_heatmap_mean_label.copy(deep=True)
    table_allocation_1 = table_allocation_1.copy(deep=True)
    table_allocation_2 = table_allocation_2.copy(deep=True)
    table_allocation_3 = table_allocation_3.copy(deep=True)
    table_allocation_4 = table_allocation_4.copy(deep=True)

    ##########
    # Heatmap Individual
    if heatmap_individual:
        figure_heatmap_individual_1 = (
            plot_heatmap_features_sets_observations_groups(
                table_signal=table_heatmap_individual_1,
                table_feature=table_allocation_1,
                index_columns=index_features,
                index_rows=index_observations,
                column_group="group",
                report=report,
        ))
        figure_heatmap_individual_2 = (
            plot_heatmap_features_sets_observations_groups(
                table_signal=table_heatmap_individual_2,
                table_feature=table_allocation_2,
                index_columns=index_features,
                index_rows=index_observations,
                column_group="group",
                report=report,
        ))
        figure_heatmap_individual_3 = (
            plot_heatmap_features_sets_observations_groups(
                table_signal=table_heatmap_individual_3,
                table_feature=table_allocation_3,
                index_columns=index_features,
                index_rows=index_observations,
                column_group="group",
                report=report,
        ))
    else:
        figure_heatmap_individual_1 = None
        figure_heatmap_individual_2 = None
        figure_heatmap_individual_3 = None
        pass

    ##########
    # Heatmap Mean
    if heatmap_mean:
        # Create heatmaps.
        figure_heatmap_mean_set = (
            plot_heatmap_features_sets_observations_labels(
                table_signal=table_heatmap_mean_set,
                table_feature=table_allocation_4,
                index_columns="group_observations",
                index_rows=index_features,
                report=False,
        ))
        figure_heatmap_mean_label = (
            plot_heatmap_features_observations_labels(
                table=table_heatmap_mean_label,
                index_columns="group_observations",
                index_rows=index_features,
                report=False,
        ))
    else:
        figure_heatmap_mean_set = None
        figure_heatmap_mean_label = None
        pass

    ##########
    # Bundle information.
    pail_write_plot_heatmap = dict()
    pail_write_plot_heatmap["heatmap_individual_1"] = (
        figure_heatmap_individual_1
    )
    pail_write_plot_heatmap["heatmap_individual_2"] = (
        figure_heatmap_individual_2
    )
    pail_write_plot_heatmap["heatmap_individual_3"] = (
        figure_heatmap_individual_3
    )
    pail_write_plot_heatmap["heatmap_mean_set"] = (
        figure_heatmap_mean_set
    )
    pail_write_plot_heatmap["heatmap_mean_label"] = (
        figure_heatmap_mean_label
    )

    ##########
    # Write product information to file.

    # Define paths to directories.
    path_directory_heatmap = os.path.join(
        path_directory_parent, "chart_heatmap",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_heatmap,
    )

    # Write figures to file.
    pplot.write_product_plots_parent_directory(
        pail_write=pail_write_plot_heatmap,
        format="jpg", # jpg, png, svg
        resolution=150,
        path_directory=path_directory_heatmap,
    )

    ##########
    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_groups_observations_sets_features.py"
        )
        print(str("module: " + module))
        print("function: manage_create_write_plot_charts()")
        putly.print_terminal_partition(level=5)
        print("files written from bundle pail:")
        print(pail_write_plot_heatmap.keys())
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    pass


################################################################################
# Procedure


##########
# Call main procedure.


def execute_procedure(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_observations=None,
    path_file_source_table_features=None,
    path_file_source_table_signals=None,
    path_file_source_table_groups_observations=None,
    path_file_source_table_sets_features=None,
    path_file_source_list_features_selection=None,
    path_file_source_list_features_cluster_one=None,
    path_file_source_list_features_cluster_two=None,
    column_identifier_observation=None,
    column_identifier_feature=None,
    column_identifier_signal=None,
    column_name_feature=None,
    transpose_table_signals=None,
    allow_replicate_observations=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review: TCW; 24 September 2025

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_file_source_table_observations (str): path to source file
        path_file_source_table_features (str): path to source file
        path_file_source_table_signals (str): path to source file
        path_file_source_table_groups_observations (str): path to source file
        path_file_source_table_sets_features (str): path to source file
        path_file_source_list_features_selection (str): path to source file
        path_file_source_list_features_cluster_one (str): path to source file
        path_file_source_list_features_cluster_two (str): path to source file
        column_identifier_observation (str): name of column in source table
        column_identifier_feature (str): name of column in source table
        column_identifier_signal (str): name of column in source table
        column_name_feature (str): name of column in source table
        transpose_table_signals (bool): whether to transpose the table of
            signals
        allow_replicate_observations (bool): whether to allow replicate
            observations or to require groups to be mutually exclusive, such
            that any individual observation can only belong to one group
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
        path_file_source_table_observations=(
            path_file_source_table_observations
        ),
        path_file_source_table_features=path_file_source_table_features,
        path_file_source_table_signals=path_file_source_table_signals,
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        path_file_source_table_sets_features=(
            path_file_source_table_sets_features
        ),
        path_file_source_list_features_selection=(
            path_file_source_list_features_selection
        ),
        path_file_source_list_features_cluster_one=(
            path_file_source_list_features_cluster_one
        ),
        path_file_source_list_features_cluster_two=(
            path_file_source_list_features_cluster_two
        ),
        column_identifier_observation=column_identifier_observation,
        column_identifier_feature=column_identifier_feature,
        column_identifier_signal=column_identifier_signal,
        column_name_feature=column_name_feature,
        transpose_table_signals=transpose_table_signals,
        allow_replicate_observations=allow_replicate_observations,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: plot_chart_volcano.py")
        print("function: execute_procedure()")
        print("system: local")
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        path_file_source_table_observations=(
            pail_parameters["path_file_source_table_observations"]
        ),
        path_file_source_table_features=(
            pail_parameters["path_file_source_table_features"]
        ),
        path_file_source_table_signals=(
            pail_parameters["path_file_source_table_signals"]
        ),
        path_file_source_table_groups_observations=(
            pail_parameters["path_file_source_table_groups_observations"]
        ),
        path_file_source_table_sets_features=(
            pail_parameters["path_file_source_table_sets_features"]
        ),
        path_file_source_list_features_selection=(
            pail_parameters["path_file_source_list_features_selection"]
        ),

        path_file_source_list_features_cluster_one=(
            pail_parameters["path_file_source_list_features_cluster_one"]
        ),
        path_file_source_list_features_cluster_two=(
            pail_parameters["path_file_source_list_features_cluster_two"]
        ),
        report=pail_parameters["report"],
    )

    # Parameters.
    pail_groups = organize_parameters_groups_observations(
        table=pail_source["table_groups"],
        column_name="abbreviation",
        report=pail_parameters["report"],
    )
    #pail_groups["table"]
    #pail_groups["names_groups_observations_sequence"]
    #pail_groups["categories_groups"]
    #pail_groups["records"]

    table_features = pail_source["table_features"]
    column_identifier_feature = pail_parameters["column_identifier_feature"]
    features_available = copy.deepcopy(
        table_features[column_identifier_feature].unique().tolist()
    )
    pail_sets = read_organize_parameters_sets_features(
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        table=pail_source["table_sets"],
        column_name="abbreviation",
        features_available=features_available,
        features_selection=pail_source["features_selection"],
        features_cluster_one=pail_source["features_cluster_one"],
        features_cluster_two=pail_source["features_cluster_two"],
        report=pail_parameters["report"],
    )
    #pail_sets["features_available"]
    #pail_sets["features_selection"]
    #pail_sets["features_cluster_one"]
    #pail_sets["features_cluster_two"]
    #pail_sets["table"]
    #pail_sets["names_sets_features_sequence"]
    #pail_sets["sets_features"]
    #pail_sets["features_sets_union"]
    #pail_sets["records"]

    # Prepare tables.
    pail_tables = manage_prepare_tables(
        table_observations=pail_source["table_observations"],
        table_features=pail_source["table_features"],
        table_signals=pail_source["table_signals"],
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_identifier_feature=(
            pail_parameters["column_identifier_feature"]
        ),
        column_identifier_signal=(
            pail_parameters["column_identifier_signal"]
        ),
        column_name_feature=(
            pail_parameters["column_name_feature"]
        ),
        instances_groups_observations=pail_groups["records"],
        key_name="abbreviation",
        names_groups_observations_sequence=(
            pail_groups["names_groups_observations_sequence"]
        ),
        categories_groups=pail_groups["categories_groups"],
        features_selection=pail_sets["features_selection"],
        features_cluster_one=pail_sets["features_cluster_one"],
        features_cluster_two=pail_sets["features_cluster_two"],
        features_sets_union=pail_sets["features_sets_union"],
        sets_features=pail_sets["sets_features"],
        names_sets_features_sequence=pail_sets["names_sets_features_sequence"],
        transpose_table_signals=(
            pail_parameters["transpose_table_signals"]
        ),
        allow_replicate_observations=(
            pail_parameters["allow_replicate_observations"]
        ),
        report=pail_parameters["report"],
    )

    # Create and write plot charts.
    manage_create_write_plot_charts(
        table_heatmap_individual_1=pail_tables["table_3"], # cluster with constraint by sets of features and by groups of observations
        table_heatmap_individual_2=pail_tables["table_4"], # cluster with constraint by groups of observations only
        table_heatmap_individual_3=pail_tables["table_5"], # cluster without constraint
        table_heatmap_mean_set=pail_tables["table_8"],
        table_heatmap_mean_label=pail_tables["table_8_translation"],
        table_allocation_1=pail_tables["table_allocation_3"],
        table_allocation_2=pail_tables["table_allocation_4"],
        table_allocation_3=pail_tables["table_allocation_5"],
        table_allocation_4=pail_tables["table_allocation_8"],
        index_features=pail_parameters["column_identifier_feature"],
        index_observations=pail_parameters["column_identifier_signal"],
        heatmap_individual=True,
        heatmap_mean=True,
        path_directory_parent=pail_parameters["path_directory_product"],
        report=report,
    )



    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_source = sys.argv[1]
    path_directory_product = sys.argv[2]
    path_directory_dock = sys.argv[3]
    path_file_source_table_observations = sys.argv[4]
    path_file_source_table_features = sys.argv[5]
    path_file_source_table_signals = sys.argv[6]
    path_file_source_table_groups_observations = sys.argv[7]
    path_file_source_table_sets_features = sys.argv[8]
    path_file_source_list_features_selection = sys.argv[9]
    path_file_source_list_features_cluster_one = sys.argv[10]
    path_file_source_list_features_cluster_two = sys.argv[11]
    column_identifier_observation = sys.argv[12]
    column_identifier_feature = sys.argv[13]
    column_identifier_signal = sys.argv[14]
    column_name_feature = sys.argv[15]
    transpose_table_signals = sys.argv[16]
    allow_replicate_observations = sys.argv[17]
    report = sys.argv[18]

    # Call function for procedure.
    execute_procedure(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_observations=(
            path_file_source_table_observations
        ),
        path_file_source_table_features=path_file_source_table_features,
        path_file_source_table_signals=path_file_source_table_signals,
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        path_file_source_table_sets_features=(
            path_file_source_table_sets_features
        ),
        path_file_source_list_features_selection=(
            path_file_source_list_features_selection
        ),
        path_file_source_list_features_cluster_one=(
            path_file_source_list_features_cluster_one
        ),
        path_file_source_list_features_cluster_two=(
            path_file_source_list_features_cluster_two
        ),
        column_identifier_observation=column_identifier_observation,
        column_identifier_feature=column_identifier_feature,
        column_identifier_signal=column_identifier_signal,
        column_name_feature=column_name_feature,
        transpose_table_signals=transpose_table_signals,
        allow_replicate_observations=allow_replicate_observations,
        report=report,
    )

    pass



#
