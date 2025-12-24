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
# Date, initialization: 26 October 2025
# Date, review or revision: 3 December 2025
# Date, review or revision: 5 November 2025
# Date, review or revision: 26 October 2025
################################################################################
# Note


##########
# Note:

# Note: TCW; 13 November 2025
# It might become practical to implement separate thresholds for the first and
# second sets of features, respectively.

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
#import matplotlib_venn
#import seaborn
#import sklearn

# Custom
import partner.utility as putly
#import partner.parallelization as prall
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
import partner.decomposition as pdcmp
import partner.plot as pplot

import utility_special as sutly
import plot_special as splot

#dir()
#importlib.reload()

###############################################################################
# Functionality


# Organize raw parameters.


def parse_text_parameters(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_features_observations=None,
    path_file_source_list_features_first=None,
    path_file_source_list_features_second=None,
    path_file_source_table_groups_observations=None,
    column_identifier_observation=None,
    column_name_observation=None,
    proportion_nonmissing_observations=None,
    type_correlation=None,
    match_pairs=None,
    prefix_match_first=None,
    suffix_match_first=None,
    prefix_match_second=None,
    suffix_match_second=None,
    intersect_features=None,
    threshold_filter_first=None,
    threshold_filter_second=None,
    threshold_type=None,
    threshold_value=None,
    threshold_proportion=None,
    sort_match_pairs_diagonal=None,
    cluster_features_first=None,
    cluster_features_second=None,
    sort_other_features=None,
    plot_scale_minimum=None,
    plot_scale_center=None,
    plot_scale_maximum=None,
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
        path_file_source_table_features_observations (str): path to source file
        path_file_source_list_features_first (str): path to source file
        path_file_source_list_features_second (str): path to source file
        path_file_source_table_groups_observations (str): path to source file
        column_identifier_observation (str): name of column in source table
        column_name_observation (str): name of column in source table
        proportion_nonmissing_observations (float): threshold by proportion of
            observations that must have nonmissing values for both features in
            each pair
        type_correlation (str): type of correlation for calculation; either
            'pearson' or 'spearman'
        ...
        cluster_features (str): how to apply the cluster operation, either
            'none', 'both', 'first', or 'second'
        sort_other_features (bool): whether to sort the set of features to
            which the cluster operation was not applied
        plot_scale_minimum (float): minimal value to set as threshold for
            representation in the plot chart
        plot_scale_maximum (float): maximal value to set as threshold for
            representation in the plot chart
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
    pail["path_file_source_table_features_observations"] = str(
        path_file_source_table_features_observations
    ).strip()
    pail["path_file_source_list_features_first"] = str(
        path_file_source_list_features_first
    ).strip()
    pail["path_file_source_list_features_second"] = str(
        path_file_source_list_features_second
    ).strip()
    pail["path_file_source_table_groups_observations"] = str(
        path_file_source_table_groups_observations
    ).strip()

    # Names of columns.
    pail["column_identifier_observation"] = str(
        column_identifier_observation
    ).strip()
    pail["column_name_observation"] = str(
        column_name_observation
    ).strip()

    # Number.
    pail["proportion_nonmissing_observations"] = float(
        str(proportion_nonmissing_observations).strip()
    )
    pail["threshold_value"] = float(
        str(threshold_value).strip()
    )
    pail["threshold_proportion"] = float(
        str(threshold_proportion).strip()
    )
    if (
        (len(str(plot_scale_minimum)) > 0) and
        (str(plot_scale_minimum) != "none")
    ):
        pail["plot_scale_minimum"] = float(
            str(plot_scale_minimum).strip()
        )
    else:
        pail["plot_scale_minimum"] = None
    if (
        (len(str(plot_scale_center).strip()) > 0) and
        (str(plot_scale_center).strip() != "none")
    ):
        pail["plot_scale_center"] = float(
            str(plot_scale_center).strip()
        )
    else:
        pail["plot_scale_center"] = None
    if (
        (len(str(plot_scale_maximum).strip()) > 0) and
        (str(plot_scale_maximum).strip() != "none")
    ):
        pail["plot_scale_maximum"] = float(
            str(plot_scale_maximum).strip()
        )
    else:
        pail["plot_scale_maximum"] = None
        pass

    # Names and categories.
    # It is problematic to pass any white space in parameters from a script in
    # Bash. Designate the hash symbol "#" as a substitute for white space.
    # It is also problematic to pass an empty string in parameters from a
    # script in Bash. Designate the word "none" as a substitute for missing or
    # empty.
    # Iterate on individual names that could be empty or missing.
    names_categories = {
        "type_correlation": type_correlation,
        "prefix_match_first": prefix_match_first,
        "suffix_match_first": suffix_match_first,
        "prefix_match_second": prefix_match_second,
        "suffix_match_second": suffix_match_second,
        "threshold_type": threshold_type,
    }
    for key_name in names_categories.keys():
        # Determine whether parameter has a valid value that is not none.
        if (
            (str(names_categories[key_name]).strip().lower() != "none")
        ):
            # Parse value.
            pail[key_name] = str(
                names_categories[key_name]
            ).strip().replace("#", " ")
        else:
            pail[key_name] = ""
            pass
        pass

    # Boolean, true or false.
    # Iterate on individual of Boolean designations.
    designations = {
        "match_pairs": match_pairs,
        "intersect_features": intersect_features,
        "threshold_filter_first": threshold_filter_first,
        "threshold_filter_second": threshold_filter_second,
        "sort_match_pairs_diagonal": sort_match_pairs_diagonal,
        "cluster_features_first": cluster_features_first,
        "cluster_features_second": cluster_features_second,
        "sort_other_features": sort_other_features,
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
            "calculate_tables_plot_charts_features_correlations.py"
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
    path_file_source_table_features_observations=None,
    path_file_source_list_features_first=None,
    path_file_source_list_features_second=None,
    path_file_source_table_groups_observations=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 26 October 2025

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
    # Tables.
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
    # Lists.
    pail["features"] = dict()
    pail["features"]["first"] = putly.read_file_text_list(
        path_file=path_file_source_list_features_first,
        delimiter="\n",
        unique=True,
    )
    pail["features"]["second"] = putly.read_file_text_list(
        path_file=path_file_source_list_features_second,
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
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def filter_organize_table_columns_rows(
    table=None,
    column_identifier_observation=None,
    column_name_observation=None,
    features_first=None,
    features_second=None,
    table_groups_observations=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 26 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    ##########
    # Copy information.
    table = table.copy(deep=True)
    table_groups_observations = table_groups_observations.copy(deep=True)
    features_first = copy.deepcopy(features_first)
    features_second = copy.deepcopy(features_second)

    ##########
    # Organize parameters for table's rows representing observations.

    pail_groups = sutly.organize_parameters_groups_observations(
        table=table_groups_observations,
        column_name="abbreviation",
        report=report,
    )
    #pail_groups["table"]
    #pail_groups["names_groups_observations_sequence"]
    #pail_groups["categories_groups"]
    #pail_groups["records"]

    pail_observations = sutly.organize_parameters_further_groups_observations(
        table_observations=table,
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        column_identifier_groups_observations=column_identifier_observation,
        instances_groups_observations=pail_groups["records"],
        key_name="abbreviation",
        names_groups_observations_sequence=(
            pail_groups["names_groups_observations_sequence"]
        ),
        report=report,
    )
    #pail_observations["table_observations_selection"]
    #pail_observations["observations_selection"]
    #pail_observations["translations_observations"]
    #pail_observations["names_groups_observations_sequence"]
    #pail_observations["groups_observations"]

    ##########
    # Filter columns and rows in table.

    # Organize preliminary parameters for table's columns representing
    # features.
    features_selection = list()
    features_selection.extend(features_first)
    features_selection.extend(features_second)
    # Copy information.
    categories_features_selection = copy.deepcopy(
        pail_groups["categories_groups"]
    )
    # Prepare inclusive list of columns.
    categories_features_selection.insert(0, column_name_observation)
    categories_features_selection.insert(0, column_identifier_observation)
    categories_features_selection.extend(features_selection)

    # Filter columns and rows in table for specific features and observations.
    table_filter = (
        porg.filter_select_table_columns_rows_by_identifiers(
            table=table,
            index_rows=column_identifier_observation,
            identifiers_columns=categories_features_selection,
            identifiers_rows=pail_observations["observations_selection"],
            report=False,
    ))

    # Filter rows in table for non-missing values across relevant columns.
    table_filter.dropna(
        axis="index",
        how="all",
        subset=features_selection,
        inplace=True,
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: filter_organize_table_columns_rows()")
        putly.print_terminal_partition(level=5)
        print("table after filters:")
        print(table_filter)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_filter


def filter_features_by_available_observations(
    table=None,
    column_identifier_observation=None,
    column_name_observation=None,
    features_first=None,
    features_second=None,
    proportion_nonmissing_observations=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 26 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    ##########
    # Copy information.
    table = table.copy(deep=True)
    features_first = copy.deepcopy(features_first)
    features_second = copy.deepcopy(features_second)

    # Extract identifiers of columns in table.
    features_available_1 = copy.deepcopy(
        table.columns.unique().tolist()
    )

    # Filter lists of features by the availability of corresponding columns in
    # table.
    features_first_filter_1 = list(filter(
        lambda item: item in features_available_1, features_first
    ))
    features_second_filter_1 = list(filter(
        lambda item: item in features_available_1, features_second
    ))
    features_first_filter_1 = putly.collect_unique_items(
        items=features_first_filter_1,
    )
    features_second_filter_1 = putly.collect_unique_items(
        items=features_second_filter_1,
    )
    features_selection_1 = list()
    features_selection_1.extend(features_first_filter_1)
    features_selection_1.extend(features_second_filter_1)

    # Filter table's columns corresponding to features by the availability of
    # observations across rows in the table.
    observations_selection = copy.deepcopy(
        table[column_identifier_observation].to_list()
    )
    table_filter = (
        porg.filter_table_columns_by_proportion_nonmissing_threshold(
            table=table,
            index_columns="features",
            index_rows=column_identifier_observation,
            columns_selection=features_selection_1,
            rows_selection=observations_selection,
            threshold_low=None,
            threshold_high=None,
            proportion=proportion_nonmissing_observations,
            report=report,
    ))

    # Extract identifiers of columns in table.
    features_available_2 = copy.deepcopy(
        table_filter.columns.unique().tolist()
    )

    # Filter lists of features by the availability of corresponding columns in
    # table.
    features_first_filter_2 = list(filter(
        lambda item: item in features_available_2, features_first
    ))
    features_second_filter_2 = list(filter(
        lambda item: item in features_available_2, features_second
    ))
    features_first_filter_2 = putly.collect_unique_items(
        items=features_first_filter_2,
    )
    features_second_filter_2 = putly.collect_unique_items(
        items=features_second_filter_2,
    )
    features_selection_2 = list()
    features_selection_2.extend(features_first_filter_1)
    features_selection_2.extend(features_second_filter_1)

    # Bundle information.
    pail = dict()
    pail["table"] = table_filter
    pail["features_first"] = features_first_filter_2
    pail["features_second"] = features_second_filter_2
    pail["features_selection"] = features_selection_2

    # Report.
    if report:
        # Organize information.
        count_first = len(features_first_filter_2)
        count_second = len(features_second_filter_2)
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: filter_features_by_available_observations()")
        putly.print_terminal_partition(level=5)
        print("table after filters:")
        print(table_filter)
        putly.print_terminal_partition(level=5)
        print("count of first features: " + str(count_first))
        print("count of second features: " + str(count_second))
        pass
    # Return information.
    return pail


def match_filter_features_pairs(
    features_first=None,
    features_second=None,
    prefix_match_first=None,
    suffix_match_first=None,
    prefix_match_second=None,
    suffix_match_second=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 11 November 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    ##########
    # Copy information.
    features_first = copy.deepcopy(features_first)
    features_second = copy.deepcopy(features_second)

    # Translate names of features to facilitate match.
    translations_first = list(map(
        lambda item: (
            str(item)
            .strip()
            .replace(prefix_match_first, "")
            .replace(suffix_match_first, "")
        ), copy.deepcopy(features_first)
    ))
    translations_second = list(map(
        lambda item: (
            str(item)
            .strip()
            .replace(prefix_match_second, "")
            .replace(suffix_match_second, "")
        ), copy.deepcopy(features_second)
    ))

    # Take union of translations.
    translations_union = putly.combine_sets_items_union_unique(
        sets_items=[
            translations_first,
            translations_second,
        ],
        report=report,
    )

    # Take intersection of translations.
    # This intersection will serve as a sort of consensus to filter both sets
    # of features.
    translations_intersection = putly.combine_sets_items_intersection_unique(
        items_first=translations_first,
        items_second=translations_second,
        report=report,
    )

    # Assemble references for pairs of features.
    # Base these references on the union of features in order to include all
    # possible pairs.
    # It will be necessary to check whether each pair actually exists.
    pairs_first_second = dict()
    pairs_second_first = dict()
    records_pairs = list()
    for item in translations_union:
        feature_first = str(prefix_match_first + item + suffix_match_first)
        feature_second = str(prefix_match_second + item + suffix_match_second)
        pairs_first_second[feature_first] = feature_second
        pairs_second_first[feature_second] = feature_first
        record = dict()
        record["first"] = feature_first
        record["second"] = feature_second
        records_pairs.append(record)
        pass

    #translations_first_second = list(map(
    #    lambda item: pairs_first_second[item], features_first_sequence
    #))
    #translations_second_first = list(map(
    #    lambda item: pairs_second_first[item], features_second_sequence
    #))

    # Filter identifiers of features by the intersection of matches between
    # first and second sets of features.
    features_first_intersection = list(filter(
        lambda item: ((
            str(item)
            .strip()
            .replace(prefix_match_first, "")
            .replace(suffix_match_first, "")
        ) in translations_intersection
        ), features_first
    ))
    features_second_intersection = list(filter(
        lambda item: ((
            str(item)
            .strip()
            .replace(prefix_match_second, "")
            .replace(suffix_match_second, "")
        ) in translations_intersection
        ), features_second
    ))


    # Bundle information.
    pail = dict()
    pail["features_first_intersection"] = features_first_intersection
    pail["features_second_intersection"] = features_second_intersection
    pail["translations_intersection"] = translations_intersection
    pail["translations_union"] = translations_union
    pail["pairs_first_second"] = pairs_first_second
    pail["pairs_second_first"] = pairs_second_first
    pail["records_pairs"] = records_pairs

    # Report.
    if report:
        # Organize information.
        count_first_original = len(features_first)
        count_second_original = len(features_second)
        count_first_novel = len(features_first_intersection)
        count_second_novel = len(features_second_intersection)
        count_intersection = len(translations_intersection)
        count_union = len(translations_union)
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: match_filter_features_pairs()")
        putly.print_terminal_partition(level=5)
        print("intersection of translations from first and second features:")
        print(translations_intersection)
        putly.print_terminal_partition(level=5)
        print("count of intersection: " + str(count_intersection))
        putly.print_terminal_partition(level=5)
        print("count of union: " + str(count_union))
        putly.print_terminal_partition(level=5)
        print("original, before intersection filter...")
        print(str(
            "count of first features: " + str(count_first_original)
        ))
        print(str(
            "count of second features: " + str(count_second_original)
        ))
        putly.print_terminal_partition(level=5)
        print("novel, after intersection filter...")
        print(str(
            "count of first features: " + str(count_first_novel)
        ))
        print(str(
            "count of second features: " + str(count_second_novel)
        ))
        pass
    # Return information.
    return pail

# TODO: TCW; 27 October 2025
# At some point, it might be convenient to reshape this table from long to wide
# so as to preserve the q-values. Refer to established functions for "stack" or
# "pivot" operations.

def calculate_correlations_populate_table_long(
    table=None,
    column_identifier_observation=None,
    features_first=None,
    features_second=None,
    proportion_nonmissing_observations=None,
    z_score=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 26 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    ##########
    # Copy information.
    table = table.copy(deep=True)
    features_first = copy.deepcopy(features_first)
    features_second = copy.deepcopy(features_second)

    # Determine whether to transform values of features to the z-score,
    # standard scale.
    if (z_score):
        columns_scale = list()
        columns_scale.extend(
            copy.deepcopy(features_first)
        )
        columns_scale.extend(
            copy.deepcopy(features_second)
        )
        table = pscl.transform_standard_z_score_by_table_columns(
                table=table,
                columns=columns_scale,
                report=report,
        )
        pass

    # Collect records for correlations across pairs of features.
    records = list()
    # Iterate on first list of features.
    for feature_first in features_first:
        # Iterate on second list of features.
        for feature_second in features_second:
            # Collect information.
            record = dict()
            record["feature_first"] = feature_first
            record["feature_second"] = feature_second
            # Copy information.
            table_pair = table.copy(deep=True)
            # Filter rows in table for non-missing values across relevant
            # columns.
            table_pair.dropna(
                axis="index",
                how="any",
                subset=[feature_first, feature_second,],
                inplace=True,
            )
            # Determine whether the pair of features has adequate availability
            # of values across observations.
            count_total = int(table.shape[0])
            count_pair = int(table_pair.shape[0])
            proportion_actual = float(count_pair / count_total)
            # Collect information.
            record["count_observations_total"] = count_total
            record["count_observations_pair"] = count_pair
            record["proportion_actual"] = proportion_actual
            record["threshold_proportion"] = proportion_nonmissing_observations
            if (
                (proportion_actual >= proportion_nonmissing_observations)
            ):
                # Calculate correlations.
                pail = pdesc.calculate_correlations_table_columns_pair(
                    table=table,
                    column_primary=feature_first,
                    column_secondary=feature_second,
                    count_minimum_observations=5,
                    report=False,
                )
                # Collect information.
                #record.update(pail)
                record["correlation_pearson"] = pail["correlation_pearson"]
                record["p_pearson"] = pail["p_pearson"]
                record["correlation_spearman"] = pail["correlation_spearman"]
                record["p_spearman"] = pail["p_spearman"]
                record["correlation_kendall"] = pail["correlation_kendall"]
                record["p_kendall"] = pail["p_kendall"]
            else:
                # Introduce missing values.
                measures = list()
                measures.append("correlation_pearson")
                measures.append("p_pearson")
                #measures.append("confidence_95_low_pearson")
                #measures.append("confidence_95_high_pearson")
                measures.append("correlation_spearman")
                measures.append("p_spearman")
                measures.append("correlation_kendall")
                measures.append("p_kendall")
                for measure in measures:
                    record[measure] = float("nan")
                    pass
                pass
            # Collect information.
            records.append(record)
            pass
        pass

    # Organize table.
    table_correlation = pandas.DataFrame(data=records)

    # Calculate Benjamini-Hochberg false discovery rate.

    # Calculate Benjamini-Hochberg q-values for False-Discovery Rate (FDR).
    # Calculate q-values across all comparisons in table.
    # FDR 5% (q <= 0.05).
    table_correlation = pdesc.calculate_table_false_discovery_rate_q_values(
        threshold=0.05, # alpha; family-wise error rate
        name_column_p_value="p_pearson",
        name_column_q_value="q_pearson",
        name_column_significance="significance_pearson",
        table=table_correlation,
    )
    table_correlation = pdesc.calculate_table_false_discovery_rate_q_values(
        threshold=0.05, # alpha; family-wise error rate
        name_column_p_value="p_spearman",
        name_column_q_value="q_spearman",
        name_column_significance="significance_spearman",
        table=table_correlation,
    )
    #table_correlation = pdesc.calculate_table_false_discovery_rate_q_values(
    #    threshold=0.05, # alpha; family-wise error rate
    #    name_column_p_value="p_kendall",
    #    name_column_q_value="q_kendall",
    #    name_column_significance="significance_kendall",
    #    table=table_correlation,
    #)

    # Filter and sort columns in table.
    columns_sequence = [
        "feature_first",
        "feature_second",
        "count_observations_total",
        "count_observations_pair",
        "proportion_actual",
        "threshold_proportion",
        "correlation_pearson", "p_pearson", "q_pearson",
        #"confidence_95_low_pearson", "confidence_95_high_pearson",
        "correlation_spearman", "p_spearman", "q_spearman",
        "correlation_kendall", "p_kendall", "q_kendall",
    ]
    table_correlation = porg.filter_sort_table_columns(
        table=table_correlation,
        columns_sequence=columns_sequence,
        report=report,
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: calculate_correlations_populate_table_long()")
        putly.print_terminal_partition(level=5)
        print("table of correlations:")
        print(table_correlation)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_correlation


def calculate_correlations_populate_table_wide(
    table=None,
    column_identifier_observation=None,
    features_first=None,
    features_second=None,
    proportion_nonmissing_observations=None,
    z_score=None,
    type_value=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 26 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    ##########
    # Copy information.
    table = table.copy(deep=True)
    features_first = copy.deepcopy(features_first)
    features_second = copy.deepcopy(features_second)

    # Determine whether to transform values of features to the z-score,
    # standard scale.
    if (z_score):
        columns_scale = list()
        columns_scale.extend(
            copy.deepcopy(features_first)
        )
        columns_scale.extend(
            copy.deepcopy(features_second)
        )
        table = pscl.transform_standard_z_score_by_table_columns(
                table=table,
                columns=columns_scale,
                report=report,
        )
        pass

    # Collect records for correlations across pairs of features.
    records = list()
    # Iterate on first list of features.
    for feature_first in features_first:
        # Collect information.
        record = dict()
        record["feature_first"] = feature_first
        # Iterate on second list of features.
        for feature_second in features_second:
            # Copy information.
            table_pair = table.copy(deep=True)
            # Filter rows in table for non-missing values across relevant
            # columns.
            table_pair.dropna(
                axis="index",
                how="any",
                subset=[feature_first, feature_second,],
                inplace=True,
            )
            # Determine whether the pair of features has adequate availability
            # of values across observations.
            count_total = int(table.shape[0])
            count_pair = int(table_pair.shape[0])
            proportion_actual = float(count_pair / count_total)
            if (
                (proportion_actual >= proportion_nonmissing_observations)
            ):
                # Calculate correlations.
                pail = pdesc.calculate_correlations_table_columns_pair(
                    table=table,
                    column_primary=feature_first,
                    column_secondary=feature_second,
                    count_minimum_observations=5,
                    report=False,
                )
                # Determine which value to collect.
                # Collect information.
                record[feature_second] = pail[type_value]
            else:
                # Introduce missing values.
                record[feature_second] = float("nan")
                pass
            pass
        # Collect information.
        records.append(record)
        pass

    # Organize table.
    table_correlation = pandas.DataFrame(data=records)

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: calculate_correlations_populate_table_wide()")
        putly.print_terminal_partition(level=5)
        print("table of correlations:")
        print(table_correlation)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_correlation


def optimize_threshold_table_wide_correlations(
    table=None,
    name_index_rows=None,
    name_index_columns=None,
    threshold_filter_rows=None,
    threshold_filter_columns=None,
    threshold_type=None,
    threshold_value=None,
    threshold_proportion=None,
    threshold_optimization=None,
    threshold_optimization_count=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 10 December 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)

    # Extract identifiers of columns and rows in table.
    columns_selection = copy.deepcopy(
        table.columns.unique().tolist()
    )
    columns_selection.remove(name_index_rows)
    rows_selection = copy.deepcopy(
        table[name_index_rows].to_list()
    )

    # Determine counts of original columns and rows in table.
    #count_columns = table.shape[1]
    #count_rows = table.shape[0]
    count_columns = int(len(columns_selection))
    count_rows = int(len(rows_selection))

    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table.columns.rename(
        name_index_columns,
        inplace=True,
    ) # single-dimensional index
    table.set_index(
        name_index_rows,
        append=False,
        drop=True,
        inplace=True
    )

    # Determine type of threshold and values.
    #if (threshold_type == "correlation"):

    # Notice
    # Parameter "threshold_proportion" designates the proportion of entities in
    # the dimension other than the dimension being filtered.
    # In contrast, parameter "threshold_optimization_count" designates the
    # optimal count of entities in the same dimension being filtered.

    # Optimize the threshold.
    if (
        (threshold_filter_columns) and
        (threshold_optimization) and
        (threshold_proportion is not None) and
        (threshold_proportion < 1.0) and
        (threshold_optimization_count is not None) and
        (threshold_optimization_count < count_columns)

    ):
        # Organize information.
        count_proportion = int(round(threshold_proportion * count_rows))
        # Collect information.
        candidates = list()
        # Iterate across columns.
        for index, series_column in table.items():
            # Extract array of values from series.
            # Array of values has length corresponding to count of rows.
            values_column = copy.deepcopy(series_column.to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            ))
            # Calculate absolute values.
            values_column = numpy.absolute(values_column)
            # Sort values in reverse sequence.
            values_column_sort = numpy.sort(
                values_column,
                axis=-1,
                kind="quicksort",
            )[::-1]
            # Determine candidate value of threshold.
            candidate = values_column_sort[count_proportion]
            # Collect information.
            candidates.append(candidate)
            pass
        # Convert list to array.
        threshold_candidates = numpy.array(
            candidates,
            dtype="float64",
            copy=True,
        )
        # Sort values in reverse sequence.
        threshold_candidates_sort = numpy.sort(
            threshold_candidates,
            axis=-1,
            kind="quicksort",
        )[::-1]
        # Determine threshold.
        threshold_selection = (
            threshold_candidates_sort[threshold_optimization_count]
        )

    elif (
        (threshold_filter_rows) and
        (threshold_optimization) and
        (threshold_proportion is not None) and
        (threshold_proportion < 1.0) and
        (threshold_optimization_count is not None) and
        (threshold_optimization_count < count_rows)
    ):
        # Organize information.
        count_proportion = int(round(threshold_proportion * count_columns))
        # Collect information.
        candidates = list()
        # Iterate across rows.
        for index, series_row in table.iterrows():
            # Extract array of values from series.
            # Array of values has length corresponding to count of columns.
            values_row = copy.deepcopy(series_row.to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            ))
            # Calculate absolute values.
            values_row = numpy.absolute(values_row)
            # Sort values in reverse sequence.
            values_row_sort = numpy.sort(
                values_row,
                axis=-1,
                kind="quicksort",
            )[::-1]
            # Determine candidate value of threshold.
            candidate = values_row_sort[count_proportion]
            # Collect information.
            candidates.append(candidate)
            pass
        # Convert list to array.
        threshold_candidates = numpy.array(
            candidates,
            dtype="float64",
            copy=True,
        )
        # Sort values in reverse sequence.
        threshold_candidates_sort = numpy.sort(
            threshold_candidates,
            axis=-1,
            kind="quicksort",
        )[::-1]
        # Determine threshold.
        threshold_selection = (
            threshold_candidates_sort[threshold_optimization_count]
        )

    else:
        threshold_selection = threshold_value

        pass

    # Bundle information.
    pail = dict()
    pail["threshold_value"] = threshold_selection

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: optimize_threshold_table_wide_correlations()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def filter_table_columns_match_rows(
    table=None,
    name_index_rows=None,
    name_index_columns=None,
    pairs_row_column=None,
    report=None,
):
    """
    Blank.

    Table orients first set of features across rows.
    Table orients second set of features across columns.

    By design this function should not filter or remove any columns or rows
    from the table.

    Review or revision: TCW; 10 December 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)
    pairs_row_column = copy.deepcopy(pairs_row_column)

    # Filter columns in table to match selection of rows in table, using a
    # prior reference of matching pairs.

    # Extract identifiers of columns and rows in table.
    columns_sequence_original = copy.deepcopy(
        table.columns.unique().tolist()
    )
    columns_sequence_original.remove(name_index_rows)
    rows_sequence = copy.deepcopy(
        table[name_index_rows].to_list()
    )

    # Determine sequence of columns to match rows.
    # Notice that reference must include all possible pairs.
    columns_sequence_novel = list(map(
        lambda item: pairs_row_column[item], rows_sequence
    ))
    # Filter columns by their availability in the table.
    columns_sequence_novel_available = list(filter(
        lambda item: (item in columns_sequence_original),
        columns_sequence_novel
    ))
    columns_sequence_novel_available.insert(0, name_index_rows,)
    columns_sequence_novel_available = putly.collect_unique_items(
        items=columns_sequence_novel_available,
    )
    # Filter and sort sequence of columns in table.
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=columns_sequence_novel_available,
        report=report,
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: filter_table_columns_match_rows()")
        putly.print_terminal_partition(level=5)
        pass

    # Return information.
    return table


def filter_table_rows_match_columns(
    table=None,
    name_index_rows=None,
    name_index_columns=None,
    pairs_column_row=None,
    report=None,
):
    """
    Blank.

    Table orients first set of features across rows.
    Table orients second set of features across columns.

    By design this function should not filter or remove any columns or rows
    from the table.

    Review or revision: TCW; 10 December 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)
    pairs_column_row = copy.deepcopy(pairs_column_row)

    # Filter rows in table to match selection of columns in table, using a
    # prior reference of matching pairs.

    # Extract identifiers of columns and rows in table.
    columns_sequence = copy.deepcopy(
        table.columns.unique().tolist()
    )
    columns_sequence.remove(name_index_rows)
    rows_sequence_original = copy.deepcopy(
        table[name_index_rows].to_list()
    )

    # Determine sequence of rows to match columns.
    # Notice that reference must include all possible pairs.
    rows_sequence_novel = list(map(
        lambda item: pairs_column_row[item], columns_sequence
    ))
    # Filter rows by their availability in the table.
    rows_sequence_novel_available = list(filter(
        lambda item: (item in rows_sequence_original),
        rows_sequence_novel
    ))
    rows_sequence_novel_available = putly.collect_unique_items(
        items=rows_sequence_novel_available,
    )
    # Filter rows in table.
    table = table.loc[
        table[name_index_rows].isin(rows_sequence_novel_available), :
    ].copy(deep=True)

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: filter_table_rows_match_columns()")
        putly.print_terminal_partition(level=5)
        pass

    # Return information.
    return table


def filter_table_wide_correlations_by_threshold(
    table_r=None,
    table_p=None,
    name_index_rows=None,
    name_index_columns=None,
    threshold_filter_rows=None,
    threshold_filter_columns=None,
    threshold_type=None,
    threshold_value=None,
    threshold_proportion=None,
    match_pairs=None,
    pairs_row_column=None,
    pairs_column_row=None,
    intersect_features=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 5 November 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table_r = table_r.copy(deep=True)
    table_r_filter = table_r.copy(deep=True)
    table_p = table_p.copy(deep=True)

    # Extract identifiers of columns and rows in table.
    columns_selection = copy.deepcopy(
        table_r.columns.unique().tolist()
    )
    columns_selection.remove(name_index_rows)
    rows_selection = copy.deepcopy(
        table_r[name_index_rows].to_list()
    )

    # Determine type of threshold and values.
    if (threshold_type == "correlation"):
        threshold_low = threshold_value
        threshold_high = None
    elif (threshold_type == "p_value"):
        threshold_low = None
        threshold_high = threshold_value
        pass

    # Filter columns and rows in table by their proportions of values that pass
    # thresholds.
    if (
        (threshold_filter_rows) and
        (threshold_type == "correlation")
    ):
        table_r_filter = (
            porg.filter_table_rows_by_proportion_nonmissing_threshold(
                table=table_r_filter,
                index_columns=name_index_columns,
                index_rows=name_index_rows,
                columns_selection=columns_selection,
                rows_selection=rows_selection,
                absolute_value=True,
                threshold_low=threshold_low,
                threshold_high=threshold_high,
                proportion=threshold_proportion,
                report=report,
        ))
        if (
            (match_pairs) and
            (intersect_features) and
            (pairs_row_column is not None)
        ):
            table_r_filter = filter_table_columns_match_rows(
                table=table_r_filter,
                name_index_rows=name_index_rows,
                name_index_columns=name_index_columns,
                pairs_row_column=pairs_row_column,
                report=report,
            )
            pass
        pass
    if (
        (threshold_filter_columns) and
        (threshold_type == "correlation")
    ):
        table_r_filter = (
            porg.filter_table_columns_by_proportion_nonmissing_threshold(
                table=table_r_filter,
                index_columns=name_index_columns,
                index_rows=name_index_rows,
                columns_selection=columns_selection,
                rows_selection=rows_selection,
                absolute_value=True,
                threshold_low=threshold_low,
                threshold_high=threshold_high,
                proportion=threshold_proportion,
                report=report,
        ))
        if (
            (match_pairs) and
            (intersect_features) and
            (pairs_column_row is not None)
        ):
            table_r_filter = filter_table_rows_match_columns(
                table=table_r_filter,
                name_index_rows=name_index_rows,
                name_index_columns=name_index_columns,
                pairs_column_row=pairs_column_row,
                report=report,
            )
            pass
        pass

    # Bundle information.
    pail = dict()
    pail["table_r_filter"] = table_r_filter

    # Report.
    if report:
        # Organize information.
        count_columns = table_r_filter.shape[1]
        count_rows = table_r_filter.shape[0]
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: filter_table_wide_correlations_by_threshold()")
        putly.print_terminal_partition(level=5)
        print("table after filters by threshold")
        print(table_r_filter)
        print("count columns: " + str(count_columns))
        print("count rows: " + str(count_rows))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def check_table_wide_correlations(
    table=None,
    count_threshold_rows=None,
    count_threshold_columns=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 10 December 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)

    # Organize information.
    count_columns = table.shape[1]
    count_rows = table.shape[0]

    # Determine whether the counts pass thresholds.
    if (
        (count_columns >= count_threshold_columns) and
        (count_rows >= count_threshold_rows)
    ):
        check = True
    else:
        check = False
        pass

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: check_table_wide_correlations()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return check


def sort_table_columns_rows_diagonal_match_pairs(
    table=None,
    name_index_rows=None,
    name_index_columns=None,
    column_sort_temporary=None,
    pairs_row_column=None,
    report=None,
):
    """
    Blank.

    Table orients first set of features across rows.
    Table orients second set of features across columns.

    By design this function should not filter or remove any columns or rows
    from the table.

    Review or revision: TCW; 12 November 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)
    pairs_row_column = copy.deepcopy(pairs_row_column)

    # Extract identifiers of columns and rows in table.
    columns_available = copy.deepcopy(
        table.columns.unique().tolist()
    )
    columns_available.remove(name_index_rows)
    rows_available = copy.deepcopy(
        table[name_index_rows].to_list()
    )

    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table.set_index(
        [name_index_rows,],
        append=False,
        drop=True,
        inplace=True,
    )

    # Collect information.
    values_pairs = list()
    # Iterate on rows in table.
    for row in rows_available:
        column = pairs_row_column[row]
        if (column in columns_available):
            value_pair = table.at[row, column]
        else:
            value_pair = float("nan")
            pass
        # Collect information
        values_pairs.append(value_pair)
        pass

    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )

    # Create new column in table.
    table[column_sort_temporary] = values_pairs

    # Sort rows in table.
    table.sort_values(
        by=[column_sort_temporary,],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Remove unnecessary columns.
    table.drop(
        labels=[column_sort_temporary,],
        axis="columns",
        inplace=True
    )
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    # Sort sequence of columns in table to match sequence of rows in table,
    # using a prior reference of matching pairs.
    table = sort_table_columns_match_rows(
        table=table,
        name_index_rows=name_index_rows,
        name_index_columns=name_index_columns,
        pairs_row_column=pairs_row_column,
        report=report,
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: sort_table_columns_rows_diagonal_match_pairs()")
        putly.print_terminal_partition(level=5)
        pass

    # Return information.
    return table


def sort_table_columns_match_rows(
    table=None,
    name_index_rows=None,
    name_index_columns=None,
    pairs_row_column=None,
    report=None,
):
    """
    Blank.

    Table orients first set of features across rows.
    Table orients second set of features across columns.

    By design this function should not filter or remove any columns or rows
    from the table.

    Review or revision: TCW; 12 November 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)
    pairs_row_column = copy.deepcopy(pairs_row_column)

    # Sort sequence of columns in table to match sequence of rows in table,
    # using a prior reference of matching pairs.

    # Extract identifiers of columns and rows in table.
    columns_sequence_original = copy.deepcopy(
        table.columns.unique().tolist()
    )
    columns_sequence_original.remove(name_index_rows)
    rows_sequence = copy.deepcopy(
        table[name_index_rows].to_list()
    )

    # Determine sequence of columns to match rows.
    # Notice that reference must include all possible pairs.
    columns_sequence_novel = list(map(
        lambda item: pairs_row_column[item], rows_sequence
    ))
    # Filter columns by their availability in the table.
    columns_sequence_novel_available = list(filter(
        lambda item: (item in columns_sequence_original),
        columns_sequence_novel
    ))
    columns_sequence_novel_available.insert(0, name_index_rows,)
    columns_sequence_novel_available = putly.collect_unique_items(
        items=columns_sequence_novel_available,
    )
    # Sort sequence of columns in table.
    # Do not remove any extra columns that do not match the rows, but sort
    # these extra columns at the end.
    table = porg.sort_table_columns_explicit_other(
        table=table,
        columns_sequence=columns_sequence_novel_available,
        sort_other=True,
        report=report,
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: sort_table_columns_match_rows()")
        putly.print_terminal_partition(level=5)
        pass

    # Return information.
    return table


def sort_table_rows_match_columns(
    table=None,
    name_index_rows=None,
    name_index_columns=None,
    pairs_column_row=None,
    report=None,
):
    """
    Blank.

    Table orients first set of features across rows.
    Table orients second set of features across columns.

    By design this function should not filter or remove any columns or rows
    from the table.

    Review or revision: TCW; 12 November 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)
    pairs_column_row = copy.deepcopy(pairs_column_row)

    # Sort sequence of rows in table to match sequence of columns in table,
    # using a prior reference of matching pairs.

    # Extract identifiers of columns and rows in table.
    columns_sequence = copy.deepcopy(
        table.columns.unique().tolist()
    )
    columns_sequence.remove(name_index_rows)
    rows_sequence_original = copy.deepcopy(
        table[name_index_rows].to_list()
    )

    # Determine sequence of rows to match columns.
    # Notice that reference must include all possible pairs.
    rows_sequence_novel = list(map(
        lambda item: pairs_column_row[item], columns_sequence
    ))
    # Filter rows by their availability in the table.
    rows_sequence_novel_available = list(filter(
        lambda item: (item in rows_sequence_original),
        rows_sequence_novel
    ))
    # Filter any rows without matching columns.
    rows_other = list(filter(
        lambda item: (item not in rows_sequence_novel_available),
        rows_sequence_original
    ))
    rows_other_sort = sorted(rows_other)
    # Organize sequence of rows.
    rows_sequence_all = list()
    rows_sequence_all.extend(rows_sequence_novel_available)
    rows_sequence_all.extend(rows_other_sort)
    rows_sequence_all = putly.collect_unique_elements(
        elements=rows_sequence_all,
    )
    # Define reference for sort.
    sequence_rows_novel = dict()
    index = 0
    for name in rows_sequence_all:
        sequence_rows_novel[name] = index
        index += 1
        pass
    # Sort rows in table.
    table = porg.sort_table_rows_by_single_column_reference(
        table=table,
        index_rows=name_index_rows,
        column_reference=name_index_rows,
        column_sort_temporary="sort_temporary",
        reference_sort=sequence_rows_novel,
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: sort_table_rows_match_columns()")
        putly.print_terminal_partition(level=5)
        pass

    # Return information.
    return table


def cluster_sort_sequence_features_correlations(
    table=None,
    name_index_rows=None,
    name_index_columns=None,
    match_pairs=None,
    pairs_first_second=None,
    pairs_second_first=None,
    sort_match_pairs_diagonal=None,
    cluster_features_first=None,
    cluster_features_second=None,
    sort_other_features=None,
    report=None,
):
    """
    Blank.

    Table orients first set of features across rows.
    Table orients second set of features across columns.

    By design this function should not filter or remove any columns or rows
    from the table.

    Review or revision: TCW; 12 November 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)

    ##########
    # Sort by self pairs along diagonal.
    # Determine whether to sort rows and columns by magnitude of correlations
    # for self pairs.
    if (
        (match_pairs) and
        (sort_match_pairs_diagonal)
    ):
        # Sort rows and columns in table by magnitude of self pairs.
        table = sort_table_columns_rows_diagonal_match_pairs(
            table=table,
            name_index_rows=name_index_rows,
            name_index_columns=name_index_columns,
            column_sort_temporary="temporary_sort_pairs_magnitude",
            pairs_row_column=pairs_first_second,
            report=report,
        )
        pass

    ##########
    # Cluster.

    # Determine whether to apply clustering to sequence of columns and rows in
    # table.
    if (
        (not sort_match_pairs_diagonal) and
        (
            (cluster_features_first) or
            (cluster_features_second)
        )
    ):
        # Organize indices in table.
        table.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table.set_index(
            [name_index_rows,],
            append=False,
            drop=True,
            inplace=True,
        )
        if (cluster_features_first):
            # Cluster sequence of first set of features.
            # Cluster rows in table.
            table = porg.cluster_table_rows(
                table=table,
            )
        if (cluster_features_second):
            # Cluster sequence of second set of features.
            # Cluster columns in table.
            table = porg.cluster_table_columns(
                table=table,
            )
            pass
        # Organize indices in table.
        table.reset_index(
            level=None,
            inplace=True,
            drop=False, # remove index; do not move to regular columns
        )
        pass

    ##########
    # Sort.

    # Determine whether to apply sorting to sequence of columns and rows in
    # table.
    if (
        (not sort_match_pairs_diagonal) and
        (cluster_features_first) and
        (not cluster_features_second) and
        (match_pairs) and
        (sort_other_features)
    ):
        # Sort sequence of second set of features to match sequence of first
        # set of features.
        # Sort sequence of columns in table to match sequence of rows in table,
        # using a prior reference of matching pairs.
        table = sort_table_columns_match_rows(
            table=table,
            name_index_rows=name_index_rows,
            name_index_columns=name_index_columns,
            pairs_row_column=pairs_first_second,
            report=report,
        )
    elif (
        (not sort_match_pairs_diagonal) and
        (cluster_features_second) and
        (not cluster_features_first) and
        (match_pairs) and
        (sort_other_features)
    ):
        # Sort sequence of first set of features to match sequence of second
        # set of features.
        # Sort sequence of rows in table to match sequence of columns in table,
        # using a prior reference of matching pairs.
        table = sort_table_rows_match_columns(
            table=table,
            name_index_rows=name_index_rows,
            name_index_columns=name_index_columns,
            pairs_column_row=pairs_second_first,
            report=report,
        )
    elif (
        (not sort_match_pairs_diagonal) and
        (not cluster_features_first) and
        (not cluster_features_second) and
        (sort_other_features)
    ):
        # Extract identifiers of rows in table.
        rows_sequence_original = copy.deepcopy(
            table[name_index_rows].to_list()
        )
        # Extract identifiers of columns in table.
        columns_sequence_original = copy.deepcopy(
            table.columns.unique().tolist()
        )
        columns_sequence_original.remove(name_index_rows)
        # Sort sequences of features from columns and rows in table.
        rows_sequence_novel = sorted(rows_sequence_original)
        columns_sequence_novel = sorted(columns_sequence_original)
        # Sort rows in table.
        table.sort_values(
            by=[name_index_rows,],
            axis="index",
            ascending=True,
            inplace=True,
        )
        # Copy information.
        columns_sequence_novel.insert(0, name_index_rows,)
        # Filter and sort columns in table.
        table = porg.filter_sort_table_columns(
            table=table,
            columns_sequence=columns_sequence_novel,
            report=report,
        )
        pass

    # Extract identifiers of rows in table.
    rows_sequence = copy.deepcopy(
        table[name_index_rows].to_list()
    )
    # Extract identifiers of columns in table.
    columns_sequence = copy.deepcopy(
        table.columns.unique().tolist()
    )

    # Bundle information.
    pail = dict()
    pail["table"] = table
    pail["rows_sequence"] = rows_sequence
    pail["columns_sequence"] = columns_sequence

    # Report.
    if report:
        # Organize information.
        count_columns = len(columns_sequence)
        count_rows = len(rows_sequence)
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: cluster_sort_sequence_features_correlations()")
        putly.print_terminal_partition(level=5)
        print("check parameters for cluster and sort...")
        print("match_pairs: " + str(match_pairs))
        print("sort_match_pairs_diagonal: " + str(sort_match_pairs_diagonal))
        print("cluster_features_first: " + str(cluster_features_first))
        print("cluster_features_second: " + str(cluster_features_second))
        print("sort_other_features: " + str(sort_other_features))
        putly.print_terminal_partition(level=5)
        print("table after cluster and sort operations:")
        print(table)
        putly.print_terminal_partition(level=5)
        print("columns:")
        print(columns_sequence)
        putly.print_terminal_partition(level=5)
        print("rows:")
        print(rows_sequence)
        putly.print_terminal_partition(level=5)
        print("count of columns: " + str(count_columns))
        print("count of rows: " + str(count_rows))
        putly.print_terminal_partition(level=5)
        pass

    # Return information.
    return pail


def manage_calculate_correlations_organize_tables_wide(
    table=None,
    column_identifier_observation=None,
    features_first=None,
    features_second=None,
    proportion_nonmissing_observations=None,
    z_score=None,
    type_correlation=None,
    match_pairs=None,
    pairs_first_second=None,
    pairs_second_first=None,
    intersect_features=None,
    threshold_filter_first=None,
    threshold_filter_second=None,
    threshold_type=None,
    threshold_value=None,
    threshold_proportion=None,
    threshold_optimization=None,
    threshold_optimization_count=None,
    sort_match_pairs_diagonal=None,
    cluster_features_first=None,
    cluster_features_second=None,
    sort_other_features=None,
    report=None,
):
    """
    Blank.

    Review or revision: TCW; 11 December 2025
    Review or revision: TCW; 13 November 2025
    Review or revision: TCW; 6 November 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Determine type of correlation.
    if (type_correlation == "pearson"):
        type_value_r = "correlation_pearson"
        type_value_p = "p_pearson"
    elif (type_correlation == "spearman"):
        type_value_r = "correlation_spearman"
        type_value_p = "p_spearman"
    elif (type_correlation == "kendall"):
        type_value_r = "correlation_kendall"
        type_value_p = "p_kendall"
        pass

    # Calculate correlations.
    table_correlation = calculate_correlations_populate_table_wide(
        table=table,
        column_identifier_observation=column_identifier_observation,
        features_first=features_first,
        features_second=features_second,
        proportion_nonmissing_observations=(
            proportion_nonmissing_observations
        ),
        z_score=z_score,
        type_value=type_value_r,
        report=report,
    )
    if True:
        table_p_value = calculate_correlations_populate_table_wide(
            table=table,
            column_identifier_observation=column_identifier_observation,
            features_first=features_first,
            features_second=features_second,
            proportion_nonmissing_observations=(
                proportion_nonmissing_observations
            ),
            z_score=z_score,
            type_value=type_value_p,
            report=report,
        )
    else:
        table_p_value = pandas.DataFrame()
        pass

    # Replace any infinite values with missing values.
    table_correlation.replace(
        value=pandas.NA,
        to_replace=[numpy.inf, -numpy.inf],
        inplace=True,
    )

    # Filter columns and rows in table by nonmissing values.
    # Even after previously filtering the first and second selections of
    # features to ensure that they have an adequate proportion of nonmissing
    # values, it is still possible to have missing values for correlations. The
    # calculation of correlation requires adequate nonmissing pairs of matching
    # values between each of the relevant features.
    table_correlation.dropna(
        axis="index",
        how="all",
        inplace=True,
    )
    table_correlation.dropna(
        axis="columns",
        how="all",
        inplace=True,
    )
    # After removing features with missing values across all correlations,
    # optionally consider filling any lingering missing values with zero.
    if False:
        # For some reason, this method (dataframe.fillna) creates problems
        # downstream.
        table_correlation.fillna(
            value=0.0,
            axis="index",
            inplace=True,
        )
        table_correlation.fillna(
            value=0.0,
            axis="columns",
            inplace=True,
        )
        pass
    table_correlation.replace(
        value=0.0,
        to_replace=pandas.NA,
        inplace=True,
    )
    table_correlation.replace(
        value=0.0,
        to_replace=numpy.nan,
        inplace=True,
    )
    table_correlation.replace(
        value=0.0,
        to_replace=math.nan,
        inplace=True,
    )
    table_correlation.replace(
        value=0.0,
        to_replace=float("nan"),
        inplace=True,
    )

    # Optionally determine optimal parameters for threshold on correlations.
    pail_optimization = optimize_threshold_table_wide_correlations(
        table=table_correlation,
        name_index_rows="feature_first",
        name_index_columns="feature_second",
        threshold_filter_rows=threshold_filter_first,
        threshold_filter_columns=threshold_filter_second,
        threshold_type=threshold_type,
        threshold_value=threshold_value,
        threshold_proportion=threshold_proportion,
        threshold_optimization=threshold_optimization,
        threshold_optimization_count=threshold_optimization_count,
        report=report,
    )
    #pail_optimization["threshold_value"]

    # Filter columns and rows in table by threshold.
    pail_threshold = filter_table_wide_correlations_by_threshold(
        table_r=table_correlation,
        table_p=table_p_value,
        name_index_rows="feature_first",
        name_index_columns="feature_second",
        threshold_filter_rows=threshold_filter_first,
        threshold_filter_columns=threshold_filter_second,
        threshold_type=threshold_type,
        threshold_value=pail_optimization["threshold_value"],
        threshold_proportion=threshold_proportion,
        match_pairs=match_pairs,
        pairs_row_column=pairs_first_second,
        pairs_column_row=pairs_second_first,
        intersect_features=intersect_features,
        report=report,
    )

    # Check that the table is not empty and otherwise suitable for cluster
    # and sort operations.
    check_table = check_table_wide_correlations(
        table=pail_threshold["table_r_filter"],
        count_threshold_rows=5,
        count_threshold_columns=5,
        report=report,
    )
    # Cluster and sort sequences of features.
    if check_table:
        pail_cluster_sort = cluster_sort_sequence_features_correlations(
            table=pail_threshold["table_r_filter"],
            name_index_rows="feature_first",
            name_index_columns="feature_second",
            match_pairs=match_pairs,
            pairs_first_second=pairs_first_second,
            pairs_second_first=pairs_second_first,
            sort_match_pairs_diagonal=sort_match_pairs_diagonal,
            cluster_features_first=cluster_features_first,
            cluster_features_second=cluster_features_second,
            sort_other_features=sort_other_features,
            report=report,
        )
        #pail["table"]
        #pail["rows_sequence"]
        #pail["columns_sequence"]
    else:
        pail_cluster_sort = dict()
        pail_cluster_sort["table"] = pail_threshold["table_r_filter"]
        pass

    ##########
    # Copy information.
    table_correlation = table_correlation.copy(deep=True)
    table_p_value = table_p_value.copy(deep=True)
    table_correlation_filter = pail_threshold["table_r_filter"].copy(deep=True)
    table_correlation_filter_cluster_sort = (
        pail_cluster_sort["table"].copy(deep=True)
    )

    # Bundle information.
    pail = dict()
    pail["table_correlation"] = table_correlation
    pail["table_p_value"] = table_p_value
    pail["table_correlation_filter"] = table_correlation_filter
    pail["table_correlation_filter_cluster_sort"] = (
        table_correlation_filter_cluster_sort
    )
    pail["threshold_value"] = pail_optimization["threshold_value"]

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: calculate_correlations_organize_tables_wide()")
        putly.print_terminal_partition(level=5)
        print("type of correlation: " + str(type_correlation))
        putly.print_terminal_partition(level=5)
        #print("table of correlations:")
        #print(table_correlation)
        #putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def prepare_text_summary_optimization_threshold(
    title=None,
    threshold_filter_first=None,
    threshold_filter_second=None,
    threshold_type=None,
    threshold_value=None,
    threshold_proportion=None,
    threshold_optimization=None,
    threshold_optimization_count=None,
    report=None,
):
    """
    Prepare text information to print to terminal or write to file.

    Review: TCW; 10 December 2025

    arguments:
        title (str): title for summary text
        ...
        report (bool): whether to print reports

    raises:

    returns:
        (str): character string representation of information

    """

    summary_text = str(
        textwrap.dedent("""\

            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------
        """) +
        str(title) +
        textwrap.dedent("""\
            --------------------------------------------------

        """) +

        str("threshold_filter_first: " + str(threshold_filter_first) + "\n") +
        str(
            "threshold_filter_second: " +
            str(threshold_filter_second) + "\n"
        ) +
        str("threshold_type: " + str(threshold_type) + "\n") +
        str("threshold_value: " + str(threshold_value) + "\n") +
        str("threshold_proportion: " + str(threshold_proportion) + "\n") +
        str("threshold_optimization: " + str(threshold_optimization) + "\n") +
        str("threshold_optimization_count: " +
        str(threshold_optimization_count) + "\n") +
        textwrap.dedent("""\

            ----------

        """) +
        textwrap.dedent("""\
            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------

        """)
    )

    # Report.
    if report:
        # Organize.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        name_module = str(
            "regression.py"
        )
        print("module: " + name_module)
        name_function = str(
            "prepare_text_summary_optimization_threshold()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("summary_text:")
        print(summary_text)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return summary_text


def create_write_plot_chart_heatmap(
    path_directory_parent=None,
    name_chart=None,
    table=None,
    name_index_columns=None,
    name_index_rows=None,
    plot_scale_minimum=None,
    plot_scale_center=None,
    plot_scale_maximum=None,
    title_chart=None,
    title_bar=None,
    title_abscissa=None,
    title_ordinate=None,
    report=None,
):
    """
    Create and plot a chart of the heatmap type.

    Original source table must not have an explicitly defined index across
    rows.

    Review or revision: TCW; 11 December 2025
    Review or revision: TCW; 21 October 2025

    arguments:
        ...
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table = table.copy(deep=True)

    # Determine whether to impose minimal and maximal thresholds on values for
    # scale representation on chart.
    if (
        ((plot_scale_minimum is None) or (plot_scale_maximum is None)) or
        (math.isnan(plot_scale_minimum) or math.isnan(plot_scale_maximum))
    ):
        table_extract = table.copy(deep=True)
        # Organize indices in table.
        table_extract.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table_extract.columns.rename(
            name_index_columns,
            inplace=True,
        ) # single-dimensional index
        table_extract.set_index(
            name_index_rows,
            append=False,
            drop=True,
            inplace=True
        )
        # Extract minimal and maximal values of signal intensity.
        matrix = numpy.copy(table_extract.to_numpy())
        #value_minimum = round((numpy.nanmin(matrix) - 0.005), 2)
        #value_maximum = round((numpy.nanmax(matrix) + 0.005), 2)
        round_offset = abs(numpy.nanmin(matrix) * 0.01)
        value_minimum = round((numpy.nanmin(matrix) - round_offset), 3)
        value_center = 0.0
        value_maximum = round((numpy.nanmax(matrix) + round_offset), 3)
    else:
        value_minimum = plot_scale_minimum
        value_center = plot_scale_center
        value_maximum = plot_scale_maximum
        pass

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = splot.plot_heatmap_signal_features_observations_labels(
        table=table,
        format_table=1, # 1: features in rows, observations or groups in columns
        index_columns=name_index_columns,
        index_rows=name_index_rows,
        transpose_table=True,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_center=value_center, # correlations, effects, fold changes, components etc
        value_maximum=value_maximum,
        title_ordinate="",
        title_abscissa="",
        title_bar=title_bar,
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        size_title_ordinate="eight",
        size_title_abscissa="eight",
        size_title_bar="twelve",
        size_label_ordinate=None, # determine automatically if "None"; "fifteen"
        size_label_abscissa=None, # determine automatically if "None"
        size_label_bar="thirteen",
        show_labels_ordinate=True,
        show_labels_abscissa=True,
        show_scale_bar=True,
        aspect="square", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )

    # Write product information to file.

    # Bundle information.
    pail_write_plot = dict()
    pail_write_plot[name_chart] = figure

    # Write figure object to file.
    pplot.write_product_plots_parent_directory(
        pail_write=pail_write_plot,
        format="jpg", # jpg, png, svg
        resolution=150, # low-resolution trial: 96 DPI; medium-resolution: 150 DPI; high-resolution print: 300 DPI
        path_directory=path_directory_parent,
    )

    # Return information.
    return figure


def manage_create_write_plot_charts(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    table=None,
    name_index_columns=None,
    name_index_rows=None,
    name_correlations=None,
    type_correlation=None,
    plot_scale_minimum=None,
    plot_scale_center=None,
    plot_scale_maximum=None,
    report=None,
):
    """
    Blank.

    Review or revision: TCW; 11 December 2025
    Review or revision: TCW; 28 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)

    # Determine title for scale bar on chart.
    if (type_correlation == "pearson"):
        title_bar = "Correlation Coefficient (Pearson)"
    elif (type_correlation == "spearman"):
        title_bar = "Correlation Coefficient (Spearman)"
        pass

    ##########
    # Define paths to directories.
    path_directory_charts = os.path.join(
        path_directory_product, "charts",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_charts,
    )

    # Check that the table is not empty and otherwise suitable for plot chart.
    check_table = check_table_wide_correlations(
        table=table,
        count_threshold_rows=3,
        count_threshold_columns=3,
        report=report,
    )

    # Create and write plot charts for the table of correlations.
    if (check_table):
        create_write_plot_chart_heatmap(
            path_directory_parent=path_directory_charts,
            name_chart="correlations",
            table=table,
            name_index_columns=name_index_columns,
            name_index_rows=name_index_rows,
            plot_scale_minimum=plot_scale_minimum,
            plot_scale_center=plot_scale_center,
            plot_scale_maximum=plot_scale_maximum,
            title_chart="",
            title_bar=title_bar,
            title_abscissa="Features First",
            title_ordinate="Features Second",
            report=report,
        )
        pass
    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
        )
        print(str("module: " + module))
        print("function: manage_create_write_plot_charts()")
        putly.print_terminal_partition(level=5)
        pass

    pass


################################################################################
# Procedure


##########
# Control procedure for a single instance of parameters.


def control_procedure(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_features_observations=None,
    path_file_source_list_features_first=None,
    path_file_source_list_features_second=None,
    path_file_source_table_groups_observations=None,
    column_identifier_observation=None,
    column_name_observation=None,
    proportion_nonmissing_observations=None,
    type_correlation=None,
    match_pairs=None,
    prefix_match_first=None,
    suffix_match_first=None,
    prefix_match_second=None,
    suffix_match_second=None,
    intersect_features=None,
    threshold_filter_first=None,
    threshold_filter_second=None,
    threshold_type=None,
    threshold_value=None,
    threshold_proportion=None,
    threshold_optimization=None,
    threshold_optimization_count=None,
    sort_match_pairs_diagonal=None,
    cluster_features_first=None,
    cluster_features_second=None,
    sort_other_features=None,
    plot_scale_minimum=None,
    plot_scale_center=None,
    plot_scale_maximum=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review or revision: TCW; 2 December 2025

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_file_source_table_features_observations (str): path to source file
        path_file_source_list_features_first (str): path to source file
        path_file_source_list_features_second (str): path to source file
        path_file_source_table_groups_observations (str): path to source file
        column_identifier_observation (str): name of column in source table
        column_name_observation (str): name of column in source table
        proportion_nonmissing_observations (float): threshold by proportion of
            observations that must have nonmissing values for both features in
            each pair
        type_correlation (str): type of correlation for calculation; either
            'pearson' or 'spearman'
        ...
        cluster_features (str): how to apply the cluster operation, either
            'none', 'both', 'first', or 'second'
        sort_other_features (bool): whether to sort the set of features to
            which the cluster operation was not applied
        plot_scale_minimum (float): minimal value to set as threshold for
            representation in the plot chart
        plot_scale_maximum (float): maximal value to set as threshold for
            representation in the plot chart
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        path_file_source_list_features_first=(
            path_file_source_list_features_first
        ),
        path_file_source_list_features_second=(
            path_file_source_list_features_second
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        report=report,
    )
    #pail_source["table_features_observations"]
    #pail_source["table_groups"]
    #pail_source["features"]["first"]
    #pail_source["features"]["second"]

    ##########
    # Report.
    if report:
        # Organize information.
        count_first = len(pail_source["features"]["first"])
        count_second = len(pail_source["features"]["second"])
        # Print information.
        putly.print_terminal_partition(level=5)
        print("count of first features: " + str(count_first))
        print("count of second features: " + str(count_second))
        putly.print_terminal_partition(level=5)
        print("features first")
        print(pail_source["features"]["first"])
        putly.print_terminal_partition(level=5)
        print("features second")
        print(pail_source["features"]["second"])
        pass

    # Filter table's columns for features and rows for observations.
    table_filter = filter_organize_table_columns_rows(
        table=pail_source["table_features_observations"],
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        features_first=pail_source["features"]["first"],
        features_second=pail_source["features"]["second"],
        table_groups_observations=pail_source["table_groups"],
        report=report,
    )

    ##########
    # Filter lists of features.

    # Filter lists of features corresponding to columns in table by their
    # availability of values across rows in table corresponding to
    # observations.
    pail_features_availability = filter_features_by_available_observations(
        table=table_filter,
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        features_first=pail_source["features"]["first"],
        features_second=pail_source["features"]["second"],
        proportion_nonmissing_observations=(
            proportion_nonmissing_observations
        ),
        report=report,
    )
    #pail_features_availability["features_first"]
    #pail_features_availability["features_second"]

    # Determie match pairs of features in first and second sets.
    if (match_pairs):
        pail_pairs = match_filter_features_pairs(
            features_first=pail_features_availability["features_first"],
            features_second=pail_features_availability["features_second"],
            prefix_match_first=prefix_match_first,
            suffix_match_first=suffix_match_first,
            prefix_match_second=prefix_match_second,
            suffix_match_second=suffix_match_second,
            report=report,
        )
    else:
        pail_pairs = dict() # need to specify appropriate missing values (None)
        pail_pairs["features_first_intersection"] = None
        pail_pairs["features_second_intersection"] = None
        pail_pairs["translations_intersection"] = None
        pail_pairs["translations_union"] = None
        pail_pairs["pairs_first_second"] = None
        pail_pairs["pairs_second_first"] = None
        pass

    # Determine whether to use features after filtering by intersection.
    if (
        (match_pairs) and
        (intersect_features)
    ):
        # Copy information.
        features_first = copy.deepcopy(
            pail_pairs["features_first_intersection"]
        )
        features_second = copy.deepcopy(
            pail_pairs["features_second_intersection"]
        )
    else:
        # Copy information.
        features_first = pail_features_availability["features_first"]
        features_second = pail_features_availability["features_second"]
        pass

    ##########
    # Calculate correlations.

    # Calculate correlations between pairs of features and organize these
    # within a table in long format.
    if False:
        table_correlation_long = calculate_correlations_populate_table_long(
            table=table_filter,
            column_identifier_observation=column_identifier_observation,
            features_first=features_first,
            features_second=features_second,
            proportion_nonmissing_observations=(
                proportion_nonmissing_observations
            ),
            z_score=True,
            report=report,
        )
    else:
        table_correlation_long = pandas.DataFrame()
        pass

    # Calculate correlations between pairs of features and organize these
    # within a table in wide format.
    pail_correlation = manage_calculate_correlations_organize_tables_wide(
        table=table_filter,
        column_identifier_observation=column_identifier_observation,
        features_first=features_first,
        features_second=features_second,
        proportion_nonmissing_observations=proportion_nonmissing_observations,
        z_score=True,
        type_correlation=type_correlation,
        match_pairs=match_pairs,
        pairs_first_second=pail_pairs["pairs_first_second"],
        pairs_second_first=pail_pairs["pairs_second_first"],
        intersect_features=intersect_features,
        threshold_filter_first=threshold_filter_first,
        threshold_filter_second=threshold_filter_second,
        threshold_type=threshold_type,
        threshold_value=threshold_value,
        threshold_proportion=threshold_proportion,
        threshold_optimization=threshold_optimization,
        threshold_optimization_count=threshold_optimization_count,
        sort_match_pairs_diagonal=sort_match_pairs_diagonal,
        cluster_features_first=cluster_features_first,
        cluster_features_second=cluster_features_second,
        sort_other_features=sort_other_features,
        report=report,
    )
    #pail_correlation["table_correlation"]
    #pail_correlation["table_p_value"]
    #pail_correlation["table_correlation_filter"]
    #pail_correlation["table_correlation_filter_cluster_sort"]
    #pail_correlation["threshold_value"]

    # Prepare text summary.
    summary_text = prepare_text_summary_optimization_threshold(
        title=str("Parameters for optimization of threshold."),
        threshold_filter_first=threshold_filter_first,
        threshold_filter_second=threshold_filter_second,
        threshold_type=threshold_type,
        threshold_value=pail_correlation["threshold_value"],
        threshold_proportion=threshold_proportion,
        threshold_optimization=threshold_optimization,
        threshold_optimization_count=threshold_optimization_count,
        report=report,
    )

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Text.
    pail_write_text = dict()
    pail_write_text[str("notes_optimization")] = summary_text
    # Lists.
    pail_write_lists = dict()
    pail_write_lists["features_first"] = features_first
    pail_write_lists["features_second"] = features_second
    # Tables.
    pail_write_tables = dict()
    pail_write_tables["table_correlation_long"] = table_correlation_long
    pail_write_tables["table_correlation_wide"] = (
        pail_correlation["table_correlation"]
    )
    pail_write_tables["table_p_value_wide"] = (
        pail_correlation["table_p_value"]
    )
    pail_write_tables["table_correlation_wide_filter"] = (
        pail_correlation["table_correlation_filter"]
    )
    pail_write_tables["table_correlation_wide_filter_cluster_sort"] = (
        pail_correlation["table_correlation_filter_cluster_sort"]
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
    # Create directories.
    putly.create_directories(
        path=path_directory_lists,
    )
    putly.create_directories(
        path=path_directory_tables,
    )
    # Text.
    putly.write_character_strings_to_file_text(
        pail_write=pail_write_text,
        path_directory=path_directory_product,
    )
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
    # Plot charts.

    # Manage the creation and write of plot charts.
    # Calculate principal components.
    manage_create_write_plot_charts(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        table=pail_correlation["table_correlation_filter_cluster_sort"],
        name_index_columns="feature_second",
        name_index_rows="feature_first",
        name_correlations="blank",
        type_correlation=type_correlation,
        plot_scale_minimum=plot_scale_minimum,
        plot_scale_center=plot_scale_center,
        plot_scale_maximum=plot_scale_maximum,
        report=report,
    )

    pass


##########
# Execute main procedure.


def execute_procedure(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_features_observations=None,
    path_file_source_list_features_first=None,
    path_file_source_list_features_second=None,
    path_file_source_table_groups_observations=None,
    column_identifier_observation=None,
    column_name_observation=None,
    proportion_nonmissing_observations=None,
    type_correlation=None,
    match_pairs=None,
    prefix_match_first=None,
    suffix_match_first=None,
    prefix_match_second=None,
    suffix_match_second=None,
    intersect_features=None,
    threshold_filter_first=None,
    threshold_filter_second=None,
    threshold_type=None,
    threshold_value=None,
    threshold_proportion=None,
    sort_match_pairs_diagonal=None,
    cluster_features_first=None,
    cluster_features_second=None,
    sort_other_features=None,
    plot_scale_minimum=None,
    plot_scale_center=None,
    plot_scale_maximum=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review or revision: TCW; 10 November 2025
    Review or revision: TCW; 6 November 2025
    Review or revision: TCW; 26 October 2025

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_file_source_table_features_observations (str): path to source file
        path_file_source_list_features_first (str): path to source file
        path_file_source_list_features_second (str): path to source file
        path_file_source_table_groups_observations (str): path to source file
        column_identifier_observation (str): name of column in source table
        column_name_observation (str): name of column in source table
        proportion_nonmissing_observations (float): threshold by proportion of
            observations that must have nonmissing values for both features in
            each pair
        type_correlation (str): type of correlation for calculation; either
            'pearson' or 'spearman'
        ...
        cluster_features (str): how to apply the cluster operation, either
            'none', 'both', 'first', or 'second'
        sort_other_features (bool): whether to sort the set of features to
            which the cluster operation was not applied
        plot_scale_minimum (float): minimal value to set as threshold for
            representation in the plot chart
        plot_scale_maximum (float): maximal value to set as threshold for
            representation in the plot chart
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
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        path_file_source_list_features_first=(
            path_file_source_list_features_first
        ),
        path_file_source_list_features_second=(
            path_file_source_list_features_second
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        proportion_nonmissing_observations=proportion_nonmissing_observations,
        type_correlation=type_correlation,
        match_pairs=match_pairs,
        prefix_match_first=prefix_match_first,
        suffix_match_first=suffix_match_first,
        prefix_match_second=prefix_match_second,
        suffix_match_second=suffix_match_second,
        intersect_features=intersect_features,
        threshold_filter_first=threshold_filter_first,
        threshold_filter_second=threshold_filter_second,
        threshold_type=threshold_type,
        threshold_value=threshold_value,
        threshold_proportion=threshold_proportion,
        sort_match_pairs_diagonal=sort_match_pairs_diagonal,
        cluster_features_first=cluster_features_first,
        cluster_features_second=cluster_features_second,
        sort_other_features=sort_other_features,
        plot_scale_minimum=plot_scale_minimum,
        plot_scale_center=plot_scale_center,
        plot_scale_maximum=plot_scale_maximum,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: execute_procedure()")
        print("system: local")
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Control procedure for current instance of parameters.
    control_procedure(
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        path_file_source_table_features_observations=(
            pail_parameters["path_file_source_table_features_observations"]
        ),
        path_file_source_list_features_first=(
            pail_parameters["path_file_source_list_features_first"]
        ),
        path_file_source_list_features_second=(
            pail_parameters["path_file_source_list_features_second"]
        ),
        path_file_source_table_groups_observations=(
            pail_parameters["path_file_source_table_groups_observations"]
        ),
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_name_observation=pail_parameters["column_name_observation"],
        proportion_nonmissing_observations=(
            pail_parameters["proportion_nonmissing_observations"]
        ),
        type_correlation=pail_parameters["type_correlation"],
        match_pairs=pail_parameters["match_pairs"],
        prefix_match_first=pail_parameters["prefix_match_first"],
        suffix_match_first=pail_parameters["suffix_match_first"],
        prefix_match_second=pail_parameters["prefix_match_second"],
        suffix_match_second=pail_parameters["suffix_match_second"],
        intersect_features=pail_parameters["intersect_features"],
        threshold_filter_first=pail_parameters["threshold_filter_first"],
        threshold_filter_second=pail_parameters["threshold_filter_second"],
        threshold_type=pail_parameters["threshold_type"],
        threshold_value=pail_parameters["threshold_value"],
        threshold_proportion=pail_parameters["threshold_proportion"],
        sort_match_pairs_diagonal=pail_parameters["sort_match_pairs_diagonal"],
        cluster_features_first=pail_parameters["cluster_features_first"],
        cluster_features_second=pail_parameters["cluster_features_second"],
        sort_other_features=pail_parameters["sort_other_features"],
        plot_scale_minimum=pail_parameters["plot_scale_minimum"],
        plot_scale_center=pail_parameters["plot_scale_center"],
        plot_scale_maximum=pail_parameters["plot_scale_maximum"],
        report=pail_parameters["report"],
    )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_source = sys.argv[1]
    path_directory_product = sys.argv[2]
    path_directory_dock = sys.argv[3]
    path_file_source_table_features_observations = sys.argv[4]
    path_file_source_list_features_first = sys.argv[5]
    path_file_source_list_features_second = sys.argv[6]
    path_file_source_table_groups_observations = sys.argv[7]
    column_identifier_observation = sys.argv[8]
    column_name_observation = sys.argv[9]
    proportion_nonmissing_observations = sys.argv[10]
    type_correlation = sys.argv[11]
    match_pairs = sys.argv[12]
    prefix_match_first = sys.argv[13]
    suffix_match_first = sys.argv[14]
    prefix_match_second = sys.argv[15]
    suffix_match_second = sys.argv[16]
    intersect_features = sys.argv[17]
    threshold_filter_first = sys.argv[18]
    threshold_filter_second = sys.argv[19]
    threshold_type = sys.argv[20]
    threshold_value = sys.argv[21]
    threshold_proportion = sys.argv[22]
    sort_match_pairs_diagonal = sys.argv[23]
    cluster_features_first = sys.argv[24]
    cluster_features_second = sys.argv[25]
    sort_other_features = sys.argv[26]
    plot_scale_minimum = sys.argv[27]
    plot_scale_center = sys.argv[28]
    plot_scale_maximum = sys.argv[29]
    report = sys.argv[30]


    # Call function for procedure.
    execute_procedure(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        path_file_source_list_features_first=(
            path_file_source_list_features_first
        ),
        path_file_source_list_features_second=(
            path_file_source_list_features_second
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        proportion_nonmissing_observations=proportion_nonmissing_observations,
        type_correlation=type_correlation,
        match_pairs=match_pairs,
        prefix_match_first=prefix_match_first,
        suffix_match_first=suffix_match_first,
        prefix_match_second=prefix_match_second,
        suffix_match_second=suffix_match_second,
        intersect_features=intersect_features,
        threshold_filter_first=threshold_filter_first,
        threshold_filter_second=threshold_filter_second,
        threshold_type=threshold_type,
        threshold_value=threshold_value,
        threshold_proportion=threshold_proportion,
        sort_match_pairs_diagonal=sort_match_pairs_diagonal,
        cluster_features_first=cluster_features_first,
        cluster_features_second=cluster_features_second,
        sort_other_features=sort_other_features,
        plot_scale_minimum=plot_scale_minimum,
        plot_scale_center=plot_scale_center,
        plot_scale_maximum=plot_scale_maximum,
        report=report,
    )

    pass



#
