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
# Date, last execution or modification: 29 September 2025
# Review: TCW; 29 September 2025
################################################################################
# Note


# The purpose of this module is to support the scripts or modules below within
# the package "partner".

# 1. plot_chart_heatmap_sets_features_groups_observations.py
# 2. calculate_principal_components_sets_features_groups_observations.py

# Recent example of usage:
# /.../pails_process/omega3/2025-09-22_heterogeneity_candidate_adipose_fibrosis

##########
# Review: TCW; 29 September 2025

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

#dir()
#importlib.reload()

###############################################################################
# Functionality


def determine_features_available(
    table_features=None,
    table_observations=None,
    table_signals=None,
    column_identifier_feature=None,
    column_name_feature=None,
    column_identifier_observation=None,
    column_identifier_signal=None,
    transpose_table_signals=None,
    report=None,
):
    """
    Determine list of features with available information in source tables.

    Available features are the union of features from "table_observations"
    and "table_signals" (dependent on orientation and transposition), but they
    must also be the intersection with features from "table_features".

    Review: TCW; 30 September 2025

    arguments:
        table_features (object): Pandas data-frame table
        table_observations (object): Pandas data-frame table
        table_signals (object): Pandas data-frame table
        column_identifier_feature (str): name of column in source table
        column_name_feature (str): name of column in source table
        column_identifier_observation (str): name of column in source table
        column_identifier_signal (str): name of column in source table
        transpose_table_signals (bool): whether to transpose the table of
            signals to match orientation in table of observations (columns:
            features; rows: observations)
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information
    """

    # Copy information.
    table_features = table_features.copy(deep=True)
    table_observations = table_observations.copy(deep=True)
    table_signals = table_signals.copy(deep=True)

    # Extract identifiers of features.
    features_features = copy.deepcopy(
        table_features[column_identifier_feature].unique().tolist()
    )
    features_observations = copy.deepcopy(
        table_features[column_identifier_feature].unique().tolist()
    )
    if (
        (table_signals is not None) and
        (transpose_table_signals is not None) and
        (not transpose_table_signals)
    ):
        features_signals = copy.deepcopy(
            table_signals.columns.to_list()
        )
        features_signals.remove(column_identifier_signal)
    elif (
        (table_signals is not None) and
        (transpose_table_signals is not None) and
        (transpose_table_signals)
    ):
        features_signals = copy.deepcopy(
            table_signals[column_identifier_feature].unique().tolist()
        )
    else:
        features_signals = list()

    # Combine identifiers of features.
    features_union = (
        putly.combine_sets_items_union_unique(
            sets_items=[features_observations, features_signals,],
            report=None,
        )
    )
    features_intersection = (
        putly.combine_sets_items_intersection_unique(
            items_first=features_union,
            items_second=features_features,
            report=None,
    ))

    # Report.
    if report:
        # Organize information.
        count_features = len(features_intersection)
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
        )
        print(str("module: " + module))
        print("function: determine_features_available()")
        putly.print_terminal_partition(level=5)
        print(str(
            "count of available features: " + str(count_features)
        ))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return features_intersection


def read_organize_parameters_sets_features(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    table=None,
    column_name=None,
    features_available=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 29 September 2025

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
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): bundle of information

    """

    # Copy information.
    table = table.copy(deep=True)
    features_available = copy.deepcopy(features_available)

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
    # Extract names of sets.
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
            "utility_special.py"
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


def organize_parameters_further_sets_features(
    table_features=None,
    column_identifier_feature=None,
    column_name_feature=None,
    features_selection=None,
    features_sets_union=None,
    sets_features=None,
    names_sets_features_sequence=None,
    prefix_name_feature=None,
    report=None,
):
    """
    Organize parameters and information about features.

    Review: TCW; 29 September 2025

    arguments:
        table_features (object): Pandas data-frame table
        column_identifier_feature (str): name of column in source table
        column_name_feature (str): name of column in source table
        features_selection (list<str>): identifiers of features for which to
            include signals across observations
        features_sets_union (list<str>): identifiers of features in union of
            all sets
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

    # Bundle information.
    pail = dict()
    pail["table_features_selection"] = table_features_selection
    pail["features_selection"] = features_selection
    pail["features_selection_translation"] = features_selection_translation
    pail["features_selection_prefix"] = features_selection_prefix
    pail["translations_features"] = translations_features
    pail["translations_features_prefix"] = translations_features_prefix
    pail["names_sets_features_sequence"] = names_sets_features_sequence
    pail["sets_features"] = sets_features

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "utility_special.py"
        )
        print(str("module: " + module))
        print("function: organize_parameters_further_sets_features()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def organize_parameters_groups_observations(
    table=None,
    column_name=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 29 September 2025

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
            "plot_chart_heatmap_sets_features_groups_observations.py"
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

    Review: TCW; 29 September 2025

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
    # Filter rows in table for selection of observations.
    if (len(observations_selection) > 0):
        table_observations_selection = table_observations.loc[
            table_observations[column_identifier_signal].isin(
                observations_selection
            ), :
        ].copy(deep=True)
        pass
    table_alias = table_observations_selection
    table_observations_selection = table_alias.loc[
        (table_alias[column_identifier_observation].str.len() > 0), :
    ].copy(deep=True)
    table_alias = table_observations_selection
    table_observations_selection = table_alias.loc[
        (table_alias[column_identifier_signal].str.len() > 0), :
    ].copy(deep=True)
    # Extract information for translation of names of columns.
    table_translations = table_observations_selection.filter(
        items=[column_identifier_signal, column_identifier_observation,],
        axis="columns",
    )
    series_translations = pandas.Series(
        table_translations[column_identifier_observation].to_list(),
        index=table_translations[column_identifier_signal],
    )
    translations_observations = series_translations.to_dict()

    # Bundle information.
    pail = dict()
    pail["table_observations_selection"] = table_observations_selection
    pail["observations_selection"] = observations_selection
    pail["translations_observations"] = translations_observations
    pail["names_groups_observations_sequence"] = (
        names_groups_observations_sequence
    )
    pail["groups_observations"] = groups_observations

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "utility_special.py"
        )
        print(str("module: " + module))
        print("function: organize_parameters_further_groups_observations()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail





################################################################################
# Procedure
# Currently, this module is not directly executable.

################################################################################
# End
