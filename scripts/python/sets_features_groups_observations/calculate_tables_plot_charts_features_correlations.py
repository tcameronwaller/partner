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
# Date, review or revision: 26 October 2025
################################################################################
# Note


##########
# Note:


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
    cluster_features=None,
    plot_threshold_minimum=None,
    plot_threshold_maximum=None,
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
        cluster_features (bool): whether to cluster the sequence of features in
            the tables and charts
        plot_threshold_minimum (float): minimal value to set as threshold for
            representation in the plot chart
        plot_threshold_maximum (float): maximal value to set as threshold for
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
    if (
        (len(str(plot_threshold_minimum)) > 0) and
        (str(plot_threshold_minimum) != "none")
    ):
        pail["plot_threshold_minimum"] = float(
            str(plot_threshold_minimum).strip()
        )
    else:
        plot_threshold_minimum = None
    if (
        (len(str(plot_threshold_maximum).strip()) > 0) and
        (str(plot_threshold_maximum).strip() != "none")
    ):
        pail["plot_threshold_maximum"] = float(
            str(plot_threshold_maximum).strip()
        )
    else:
        plot_threshold_maximum = None
        pass

    # Category.
    pail["type_correlation"] = str(type_correlation).strip()

    # Boolean, true or false.
    if (
        (cluster_features is not None) and
        (str(cluster_features) != "") and
        (str(cluster_features) != "none") and
        (str(cluster_features) == "true")
    ):
        pail["cluster_features"] = True
    else:
        pail["cluster_features"] = False
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
            "calculate_tables_plot_charts_features_correlations.py"
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
        column_identifier_signal=column_identifier_observation,
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

    # Filter table's columns corresponding to features by the available of
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
            record["proportion_threshold"] = proportion_nonmissing_observations
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
            else:
                # Introduce missing values.
                measures = list()
                measures.append("correlation_pearson")
                measures.append("p_pearson")
                #measures.append("confidence_95_low_pearson")
                #measures.append("confidence_95_high_pearson")
                measures.append("correlation_spearman")
                measures.append("p_spearman")
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

    # Filter and sort columns in table.
    columns_sequence = [
        "feature_first",
        "feature_second",
        "count_observations_total",
        "count_observations_pair",
        "proportion_actual",
        "proportion_threshold",
        "correlation_pearson", "p_pearson", "q_pearson",
        #"confidence_95_low_pearson", "confidence_95_high_pearson",
        "correlation_spearman", "p_spearman", "q_spearman",
        #"correlation_kendall", "p_kendall", "q_kendall",
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


def calculate_correlations_organize_tables_wide(
    table=None,
    column_identifier_observation=None,
    features_first=None,
    features_second=None,
    proportion_nonmissing_observations=None,
    z_score=None,
    type_correlation=None,
    cluster_features=None,
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

    # Determine type of correlation.
    if (type_correlation == "pearson"):
        type_value_r = "correlation_pearson"
        type_value_p = "p_pearson"
    elif (type_correlation == "spearman"):
        type_value_r = "correlation_spearman"
        type_value_p = "p_spearman"
        pass

    # Calculate correlations.
    table_correlation_r = calculate_correlations_populate_table_wide(
        table=table,
        column_identifier_observation=column_identifier_observation,
        features_first=features_first,
        features_second=features_second,
        proportion_nonmissing_observations=(
            proportion_nonmissing_observations
        ),
        z_score=z_score,
        type_value=type_value_r,
        report=True,
    )
    table_correlation_p = calculate_correlations_populate_table_wide(
        table=table,
        column_identifier_observation=column_identifier_observation,
        features_first=features_first,
        features_second=features_second,
        proportion_nonmissing_observations=(
            proportion_nonmissing_observations
        ),
        z_score=z_score,
        type_value=type_value_p,
        report=True,
    )

    # Cluster sequence of features by their correlations.
    # Copy information.
    table_correlation_r_cluster = table_correlation_r.copy(deep=True)
    table_correlation_p_cluster = table_correlation_p.copy(deep=True)
    # Organize indices in table.
    table_correlation_r_cluster.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_correlation_r_cluster.set_index(
        ["feature_first"],
        append=False,
        drop=True,
        inplace=True,
    )
    # Cluster columns in table.
    table_correlation_r_cluster = porg.cluster_table_columns(
        table=table_correlation_r_cluster,
    )
    # Cluster rows in table.
    table_correlation_r_cluster = porg.cluster_table_rows(
        table=table_correlation_r_cluster,
    )
    # Organize indices in table.
    table_correlation_r_cluster.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )

    # Cluster columns and rows in table of p-values to match the sequence in
    # the table of correlations.
    # TODO: filter by indices of 'features_first' and 'features_second'

    # Bundle information.
    pail = dict()
    pail["table_correlation"] = table_correlation_r
    pail["table_p_value"] = table_correlation_p
    pail["table_correlation_cluster"] = table_correlation_r_cluster
    #pail["table_p_value_cluster"] = table_p_value_cluster

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
        print("table of correlations:")
        print(table_correlation_r)
        putly.print_terminal_partition(level=5)
        print("table of p-values:")
        print(table_correlation_p)
        putly.print_terminal_partition(level=5)
        print("clustered table of correlations:")
        print(table_correlation_r_cluster)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def create_write_plot_chart_heatmap(
    path_directory_parent=None,
    name_chart=None,
    table=None,
    name_index_columns=None,
    name_index_rows=None,
    threshold_minimum=None,
    threshold_maximum=None,
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

    Review: TCW; 21 October 2025

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
        ((threshold_minimum is None) or (threshold_maximum is None)) or
        (math.isnan(threshold_minimum) or math.isnan(threshold_maximum))
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
        value_maximum = round((numpy.nanmax(matrix) + round_offset), 3)
    else:
        value_minimum = threshold_minimum
        value_maximum = threshold_maximum
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
        resolution=150,
        path_directory=path_directory_parent,
    )

    # Return information.
    return figure


def manage_create_write_plot_charts(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    table_correlation=None,
    table_p_value=None,
    table_correlation_cluster=None,
    column_name_feature=None,
    name_correlations=None,
    type_correlation=None,
    threshold_minimum=None,
    threshold_maximum=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 28 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table_correlation = table_correlation.copy(deep=True)
    table_p_value = table_p_value.copy(deep=True)
    table_correlation_cluster = table_correlation_cluster.copy(deep=True)

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

    # Create and write plot charts for the table of correlations.
    create_write_plot_chart_heatmap(
        path_directory_parent=path_directory_charts,
        name_chart="correlations_sort",
        table=table_correlation,
        name_index_columns="feature_second",
        name_index_rows="feature_first",
        threshold_minimum=threshold_minimum,
        threshold_maximum=threshold_maximum,
        title_chart="",
        title_bar=title_bar,
        title_abscissa="Features First",
        title_ordinate="Features Second",
        report=report,
    )
    create_write_plot_chart_heatmap(
        path_directory_parent=path_directory_charts,
        name_chart="correlations_cluster",
        table=table_correlation_cluster,
        name_index_columns="feature_second",
        name_index_rows="feature_first",
        threshold_minimum=threshold_minimum,
        threshold_maximum=threshold_maximum,
        title_chart="",
        title_bar=title_bar,
        title_abscissa="Features First",
        title_ordinate="Features Second",
        report=report,
    )
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
# Call main procedure.


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
    cluster_features=None,
    plot_threshold_minimum=None,
    plot_threshold_maximum=None,
    report=None,
):
    """
    Function to execute module's main behavior.

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
        cluster_features (bool): whether to cluster the sequence of features in
            the tables and charts
        plot_threshold_minimum (float): minimal value to set as threshold for
            representation in the plot chart
        plot_threshold_maximum (float): maximal value to set as threshold for
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
        cluster_features=cluster_features,
        plot_threshold_minimum=plot_threshold_minimum,
        plot_threshold_maximum=plot_threshold_maximum,
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
    # Read source information from file.
    pail_source = read_source(
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
        report=pail_parameters["report"],
    )
    #pail_source["table_features_observations"]
    #pail_source["table_groups"]
    #pail_source["features"]["first"]
    #pail_source["features"]["second"]

    ##########
    # Report.
    if pail_parameters["report"]:
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
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_name_observation=(
            pail_parameters["column_name_observation"]
        ),
        features_first=pail_source["features"]["first"],
        features_second=pail_source["features"]["second"],
        table_groups_observations=pail_source["table_groups"],
        report=pail_parameters["report"],
    )

    ##########
    # Filter lists of features corresponding to columns in table by their
    # availability of values across rows in table corresponding to
    # observations.

    pail_features = filter_features_by_available_observations(
        table=table_filter,
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_name_observation=(
            pail_parameters["column_name_observation"]
        ),
        features_first=pail_source["features"]["first"],
        features_second=pail_source["features"]["second"],
        proportion_nonmissing_observations=(
            pail_parameters["proportion_nonmissing_observations"]
        ),
        report=pail_parameters["report"],
    )
    #pail_features["features_first"]
    #pail_features["features_second"]

    ##########
    # Calculate correlations.

    # Calculate correlations between pairs of features and organize these
    # within a table in long format.
    table_correlation_long = calculate_correlations_populate_table_long(
        table=table_filter,
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        features_first=pail_features["features_first"],
        features_second=pail_features["features_second"],
        proportion_nonmissing_observations=(
            pail_parameters["proportion_nonmissing_observations"]
        ),
        z_score=True,
        report=pail_parameters["report"],
    )

    # Calculate correlations between pairs of features and organize these
    # within a table in wide format.
    pail_correlation = calculate_correlations_organize_tables_wide(
        table=table_filter,
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        features_first=pail_features["features_first"],
        features_second=pail_features["features_second"],
        proportion_nonmissing_observations=(
            pail_parameters["proportion_nonmissing_observations"]
        ),
        z_score=True,
        type_correlation=pail_parameters["type_correlation"],
        cluster_features=pail_parameters["cluster_features"],
        report=pail_parameters["report"],
    )
    #pail["table_correlation"]
    #pail["table_p_value"]
    #pail["table_correlation_cluster"]

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Lists.
    pail_write_lists = dict()
    pail_write_lists["features_first"] = pail_features["features_first"]
    pail_write_lists["features_second"] = pail_features["features_second"]
    # Tables.
    pail_write_tables = dict()
    pail_write_tables["table_correlation_long"] = table_correlation_long
    pail_write_tables["table_correlation_wide"] = (
        pail_correlation["table_correlation"]
    )
    pail_write_tables["table_p_value_wide"] = pail_correlation["table_p_value"]
    pail_write_tables["table_correlation_wide_cluster"] = (
        pail_correlation["table_correlation_cluster"]
    )

    ##########
    # Write product information to file.

    # Define paths to directories.
    path_directory_lists = os.path.join(
        pail_parameters["path_directory_product"], "lists",
    )
    path_directory_tables = os.path.join(
        pail_parameters["path_directory_product"], "tables",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_lists,
    )
    putly.create_directories(
        path=path_directory_tables,
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
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        table_correlation=pail_correlation["table_correlation"],
        table_p_value=pail_correlation["table_p_value"],
        table_correlation_cluster=(
            pail_correlation["table_correlation_cluster"]
        ),
        column_name_feature="feature_first",
        name_correlations="blank",
        type_correlation=pail_parameters["type_correlation"],
        threshold_minimum=pail_parameters["plot_threshold_minimum"],
        threshold_maximum=pail_parameters["plot_threshold_maximum"],
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
    cluster_features = sys.argv[12]
    plot_threshold_minimum=sys.argv[13]
    plot_threshold_maximum=sys.argv[14]
    report = sys.argv[15]

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
        cluster_features=cluster_features,
        plot_threshold_minimum=plot_threshold_minimum,
        plot_threshold_maximum=plot_threshold_maximum,
        report=report,
    )

    pass



#
