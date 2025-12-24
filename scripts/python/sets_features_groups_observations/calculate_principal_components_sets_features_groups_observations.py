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
# Date, initialization: 30 September 2025
# Date, restoration (update, modification): 20 October 2025
# Date, restoration (update, modification): 30 September 2025
# Review: TCW; 20 October 2025
# Review: TCW; 30 September 2025
################################################################################
# Note

# The specialty of this Python script is to manage a general and versatile
# design of principal components analysis.



# TODO: TCW; 1 October 2025
# Implement heatmap (filter to < 100 features and PCs) for PC loadings

# TODO: TCW; 20 October 2025
# For all distinct sets of features, collect principal components within the
# main table.


# TODO: TCW; 14 November 2025
# Make the "table_features" optional. It only exists to provide name translations
# for the features, I think...


##########
# Note:

# Note: TCW; 30 September 2025
# Features of interest may optionally exist within the "table_observations" or
# within the "table_signals", with or without the need for a transposition
# before merging with "table_observations".


# Note: TCW; 11 December 2025
# Consider simplifying this procedure by removing all accommodation for
# supplemental features that are not already in the main table of features and
# observations. Leave that job to the separate procedure,
# "merge_features_main_supplement.py".


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
    path_file_source_table_features=None,
    path_file_source_table_sets_features=None,
    path_file_source_table_groups_observations=None,
    column_identifier_observation=None,
    column_name_observation=None,
    column_identifier_feature=None,
    column_name_feature=None,
    proportion_nonmissing_observations=None,
    allow_replicate_observations=None,
    count_components=None,
    variance_components=None,
    identifiers_emphasis=None,
    features_response_quantitative=None,
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
    pail["path_file_source_table_features_observations"] = str(
        path_file_source_table_features_observations
    ).strip()
    pail["path_file_source_table_features"] = str(
        path_file_source_table_features
    ).strip()
    pail["path_file_source_table_sets_features"] = str(
        path_file_source_table_sets_features
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
    pail["column_identifier_feature"] = str(
        column_identifier_feature
    ).strip()
    pail["column_name_feature"] = str(
        column_name_feature
    ).strip()

    # Number.
    # Iterate on individual of numeric parameters.
    numbers = {
        "proportion_nonmissing_observations": (
            proportion_nonmissing_observations
        ),
        "count_components": count_components,
        "variance_components": variance_components,
    }
    nulls = [
        "na", "nan", "NA", "NAN",
    ]
    for key_number in numbers.keys():
        # Determine whether parameter has a valid value.
        if (
            (numbers[key_number] is not None) and
            (len(str(numbers[key_number])) > 0) and
            (str(numbers[key_number]) != "") and
            (str(numbers[key_number]).strip().lower() != "none") and
            (str(numbers[key_number]).strip().lower() not in nulls)
        ):
            # Number is valid.
            pail[key_number] = float(str(numbers[key_number]).strip())
        else:
            # Number is missing or null.
            pail[key_number] = None
            pass
        pass

    # List.
    pail["identifiers_emphasis"] = putly.parse_text_list_values(
        text=identifiers_emphasis,
        delimiter=",",
    )
    pail["identifiers_emphasis"] = putly.collect_unique_items(
        items=pail["identifiers_emphasis"],
    )
    pail["features_response_quantitative"] = putly.parse_text_list_values(
        text=features_response_quantitative,
        delimiter=",",
    )
    pail["features_response_quantitative"] = putly.collect_unique_items(
        items=pail["features_response_quantitative"],
    )

    # Boolean, true or false.
    # Iterate on individual of Boolean designations.
    designations = {
        "allow_replicate_observations": allow_replicate_observations,
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
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
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
    path_file_source_table_features=None,
    path_file_source_table_sets_features=None,
    path_file_source_table_groups_observations=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review or revision: TCW; 16 December 2025
    Review or revision: TCW; 29 September 2025

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

    # Determine whether paths point to a directory or file that exist.
    #existence_directory = os.path.exists(path_directory)
    existence_file_features = os.path.exists(path_file_source_table_features)

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
    if (existence_file_features):
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
    else:
        pail["table_features"] = None
        pass
    pail["table_sets_features"] = pandas.read_csv(
        path_file_source_table_sets_features,
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
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def filter_organize_table_features_observations(
    table=None,
    column_identifier_observation=None,
    column_name_observation=None,
    column_identifier_groups_observations=None,
    columns_categories=None,
    features_selection=None,
    features_response_quantitative=None,
    observations_selection=None,
    filter_table_main=None,
    report=None,
):
    """
    Blank.

    Review or revision: TCW; 16 December 2025
    Review or revision: TCW; 2 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)
    columns_categories = copy.deepcopy(columns_categories)
    features_selection = copy.deepcopy(features_selection)
    features_response_quantitative = copy.deepcopy(
        features_response_quantitative
    )
    observations_selection = copy.deepcopy(observations_selection)

    # Filter columns and rows in main table for specific features and
    # observations.
    # Notice that the function below already includes the name of the column
    # for the index across rows.
    if (filter_table_main):
        # Copy information.
        categories_features_selection = copy.deepcopy(columns_categories)
        # Prepare inclusive list of columns.
        categories_features_selection.insert(0, column_name_observation)
        categories_features_selection.insert(0, column_identifier_observation)
        categories_features_selection.insert(
            0, column_identifier_groups_observations
        )
        categories_features_selection.extend(features_selection)
        categories_features_selection.extend(features_response_quantitative)
        # Filter columns and rows in table.
        table_selection = (
            porg.filter_select_table_columns_rows_by_identifiers(
                table=table,
                index_rows=column_identifier_groups_observations,
                identifiers_columns=categories_features_selection,
                identifiers_rows=observations_selection,
                report=False,
        ))
        # Filter rows in table for non-missing values across relevant columns.
        table_filter.dropna(
            axis="index",
            how="all",
            subset=features_selection,
            inplace=True,
        )
    else:
        # Copy information.
        table_selection = table.copy(deep=True)
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
        print("function: filter_organize_table_features_observations()")
        putly.print_terminal_partition(level=5)
        print("table before filtering columns and rows")
        print(table)
        print("table after filtering columns and rows")
        print(table_selection)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_selection


def filter_features_by_variance_available_observations(
    table=None,
    column_identifier_observation=None,
    column_name_observation=None,
    features_selection=None,
    features_response_quantitative=None,
    sets_features=None,
    observations_selection=None,
    proportion_nonmissing_observations=None,
    report=None,
):
    """
    Blank.

    Review or revision: TCW; 16 December 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    ##########
    # Copy information.
    table = table.copy(deep=True)
    features_selection = copy.deepcopy(features_selection)
    features_response_quantitative = copy.deepcopy(
        features_response_quantitative
    )
    sets_features = copy.deepcopy(sets_features)
    observations_selection = copy.deepcopy(observations_selection)

    # Assemble list of features.
    features_selection_total = list()
    features_selection_total.extend(features_selection)
    features_selection_total.extend(features_response_quantitative)
    features_selection_total = putly.collect_unique_items(
        items=features_selection_total,
    )

    # Filter table's columns corresponding to features by the availability of
    # observations across rows in the table.
    table_filter = (
        porg.filter_table_columns_by_proportion_nonmissing_threshold(
            table=table,
            index_columns="features",
            index_rows=column_identifier_observation,
            columns_selection=features_selection_total,
            rows_selection=observations_selection,
            threshold_low=None,
            threshold_high=None,
            proportion=proportion_nonmissing_observations,
            report=report,
    ))

    # Filter table's columns corresponding to features by the corresponding
    # variance across observations in rows of the table.
    table_filter = (
        porg.filter_table_columns_by_variance_threshold(
            table=table_filter,
            index_columns="features",
            index_rows=column_identifier_observation,
            columns_selection=features_selection_total,
            rows_selection=observations_selection,
            measure_variance="coefficient_variation",
            threshold_variance=0.001,
            report=report,
    ))

    # Extract identifiers of columns in table.
    features_available = copy.deepcopy(
        table_filter.columns.unique().tolist()
    )

    # Filter sets of features by the availability of nonmissing observations
    # across rows for corresponding columns in table.
    sets_features_filter = dict()
    features_selection_filter = list()
    for name_set in sets_features.keys():
        features_set = copy.deepcopy(sets_features[name_set])
        sets_features_filter[name_set] = list(filter(
            lambda feature: (feature in features_available),
            features_set
        ))
        features_selection_filter.extend(sets_features_filter[name_set])
        pass
    # Filter features for quantitative responses.
    features_response_quantitative_filter = list(filter(
        lambda feature: (feature in features_available),
        features_response_quantitative
    ))
    # Collect unique features.
    features_selection_filter = putly.collect_unique_items(
        items=features_selection_filter,
    )

    # Bundle information.
    pail = dict()
    pail["features_selection"] = features_selection_filter
    pail["features_response_quantitative"] = (
        features_response_quantitative_filter
    )
    pail["sets_features"] = sets_features_filter

    # Report.
    if report:
        # Organize information.
        count = len(features_selection_filter)
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_tables_plot_charts_features_correlations.py"
        )
        print(str("module: " + module))
        print("function: filter_features_by_available_observations()")
        putly.print_terminal_partition(level=5)
        print("count of selection features: " + str(count))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def create_write_plot_chart_heatmap_loadings(
    path_directory_parent=None,
    name_chart=None,
    table=None,
    name_index_columns=None,
    name_index_rows=None,
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
        transpose_table=False,
        fill_missing=True,
        value_missing_fill=0.0,
        constrain_signal_values=True,
        value_minimum=value_minimum,
        value_center=0.0, # correlations, effects, fold changes, components etc
        value_maximum=value_maximum,
        title_ordinate="",
        title_abscissa="",
        title_bar=title_bar,
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        size_title_ordinate="eight",
        size_title_abscissa="eight",
        size_title_bar="ten",
        size_label_ordinate=None, # determine automatically if "None"; "fifteen"
        size_label_abscissa=None, # determine automatically if "None"
        size_label_bar="twelve",
        show_labels_ordinate=True,
        show_labels_abscissa=True,
        show_scale_bar=True,
        aspect="portrait", # square, portrait, landscape, ...
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


def create_write_plot_charts_principal_component_loadings(
    path_directory_parent=None,
    table_loadings=None,
    column_identifier_feature=None,
    column_name_feature=None,
    columns_components=None,
    name_set_features=None,
    cluster_features=None,
    write_charts=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 20 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table_loadings = table_loadings.copy(deep=True)
    columns_components = copy.deepcopy(columns_components)

    # Define prefix for names of principal components.
    prefix_name_components = str(
        name_set_features + "_pc_"
    )

    # Determine translations for names of principal components.
    translations_components = dict()
    for column_component in columns_components:
        count = str(column_component).replace(prefix_name_components, "")
        translation = str("PC-" + str(count))
        translations_components[column_component] = translation
        pass

    # Translate names of columns in table.
    table_loadings_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=table_loadings,
            index_rows=column_identifier_feature,
            translations_columns=translations_components,
            translations_rows=None,
            remove_redundancy=False,
            report=False,
    ))

    # Cluster rows in table.
    if (cluster_features):
        # Organize indices in table.
        table_loadings_translation.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table_loadings_translation.set_index(
            [column_identifier_feature],
            append=False,
            drop=True,
            inplace=True,
        )
        # Cluster rows in table.
        table_loadings_translation = porg.cluster_table_rows(
            table=table_loadings_translation,
        )
        # Organize indices in table.
        table_loadings_translation.reset_index(
            level=None,
            inplace=True,
            drop=False, # remove index; do not move to regular columns
        )
        pass

    # Determine name for chart.
    name_chart = str("chart_" + name_set_features + "_loadings")

    # Determine titles for abscissa and ordinate axes.
    title_abscissa = str("Principal Components (" + name_set_features + ")")
    title_ordinate = str("Features")

    ##########
    # Plot chart: scatter with colors for categorical groups of
    # observations.
    # Define paths to directories.
    path_directory_chart = os.path.join(
        path_directory_parent, "heatmap",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_chart,
    )
    # Create plot chart and write to file.
    create_write_plot_chart_heatmap_loadings(
        path_directory_parent=path_directory_chart,
        name_chart=name_chart,
        table=table_loadings_translation,
        name_index_columns="components",
        name_index_rows=column_identifier_feature,
        title_chart="",
        title_bar="Loading Weight",
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
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
        function = str(
            "create_write_plot_charts_principal_component_loadings()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass

    pass


def create_write_plot_charts_principal_component_scores(
    path_directory_parent=None,
    table_features_observations=None,
    table_scores=None,
    column_identifier_observation=None,
    column_name_observation=None,
    columns_categories=None,
    columns_components=None,
    name_set_features=None,
    features_response_quantitative=None,
    observations_selection=None,
    groups_observations=None,
    names_groups_observations_sequence=None,
    translations_observations=None,
    identifiers_emphasis=None,
    write_charts=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 1 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table_features_observations = table_features_observations.copy(deep=True)
    table_scores = table_scores.copy(deep=True)
    columns_categories = copy.deepcopy(columns_categories)
    columns_components = copy.deepcopy(columns_components)
    features_response_quantitative = copy.deepcopy(
        features_response_quantitative
    )
    observations_selection = copy.deepcopy(observations_selection)
    groups_observations = copy.deepcopy(groups_observations)
    names_groups_observations_sequence = copy.deepcopy(
        names_groups_observations_sequence
    )
    translations_observations = copy.deepcopy(translations_observations)
    identifiers_emphasis = copy.deepcopy(identifiers_emphasis)

    # Determine translations for names of principal components.
    translations_components = dict()
    for column_component in columns_components:
        name_prefix = str(name_set_features + "_")
        translation = str(column_component).replace(name_prefix, "")
        translations_components[column_component] = translation
        pass

    # Prepare reference to sort groups of rows in table.
    sequence_groups_observations = dict()
    index = 0
    for name in names_groups_observations_sequence:
        sequence_groups_observations[name] = index
        index += 1
        pass

    # Filter columns and rows in table.
    columns_sequence = [
        column_identifier_observation,
        column_name_observation,
    ]
    columns_sequence.extend(columns_categories)
    columns_sequence.extend(columns_components)
    columns_sequence.extend(features_response_quantitative)
    table_selection = porg.filter_select_table_columns_rows_by_identifiers(
        table=table_features_observations,
        index_rows=column_identifier_observation,
        identifiers_columns=columns_sequence,
        identifiers_rows=observations_selection,
        report=False,
    )

    # Filter rows in table for non-missing values across relevant columns.
    table_selection.dropna(
        axis="index",
        how="any",
        subset=columns_components,
        inplace=True,
    )

    # Determine and fill groups of observations.
    # These next two functions create a new column named "group" and assign
    # categorical names corresponding to specific groups of observations.
    # Each observation can only belong to a single group.
    table_group = porg.determine_fill_table_groups_rows(
        table=table_selection,
        index_rows=column_identifier_observation,
        column_group="group",
        groups_rows=groups_observations,
        report=False,
    )
    # Sort rows in table by groups.
    table_group = porg.sort_table_rows_by_single_column_reference(
        table=table_group,
        index_rows=column_identifier_observation,
        column_reference="group",
        column_sort_temporary="sort_temporary",
        reference_sort=sequence_groups_observations,
    )

    # Define comparisons between principal components.
    comparisons = [
        ["1", "2",],
        ["1", "3",],
        ["1", "4",],
        ["1", "5",],
        ["2", "3",],
        ["2", "4",],
        ["2", "5",],
        ["3", "4",],
        ["3", "5",],
    ]
    # Iterate on comparisons.
    for comparison in comparisons:
        # Determine name for chart.
        name_chart = str(
            "chart_" +
            name_set_features +
            "_pc_" + comparison[0] + "_" + comparison[1]
        )
        # Determine names of columns for values to represent on abscissa and
        # ordinate axes.
        column_abscissa = str(name_set_features + "_pc_" + comparison[0])
        column_ordinate = str(name_set_features + "_pc_" + comparison[1])
        title_abscissa = str("PC-" + comparison[0])
        title_ordinate = str("PC-" + comparison[1])

        ##########
        # Plot chart: scatter with colors for categorical groups of
        # observations.
        # Define paths to directories.
        path_directory_chart = os.path.join(
            path_directory_parent, "scatter_categories",
        )
        # Create directories.
        putly.create_directories(
            path=path_directory_chart,
        )
        # Create plot chart and write to file.
        splot.create_write_plot_chart_scatter_point_response(
            path_directory_parent=path_directory_chart,
            name_chart=name_chart,
            table=table_group,
            column_identifier=column_identifier_observation,
            column_name=column_name_observation,
            column_response_markers="group",
            column_group_ellipses="group",
            column_abscissa=column_abscissa,
            column_ordinate=column_ordinate,
            type_response="category",
            title_chart="",
            title_response="group",
            title_abscissa=title_abscissa,
            title_ordinate=title_ordinate,
            #identifiers_emphasis=list(),
            identifiers_emphasis=identifiers_emphasis,
            size_marker=15,
            size_edge_marker=1.0,
            size_edge_ellipse=1.0,
            factor_confidence_ellipse=2.0,
            colors_fill_markers=None,
            colors_fill_ellipses=(0.500,0.500,0.500,1.0),
            color_edge_markers="black",
            color_edge_ellipses="black",
            color_emphasis="orange",
            show_confidence_ellipse=True,
            show_emphasis_marker=True,
            show_emphasis_label=True,
            show_legend_bar=True,
            report=report,
        )

        ##########
        # Plot chart: scatter with color gradient for response feature on a
        # quantitative scale.
        for feature in features_response_quantitative:
            # Define paths to directories.
            path_directory_chart = os.path.join(
                path_directory_parent, "scatter_response", feature,
            )
            # Create directories.
            putly.create_directories(
                path=path_directory_chart,
            )
            # Create plot chart and write to file.
            splot.create_write_plot_chart_scatter_point_response(
                path_directory_parent=path_directory_chart,
                name_chart=name_chart,
                table=table_group,
                column_identifier=column_identifier_observation,
                column_name=column_name_observation,
                column_response_markers=feature,
                column_group_ellipses="group",
                column_abscissa=column_abscissa,
                column_ordinate=column_ordinate,
                type_response="continuity",
                title_chart="",
                title_response=feature,
                title_abscissa=title_abscissa,
                title_ordinate=title_ordinate,
                #identifiers_emphasis=list(),
                identifiers_emphasis=identifiers_emphasis,
                size_marker=15,
                size_edge_marker=1.0,
                size_edge_ellipse=1.0,
                factor_confidence_ellipse=2.0,
                colors_fill_markers=None,
                colors_fill_ellipses=(0.500,0.500,0.500,1.0),
                color_edge_markers="black",
                color_edge_ellipses="black",
                color_emphasis="orange",
                show_confidence_ellipse=True,
                show_emphasis_marker=True,
                show_emphasis_label=True,
                show_legend_bar=True,
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
        print("function: create_write_plot_charts()")
        putly.print_terminal_partition(level=5)
        pass

    pass


def manage_create_write_plot_charts(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    table_features_observations=None,
    table_loadings=None,
    table_variances=None,
    table_scores=None,
    column_identifier_observation=None,
    column_name_observation=None,
    columns_categories=None,
    columns_components=None,
    name_set_features=None,
    features_selection=None,
    features_response_quantitative=None,
    observations_selection=None,
    groups_observations=None,
    names_groups_observations_sequence=None,
    translations_features=None,
    translations_observations=None,
    identifiers_emphasis=None,
    write_charts=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 20 October 2025
    Review: TCW; 1 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table_features_observations = table_features_observations.copy(deep=True)
    table_loadings = table_loadings.copy(deep=True)
    table_variances = table_variances.copy(deep=True)
    table_scores = table_scores.copy(deep=True)
    columns_categories = copy.deepcopy(columns_categories)
    columns_components = copy.deepcopy(columns_components)
    features_selection = copy.deepcopy(features_selection)
    features_response_quantitative = copy.deepcopy(
        features_response_quantitative
    )
    observations_selection = copy.deepcopy(observations_selection)
    groups_observations = copy.deepcopy(groups_observations)
    names_groups_observations_sequence = copy.deepcopy(
        names_groups_observations_sequence
    )
    translations_features = copy.deepcopy(translations_features)
    translations_observations = copy.deepcopy(translations_observations)
    identifiers_emphasis = copy.deepcopy(identifiers_emphasis)

    ##########
    # Write product information to file.
    # Define paths to directories.
    path_directory_instance = os.path.join(
        path_directory_product, name_set_features,
    )
    path_directory_instance_charts = os.path.join(
        path_directory_instance, "charts",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_instance_charts,
    )

    # Create and write plot charts for the main table of features and
    # observations, including principal components as features.
    create_write_plot_charts_principal_component_scores(
        path_directory_parent=path_directory_instance_charts,
        table_features_observations=table_features_observations,
        table_scores=table_scores,
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        columns_categories=columns_categories,
        columns_components=columns_components,
        name_set_features=name_set_features,
        features_response_quantitative=features_response_quantitative,
        observations_selection=observations_selection,
        groups_observations=groups_observations,
        names_groups_observations_sequence=names_groups_observations_sequence,
        translations_observations=translations_observations,
        identifiers_emphasis=identifiers_emphasis,
        write_charts=write_charts,
        report=report,
    )

    # Create and write plot charts for the table of principal component
    # loadings.
    create_write_plot_charts_principal_component_loadings(
        path_directory_parent=path_directory_instance_charts,
        table_loadings=table_loadings,
        column_identifier_feature="features",
        column_name_feature="features",
        columns_components=columns_components,
        name_set_features=name_set_features,
        cluster_features=True,
        write_charts=write_charts,
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


def manage_components_set_features_groups_observations(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    table=None,
    column_identifier_observation=None,
    column_name_observation=None,
    columns_categories=None,
    name_set_features=None,
    features_selection=None,
    observations_selection=None,
    groups_observations=None,
    names_groups_observations_sequence=None,
    translations_features=None,
    translations_observations=None,
    count_components=None,
    variance_components=None,
    identifiers_emphasis=None,
    features_response_quantitative=None,
    allow_replicate_observations=None,
    write_lists=None,
    write_tables=None,
    write_charts=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 20 October 2025

    arguments:

    TODO: update documentation

        count_components (int): count of principal components to keep
        variance_components (float): threshold on the cumulative proportion of
            total variance to determine how many principal components to keep


    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)
    columns_categories = copy.deepcopy(columns_categories)
    features_selection = copy.deepcopy(features_selection)
    observations_selection = copy.deepcopy(observations_selection)
    groups_observations = copy.deepcopy(groups_observations)
    names_groups_observations_sequence = copy.deepcopy(
        names_groups_observations_sequence
    )
    translations_features = copy.deepcopy(translations_features)
    translations_observations = copy.deepcopy(translations_observations)

    # Define prefix for names of principal components.
    prefix_name_components = str(
        name_set_features + "_" + "pc"
    )

    # Calculate principal components and combine with table of main features
    # and observations.
    pail_components = (
        pdcmp.calculate_principal_components_table_features_observations(
            table=table,
            name_index_columns="features",
            name_index_rows=column_identifier_observation,
            columns_selection=features_selection,
            rows_selection=observations_selection,
            prefix=prefix_name_components,
            separator="_",
            count_components=count_components,
            variance_components=variance_components, # cumulative proportion of variance; float or None to keep all
            explicate_indices=False,
            report=report,
    ))
    #pail_components["table_merge"] = table_merge
    #pail_components["table_loadings"] = table_loadings
    #pail_components["table_variances"] = table_variances
    #pail_components["table_scores"] = table_scores
    #pail_components["columns_scores"] = columns_scores

    # Simplify main table of features and observations.
    columns_all = copy.deepcopy(
        pail_components["table_merge"].columns.to_list()
    )
    columns_sequence = list(filter(
        lambda item: item not in features_selection, columns_all
    ))
    # Filter and sort columns in table.
    if False:
        table_features_observations = porg.filter_sort_table_columns(
            table=pail_components["table_merge"],
            columns_sequence=columns_sequence,
            report=report,
        )
    else:
        table_features_observations = pail_components["table_merge"].copy(
            deep=True,
        )
        pass

    # Translate names of features and observations.
    table_scores_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=pail_components["table_scores"],
            index_rows=column_identifier_observation,
            translations_columns=None,
            translations_rows=translations_observations,
            remove_redundancy=False,
            report=False,
    ))
    table_loadings_translation = (
        porg.translate_identifiers_table_indices_columns_rows(
            table=pail_components["table_loadings"],
            index_rows="features",
            translations_columns=None,
            translations_rows=translations_features,
            remove_redundancy=False,
            report=False,
    ))

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Lists.
    pail_write_lists = dict()
    pail_write_lists["columns_scores"] = (
        pail_components["columns_scores"]
    )
    # Tables.
    pail_write_tables = dict()
    pail_write_tables["table_loadings"] = table_loadings_translation
    pail_write_tables["table_variances"] = (
        pail_components["table_variances"]
    )
    pail_write_tables["table_scores"] = table_scores_translation

    ##########
    # Write product information to file.

    # Define paths to directories.
    path_directory_instance = os.path.join(
        path_directory_product, name_set_features,
    )
    path_directory_instance_lists = os.path.join(
        path_directory_instance, "lists",
    )
    path_directory_instance_tables = os.path.join(
        path_directory_instance, "tables",
    )
    path_directory_instance_charts = os.path.join(
        path_directory_instance, "charts",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_instance_lists,
    )
    putly.create_directories(
        path=path_directory_instance_tables,
    )
    putly.create_directories(
        path=path_directory_instance_charts,
    )
    # Lists.
    if (write_lists):
        putly.write_lists_to_file_text(
            pail_write=pail_write_lists,
            path_directory=path_directory_instance_lists,
            delimiter="\n",
        )
        pass
    # Tables.
    if (write_tables):
        putly.write_tables_to_file(
            pail_write=pail_write_tables,
            path_directory=path_directory_instance_tables,
            reset_index_rows=False,
            write_index_rows=False,
            write_index_columns=True,
            type="text",
            delimiter="\t",
            suffix=".tsv",
        )
        pass

    # Manage the creation and write of plot charts.
    # Calculate principal components.
    manage_create_write_plot_charts(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        table_features_observations=table_features_observations,
        table_loadings=table_loadings_translation,
        table_variances=pail_components["table_variances"],
        table_scores=pail_components["table_scores"],
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        columns_categories=columns_categories,
        columns_components=pail_components["columns_scores"],
        name_set_features=name_set_features,
        features_selection=features_selection,
        features_response_quantitative=features_response_quantitative,
        observations_selection=observations_selection,
        groups_observations=groups_observations,
        names_groups_observations_sequence=names_groups_observations_sequence,
        translations_features=translations_features,
        translations_observations=translations_observations,
        identifiers_emphasis=identifiers_emphasis,
        write_charts=write_charts,
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
        print("function: manage_components_set_features_groups_observations()")
        putly.print_terminal_partition(level=5)
        pass


    # Bundle information.
    pail_return = dict()
    pail_return["table"] = table_features_observations
    # Return information.
    return pail_return


################################################################################
# Procedure


##########
# Call main procedure.


def execute_procedure(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_features_observations=None,
    path_file_source_table_features=None,
    path_file_source_table_sets_features=None,
    path_file_source_table_groups_observations=None,
    column_identifier_observation=None,
    column_name_observation=None,
    column_identifier_feature=None,
    column_name_feature=None,
    proportion_nonmissing_observations=None,
    allow_replicate_observations=None,
    count_components=None,
    variance_components=None,
    identifiers_emphasis=None,
    features_response_quantitative=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review or revision: TCW; 15 December 2025
    Review or revision: TCW; 2 October 2025

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files

        path_file_source_table_features_observations (str): path to source file
        path_file_source_table_features (str): path to source file
        path_file_source_table_sets_features (str): path to source file
        path_file_source_table_groups_observations (str): path to source file
        column_identifier_observation (str): name of column in source table
        column_name_observation (str): name of column in source table
        column_identifier_feature (str): name of column in source table
        column_name_feature (str): name of column in source table
        proportion_nonmissing_observations (float): threshold by proportion of
            observations that must have nonmissing values for each feature
        allow_replicate_observations (bool): whether to allow replicate
            observations or to require groups to be mutually exclusive, such
            that any individual observation can only belong to one group
        count_components (int): count of principal components to keep for each
            set of features
        variance_components (float):
        identifiers_emphasis (list<str>):
        features_response_quantitative (list<str>):
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
        path_file_source_table_features=path_file_source_table_features,
        path_file_source_table_sets_features=(
            path_file_source_table_sets_features
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        column_identifier_feature=column_identifier_feature,
        column_name_feature=column_name_feature,
        proportion_nonmissing_observations=proportion_nonmissing_observations,
        allow_replicate_observations=allow_replicate_observations,
        count_components=count_components,
        variance_components=variance_components,
        identifiers_emphasis=identifiers_emphasis,
        features_response_quantitative=features_response_quantitative,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
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
        path_file_source_table_features=(
            pail_parameters["path_file_source_table_features"]
        ),
        path_file_source_table_sets_features=(
            pail_parameters["path_file_source_table_sets_features"]
        ),
        path_file_source_table_groups_observations=(
            pail_parameters["path_file_source_table_groups_observations"]
        ),
        report=pail_parameters["report"],
    )

    # Parameters.

    # Determine preliminary selection of features with available observations.
    # This preliminary selection of features does not consider their respective
    # proportions of observations with nonmissing values.
    features_available = sutly.determine_features_available(
        table_features=pail_source["table_features"],
        table_observations=pail_source["table_features_observations"],
        table_signals=None,
        column_identifier_feature=pail_parameters["column_identifier_feature"],
        column_name_feature=pail_parameters["column_name_feature"],
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_name_observation=pail_parameters["column_name_observation"],
        column_identifier_signal=None,
        transpose_table_signals=None,
        report=pail_parameters["report"],
    )

    # Read and organize parameters about sets of features.
    pail_sets = sutly.read_organize_parameters_sets_features(
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        table=pail_source["table_sets_features"],
        column_name="abbreviation",
        features_available=features_available,
        report=pail_parameters["report"],
    )
    #pail_sets["table"]
    #pail_sets["names_sets_features_sequence"]
    #pail_sets["sets_features"]
    #pail_sets["features_sets_union"]
    #pail_sets["records"]
    pail_features = sutly.organize_parameters_further_sets_features(
        table_features=pail_source["table_features"],
        column_identifier_feature=(
            pail_parameters["column_identifier_feature"]
        ),
        column_name_feature=(
            pail_parameters["column_name_feature"]
        ),
        features_selection=pail_sets["features_sets_union"],
        features_sets_union=pail_sets["features_sets_union"],
        sets_features=pail_sets["sets_features"],
        names_sets_features_sequence=pail_sets["names_sets_features_sequence"],
        prefix_name_feature="",
        report=pail_parameters["report"],
    )
    #pail_features["table_features_selection"]
    #pail_features["features_selection"]
    #pail_features["features_selection_translation"]
    #pail_features["features_selection_prefix"]
    #pail_features["translations_features"]
    #pail_features["names_sets_features_sequence"]
    #pail_features["sets_features"]

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
        column_identifier_groups_observations=(
            pail_parameters["column_identifier_observation"]
        ),
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

    # Filter sets of features corresponding to columns in table by their
    # availability of values across rows in table corresponding to
    # observations.
    pail_features_availability = (
        filter_features_by_variance_available_observations(
            table=pail_source["table_features_observations"],
            column_identifier_observation=column_identifier_observation,
            column_name_observation=column_name_observation,
            features_selection=pail_features["features_selection"],
            features_response_quantitative=(
                pail_parameters["features_response_quantitative"]
            ),
            sets_features=pail_features["sets_features"],
            observations_selection=pail_observations["observations_selection"],
            proportion_nonmissing_observations=(
                pail_parameters["proportion_nonmissing_observations"]
            ),
            report=report,
    ))
    #pail_features_availability["features_selection"]
    #pail_features_availability["features_response_quantitative"]
    #pail_features_availability["sets_features"]

    # Filter features and observations in main table.
    table_features_observations = filter_organize_table_features_observations(
        table=pail_source["table_features_observations"],
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_name_observation=pail_parameters["column_name_observation"],
        column_identifier_groups_observations=(
            pail_parameters["column_identifier_observation"]
        ),
        columns_categories=pail_groups["categories_groups"],
        features_selection=pail_features_availability["features_selection"],
        features_response_quantitative=(
            pail_features_availability["features_response_quantitative"]
        ),
        observations_selection=pail_observations["observations_selection"],
        filter_table_main=False,
        report=pail_parameters["report"],
    )

    # Iterate on sets of features.
    for name_set_features in pail_features["names_sets_features_sequence"]:

        # Access set of features for current instance.
        features_set = copy.deepcopy(
            pail_features_availability["sets_features"][name_set_features]
        )

        # Calculate principal components.
        pail_instance = manage_components_set_features_groups_observations(
            path_directory_source=pail_parameters["path_directory_source"],
            path_directory_product=pail_parameters["path_directory_product"],
            path_directory_dock=pail_parameters["path_directory_dock"],
            table=table_features_observations,
            column_identifier_observation=(
                pail_parameters["column_identifier_observation"]
            ),
            column_name_observation=(
                pail_parameters["column_name_observation"]
            ),
            columns_categories=pail_groups["categories_groups"],
            name_set_features=name_set_features,
            features_selection=features_set,
            observations_selection=pail_observations["observations_selection"],
            groups_observations=pail_observations["groups_observations"],
            names_groups_observations_sequence=(
                pail_observations["names_groups_observations_sequence"]
            ),
            translations_features=(
                pail_features["translations_features"]
            ),
            translations_observations=(
                pail_observations["translations_observations"]
            ),
            count_components=(pail_parameters["count_components"]),
            variance_components=(
                pail_parameters["variance_components"]
            ),
            identifiers_emphasis=(
                pail_parameters["identifiers_emphasis"]
            ),
            features_response_quantitative=(
                pail_features_availability["features_response_quantitative"]
            ),
            allow_replicate_observations=(
                pail_parameters["allow_replicate_observations"]
            ),
            write_lists=True,
            write_tables=True,
            write_charts=True,
            report=pail_parameters["report"],
        )
        # Update reference to table of features and observations.
        # Collect principal components for all sets of features.
        # Copy information.
        table_features_observations = pail_instance["table"].copy(deep=True)
        pass

    # Review and write the table after introduction of all relevant PCs.

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Lists.
    pail_write_lists = dict()
    # Tables.
    pail_write_tables = dict()
    pail_write_tables["table_merge"] = table_features_observations

    ##########
    # Write product information to file.

    # Define paths to directories.
    path_directory_tables = os.path.join(
        path_directory_product, "tables_combination",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_tables,
    )
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
    path_directory_source = sys.argv[1]
    path_directory_product = sys.argv[2]
    path_directory_dock = sys.argv[3]
    path_file_source_table_features_observations = sys.argv[4]
    path_file_source_table_features = sys.argv[5]
    path_file_source_table_sets_features = sys.argv[6]
    path_file_source_table_groups_observations = sys.argv[7]
    column_identifier_observation = sys.argv[8]
    column_name_observation = sys.argv[9]
    column_identifier_feature = sys.argv[10]
    column_name_feature = sys.argv[11]
    proportion_nonmissing_observations = sys.argv[12]
    allow_replicate_observations = sys.argv[13]
    count_components = sys.argv[14]
    variance_components = sys.argv[15]
    identifiers_emphasis = sys.argv[16]
    features_response_quantitative = sys.argv[17]
    report = sys.argv[18]

    # Call function for procedure.
    execute_procedure(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        path_file_source_table_features=path_file_source_table_features,
        path_file_source_table_sets_features=(
            path_file_source_table_sets_features
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        column_identifier_feature=column_identifier_feature,
        column_name_feature=column_name_feature,
        proportion_nonmissing_observations=proportion_nonmissing_observations,
        allow_replicate_observations=allow_replicate_observations,
        count_components=count_components,
        variance_components=variance_components,
        identifiers_emphasis=identifiers_emphasis,
        features_response_quantitative=features_response_quantitative,
        report=report,
    )

    pass



#
