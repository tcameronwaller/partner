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
# Date, first execution: 2 May 2025
# Date, last execution or modification: 23 July 2025
# Review: TCW; 23 July 2025
################################################################################
# Note

# The specialty of this Python script is to drive the creation of a single plot
# chart for visual representation of information such as coefficients from
# regression analyses. This script calls versatile functionality from the
# "plot.py" module within the "partner" Python package.

# Find a table of artificial data for demonstration in the file path below.
# /.../partner/repository/demonstration/table_dot_plot_forest.tsv

##########
# Review: TCW;

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


##########
# Read source information from file.


def parse_values_intervals(
    text=None,
):
    """
    Parse information from text.

    arguments:
        text (str): text parameters

    raises:

    returns:
        (dict): collection of information
    """

    # Parse and extract information from text.
    if (
        (text is not None) and
        (len(str(text)) > 0) and
        (str(text) != "") and
        (str(text).strip().lower() != "none")
    ):
        extraction = (
            putly.parse_extract_text_keys_values_semicolon_colon_comma(
                text=text,
            )
        )["features_values"]
        values_intervals_parse = dict()
        for key in extraction.keys():
            values_intervals_parse[key] = copy.deepcopy(
                extraction[key][0]
            )
            pass
    else:
        values_intervals_parse = None
        pass
    # Return information.
    return values_intervals_parse


def parse_text_parameters(
    path_file_source_table=None,
    name_file_product_chart=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    feature=None,
    features=None,
    translation_features=None,
    values_intervals_primary=None,
    values_intervals_secondary=None,
    values_intervals_tertiary=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    legend_series_primary=None,
    legend_series_secondary=None,
    legend_series_tertiary=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_file_source_table (str): path to source file in text format as a
            table with tab delimiters between columns and newline delimiters
            between rows, with data for creation of plot chart
        name_file_product_chart (str): name of product file
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        feature (str): parameters for extraction of features
        features (str): parameters for extraction of features
        translation_features (str): parameters for extraction of features
        values_intervals_primary (str): parameters for extraction of values
            and intervals
        values_intervals_secondary (str): parameters for extraction of values
            and intervals
        values_intervals_tertiary (str): parameters for extraction of values
            and intervals
        minimum_abscissa (float): value for minimal limit to represent on the
            abscissa horizontal axis
        maximum_abscissa (float): value for maximal limit to represent on the
            abscissa horizontal axis
        title_chart (str): title for chart figure
        title_abscissa (str): title for abscissa horizontal axis
        title_ordinate (str): title for ordinate vertical axis
        legend_series_primary (str): description in legend for primary series
        legend_series_secondary (str): description in legend for secondary
            series
        legend_series_tertiary (str): description in legend for tertiary series
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information
            columns_source (list<str>): names of relevant columns in table
            columns_continuity_source (list<str>): names of relevant columns in
                table
            columns_product (list<str>): names of relevant columns in table
            columns_continuity_product (list<str>): names of relevant columns
                in table
            translation_columns (dict<str>): translations of names of columns
                in table
            translation_features (dict<str>): translations of names or
                identifiers of features in rows of table
            feature (dict<str>): translation of column in table for names or
                identifiers of features
            feature_source (str): original, source name of column in table for
                names or identifiers of features
            feature_product (str): novel, product name of column in table for
                names or identifiers of features
            features (list<str>): names or identifiers of features to include
                for visual representation on a plot chart
            sequence_features (dict<int>): indices for sort that match
                categorical values in column of table
            values_intervals_primary (dict<str>): translations of columns in
                table for values and their intervals
            values_intervals_secondary (dict<str>): translations of columns in
                table for values and their intervals
            values_intervals_tertiary (dict<str>): translations of columns in
                table for values and their intervals
            report (bool): whether to print reports
    """

    # Bundle information.
    pail = dict()

    # Parse and extract information from text.

    ##########
    # Paths to directories and files.
    pail["path_file_source_table"] = str(
        path_file_source_table
    ).strip()
    pail["name_file_product_chart"] = str(
        name_file_product_chart
    ).strip()
    pail["path_directory_source"] = str(path_directory_source).strip()
    pail["path_directory_product"] = str(path_directory_product).strip()
    pail["path_directory_dock"] = str(path_directory_dock).strip()

    ##########
    # Features.
    # feature
    if (
        (feature is not None) and
        (len(str(feature)) > 0) and
        (str(feature) != "") and
        (str(feature).strip().lower() != "none")
    ):
        extraction = (
            putly.parse_extract_text_keys_values_semicolon_colon_comma(
                text=feature,
            )
        )["features_values"]
        pail["feature"] = dict()
        for key in extraction.keys():
            pail["feature"][key] = copy.deepcopy(extraction[key][0])
            pass
        pail["feature_product"] = list(pail["feature"].keys())[0]
        pail["feature_source"] = pail["feature"][pail["feature_product"]]
    else:
        pail["feature"] = None
        pail["feature_product"] = None
        pail["feature_source"] = None
        pass
    # features
    pail["features"] = putly.parse_text_list_values(
        text=features,
        delimiter=",",
    )
    # translation_features
    translation_features = str(translation_features).replace("#", " ")
    if (
        (translation_features is not None) and
        (len(str(translation_features)) > 0) and
        (str(translation_features) != "") and
        (str(translation_features).strip().lower() != "none")
    ):
        extraction = (
            putly.parse_extract_text_keys_values_semicolon_colon_comma(
                text=translation_features,
            )
        )["features_values"]
        pail["translation_features"] = dict()
        for key in extraction.keys():
            pail["translation_features"][key] = copy.deepcopy(
                extraction[key][0]
            )
            pass
    else:
        pail["translation_features"] = None
        pass
    # Prepare reference to sort rows by group in a subsequent table.
    if (pail["features"] is not None):
        pail["sequence_features"] = dict(zip(
            pail["features"],
            range(len(pail["features"]))
        ))
    else:
        pail["sequence_features"] = None
        pass

    ##########
    # Values and confidence intervals.
    pail["values_intervals_primary"] = parse_values_intervals(
        text=values_intervals_primary,
    )
    pail["values_intervals_secondary"] = parse_values_intervals(
        text=values_intervals_secondary,
    )
    pail["values_intervals_tertiary"] = parse_values_intervals(
        text=values_intervals_tertiary,
    )

    ##########
    # Numerical values.
    pail["minimum_abscissa"] = float(str(minimum_abscissa).strip())
    pail["maximum_abscissa"] = float(str(maximum_abscissa).strip())


    ##########
    # Titles.
    title_chart_parse = str(title_chart).strip().replace("#", " ")
    if (title_chart_parse != "none"):
        pail["title_chart"] = title_chart_parse
    else:
        pail["title_chart"] = ""
        pass
    title_abscissa_parse = str(title_abscissa).strip().replace("#", " ")
    if (title_chart_parse != "none"):
        pail["title_abscissa"] = title_abscissa_parse
    else:
        pail["title_abscissa"] = ""
        pass
    title_ordinate_parse = str(title_ordinate).strip().replace("#", " ")
    if (title_ordinate_parse != "none"):
        pail["title_ordinate"] = title_ordinate_parse
    else:
        pail["title_ordinate"] = ""
        pass

    ##########
    # Legends.
    # legend_series_primary
    pail["legend_series_primary"] = str(
        legend_series_primary
    ).replace("#", " ")
    # legend_series_secondary
    pail["legend_series_secondary"] = str(
        legend_series_secondary
    ).replace("#", " ")
    # legend_series_tertiary
    pail["legend_series_tertiary"] = str(
        legend_series_tertiary
    ).replace("#", " ")

    # report
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

    # Collect names of relevant columns in original source table.
    pail["columns_source"] = list()
    pail["columns_continuity_source"] = list()
    if (pail["feature_source"] is not None):
        pail["columns_source"].append(pail["feature_source"])
    if (pail["values_intervals_primary"] is not None):
        for key in pail["values_intervals_primary"].keys():
            pail["columns_continuity_source"].append(
                pail["values_intervals_primary"][key]
            )
    if (pail["values_intervals_secondary"] is not None):
        for key in pail["values_intervals_secondary"].keys():
            pail["columns_continuity_source"].append(
                pail["values_intervals_secondary"][key]
            )
    if (pail["values_intervals_tertiary"] is not None):
        for key in pail["values_intervals_tertiary"].keys():
            pail["columns_continuity_source"].append(
                pail["values_intervals_tertiary"][key]
            )
            pass
        pass
    pail["columns_source"].extend(copy.deepcopy(
        pail["columns_continuity_source"]
    ))

    # For convenience in translation of names of columns in table, invert keys
    # and values in dictionaries.
    pail["translation_columns"] = dict()
    extractions = [
        pail["feature"],
        pail["values_intervals_primary"],
        pail["values_intervals_secondary"],
        pail["values_intervals_tertiary"],
    ]
    for extraction in extractions:
        if (extraction is not None):
            for key_original in extraction.keys():
                key_novel = copy.deepcopy(extraction[key_original])
                pail["translation_columns"][key_novel] = key_original
                pass
            pass
        pass

    # Combine names of columns.
    pail["columns_product"] = list()
    pail["columns_continuity_product"] = list()
    pail["columns_product"].append(pail["feature_product"])
    extractions = [
        pail["values_intervals_primary"],
        pail["values_intervals_secondary"],
        pail["values_intervals_tertiary"],
    ]
    for extraction in extractions:
        if (extraction is not None):
            pail["columns_continuity_product"].extend(copy.deepcopy(list(
                extraction.keys()
            )))
            pass
        pass
    pail["columns_product"].extend(copy.deepcopy(
        pail["columns_continuity_product"]
    ))

    # Report.
    if pail["report"]:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_plot_dot_forest_from_table_data.py")
        print("function: parse_text_parameters()")
        putly.print_terminal_partition(level=5)
        print("parameters:")
        print(pail)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def define_column_types_table_data(
    feature=None,
    columns_continuity_source=None,
    columns_source=None,
):
    """
    Defines the types of variables for columns in table of parameters.

    Review: TCW; 5 May 2025

    arguments:
        feature (str): name of column in table for names or identifiers of
            features
        columns_continuity_source (list<str>): names of relevant columns in
            table
        columns_source (list<str>): names of relevant columns in table

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify types of variables in columns of table.
    types_columns = dict()
    types_columns[feature] = "string"
    for column in columns_continuity_source:
        types_columns[column] = "float32"
        pass
    # Return information.
    return types_columns


def read_organize_source_table_data(
    path_file_source_table=None,
    columns_source=None,
    columns_continuity_source=None,
    columns_product=None,
    columns_continuity_product=None,
    translation_columns=None,
    translation_features=None,
    feature=None,
    feature_source=None,
    feature_product=None,
    features=None,
    sequence_features=None,
    values_intervals_primary=None,
    values_intervals_secondary=None,
    values_intervals_tertiary=None,
    path_directory_dock=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 13 May 2025

    arguments:
        path_file_source_table (str): path to source file in text format as a
            table with tab delimiters between columns and newline delimiters
            between rows, with data for creation of plot chart
        columns_source (list<str>): names of relevant columns in table
        columns_continuity_source (list<str>): names of relevant columns in
            table
        columns_product (list<str>): names of relevant columns in table
        columns_continuity_product (list<str>): names of relevant columns in
            table
        translation_columns (dict<str>): translations of names of columns in
            table
        translation_features (dict<str>): translations of names or identifiers
            of features in rows of table
        feature (dict<str>): translation of column in table for names or
            identifiers of features
        feature_source (str): original, source name of column in table for
            names or identifiers of features
        feature_product (str): novel, product name of column in table for names
            or identifiers of features
        features (list<str>): names or identifiers of features to include for
            visual representation on a plot chart
        sequence_features (dict<int>): indices for sort that match categorical
            values in column of table
        values_intervals_primary (dict<str>): translations of columns in table
            for values and their intervals
        values_intervals_secondary (dict<str>): translations of columns in
            table for values and their intervals
        values_intervals_tertiary (dict<str>): translations of columns in
            table for values and their intervals
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Read information from file.
    types_columns = define_column_types_table_data(
        feature=feature_source,
        columns_continuity_source=columns_continuity_source,
        columns_source=columns_source,
    )
    table = pandas.read_csv(
        path_file_source_table,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Translate names of columns.
    table.rename(
        columns=translation_columns,
        inplace=True,
    )
    # Filter columns in table.
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=columns_product,
        report=report,
    )
    # Filter rows in table for observations.
    if (
        (features is not None) and
        (len(features) > 0)
    ):
        selection = dict()
        selection[feature_product] = copy.deepcopy(features)
        table = porg.filter_select_table_rows_by_columns_categories(
            table=table,
            columns_categories=selection,
            report=report,
        )
        pass
    # Sort rows in table by names or identifiers of features.
    if (sequence_features is not None):
        table = porg.sort_table_rows_by_single_column_reference(
            table=table,
            index_rows=feature_product,
            column_reference=feature_product,
            column_sort_temporary="sort_temporary_312",
            reference_sort=sequence_features,
        )
    # Translate names or identifiers of features.
    table["feature_translation"] = table[feature_product].copy(deep=True)
    table["feature_translation"] = table["feature_translation"].replace(
        translation_features
    )
    #table["feature_translation"] = table.apply(
    #    lambda row: (
    #        translation_features[row[feature_product]]
    #        ) if (
    #            row[feature_product] in translation_features.keys()
    #        ) else row[feature_product],
    #    axis="columns", # apply function to each row
    #)

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_plot_dot_forest_from_table_data.py")
        print("function: read_source_table_data()")
        putly.print_terminal_partition(level=5)
        print("table of data:")
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


##########
# Create plot chart and write to file.


def plot_dot_forest_category_ordinate_three_series(
    table=None,
    column_feature=None,
    column_feature_name=None,
    column_value_primary=None,
    column_interval_low_primary=None,
    column_interval_high_primary=None,
    column_value_secondary=None,
    column_interval_low_secondary=None,
    column_interval_high_secondary=None,
    column_value_tertiary=None,
    column_interval_low_tertiary=None,
    column_interval_high_tertiary=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    show_legend=None,
    legend_series_primary=None,
    legend_series_secondary=None,
    legend_series_tertiary=None,
    size_title_chart=None,
    size_title_abscissa=None,
    size_title_ordinate=None,
    size_label_abscissa=None,
    size_label_ordinate=None,
    size_label_legend=None,
    aspect=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    factor_space_series=None,
    space_between_series=None,
    position_line_origin=None,
    size_marker_primary=None,
    size_marker_secondary=None,
    size_marker_tertiary=None,
    size_line_origin=None,
    size_line_interval=None,
    color_marker_primary=None,
    color_marker_secondary=None,
    color_marker_tertiary=None,
    color_interval_primary=None,
    color_interval_secondary=None,
    color_interval_tertiary=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Dependency:
    This function is a dependency of the functions below.
    1.
    package: partner
    module or script: drive_plot_dot_forest_from_table_data.py
    function: create_write_plot_dot_forest()

    Create a plot chart of type scatter with discrete, categorical, text labels
    corresponding to features along the vertical ordinate axis and with
    points and error bars corresponding to quantitative, continuous values
    along the horizontal abscissa axis.

    A common name of this chart design is "Forest Plot".

    This function accommodates exactly two groups of values (series), along
    with their respective intervals.

    This plot's ordinate (vertical, y axis) fits best for a variable with
    discrete categorical values that serve as names.

    This plot's abscissa (horizontal, x axis) fits best for a variable with
    quantitative, continuous values on an interval or ratio scale of
    measurement that is centered on zero.

    MatPlotLib accepts intervals, not ranges for error bars. The function does
    the arithmetic to calculate ranges below and above the central value.

    Review: 13 May 2025

    arguments:
        table (object): Pandas data-frame table of features across rows with
            statistical measures such as correlation coefficients or regression
            coefficients and their confidence intervals across columns
        column_feature (str): name of column in table corresponding to unique
            names or identifiers of features in their proper sequence across
            rows
        column_feature_name (str): name of column in table corresponding to
            readable names of features for visual representation on labels
        column_value_primary (str): name of column in table corresponding to
            values
        column_interval_low_primary (str): name of column in table
            corresponding to intervals below values
        column_interval_high_primary (str): name of column in table
            corresponding to intervals above values
        column_value_secondary (str): name of column in table corresponding to
            values
        column_interval_low_secondary (str): name of column in table
            corresponding to intervals below values
        column_interval_high_secondary (str): name of column in table
            corresponding to intervals above values
        column_value_tertiary (str): name of column in table corresponding to
            values
        column_interval_low_tertiary (str): name of column in table
            corresponding to intervals below values
        column_interval_high_tertiary (str): name of column in table
            corresponding to intervals above values
        title_chart (str): title for plot chart as a whole
        title_abscissa (str): title for abscissa or horizontal axis
        title_ordinate (str): title for ordinate or vertical axis
        show_legend (bool): whether to create legend on plot chart for
            explanation of the series
        legend_series_primary (str): description in legend for primary series
        legend_series_secondary (str): description in legend for secondary
            series
        legend_series_tertiary (str): description in legend for tertiary
            series
        size_title_chart (str): font size for title of plot chart as a whole
        size_title_abscissa (str): font size for title of abscissa horizontal
            axis
        size_title_ordinate (str): font size for title of ordinate vertical
            axis
        size_label_abscissa (str): font size for labels of abscissa horizontal
            axis
        size_label_ordinate (str): font size for labels of ordinate vertical
            axis
        size_label_legend (str): font size for labels in legend
        aspect (str): aspect ratio for MatPlotLib chart figure
        minimum_abscissa (float): minimal value for range of abscissa axis
        maximum_abscissa (float): maximal value for range of abscissa axis
        factor_space_series (float): factor for automatic calculation of
            space between series on the basis of counts of features
        space_between_series (float): vertical space between markers for series
        position_line_origin (float): value of intercept on abscissa horizontal
            axis at which for a vertical line to cross, representing the origin
            or visual reference
        size_marker_primary (int): size of markers for primary series of values
        size_marker_secondary (int): size of markers for secondary series of
            values
        size_marker_tertiary (int): size of markers for tertiary series of
            values
        size_line_origin (int): size, width, or thickness of line for origin
        size_line_interval (int): size, width, or thickness of lines for
            representation of the interval of error or confidence
        color_marker_primary (str): color red, green, blue, alpha parameters
            for markers of primary series of values
        color_marker_secondary (str): color red, green, blue, alpha parameters
            for markers of secondary series of values
        color_marker_tertiary (str): color red, green, blue, alpha parameters
            for markers of tertiary series of values
        color_interval_primary (str): color red, green, blue, alpha parameters
            for markers of primary intervals
        color_interval_secondary (str): color red, green, blue, alpha
            parameters for markers of secondary intervals
        color_interval_tertiary (str): color red, green, blue, alpha
            parameters for markers of tertiary intervals
        fonts (dict<object>): definitions of font properties
        colors (dict<tuple>): definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object

    """

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.plot.py")
        print("function: plot_dot_forest_category_ordinate_two_series()")
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Organize information for chart.

    # Copy information.
    table = table.copy(deep=True)
    columns_available = table.columns.to_list()

    # Determine whether to represent values for one, two, or three series.
    # Assume that there is definition of parameters for the primary, secondary,
    # and tertiary series in sequential order.
    count_series = 0
    columns_series = [
        column_value_primary,
        column_value_secondary,
        column_value_tertiary,
    ]
    for column in columns_series:
        if (
            (column is not None) and
            (len(str(column)) > 0) and
            (str(column) != "") and
            (str(column).strip().lower() != "none") and
            (column in columns_available)
        ):
            count_series += 1
            pass
        pass

    # Organize information for categorical labels and positions of those labels
    # and their corresponding dot points along the ordinate vertical axis.
    # Assign positions for primary series to be above center point.
    # Assign positions for secondary series to be below center point.
    labels_ordinate = table[column_feature_name].to_list()
    count_ordinate = len(labels_ordinate)
    positions_ordinate_center = list(reversed(list(map(
        lambda position: (position + 1),
        range(count_ordinate)
    ))))
    if (
        (factor_space_series is not None) and
        (factor_space_series > 0)
    ):
        space_between_series = float(factor_space_series / count_ordinate)
        pass
    if (count_series == 1):
        positions_ordinate_primary = positions_ordinate_center
    elif (count_series == 2):
        positions_ordinate_primary = list(map(
            lambda position: (position + (2 * (space_between_series / 3))),
            positions_ordinate_center
        ))
        positions_ordinate_secondary = list(map(
            lambda position: (position - (2 * (space_between_series / 3))),
            positions_ordinate_center
        ))
    elif (count_series == 3):
        positions_ordinate_primary = list(map(
            lambda position: (position + (space_between_series / 2)),
            positions_ordinate_center
        ))
        positions_ordinate_secondary = positions_ordinate_center
        positions_ordinate_tertiary = list(map(
            lambda position: (position - (space_between_series / 2)),
            positions_ordinate_center
        ))
        pass

    # Organize information for labels, dot points, and error bars along the
    # abscissa horizontal axis.
    # Shape (n, 2)
    #errors_ordinate = numpy.array(list(zip(
    #    errors_ordinate_low, errors_ordinate_high
    #)))
    # Shape (2, n)
    if (
        (count_series == 1) or
        (count_series == 2) or
        (count_series == 3)
    ):
        positions_abscissa_primary = table[column_value_primary].to_numpy()
        intervals_abscissa_low_primary = (
            table[column_interval_low_primary].to_numpy()
        )
        intervals_abscissa_high_primary = (
            table[column_interval_high_primary].to_numpy()
        )
        intervals_abscissa_primary = numpy.array([
            intervals_abscissa_low_primary,
            intervals_abscissa_high_primary
        ])
    if (
        (count_series == 2) or
        (count_series == 3)
    ):
        positions_abscissa_secondary = table[column_value_secondary].to_numpy()
        intervals_abscissa_low_secondary = (
            table[column_interval_low_secondary].to_numpy()
        )
        intervals_abscissa_high_secondary = (
            table[column_interval_high_secondary].to_numpy()
        )
        intervals_abscissa_secondary = numpy.array([
            intervals_abscissa_low_secondary,
            intervals_abscissa_high_secondary
        ])
    if (count_series == 3):
        positions_abscissa_tertiary = table[column_value_tertiary].to_numpy()
        intervals_abscissa_low_tertiary = (
            table[column_interval_low_tertiary].to_numpy()
        )
        intervals_abscissa_high_tertiary = (
            table[column_interval_high_tertiary].to_numpy()
        )
        intervals_abscissa_tertiary = numpy.array([
            intervals_abscissa_low_tertiary,
            intervals_abscissa_high_tertiary
        ])
        pass

    ##########
    # Create and initialize figure chart object.

    # Create figure.
    figure = pplot.initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Create axes.
    #axes = matplotlib.pyplot.axes()
    axes = figure.add_subplot(111)
    # Define limits for axes.
    if (minimum_abscissa is not None):
        axes.set_xlim(xmin=minimum_abscissa)
    if (maximum_abscissa is not None):
        axes.set_xlim(xmax=maximum_abscissa)
    axes.set_ylim(ymin=float(min(positions_ordinate_center) - 1))
    axes.set_ylim(ymax=float(max(positions_ordinate_center) + 1))

    # Include title label on chart.
    if len(title_chart) > 0:
        axes.set_title(
            title_chart,
            fontproperties=fonts["properties"][size_title_chart],
            loc="right",
            horizontalalignment="right",
            verticalalignment="top",
            pad=5,
        )
        pass
    # Set titles for axes.
    if (len(title_abscissa) > 0):
        axes.set_xlabel(
            xlabel=title_abscissa,
            labelpad=15,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
    if (len(title_ordinate) > 0):
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=15,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_ordinate]
        )
    # Set parameters for tick labels on axes.
    axes.tick_params(
        axis="both", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        direction="out",
        top=False,
        labeltop=False,
        bottom=True,
        labelbottom=True,
        left=True,
        labelleft=True,
        right=False,
        labelright=False,
        length=7.5, # 5.0
        width=5.0, # 3.0, 5.0
        pad=10.0, # 5.0, 7.5
        color=colors["black"],
        labelcolor=colors["black"],
    )
    axes.tick_params(
        axis="x",
        which="both",
        labelsize=fonts["values"][size_label_abscissa]["size"],
    )
    axes.tick_params(
        axis="y",
        which="both",
        labelsize=fonts["values"][size_label_ordinate]["size"],
    )
    # Keep axes, ticks, and labels, but remove border.
    # ["left", "top", "right", "bottom",]
    for position in ["top", "right",]:
        matplotlib.pyplot.gca().spines[position].set_visible(False)

    # Set explicit tick positions and labels on vertical ordinate axis.
    # (https://matplotlib.org/3.5.1/api/_as_gen/
    # matplotlib.axes.Axes.set_yticks.html)
    axes.set_xticks(
        [round(minimum_abscissa, 1), 0.0, round(maximum_abscissa, 1)],
        labels=None,
        minor=False,
    )
    axes.set_yticks(
        positions_ordinate_center, # center positions with even spacing
        labels=labels_ordinate, # place labels at center positions
        minor=False,
    )
    # Plot dashed line at origin.
    # Consider the "dashes" argument for fine control of line dash pattern.
    axes.axvline(
        x=position_line_origin,
        ymin=0,
        ymax=1,
        alpha=1.0,
        color=colors["gray"],
        linestyle="--",
        linewidth=size_line_origin,
    )

    ##########
    # Represent main information on the chart figure object.

    # Plot points and error bars for values and intervals from each series.
    # First plot markers for group two so that these are below.
    # Second plot markers for group one so that these are above.
    # (https://matplotlib.org/3.5.1/api/_as_gen/
    # matplotlib.axes.Axes.errorbar.html)
    if (
        (count_series == 1) or
        (count_series == 2) or
        (count_series == 3)
    ):
        handle_primary = axes.errorbar(
            positions_abscissa_primary,
            positions_ordinate_primary,
            yerr=None,
            xerr=intervals_abscissa_primary,
            ecolor=color_interval_primary,
            elinewidth=size_line_interval, # 7.5
            barsabove=False, # whether to print error bars in layer above points
            linestyle="",
            marker="o", # marker shape: circle
            markersize=size_marker_primary, # 5, 15, 50, 70
            markeredgecolor=color_marker_primary, # colors["purple"],
            markerfacecolor=color_marker_primary, # colors["purple"],
        )
        pass
    if (
        (count_series == 2) or
        (count_series == 3)
    ):
        handle_secondary = axes.errorbar(
            positions_abscissa_secondary,
            positions_ordinate_secondary,
            yerr=None,
            xerr=intervals_abscissa_secondary,
            ecolor=color_interval_secondary,
            elinewidth=size_line_interval, # 7.5
            barsabove=False, # whether to print error bars in layer above points,
            linestyle="",
            marker="D", # "D" marker shape: diamond
            markersize=size_marker_secondary, # 5, 15
            markeredgecolor=color_marker_secondary, # colors["green"],
            markerfacecolor=color_marker_secondary, # colors["green"],
        )
        pass
    if (count_series == 3):
        handle_tertiary = axes.errorbar(
            positions_abscissa_tertiary,
            positions_ordinate_tertiary,
            yerr=None,
            xerr=intervals_abscissa_tertiary,
            ecolor=color_interval_tertiary,
            elinewidth=size_line_interval, # 7.5
            barsabove=False, # whether to print error bars in layer above points,
            linestyle="",
            marker="s", # "s" marker shape: square
            markersize=size_marker_tertiary, # 5, 15
            markeredgecolor=color_marker_tertiary, # colors["green"],
            markerfacecolor=color_marker_tertiary, # colors["green"],
        )
        pass

    # Create legend.
    # Create custom elements for the legend.
    if (show_legend):
        if (
            (count_series == 1) or
            (count_series == 2) or
            (count_series == 3)
        ):
            handle_primary = axes.errorbar(
                [0],
                [-1], # create outside of visible portion of axes
                yerr=None,
                xerr=5.0,
                label=legend_series_primary,
                ecolor=color_interval_primary,
                elinewidth=(size_line_interval/1.5),
                barsabove=False, # whether to print error bars in layer above points
                linestyle="",
                marker="o", # "o" marker shape: circle
                markersize=(size_marker_primary/1.5), # 5, 15, 50, 70
                markeredgecolor=color_marker_primary, # colors["purple"],
                markerfacecolor=color_marker_primary, # colors["purple"],
            )
            pass
        if (
            (count_series == 2) or
            (count_series == 3)
        ):
            handle_secondary = axes.errorbar(
                [0],
                [-1], # create outside of visible portion of axes
                yerr=None,
                xerr=5.0,
                label=legend_series_secondary,
                ecolor=color_interval_secondary,
                elinewidth=(size_line_interval/1.5),
                barsabove=False, # whether to print error bars in layer above points
                linestyle="",
                marker="D", # "D" marker shape: diamond
                markersize=(size_marker_secondary/1.5), # 5, 15
                markeredgecolor=color_marker_secondary, # colors["green"],
                markerfacecolor=color_marker_secondary, # colors["green"],
            )
            pass
        if (count_series == 3):
            handle_tertiary = axes.errorbar(
                [0],
                [-1], # create outside of visible portion of axes
                yerr=None,
                xerr=5.0,
                label=legend_series_tertiary,
                ecolor=color_interval_tertiary,
                elinewidth=(size_line_interval/1.5),
                barsabove=False, # whether to print error bars in layer above points
                linestyle="",
                marker="s", # "s" marker shape: square
                markersize=(size_marker_tertiary/1.5), # 5, 15
                markeredgecolor=color_marker_tertiary, # colors["green"],
                markerfacecolor=color_marker_tertiary, # colors["green"],
            )
            pass
        if (count_series == 1):
            axes.legend(
                handles=[
                    handle_primary,
                ],
                loc="upper right",
                prop=fonts["properties"][size_label_legend],
                title="",
                title_fontsize=fonts["values"][size_label_legend]["size"]
            )
            pass
        elif (count_series == 2):
            axes.legend(
                handles=[
                    handle_primary, handle_secondary,
                ],
                loc="upper right",
                prop=fonts["properties"][size_label_legend],
                title="",
                title_fontsize=fonts["values"][size_label_legend]["size"]
            )
            pass
        elif (count_series == 3):
            axes.legend(
                handles=[
                    handle_primary, handle_secondary, handle_tertiary,
                ],
                loc="upper right",
                prop=fonts["properties"][size_label_legend],
                title="",
                title_fontsize=fonts["values"][size_label_legend]["size"]
            )
            pass
        pass

    # Return figure.
    return figure


def create_write_plot_chart_dot_forest(
    table=None,
    column_feature=None,
    column_feature_name=None,
    column_value_primary=None,
    column_interval_low_primary=None,
    column_interval_high_primary=None,
    column_value_secondary=None,
    column_interval_low_secondary=None,
    column_interval_high_secondary=None,
    column_value_tertiary=None,
    column_interval_low_tertiary=None,
    column_interval_high_tertiary=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    legend_series_primary=None,
    legend_series_secondary=None,
    legend_series_tertiary=None,
    name_file_product_chart=None,
    path_directory_product=None,
    report=None,
):
    """
    Create and write to file a plot chart of type dot or forest.

    The sequence of features along the vertical, ordinate axis will correspond
    to their sequence in the table.
    The values designated as primary will appear on top with more prominent
    markings, whereas the values designated as secondary will appear on bottom
    with less prominent markings.

    Review: TCW; 24 July 2025

    arguments:
        table (object): Pandas data-frame table of features across rows with
            statistical measures such as correlation coefficients or regression
            coefficients and their confidence intervals across columns
        column_feature (str): name of column in table corresponding to unique
            names or identifiers of features in their proper sequence across
            rows
        column_feature_name (str): name of column in table corresponding to
            readable names of features for visual representation on labels
        column_value_primary (str): name of column in table corresponding to
            values
        column_interval_low_primary (str): name of column in table
            corresponding to intervals below values
        column_interval_high_primary (str): name of column in table
            corresponding to intervals above values
        column_value_secondary (str): name of column in table corresponding to
            values
        column_interval_low_secondary (str): name of column in table
            corresponding to intervals below values
        column_interval_high_secondary (str): name of column in table
            corresponding to intervals above values
        column_value_tertiary (str): name of column in table corresponding to
            values
        column_interval_low_tertiary (str): name of column in table
            corresponding to intervals below values
        column_interval_high_tertiary (str): name of column in table
            corresponding to intervals above values
        minimum_abscissa (float): value for minimal limit to represent on the
            abscissa horizontal axis
        maximum_abscissa (float): value for maximal limit to represent on the
            abscissa horizontal axis
        title_chart (str): title for chart figure
        title_abscissa (str): title for abscissa horizontal axis
        title_ordinate (str): title for ordinate vertical axis
        legend_series_primary (str): description in legend for primary series
        legend_series_secondary (str): description in legend for secondary
            series
        legend_series_tertiary (str): description in legend for tertiary
            series
        name_file_product_chart (str): name of product file
        path_directory_product (str): path to directory for procedure's product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()

    # Create figure.
    figure = pplot.plot_dot_forest_category_ordinate_three_series(
        table=table,
        column_feature=column_feature,
        column_feature_name=column_feature_name,
        column_value_primary=column_value_primary,
        column_interval_low_primary=column_interval_low_primary,
        column_interval_high_primary=column_interval_high_primary,
        column_value_secondary=column_value_secondary,
        column_interval_low_secondary=column_interval_low_secondary,
        column_interval_high_secondary=column_interval_high_secondary,
        column_value_tertiary=column_value_tertiary,
        column_interval_low_tertiary=column_interval_low_tertiary,
        column_interval_high_tertiary=column_interval_high_tertiary,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        show_legend=True,
        legend_series_primary=legend_series_primary,
        legend_series_secondary=legend_series_secondary,
        legend_series_tertiary=legend_series_tertiary,
        size_title_chart="ten",
        size_title_abscissa="ten",
        size_title_ordinate="ten",
        size_label_abscissa="nine",
        size_label_ordinate="nine",
        size_label_legend="thirteen",
        aspect="portrait",
        minimum_abscissa=minimum_abscissa, # -2.5
        maximum_abscissa=maximum_abscissa, # 2.5
        factor_space_series=None, # if not zero or None, overrides space
        space_between_series=0.33,
        position_line_origin=0.0,
        size_marker_primary=15, # 25
        size_marker_secondary=13, # 20
        size_marker_tertiary=9,
        size_line_origin=3.5,
        size_line_interval=2, # 3.5
        color_marker_primary=colors["blue_navy"],
        color_marker_secondary=colors["orange_burnt"],
        color_marker_tertiary=colors["red_burgundy"],
        color_interval_primary=colors["black"],
        color_interval_secondary=colors["gray_dark"],
        color_interval_tertiary=colors["gray_dark"],
        fonts=fonts,
        colors=colors,
        report=report,
    )

    ##########
    # Bundle information.
    pail_write_plot = dict()
    pail_write_plot["forest_plot"] = figure

    ##########
    # Write product information to file.

    # Define paths to directories.
    path_directory = os.path.join(
        path_directory_product, "forest_plot",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory,
    )
    # Write figures to file.
    #pplot.write_product_plots_parent_directory(
    #    pail_write=pail_write_plot_box,
    #    format="svg", # jpg, png, svg
    #    resolution=300,
    #    path_directory=path_directory_box,
    #)
    pplot.write_product_plots_parent_directory(
        pail_write=pail_write_plot,
        format="jpg", # jpg, png, svg
        resolution=150,
        path_directory=path_directory,
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_plot_dot_forest_from_table_data.py")
        print("function: create_write_plot_dot_forest()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    pass


################################################################################
# Procedure


##########
# Call main procedure.


def execute_procedure(
    path_file_source_table=None,
    name_file_product_chart=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    feature=None,
    features=None,
    translation_features=None,
    values_intervals_primary=None,
    values_intervals_secondary=None,
    values_intervals_tertiary=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    legend_series_primary=None,
    legend_series_secondary=None,
    legend_series_tertiary=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_file_source_table (str): path to source file in text format as a
            table with tab delimiters between columns and newline delimiters
            between rows, with data for creation of plot chart
        name_file_product_chart (str): name of product file
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        feature (str): parameters for extraction of features
        features (str): parameters for extraction of features
        translation_features (str): parameters for extraction of features
        values_intervals_primary (str): parameters for extraction of values
            and intervals
        values_intervals_secondary (str): parameters for extraction of values
            and intervals
        values_intervals_tertiary (str): parameters for extraction of values
            and intervals
        minimum_abscissa (float): value for minimal limit to represent on the
            abscissa horizontal axis
        maximum_abscissa (float): value for maximal limit to represent on the
            abscissa horizontal axis
        title_chart (str): title for chart figure
        title_abscissa (str): title for abscissa horizontal axis
        title_ordinate (str): title for ordinate vertical axis
        legend_series_primary (str): description in legend for primary series
        legend_series_secondary (str): description in legend for secondary
            series
        legend_series_tertiary (str): description in legend for tertiary series
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parameters.
    pail_parameters = parse_text_parameters(
        path_file_source_table=path_file_source_table,
        name_file_product_chart=name_file_product_chart,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        feature=feature,
        features=features,
        translation_features=translation_features,
        values_intervals_primary=values_intervals_primary,
        values_intervals_secondary=values_intervals_secondary,
        values_intervals_tertiary=values_intervals_tertiary,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        legend_series_primary=legend_series_primary,
        legend_series_secondary=legend_series_secondary,
        legend_series_tertiary=legend_series_tertiary,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_plot_dot_forest_from_table_data.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print(
            "path_file_source_table: " +
            str(pail_parameters["path_file_source_table"])
        )
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    table = read_organize_source_table_data(
        path_file_source_table=pail_parameters["path_file_source_table"],
        columns_source=pail_parameters["columns_source"],
        columns_continuity_source=pail_parameters["columns_continuity_source"],
        columns_product=pail_parameters["columns_product"],
        columns_continuity_product=(
            pail_parameters["columns_continuity_product"]
        ),
        translation_columns=pail_parameters["translation_columns"],
        translation_features=pail_parameters["translation_features"],
        feature=pail_parameters["feature"],
        feature_source=pail_parameters["feature_source"],
        feature_product=pail_parameters["feature_product"],
        features=pail_parameters["features"],
        sequence_features=pail_parameters["sequence_features"],
        values_intervals_primary=pail_parameters["values_intervals_primary"],
        values_intervals_secondary=(
            pail_parameters["values_intervals_secondary"]
        ),
        values_intervals_tertiary=(
            pail_parameters["values_intervals_tertiary"]
        ),
        path_directory_dock=pail_parameters["path_directory_dock"],
        report=pail_parameters["report"],
    )

    ##########
    # Create visual representation in plot chart and write to file.
    create_write_plot_chart_dot_forest(
        table=table,
        column_feature=pail_parameters["feature_product"],
        column_feature_name="feature_translation",
        column_value_primary="value_primary",
        column_interval_low_primary="interval_low_primary",
        column_interval_high_primary="interval_high_primary",
        column_value_secondary="value_secondary",
        column_interval_low_secondary="interval_low_secondary",
        column_interval_high_secondary="interval_high_secondary",
        column_value_tertiary="value_tertiary",
        column_interval_low_tertiary="interval_low_tertiary",
        column_interval_high_tertiary="interval_high_tertiary",
        minimum_abscissa=pail_parameters["minimum_abscissa"],
        maximum_abscissa=pail_parameters["maximum_abscissa"],
        title_chart=pail_parameters["title_chart"],
        title_abscissa=pail_parameters["title_abscissa"],
        title_ordinate=pail_parameters["title_ordinate"],
        legend_series_primary=pail_parameters["legend_series_primary"],
        legend_series_secondary=pail_parameters["legend_series_secondary"],
        legend_series_tertiary=pail_parameters["legend_series_tertiary"],
        name_file_product_chart=pail_parameters["name_file_product_chart"],
        path_directory_product=pail_parameters["path_directory_product"],
        report=pail_parameters["report"],
    )

    pass


# Execute program process in Python.


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_source_table = sys.argv[1]
    name_file_product_chart = sys.argv[2]
    path_directory_source = sys.argv[3]
    path_directory_product = sys.argv[4]
    path_directory_dock = sys.argv[5]
    feature = sys.argv[6]
    features = sys.argv[7]
    translation_features = sys.argv[8]
    values_intervals_primary = sys.argv[9]
    values_intervals_secondary = sys.argv[10]
    values_intervals_tertiary = sys.argv[11]
    minimum_abscissa = sys.argv[12]
    maximum_abscissa = sys.argv[13]
    title_chart = sys.argv[14]
    title_abscissa = sys.argv[15]
    title_ordinate = sys.argv[16]
    legend_series_primary = sys.argv[17]
    legend_series_secondary = sys.argv[18]
    legend_series_tertiary = sys.argv[19]
    report = sys.argv[20]

    # Call function for procedure.
    execute_procedure(
        path_file_source_table=path_file_source_table,
        name_file_product_chart=name_file_product_chart,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        feature=feature,
        features=features,
        translation_features=translation_features,
        values_intervals_primary=values_intervals_primary,
        values_intervals_secondary=values_intervals_secondary,
        values_intervals_tertiary=values_intervals_tertiary,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        legend_series_primary=legend_series_primary,
        legend_series_secondary=legend_series_secondary,
        legend_series_tertiary=legend_series_tertiary,
        report=report,
    )

    pass



#
