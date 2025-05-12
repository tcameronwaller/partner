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
# Date, last execution or modification: 2 May 2025
# Review: TCW; 2 May 2025
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


# TODO: TCW; 9 May 2025
# pass a new parameter for the name of the file for the plot chart


##########
# Read source information from file.


def parse_text_parameters(
    title=None,
    feature=None,
    features=None,
    translation_features=None,
    legend_series_primary=None,
    legend_series_secondary=None,
    values_intervals_primary=None,
    values_intervals_secondary=None,
    report=None,
):
    """
    Parse parameters from text.

    arguments:
        title (str): title for top center of plot chart
        feature (str): parameters for extraction of features
        features (str): parameters for extraction of features
        translation_features (str): parameters for extraction of features
        legend_series_primary (str): description in legend for primary series
        legend_series_secondary (str): description in legend for secondary
            series
        values_intervals_primary (str): parameters for extraction of values
            and intervals
        values_intervals_secondary (str): parameters for extraction of values
            and intervals
        report (str): whether to print reports

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
            report (bool): whether to print reports
    """

    # Parse and extract information from text.
    # title.
    title_parse = str(title).replace("#", " ")
    # feature.
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
        feature_parse = dict()
        for key in extraction.keys():
            feature_parse[key] = copy.deepcopy(extraction[key][0])
            pass
        feature_product = list(feature_parse.keys())[0]
        feature_source = feature_parse[feature_product]
    else:
        feature_parse = None
        feature_product = None
        feature_source = None
        pass
    # features.
    features_parse = putly.parse_text_list_values(
        text=features,
        delimiter=",",
    )

    # translation_features.
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
        translation_features_parse = dict()
        for key in extraction.keys():
            translation_features_parse[key] = copy.deepcopy(extraction[key][0])
            pass
    else:
        translation_features_parse = None
        pass

    # legend_series_primary.
    legend_series_primary_parse = str(
        legend_series_primary
    ).replace("#", " ")
    # legend_series_secondary.
    legend_series_secondary_parse = str(
        legend_series_secondary
    ).replace("#", " ")
    # values_intervals_primary.
    if (
        (values_intervals_primary is not None) and
        (len(str(values_intervals_primary)) > 0) and
        (str(values_intervals_primary) != "") and
        (str(values_intervals_primary).strip().lower() != "none")
    ):
        extraction = (
            putly.parse_extract_text_keys_values_semicolon_colon_comma(
                text=values_intervals_primary,
            )
        )["features_values"]
        values_intervals_primary_parse = dict()
        for key in extraction.keys():
            values_intervals_primary_parse[key] = copy.deepcopy(
                extraction[key][0]
            )
            pass
    else:
        values_intervals_primary_parse = None
        pass
    # values_intervals_secondary.
    if (
        (values_intervals_secondary is not None) and
        (len(str(values_intervals_secondary)) > 0) and
        (str(values_intervals_secondary) != "") and
        (str(values_intervals_secondary).strip().lower() != "none")
    ):
        extraction = (
            putly.parse_extract_text_keys_values_semicolon_colon_comma(
                text=values_intervals_secondary,
            )
        )["features_values"]
        values_intervals_secondary_parse = dict()
        for key in extraction.keys():
            values_intervals_secondary_parse[key] = copy.deepcopy(
                extraction[key][0]
            )
            pass
    else:
        values_intervals_secondary_parse = None
        pass
    # report.
    if (
        (report is not None) and
        (str(report) != "") and
        (str(report) != "none") and
        (str(report) == "true")
    ):
        report_parse = True
    else:
        report_parse = False
        pass

    # Collect names of relevant columns in original source table.
    columns_source = list()
    columns_continuity_source = list()
    if (feature_source is not None):
        columns_source.append(feature_source)
    if (values_intervals_primary_parse is not None):
        for key in values_intervals_primary_parse.keys():
            columns_continuity_source.append(
                values_intervals_primary_parse[key]
            )
    if (values_intervals_secondary_parse is not None):
        for key in values_intervals_secondary_parse.keys():
            columns_continuity_source.append(
                values_intervals_secondary_parse[key]
            )
            pass
        pass
    columns_source.extend(copy.deepcopy(columns_continuity_source))

    # For convenience in translation of names of columns in table, invert keys
    # and values in dictionaries.
    translation_columns = dict()
    extractions = [
        feature_parse,
        values_intervals_primary_parse,
        values_intervals_secondary_parse,
    ]
    for extraction in extractions:
        if (extraction is not None):
            for key_original in extraction.keys():
                key_novel = copy.deepcopy(extraction[key_original])
                translation_columns[key_novel] = key_original
                pass
            pass
        pass

    # Combine names of columns.
    columns_product = list()
    columns_continuity_product = list()
    columns_product.append(feature_product)
    columns_continuity_product.extend(copy.deepcopy(list(
        values_intervals_primary_parse.keys()
    )))
    columns_continuity_product.extend(copy.deepcopy(list(
        values_intervals_secondary_parse.keys()
    )))
    columns_product.extend(copy.deepcopy(columns_continuity_product))

    # Prepare reference to sort rows by group in a subsequent table.
    sequence_features = dict(zip(
        features_parse,
        range(len(features_parse))
    ))

    # Collect information.
    pail = dict()
    pail["columns_source"] = columns_source
    pail["columns_continuity_source"] = columns_continuity_source
    pail["columns_product"] = columns_product
    pail["columns_continuity_product"] = columns_continuity_product
    pail["translation_columns"] = translation_columns
    pail["translation_features"] = translation_features_parse
    pail["title"] = title_parse
    pail["feature"] = feature_parse
    pail["feature_source"] = feature_source
    pail["feature_product"] = feature_product
    pail["features"] = features_parse
    pail["sequence_features"] = sequence_features
    pail["legend_series_primary"] = legend_series_primary_parse
    pail["legend_series_secondary"] = legend_series_secondary_parse
    pail["values_intervals_primary"] = values_intervals_primary_parse
    pail["values_intervals_secondary"] = values_intervals_secondary_parse
    pail["report"] = report_parse

    # Report.
    if report_parse:
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
    path_file_table_data=None,
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
    path_directory_dock=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 5 May 2025

    arguments:
        path_file_table_data (str): path to source file in text format as a
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
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    if False:
        # Fill missing values.
        table["value_primary"].fillna(
            value=-3,
            inplace=True,
        )
        table["interval_low_primary"].fillna(
            value=0,
            inplace=True,
        )
        pass

    # Read information from file.
    types_columns = define_column_types_table_data(
        feature=feature_source,
        columns_continuity_source=columns_continuity_source,
        columns_source=columns_source,
    )
    table = pandas.read_csv(
        path_file_table_data,
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
    selection = dict()
    selection[feature_product] = copy.deepcopy(features)
    table = porg.filter_select_table_rows_by_columns_categories(
        table=table,
        columns_categories=selection,
        report=report,
    )
    # Sort rows in table by names or identifiers of features.
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


def create_write_plot_dot_forest(
    title=None,
    legend_series_primary=None,
    legend_series_secondary=None,
    table=None,
    column_feature=None,
    column_feature_name=None,
    column_value_primary=None,
    column_interval_low_primary=None,
    column_interval_high_primary=None,
    column_value_secondary=None,
    column_interval_low_secondary=None,
    column_interval_high_secondary=None,
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

    Review: TCW; 12 May 2025

    arguments:
        title (str): title for top center of plot chart
        legend_series_primary (str): description in legend for primary series
        legend_series_secondary (str): description in legend for secondary
            series
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
    figure = pplot.plot_dot_forest_category_ordinate_two_series(
        table=table,
        column_feature=column_feature,
        column_feature_name=column_feature_name,
        column_value_primary=column_value_primary,
        column_interval_low_primary=column_interval_low_primary,
        column_interval_high_primary=column_interval_high_primary,
        column_value_secondary=column_value_secondary,
        column_interval_low_secondary=column_interval_low_secondary,
        column_interval_high_secondary=column_interval_high_secondary,
        title_chart=title,
        title_abscissa="Regression Coefficient (95% CI)",
        title_ordinate="",
        show_legend=True,
        legend_series_primary=legend_series_primary,
        legend_series_secondary=legend_series_secondary,
        size_title_chart="ten",
        size_title_abscissa="ten",
        size_title_ordinate="ten",
        size_label_abscissa="ten",
        size_label_ordinate="ten",
        size_label_legend="fifteen",
        aspect="portrait",
        minimum_abscissa=-5.0,
        maximum_abscissa=5.0,
        factor_space_series=100.0, # if not zero or None, overrides space
        space_between_series=0.15,
        position_line_origin=0.0,
        size_marker_primary=15,
        size_marker_secondary=10,
        size_line_origin=3,
        size_line_interval=3,
        color_marker_primary=colors["blue_navy"],
        color_marker_secondary=colors["orange_burnt"],
        color_interval_primary=colors["black"],
        color_interval_secondary=colors["gray_dark"],
        fonts=fonts,
        colors=colors,
        report=report,
    )

    ##########
    # Collect information.
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
        resolution=300,
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
    path_file_table_data=None,
    title=None,
    feature=None,
    features=None,
    translation_features=None,
    legend_series_primary=None,
    legend_series_secondary=None,
    values_intervals_primary=None,
    values_intervals_secondary=None,
    path_directory_product=None,
    path_directory_dock=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_file_table_data (str): path to source file in text format as a
            table with tab delimiters between columns and newline delimiters
            between rows, with data for creation of plot chart
        title (str): title for top center of plot chart
        feature (str): parameters for extraction of features
        features (str): parameters for extraction of features
        translation_features (str): parameters for extraction of features
        legend_series_primary (str): description in legend for primary series
        legend_series_secondary (str): description in legend for secondary
            series
        values_intervals_primary (str): parameters for extraction of values
            and intervals
        values_intervals_secondary (str): parameters for extraction of values
            and intervals
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
    pail_parameters = parse_text_parameters(
        title=title,
        feature=feature,
        features=features,
        translation_features=translation_features,
        legend_series_primary=legend_series_primary,
        legend_series_secondary=legend_series_secondary,
        values_intervals_primary=values_intervals_primary,
        values_intervals_secondary=values_intervals_secondary,
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
        print("path_file_table_data: " + str(path_file_table_data))
        print("path_directory_product: " + str(path_directory_product))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    table = read_organize_source_table_data(
        path_file_table_data=path_file_table_data,
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
        path_directory_dock=path_directory_dock,
        report=pail_parameters["report"],
    )

    ##########
    # Create visual representation in plot chart and write to file.
    create_write_plot_dot_forest(
        title=pail_parameters["title"],
        legend_series_primary=pail_parameters["legend_series_primary"],
        legend_series_secondary=pail_parameters["legend_series_secondary"],
        table=table,
        column_feature=pail_parameters["feature_product"],
        column_feature_name="feature_translation",
        column_value_primary="value_primary",
        column_interval_low_primary="interval_low_primary",
        column_interval_high_primary="interval_high_primary",
        column_value_secondary="value_secondary",
        column_interval_low_secondary="interval_low_secondary",
        column_interval_high_secondary="interval_high_secondary",
        path_directory_product=path_directory_product,
        report=report,
    )

    pass


# Execute program process in Python.


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_table_data = sys.argv[1]
    title=sys.argv[2]
    feature = sys.argv[3]
    features = sys.argv[4]
    translation_features = sys.argv[5]
    legend_series_primary = sys.argv[6]
    legend_series_secondary = sys.argv[7]
    values_intervals_primary = sys.argv[8]
    values_intervals_secondary = sys.argv[9]
    path_directory_product = sys.argv[10]
    path_directory_dock = sys.argv[11]
    report = sys.argv[12]

    # Call function for procedure.
    execute_procedure(
        path_file_table_data=path_file_table_data,
        title=title,
        feature=feature,
        features=features,
        translation_features=translation_features,
        legend_series_primary=legend_series_primary,
        legend_series_secondary=legend_series_secondary,
        values_intervals_primary=values_intervals_primary,
        values_intervals_secondary=values_intervals_secondary,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass



#
