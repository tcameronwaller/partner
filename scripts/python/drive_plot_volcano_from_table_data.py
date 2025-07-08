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
# Date, first execution: 8 July 2025
# Date, last execution or modification: 8 July 2025
# Review: TCW; 8 July 2025
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

# TODO: TCW; 12 June 2025
# pass a new parameter for the limits of the horizontal axis
# maybe?


# TODO: TCW; 2 July 2025
# In text titles, replace "#" with " " as I did in the program for forest plots.


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
    legend_series_tertiary=None,
    values_intervals_primary=None,
    values_intervals_secondary=None,
    values_intervals_tertiary=None,
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
        legend_series_tertiary (str): description in legend for tertiary series
        values_intervals_primary (str): parameters for extraction of values
            and intervals
        values_intervals_secondary (str): parameters for extraction of values
            and intervals
        values_intervals_tertiary (str): parameters for extraction of values
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
        legend_series_tertiary=legend_series_tertiary,
        values_intervals_primary=values_intervals_primary,
        values_intervals_secondary=values_intervals_secondary,
        values_intervals_tertiary=values_intervals_tertiary,
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
        values_intervals_tertiary=(
            pail_parameters["values_intervals_tertiary"]
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
        legend_series_tertiary=pail_parameters["legend_series_tertiary"],
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
        path_directory_product=path_directory_product,
        report=report,
    )

    pass


# Execute program process in Python.


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_table_source = sys.argv[1]
    title_chart=sys.argv[2]
    title_axis_abscissa=sys.argv[2]
    title_axis_ordinate=sys.argv[2]
    feature = sys.argv[3]
    features = sys.argv[4]
    translation_features = sys.argv[5]
    legend_series_primary = sys.argv[6]
    legend_series_secondary = sys.argv[7]
    legend_series_tertiary = sys.argv[8]
    values_intervals_primary = sys.argv[9]
    values_intervals_secondary = sys.argv[10]
    values_intervals_tertiary = sys.argv[11]
    path_file_figure_product = sys.argv[12]
    path_directory_product = sys.argv[13]
    path_directory_dock = sys.argv[14]
    report = sys.argv[15]

    # Call function for procedure.
    execute_procedure(
        path_file_table_data=path_file_table_data,
        title=title,
        feature=feature,
        features=features,
        translation_features=translation_features,
        legend_series_primary=legend_series_primary,
        legend_series_secondary=legend_series_secondary,
        legend_series_tertiary=legend_series_tertiary,
        values_intervals_primary=values_intervals_primary,
        values_intervals_secondary=values_intervals_secondary,
        values_intervals_tertiary=values_intervals_tertiary,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass



#
