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
# Date, first execution: 2 July 2025
# Date, last execution or modification: 16 July 2025
# Review: TCW; 16 July 2025
################################################################################
# Note


##########
# Review: TCW;

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
import partner.plot as pplot


#dir()
#importlib.reload()

###############################################################################
# Functionality


def parse_text_parameters(
    path_file_source_table_effects=None,
    path_file_source_list_features_emphasis=None,
    name_file_product_chart=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    column_identifier=None,
    column_name=None,
    column_effect=None,
    column_p_value=None,
    column_significance=None,
    threshold_effect=None,
    threshold_significance=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    lines_thresholds=None,
    label_emphasis=None,
    label_count=None,
    report=None,
):
    """
    Parse parameters from text.

    arguments:
        path_file_source_table_effects (str): path to source file
        path_file_source_list_features_emphasis (str): path to source file
        name_file_product_chart (str): name of product file
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        column_identifier (str): name of column in source table
        column_name (str): name of column in source table
        column_effect (str): name of column in source table
        column_p_value (str): name of column in source table
        column_significance (str): name of column in source table
        threshold_effect (str): value for threshold on magnitude of effect
            (|effect| > threshold) that is on the same scale, such as
            base-two logarithm, as the actual values themselves
        threshold_significance (str): value for threshold on p-values or
            q-values (p-value or q-value < threshold) that is not on a scale of
            the negative logarithm
        minimum_abscissa (str): value for minimal limit to represent on the
            abscissa horizontal axis
        maximum_abscissa (str): value for maximal limit to represent on the
            abscissa horizontal axis
        minimum_ordinate (str): value for minimal limit to represent on the
            ordinate horizontal axis
        maximum_ordinate (str): value for maximal limit to represent on the
            ordinate horizontal axis
        title_chart (str): title for chart
        title_abscissa (str): title for abscissa horizontal axis
        title_ordinate (str): title for ordinate vertical axis
        lines_thresholds (str): whether to draw lines to represent thresholds
        label_emphasis (str): whether to create text labels adjacent to
            the special selection of points for special emphasis
        label_count (str): whether to create text labels on chart to report
            counts of fold changes that pass thresholds
        report (str): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information
    """

    # Bundle information.
    pail = dict()
    # Parse information.
    pail["path_file_source_table_effects"] = str(
        path_file_source_table_effects
    ).strip()
    pail["path_file_source_list_features_emphasis"] = str(
        path_file_source_list_features_emphasis
    ).strip()
    pail["name_file_product_chart"] = str(
        name_file_product_chart
    ).strip()

    pail["path_directory_source"] = str(path_directory_source).strip()
    pail["path_directory_product"] = str(path_directory_product).strip()
    pail["path_directory_dock"] = str(path_directory_dock).strip()

    pail["column_identifier"] = str(column_identifier).strip()
    pail["column_name"] = str(column_name).strip()
    pail["column_effect"] = str(column_effect).strip()
    pail["column_p_value"] = str(column_p_value).strip()
    pail["column_significance"] = str(column_significance).strip()

    pail["threshold_effect"] = float(str(threshold_effect).strip())
    pail["threshold_significance"] = float(str(threshold_significance).strip())
    pail["minimum_abscissa"] = float(str(minimum_abscissa).strip())
    pail["maximum_abscissa"] = float(str(maximum_abscissa).strip())
    pail["minimum_ordinate"] = float(str(minimum_ordinate).strip())
    pail["maximum_ordinate"] = float(str(maximum_ordinate).strip())

    pail["title_chart"] = str(title_chart).strip().replace("#", " ")
    pail["title_abscissa"] = str(title_abscissa).strip().replace("#", " ")
    pail["title_ordinate"] = str(title_ordinate).strip().replace("#", " ")

    if (
        (lines_thresholds is not None) and
        (str(lines_thresholds) != "") and
        (str(lines_thresholds) != "none") and
        (str(lines_thresholds) == "true")
    ):
        pail["lines_thresholds"] = True
    else:
        pail["lines_thresholds"] = False
        pass
    if (
        (label_emphasis is not None) and
        (str(label_emphasis) != "") and
        (str(label_emphasis) != "none") and
        (str(label_emphasis) == "true")
    ):
        pail["label_emphasis"] = True
    else:
        pail["label_emphasis"] = False
        pass
    if (
        (label_count is not None) and
        (str(label_count) != "") and
        (str(label_count) != "none") and
        (str(label_count) == "true")
    ):
        pail["label_count"] = True
    else:
        pail["label_count"] = False
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

    # Collect names of relevant columns in original source table.
    pail["columns_source_text"] = list()
    pail["columns_source_number"] = list()
    pail["columns_source_text"].append(pail["column_identifier"])
    pail["columns_source_text"].append(pail["column_name"])
    pail["columns_source_number"].append(pail["column_effect"])
    pail["columns_source_number"].append(pail["column_p_value"])
    pail["columns_source_number"].append(pail["column_significance"])

    # Report.
    if pail["report"]:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: plot_chart_volcano.py")
        print("function: parse_text_parameters()")
        putly.print_terminal_partition(level=5)
        print("parameters:")
        print(pail)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


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
    path_file_source_table_effects=None,
    path_file_source_list_features_emphasis=None,
    path_directory_dock=None,
    columns_source_text=None,
    columns_source_number=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 2 July 2025

    arguments:
        path_file_source_table_effects (str): path to source file
        path_file_source_list_features_emphasis (str): path to source file
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        columns_source_text (list<str>): names of relevant columns in table
        columns_source_number (list<str>): names of relevant columns in table
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Define types of variables in columns of table.
    types_columns = define_column_types_table_source(
        columns_source_text=columns_source_text,
        columns_source_number=columns_source_number,
    )
    # Determine paths point to files that exist.
    existence_effects = os.path.exists(path_file_source_table_effects)
    existence_features = os.path.exists(
        path_file_source_list_features_emphasis
    )
    # Read information from file.
    if (existence_effects):
        table_effects = pandas.read_csv(
            path_file_source_table_effects,
            sep="\t",
            header=0,
            dtype=types_columns,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        table_effects = None
        pass
    if (existence_features):
        # Read information from file.
        features_emphasis = putly.read_file_text_list(
            path_file=path_file_source_list_features_emphasis,
            delimiter="\n",
            unique=True,
        )
    else:
        features_emphasis = list()
        pass

    # Bundle information.
    pail = dict()
    pail["table_effects"] = table_effects
    pail["features_emphasis"] = features_emphasis

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: plot_chart_volcano.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("table of effects:")
        print(pail["table_effects"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def create_write_plot_chart_volcano(
    table=None,
    column_identifier=None,
    column_name=None,
    column_plot_effect=None,
    column_plot_p=None,
    column_effect=None,
    column_significance=None,
    threshold_effect=None,
    threshold_significance=None,
    identifiers_emphasis=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    lines_thresholds=None,
    label_emphasis=None,
    label_count=None,
    name_file_product_chart=None,
    path_directory_product=None,
    report=None,
):
    """
    Create and write to file a plot chart of type volcano.

    Review: TCW; 16 July 2025

    arguments:
        table (object): Pandas data-frame table of features across columns and
            values for observations across rows
        column_identifier (str): name of column in table for the unique
            identifier corresponding to the fold change
        column_name (str): name of column in table for the name corresponding
            to the fold change
        column_plot_effect (str): name of column in table for values of fold
            change for representation on the plot
        column_plot_p (str): name of column in table for values of p-value or
            representation of significance for representation on the plot
        column_effect (str): name of column in table on which to apply the
            threshold for the effect magnitude and direction
        column_significance (str): name of column in table on which to apply
            the threshold for significance, corresponding to the p-value or
            q-value corresponding to the estimate of effect
        threshold_effect (float): value for threshold on magnitude of effect
            (|effect| > threshold) that is on the same scale, such as
            base-two logarithm, as the actual values themselves
        threshold_significance (float): value for threshold on p-values or
            q-values (p-value or q-value < threshold) that is not on a scale of
            the negative logarithm
        identifiers_emphasis (list<str>): identifiers corresponding to a
            special selection of effects for which to emphasize points on chart
            and for which to create text labels adjacent to the points on the
            chart
        minimum_abscissa (float): value for minimal limit to represent on the
            abscissa horizontal axis
        maximum_abscissa (float): value for maximal limit to represent on the
            abscissa horizontal axis
        minimum_ordinate (float): value for minimal limit to represent on the
            ordinate horizontal axis
        maximum_ordinate (float): value for maximal limit to represent on the
            ordinate horizontal axis
        title_chart (str): title for chart figure
        title_abscissa (str): title for abscissa horizontal axis
        title_ordinate (str): title for ordinate vertical axis
        lines_thresholds (bool): whether to draw lines to represent thresholds
        label_emphasis (bool): whether to create text labels adjacent to
            the special selection of points for special emphasis
        label_count (bool): whether to create text labels on chart to report
            counts of fold changes that pass thresholds
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

    # Create plot chart figure.
    figure = pplot.plot_scatter_fold_change_volcano(
        table=table,
        column_identifier=column_identifier,
        column_name=column_name,
        column_plot_effect=column_plot_effect,
        column_plot_p=column_plot_p,
        column_effect=column_effect,
        column_significance=column_significance,
        threshold_effect=threshold_effect,
        threshold_significance=threshold_significance,
        identifiers_emphasis=identifiers_emphasis,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        minimum_ordinate=minimum_ordinate,
        maximum_ordinate=maximum_ordinate,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        lines_thresholds=lines_thresholds,
        label_emphasis=label_emphasis,
        label_count=label_count,
        size_title_abscissa="seven", # ten
        size_title_ordinate="seven", # ten
        size_label_abscissa="eleven", # multi-panel: ten; individual: twelve
        size_label_ordinate="eleven", # multi-panel: ten; individual: twelve
        size_label_emphasis="thirteen", # eleven
        size_label_count="seven",
        aspect="landscape", # square, portrait, landscape, ...
        fonts=fonts,
        colors=colors,
        report=report,
    )

    ##########
    # Collect information.
    pail_write_plot = dict()
    pail_write_plot[name_file_product_chart] = figure

    ##########
    # Write product information to file.

    # Define paths to directories.
    #path_directory = os.path.join(
    #    path_directory_product, "forest_plot",
    #)
    # Create directories.
    putly.create_directories(
        path=path_directory_product,
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
        resolution=100,
        path_directory=path_directory_product,
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: plot_chart_volcano.py")
        print("function: create_write_plot_chart_volcano()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    pass



################################################################################
# Procedure


##########
# Call main procedure.


def execute_procedure(
    path_file_source_table_effects=None,
    path_file_source_list_features_emphasis=None,
    name_file_product_chart=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    column_identifier=None,
    column_name=None,
    column_effect=None,
    column_p_value=None,
    column_significance=None,
    threshold_effect=None,
    threshold_significance=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    lines_thresholds=None,
    label_emphasis=None,
    label_count=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review: TCW; 16 July 2025

    arguments:
        path_file_source_table_effects (str): path to source file
        path_file_source_list_features_emphasis (str): path to source file
        name_file_product_chart (str): name of product file
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        column_identifier (str): name of column in source table
        column_name (str): name of column in source table
        column_effect (str): name of column in table on which to apply the
            threshold for the effect magnitude and direction
        column_p_value (str): name of column in source table
        column_significance (str): name of column in table on which to apply
            the threshold for significance, corresponding to the p-value or
            q-value corresponding to the estimate of effect
        threshold_effect (float): value for threshold on effect
            (|effect| >= threshold) that is on the same scale, such as base-two
            logarithm, as the actual values themselves
        threshold_significance (float): value for threshold on p-values or
            q-values (p-value or q-value < threshold) that is not on a scale of
            the negative logarithm
        minimum_abscissa (float): value for minimal limit to represent on the
            abscissa horizontal axis
        maximum_abscissa (float): value for maximal limit to represent on the
            abscissa horizontal axis
        minimum_ordinate (float): value for minimal limit to represent on the
            ordinate horizontal axis
        maximum_ordinate (float): value for maximal limit to represent on the
            ordinate horizontal axis
        title_chart (str): title for chart figure
        title_abscissa (str): title for abscissa horizontal axis
        title_ordinate (str): title for ordinate vertical axis
        lines_thresholds (bool): whether to draw lines to represent thresholds
        label_emphasis (bool): whether to create text labels adjacent to
            the special selection of points for special emphasis
        label_count (bool): whether to create text labels on chart to report
            counts of fold changes that pass thresholds
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parse parameters.
    pail_parameters = parse_text_parameters(
        path_file_source_table_effects=path_file_source_table_effects,
        path_file_source_list_features_emphasis=(
            path_file_source_list_features_emphasis
        ),
        name_file_product_chart=name_file_product_chart,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        column_identifier=column_identifier,
        column_name=column_name,
        column_effect=column_effect,
        column_p_value=column_p_value,
        column_significance=column_significance,
        threshold_effect=threshold_effect,
        threshold_significance=threshold_significance,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        minimum_ordinate=minimum_ordinate,
        maximum_ordinate=maximum_ordinate,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        lines_thresholds=lines_thresholds,
        label_emphasis=label_emphasis,
        label_count=label_count,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: plot_chart_volcano.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("path_directory_product: " + str(path_directory_product))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_file_source_table_effects=(
            pail_parameters["path_file_source_table_effects"]
        ),
        path_file_source_list_features_emphasis=(
            pail_parameters["path_file_source_list_features_emphasis"]
        ),
        path_directory_dock=path_directory_dock,
        columns_source_text=pail_parameters["columns_source_text"],
        columns_source_number=pail_parameters["columns_source_number"],
        report=pail_parameters["report"],
    )
    #pail_source["table_effects"]
    #pail_source["features_emphasis"]

    ##########
    # Create and write to file plot chart volcano.
    create_write_plot_chart_volcano(
        table=pail_source["table_effects"],
        column_identifier=pail_parameters["column_identifier"],
        column_name=pail_parameters["column_name"],
        column_plot_effect=pail_parameters["column_effect"],
        column_plot_p=pail_parameters["column_p_value"],
        column_effect=pail_parameters["column_effect"],
        column_significance=pail_parameters["column_significance"],
        threshold_effect=pail_parameters["threshold_effect"],
        threshold_significance=pail_parameters["threshold_significance"],
        identifiers_emphasis=pail_source["features_emphasis"],
        minimum_abscissa=pail_parameters["minimum_abscissa"],
        maximum_abscissa=pail_parameters["maximum_abscissa"],
        minimum_ordinate=pail_parameters["minimum_ordinate"],
        maximum_ordinate=pail_parameters["maximum_ordinate"],
        title_chart=pail_parameters["title_chart"],
        title_abscissa=pail_parameters["title_abscissa"],
        title_ordinate=pail_parameters["title_ordinate"],
        lines_thresholds=pail_parameters["lines_thresholds"],
        label_emphasis=pail_parameters["label_emphasis"],
        label_count=pail_parameters["label_count"],
        name_file_product_chart=pail_parameters["name_file_product_chart"],
        path_directory_product=pail_parameters["path_directory_product"],
        report=pail_parameters["report"],
    )


    pass


# Execute program process in Python.



if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_source_table_effects = sys.argv[1]
    path_file_source_list_features_emphasis = sys.argv[2]
    name_file_product_chart = sys.argv[3]
    path_directory_source = sys.argv[4]
    path_directory_product = sys.argv[5]
    path_directory_dock = sys.argv[6]
    column_identifier = sys.argv[7]
    column_name = sys.argv[8]
    column_effect = sys.argv[9]
    column_p_value = sys.argv[10]
    column_significance = sys.argv[11]
    threshold_effect = sys.argv[12]
    threshold_significance = sys.argv[13]
    minimum_abscissa = sys.argv[14]
    maximum_abscissa = sys.argv[15]
    minimum_ordinate = sys.argv[16]
    maximum_ordinate = sys.argv[17]
    title_chart = sys.argv[18]
    title_abscissa = sys.argv[19]
    title_ordinate = sys.argv[20]
    lines_thresholds = sys.argv[21]
    label_emphasis = sys.argv[22]
    label_count = sys.argv[23]
    report = sys.argv[24]

    # Call function for procedure.
    execute_procedure(
        path_file_source_table_effects=path_file_source_table_effects,
        path_file_source_list_features_emphasis=(
            path_file_source_list_features_emphasis
        ),
        name_file_product_chart=name_file_product_chart,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        column_identifier=column_identifier,
        column_name=column_name,
        column_effect=column_effect,
        column_p_value=column_p_value,
        column_significance=column_significance,
        threshold_effect=threshold_effect,
        threshold_significance=threshold_significance,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        minimum_ordinate=minimum_ordinate,
        maximum_ordinate=maximum_ordinate,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        lines_thresholds=lines_thresholds,
        label_emphasis=label_emphasis,
        label_count=label_count,
        report=report,
    )

    pass



#
