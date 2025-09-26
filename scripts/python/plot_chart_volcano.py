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
    path_file_source_list_emphasis=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    name_file_product_chart=None,
    suffix_file_product_plot_chart=None,
    column_effect_identifier=None,
    column_effect_estimate=None,
    column_effect_error=None,
    column_effect_p=None,
    column_effect_q=None,
    column_feature_identifier=None,
    column_feature_name=None,
    column_plot_spread=None,
    column_plot_rise=None,
    column_significance=None,
    identifiers_emphasis=None,
    threshold_spread=None,
    threshold_rise=None,
    threshold_significance=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    title_abscissa=None,
    title_ordinate=None,
    size_mark_point=None,
    size_mark_point_emphasis=None,
    size_font_title_abscissa=None,
    size_font_title_ordinate=None,
    size_font_label_abscissa=None,
    size_font_label_ordinate=None,
    size_font_label_emphasis=None,
    size_font_label_count=None,
    line_threshold_spread=None,
    line_threshold_rise=None,
    label_emphasis=None,
    label_count=None,
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

    # Paths to files and directories.
    pail["path_file_source_table_effects"] = str(
        path_file_source_table_effects
    ).strip()
    pail["path_file_source_list_emphasis"] = str(
        path_file_source_list_emphasis
    ).strip()
    pail["path_directory_source"] = str(path_directory_source).strip()
    pail["path_directory_product"] = str(path_directory_product).strip()
    pail["path_directory_dock"] = str(path_directory_dock).strip()

    # Other.

    pail["name_file_product_chart"] = str(name_file_product_chart).strip()
    pail["suffix_file_product_plot_chart"] = str(
        suffix_file_product_plot_chart
    ).strip()

    pail["column_effect_identifier"] = str(column_effect_identifier).strip()
    pail["column_effect_estimate"] = str(column_effect_estimate).strip()
    pail["column_effect_error"] = str(column_effect_error).strip()
    pail["column_effect_p"] = str(column_effect_p).strip()
    pail["column_effect_q"] = str(column_effect_q).strip()
    pail["column_feature_identifier"] = str(column_feature_identifier).strip()
    pail["column_feature_name"] = str(column_feature_name).strip()
    pail["column_plot_spread"] = str(column_plot_spread).strip()
    pail["column_plot_rise"] = str(column_plot_rise).strip()
    pail["column_significance"] = str(column_significance).strip()

    pail["identifiers_emphasis"] = putly.parse_text_list_values(
        text=identifiers_emphasis,
        delimiter=",",
    )
    pail["identifiers_emphasis"] = putly.collect_unique_items(
        items=pail["identifiers_emphasis"],
    )

    pail["threshold_spread"] = float(str(threshold_spread).strip())
    pail["threshold_rise"] = float(str(threshold_rise).strip())
    pail["threshold_significance"] = float(str(threshold_significance).strip())
    pail["minimum_abscissa"] = float(str(minimum_abscissa).strip())
    pail["maximum_abscissa"] = float(str(maximum_abscissa).strip())
    pail["minimum_ordinate"] = float(str(minimum_ordinate).strip())
    pail["maximum_ordinate"] = float(str(maximum_ordinate).strip())

    #pail["title_chart"] = str(title_chart).strip().replace("#", " ")
    pail["title_abscissa"] = str(title_abscissa).strip().replace("#", " ")
    pail["title_ordinate"] = str(title_ordinate).strip().replace("#", " ")

    pail["size_mark_point"] = float(str(size_mark_point).strip())
    pail["size_mark_point_emphasis"] = float(str(
        size_mark_point_emphasis
    ).strip())
    pail["size_font_title_abscissa"] = str(size_font_title_abscissa).strip()
    pail["size_font_title_ordinate"] = str(size_font_title_ordinate).strip()
    pail["size_font_label_abscissa"] = str(size_font_label_abscissa).strip()
    pail["size_font_label_ordinate"] = str(size_font_label_ordinate).strip()
    pail["size_font_label_emphasis"] = str(size_font_label_emphasis).strip()
    pail["size_font_label_count"] = str(size_font_label_count).strip()

    # Boolean, true or false.
    if (
        (line_threshold_spread is not None) and
        (str(line_threshold_spread) != "") and
        (str(line_threshold_spread) != "none") and
        (str(line_threshold_spread) == "true")
    ):
        pail["line_threshold_spread"] = True
    else:
        pail["line_threshold_spread"] = False
        pass
    if (
        (line_threshold_rise is not None) and
        (str(line_threshold_rise) != "") and
        (str(line_threshold_rise) != "none") and
        (str(line_threshold_rise) == "true")
    ):
        pail["line_threshold_rise"] = True
    else:
        pail["line_threshold_rise"] = False
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
    pail["columns_source_text"].append(pail["column_effect_identifier"])
    pail["columns_source_text"].append(pail["column_feature_identifier"])
    pail["columns_source_text"].append(pail["column_feature_name"])
    pail["columns_source_text"] = putly.collect_unique_items(
        items=pail["columns_source_text"],
    )
    pail["columns_source_number"].append(pail["column_effect_estimate"])
    pail["columns_source_number"].append(pail["column_effect_error"])
    pail["columns_source_number"].append(pail["column_effect_p"])
    pail["columns_source_number"].append(pail["column_effect_q"])
    pail["columns_source_number"].append(pail["column_plot_spread"])
    pail["columns_source_number"].append(pail["column_plot_rise"])
    pail["columns_source_number"].append(pail["column_significance"])
    pail["columns_source_number"] = putly.collect_unique_items(
        items=pail["columns_source_number"],
    )

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
    path_file_source_list_emphasis=None,
    path_directory_source=None,
    path_directory_product=None,
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
        path_file_source_list_emphasis (str): path to source file
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        columns_source_text (list<str>): names of columns in table to read as
            text character strings
        columns_source_number (list<str>): names of columns in table to read as
            floating point numbers
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
        path_file_source_list_emphasis
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
            path_file=path_file_source_list_emphasis,
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
    path_file_source_list_emphasis=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    name_file_product_chart=None,
    suffix_file_product_plot_chart=None,
    column_effect_identifier=None,
    column_effect_estimate=None,
    column_effect_error=None,
    column_effect_p=None,
    column_effect_q=None,
    column_feature_identifier=None,
    column_feature_name=None,
    column_plot_spread=None,
    column_plot_rise=None,
    column_significance=None,
    identifiers_emphasis=None,
    threshold_spread=None,
    threshold_rise=None,
    threshold_significance=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    title_abscissa=None,
    title_ordinate=None,
    size_mark_point=None,
    size_mark_point_emphasis=None,
    size_font_title_abscissa=None,
    size_font_title_ordinate=None,
    size_font_label_abscissa=None,
    size_font_label_ordinate=None,
    size_font_label_emphasis=None,
    size_font_label_count=None,
    line_threshold_spread=None,
    line_threshold_rise=None,
    label_emphasis=None,
    label_count=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review: TCW; 19 September 2025

    arguments:
        path_file_source_table_effects (str): path to source file
        path_file_source_list_emphasis (str): path to source file
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        name_file_product_chart (str): name of product file
        suffix_file_product_plot_chart (str):

        TODO: TCW; 19 September 2025
        # need to update the documentation

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
        path_file_source_list_emphasis=path_file_source_list_emphasis,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        name_file_product_chart=name_file_product_chart,
        suffix_file_product_plot_chart=Nsuffix_file_product_plot_chartone,
        column_effect_identifier=column_effect_identifier,
        column_effect_estimate=column_effect_estimate,
        column_effect_error=column_effect_error,
        column_effect_p=column_effect_p,
        column_effect_q=column_effect_q,
        column_feature_identifier=column_feature_identifier,
        column_feature_name=column_feature_name,
        column_plot_spread=column_plot_spread,
        column_plot_rise=column_plot_rise,
        column_significance=column_significance,
        identifiers_emphasis=identifiers_emphasis,
        threshold_spread=threshold_spread,
        threshold_rise=threshold_rise,
        threshold_significance=threshold_significance,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        minimum_ordinate=minimum_ordinate,
        maximum_ordinate=maximum_ordinate,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        size_mark_point=size_mark_point,
        size_mark_point_emphasis=size_mark_point_emphasis,
        size_font_title_abscissa=size_font_title_abscissa,
        size_font_title_ordinate=size_font_title_ordinate,
        size_font_label_abscissa=size_font_label_abscissa,
        size_font_label_ordinate=size_font_label_ordinate,
        size_font_label_emphasis=size_font_label_emphasis,
        size_font_label_count=size_font_label_count,
        line_threshold_spread=line_threshold_spread,
        line_threshold_rise=line_threshold_rise,
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
        path_file_source_list_emphasis=(
            pail_parameters["path_file_source_list_emphasis"]
        ),
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        columns_source_text=pail_parameters["columns_source_text"],
        columns_source_number=pail_parameters["columns_source_number"],
        report=pail_parameters["report"],
    )
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1")
    print(pail["table_effects"])
    print(pail["features_emphasis"])

    sys.exit()

    # Combine identifiers of features for emphasis.

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
    path_file_source_list_emphasis = sys.argv[2]
    path_directory_source = sys.argv[3]
    path_directory_product = sys.argv[4]
    path_directory_dock = sys.argv[5]
    name_file_product_chart = sys.argv[6]
    suffix_file_product_plot_chart = sys.argv[7]
    column_effect_identifier = sys.argv[8]
    column_effect_estimate = sys.argv[9]
    column_effect_error = sys.argv[10]
    column_effect_p = sys.argv[11]
    column_effect_q = sys.argv[12]
    column_feature_identifier = sys.argv[13]
    column_feature_name = sys.argv[14]
    column_plot_spread = sys.argv[15]
    column_plot_rise = sys.argv[16]
    column_significance = sys.argv[17]
    identifiers_emphasis = sys.argv[18]
    threshold_spread = sys.argv[19]
    threshold_rise = sys.argv[20]
    threshold_significance = sys.argv[21]
    minimum_abscissa = sys.argv[22]
    maximum_abscissa = sys.argv[23]
    minimum_ordinate = sys.argv[24]
    maximum_ordinate = sys.argv[25]
    title_abscissa = sys.argv[26]
    title_ordinate = sys.argv[27]
    size_mark_point = sys.argv[28]
    size_mark_point_emphasis = sys.argv[29]
    size_font_title_abscissa = sys.argv[30]
    size_font_title_ordinate = sys.argv[31]
    size_font_label_abscissa = sys.argv[32]
    size_font_label_ordinate = sys.argv[33]
    size_font_label_emphasis = sys.argv[34]
    size_font_label_count = sys.argv[35]
    line_threshold_spread = sys.argv[36]
    line_threshold_rise = sys.argv[37]
    label_emphasis = sys.argv[38]
    label_count = sys.argv[39]
    report = sys.argv[40]

    # Call function for procedure.
    execute_procedure(
        path_file_source_table_effects=path_file_source_table_effects,
        path_file_source_list_emphasis=path_file_source_list_emphasis,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        name_file_product_chart=name_file_product_chart,
        suffix_file_product_plot_chart=Nsuffix_file_product_plot_chartone,
        column_effect_identifier=column_effect_identifier,
        column_effect_estimate=column_effect_estimate,
        column_effect_error=column_effect_error,
        column_effect_p=column_effect_p,
        column_effect_q=column_effect_q,
        column_feature_identifier=column_feature_identifier,
        column_feature_name=column_feature_name,
        column_plot_spread=column_plot_spread,
        column_plot_rise=column_plot_rise,
        column_significance=column_significance,
        identifiers_emphasis=identifiers_emphasis,
        threshold_spread=threshold_spread,
        threshold_rise=threshold_rise,
        threshold_significance=threshold_significance,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        minimum_ordinate=minimum_ordinate,
        maximum_ordinate=maximum_ordinate,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        size_mark_point=size_mark_point,
        size_mark_point_emphasis=size_mark_point_emphasis,
        size_font_title_abscissa=size_font_title_abscissa,
        size_font_title_ordinate=size_font_title_ordinate,
        size_font_label_abscissa=size_font_label_abscissa,
        size_font_label_ordinate=size_font_label_ordinate,
        size_font_label_emphasis=size_font_label_emphasis,
        size_font_label_count=size_font_label_count,
        line_threshold_spread=line_threshold_spread,
        line_threshold_rise=line_threshold_rise,
        label_emphasis=label_emphasis,
        label_count=label_count,
        report=report,
    )

    pass



#
