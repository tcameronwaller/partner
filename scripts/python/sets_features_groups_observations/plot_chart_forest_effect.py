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
# Date, initialization: 31 October 2025
# Date, review or revision: 25 February 2026
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
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_features_observations=None,
    column_filter_binary=None,
    column_category=None,
    column_response_identifier=None,
    column_response_name=None,
    column_effect_primary=None,
    column_interval_low_primary=None,
    column_interval_high_primary=None,
    column_effect_secondary=None,
    column_interval_low_secondary=None,
    column_interval_high_secondary=None,
    column_effect_tertiary=None,
    column_interval_low_tertiary=None,
    column_interval_high_tertiary=None,
    name_chart=None,
    title_chart=None,
    title_legend=None,
    title_abscissa=None,
    title_ordinate=None,
    label_effect_primary=None,
    label_effect_secondary=None,
    label_effect_tertiary=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    size_marker=None,
    size_edge_marker=None,
    colors_fill_markers=None,
    color_edge_markers=None,
    color_edge_intervals=None,
    show_legend=None,
    name_batch=None,
    categories_batch=None,
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
    pail["path_directory_dock"] = str(path_directory_dock).strip()
    pail["path_directory_dock_pail"] = str(path_directory_dock).strip()
    pail["path_directory_source"] = str(path_directory_source).strip()
    pail["path_directory_product"] = str(path_directory_product).strip()
    # Paths to files.
    pail["path_file_source_table_features_observations"] = str(
        path_file_source_table_features_observations
    ).strip()

    # Names and categories.
    # It is problematic to pass any white space in parameters from a script in
    # Bash. Designate the hash symbol "#" as a substitute for white space.
    # It is also problematic to pass an empty string in parameters from a
    # script in Bash. Designate the word "none" as a substitute for missing or
    # empty.
    # Iterate on individual names that could be empty or missing.
    names_categories = {
        "column_filter_binary": column_filter_binary,
        "column_category": column_category,
        "column_response_identifier": column_response_identifier,
        "column_response_name": column_response_name,
        "column_effect_primary": column_effect_primary,
        "column_interval_low_primary": column_interval_low_primary,
        "column_interval_high_primary": column_interval_high_primary,
        "column_effect_secondary": column_effect_secondary,
        "column_interval_low_secondary": column_interval_low_secondary,
        "column_interval_high_secondary": column_interval_high_secondary,
        "column_effect_tertiary": column_effect_tertiary,
        "column_interval_low_tertiary": column_interval_low_tertiary,
        "column_interval_high_tertiary": column_interval_high_tertiary,
        "name_chart": name_chart,
        "title_chart": title_chart,
        "title_legend": title_legend,
        "title_abscissa": title_abscissa,
        "title_ordinate": title_ordinate,
        "label_effect_primary": label_effect_primary,
        "label_effect_secondary": label_effect_secondary,
        "label_effect_tertiary": label_effect_tertiary,
        "name_batch": name_batch,
        "categories_batch": categories_batch,
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

    # Number.
    pail["minimum_abscissa"] = float(str(minimum_abscissa).strip())
    pail["maximum_abscissa"] = float(str(maximum_abscissa).strip())
    pail["size_marker"] = float(str(size_marker).strip())
    pail["size_edge_marker"] = float(str(size_edge_marker).strip())

    # Lists, simple text.
    # Iterate on individual parameters for names and categories.
    lists = {
        "categories_batch": categories_batch,
    }
    for key_list in lists.keys():
        # Parse information.
        pail[key_list] = putly.parse_text_list_values(
            text=lists[key_list],
            delimiter=",",
        )
        pail[key_list] = putly.collect_unique_items(
            items=pail[key_list],
        )
        pass

    # List, objects.

    # Iterate on lists of colors.
    colors_lists = {
        "colors_fill_markers": colors_fill_markers,
    }
    for key_colors in colors_lists.keys():
        # Determine whether parameter has a valid value.
        if (
            (colors_lists[key_colors] != "none")
        ):
            # Determine type of definitions of colors.
            colors_split_simple = putly.parse_text_list_values(
                text=colors_lists[key_colors],
                delimiter=";",
            )
            if (
                ("(" not in colors_split_simple[0]) and
                (")" not in colors_split_simple[0])
            ):
                pail[key_colors] = colors_split_simple
            else:
                # Collect color items in a list.
                pail[key_colors] = list()
                for color_item in colors_split_simple:
                    color_tuple = putly.parse_extract_text_tuple(
                        text=color_item,
                        delimiter=",",
                        type_value="float",
                    )
                    pail[key_colors].append(color_tuple)
                    pass
                pass
        else:
            pail[key_colors] = None
            pass

    # Iterate on individual colors.
    colors = {
        "color_edge_markers": color_edge_markers,
        "color_edge_intervals": color_edge_intervals,
    }
    for key_color in colors.keys():
        # Determine whether parameter has a valid value.
        if (
            (colors[key_color] != "none")
        ):
            # Determine type of definitions of colors.
            if (
                ("(" not in colors[key_color]) and
                (")" not in colors[key_color])
            ):
                pail[key_color] = colors[key_color]
            else:
                # Collect color items in a list.
                pail[key_color] = putly.parse_extract_text_tuple(
                        text=colors[key_color],
                        delimiter=",",
                        type_value="float",
                    )
                pass
        else:
            pail[key_color] = None
            pass

    # Boolean, true or false.
    # Iterate on individual of Boolean designations.
    designations = {
        "show_legend": show_legend,
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
            "plot_chart_scatter_response.py"
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


def read_source(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_features_observations=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Date, revision or review: 25 February 2026

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_file_source_table_features_observations (str): path to source file
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Determine whether paths point to a directory or file that exist.
    #existence_directory = os.path.exists(path_directory)
    existence_file = os.path.exists(
        path_file_source_table_features_observations
    )

    # Bundle information.
    pail = dict()

    # Read information from file.
    if (existence_file):
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
    else:
        pail["table_features_observations"] = None
        pass

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_scatter_response.py"
        )
        print(str("module: " + module))
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def organize_table_features_observations(
    table=None,
    column_filter_binary=None,
    column_category=None,
    name_batch=None,
    categories_batch=None,
    filter_rows=None,
    report=None,
):
    """
    Organize and filter information in a table of parameters.

    Date, revision or review: 25 February 2026

    arguments:
        table (object): Pandas data-frame table
        column_filter_binary (str): name of column in source table
        column_category (str): name of column in source table
        name_batch (str): name for a set or group of categories that designate
            instances of parameters in a batch for execution
        categories_batch (list<str>): names of categories that designate sets
            or groups of instances of parameters in a batch for execution
        filter_rows (bool): whether to filter rows in table
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Copy information.
    table = table.copy(deep=True)

    # Organize information.
    table[column_filter_binary] = pandas.to_numeric(
        table[column_filter_binary],
        downcast="integer",
        errors="coerce",
    )

    # Filter rows in table by names of categories.
    if filter_rows:
        table = table.loc[(
            (table[column_filter_binary] == 1) &
            (table["category"].isin(categories_batch))
        ), :].copy(deep=True)
        pass

    # Bundle information.
    pail = dict()
    pail["table"] = table

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: drive_correlations_from_table_parameters.py")
        print("function: organize_filter_table_parameters()")
        putly.print_terminal_partition(level=5)
        print(str("name_batch: " + name_batch))
        print("categories_batch:")
        print(categories_batch)
        putly.print_terminal_partition(level=5)
        print("parameter table:")
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


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
    title_legend=None,
    legend_series_primary=None,
    legend_series_secondary=None,
    legend_series_tertiary=None,
    aspect=None,
    minimum_abscissa=None,
    center_abscissa=None,
    maximum_abscissa=None,
    position_line_origin=None,
    factor_space_series=None,
    space_between_series=None,
    size_title_chart=None,
    size_title_legend=None,
    size_title_abscissa=None,
    size_title_ordinate=None,
    size_label_abscissa=None,
    size_label_ordinate=None,
    size_label_legend=None,
    size_marker_primary=None,
    size_marker_secondary=None,
    size_marker_tertiary=None,
    size_edge_marker=None,
    size_line_origin=None,
    size_line_interval=None,
    colors_fill_markers=None,
    color_edge_markers=None,
    color_edge_intervals=None,
    fonts=None,
    show_legend=None,
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

    Date, revision or review: 25 February 2026
    Date, revision or review: 13 May 2025

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

        colors_fill_markers (list<tuple>): definitions of color properties
        color_edge_markers (tuple): definition of color properties
        color_edge_intervals (tuple): definition of color properties

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
        print("module: plot_chart_forest_effect.py")
        print("function: plot_dot_forest_category_ordinate_three_series()")
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Organize information for chart.

    # Copy information.
    table = table.copy(deep=True)
    columns_available = copy.deepcopy(table.columns.to_list())

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
    # Organize parameters for visual representation.

    # Colors.
    # Determine whether the parameters define adequate custom colors.
    # Lists of colors.
    if (
        (colors_fill_markers is None) or
        (len(colors_fill_markers) < count_series)
    ):
        # Define colors.
        # Create discrete color map for categorical values.
        # color_map = matplotlib.pyplot.get_cmap()
        color_map_fill_markers = matplotlib.colormaps.get_cmap("tab10") # "Set1", "Set2", "Dark2", "tab10",
    else:
        color_map_fill_markers = matplotlib.colors.ListedColormap(
            colors_fill_markers
        )
    # Individual colors.
    # Bundle information.
    pail_colors = dict()
    # Iterate on individual parameters for red, blue, green, and alpha channels
    # in colors.
    colors = {
        "edge_markers": color_edge_markers,
        "edge_intervals": color_edge_intervals,
    }
    for key_color in colors.keys():
        # Parse information.
        if (colors[key_color] is None):
            pail_colors[key_color] = matplotlib.colors.to_rgba("black", 1.0)
        elif (not isinstance(colors[key_color], tuple)):
            pail_colors[key_color] = matplotlib.colors.to_rgba(
                colors[key_color], 1.0
            )
        else:
            pail_colors[key_color] = colors[key_color]
            pass
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
            backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
            color=matplotlib.colors.to_rgba("black", 1.0),
            fontproperties=fonts["properties"][size_title_abscissa]
        )
    if (len(title_ordinate) > 0):
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=15,
            alpha=1.0,
            backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
            color=matplotlib.colors.to_rgba("black", 1.0),
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
        length=5.0, # 5.0
        width=2.5, # 3.0, 5.0
        pad=10.0, # 5.0, 7.5
        color=matplotlib.colors.to_rgba("black", 1.0),
        labelcolor=matplotlib.colors.to_rgba("black", 1.0),
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
        [
            round(minimum_abscissa, 1),
            round(center_abscissa, 1),
            round(maximum_abscissa, 1),
        ],
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
        color=matplotlib.colors.to_rgba("black", 1.0),
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
            ecolor=pail_colors["edge_intervals"],
            elinewidth=size_line_interval, # 7.5
            barsabove=False, # whether to print error bars in layer above points
            linestyle="",
            linewidth=size_edge_marker,
            marker="o", # marker shape: circle
            markersize=size_marker_primary, # 5, 15, 50, 70
            markeredgecolor=pail_colors["edge_markers"], # colors["purple"],
            markerfacecolor=color_map_fill_markers(0), # colors["purple"],
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
            ecolor=pail_colors["edge_intervals"],
            elinewidth=size_line_interval, # 7.5
            barsabove=False, # whether to print error bars in layer above points,
            linestyle="",
            linewidth=size_edge_marker,
            marker="D", # "D" marker shape: diamond
            markersize=size_marker_secondary, # 5, 15
            markeredgecolor=pail_colors["edge_markers"], # colors["green"],
            markerfacecolor=color_map_fill_markers(1), # colors["green"],
        )
        pass
    if (count_series == 3):
        handle_tertiary = axes.errorbar(
            positions_abscissa_tertiary,
            positions_ordinate_tertiary,
            yerr=None,
            xerr=intervals_abscissa_tertiary,
            ecolor=pail_colors["edge_intervals"],
            elinewidth=size_line_interval, # 7.5
            barsabove=False, # whether to print error bars in layer above points,
            linestyle="",
            linewidth=size_edge_marker,
            marker="s", # "s" marker shape: square
            markersize=size_marker_tertiary, # 5, 15
            markeredgecolor=pail_colors["edge_markers"], # colors["green"],
            markerfacecolor=color_map_fill_markers(2), # colors["green"],
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
                ecolor=pail_colors["edge_intervals"],
                elinewidth=(size_line_interval/1.3),
                barsabove=False, # whether to print error bars in layer above points
                linestyle="",
                linewidth=(size_edge_marker/1.3),
                marker="o", # "o" marker shape: circle
                markersize=(size_marker_primary/1.3), # 5, 15, 50, 70
                markeredgecolor=pail_colors["edge_markers"], # colors["purple"],
                markerfacecolor=color_map_fill_markers(0), # colors["purple"],
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
                ecolor=pail_colors["edge_intervals"],
                elinewidth=(size_line_interval/1.3),
                barsabove=False, # whether to print error bars in layer above points
                linestyle="",
                linewidth=(size_edge_marker/1.3),
                marker="D", # "D" marker shape: diamond
                markersize=(size_marker_secondary/1.3), # 5, 15
                markeredgecolor=pail_colors["edge_markers"], # colors["green"],
                markerfacecolor=color_map_fill_markers(1), # colors["green"],
            )
            pass
        if (count_series == 3):
            handle_tertiary = axes.errorbar(
                [0],
                [-1], # create outside of visible portion of axes
                yerr=None,
                xerr=5.0,
                label=legend_series_tertiary,
                ecolor=pail_colors["edge_intervals"],
                elinewidth=(size_line_interval/1.3),
                barsabove=False, # whether to print error bars in layer above points
                linestyle="",
                linewidth=(size_edge_marker/1.3),
                marker="s", # "s" marker shape: square
                markersize=(size_marker_tertiary/1.3), # 5, 15
                markeredgecolor=pail_colors["edge_markers"], # colors["green"],
                markerfacecolor=color_map_fill_markers(2), # colors["green"],
            )
            pass
        if (
            (show_legend) and
            (count_series == 1)
        ):
            axes.legend(
                handles=[
                    handle_primary,
                ],
                loc="upper right",
                prop=fonts["properties"][size_label_legend],
                title=title_legend,
                title_fontsize=fonts["values"][size_title_legend]["size"]
            )
            pass
        elif (
            (show_legend) and
            (count_series == 2)
        ):
            axes.legend(
                handles=[
                    handle_primary, handle_secondary,
                ],
                loc="upper right",
                prop=fonts["properties"][size_label_legend],
                title=title_legend,
                title_fontsize=fonts["values"][size_title_legend]["size"]
            )
            pass
        elif (
            (show_legend) and
            (count_series == 3)
        ):
            axes.legend(
                handles=[
                    handle_primary, handle_secondary, handle_tertiary,
                ],
                loc="upper right",
                prop=fonts["properties"][size_label_legend],
                title=title_legend,
                title_fontsize=fonts["values"][size_title_legend]["size"]
            )
            pass
        pass

    # Return figure.
    return figure


def create_write_plot_chart_dot_forest(
    path_directory_parent=None,
    name_chart=None,
    table=None,
    column_response_identifier=None,
    column_response_name=None,
    column_effect_primary=None,
    column_interval_low_primary=None,
    column_interval_high_primary=None,
    column_effect_secondary=None,
    column_interval_low_secondary=None,
    column_interval_high_secondary=None,
    column_effect_tertiary=None,
    column_interval_low_tertiary=None,
    column_interval_high_tertiary=None,
    title_chart=None,
    title_legend=None,
    title_abscissa=None,
    title_ordinate=None,
    label_effect_primary=None,
    label_effect_secondary=None,
    label_effect_tertiary=None,
    minimum_abscissa=None,
    center_abscissa=None,
    maximum_abscissa=None,
    size_marker=None,
    size_edge_marker=None,
    colors_fill_markers=None,
    color_edge_markers=None,
    color_edge_intervals=None,
    show_legend=None,
    report=None,
):
    """
    Create, plot, and write to file a chart of the type dot forest.

    Date, revision or review: 25 February 2026

    arguments:
        path_directory_parent (str): path to parent directory for procedure's
            product directories and files
        name_chart (str): name for writing figure object to file
        table (object): Pandas data-frame table of features across columns and
            observations across rows with values on quantitative, continuous
            interval or ratio scales of measurement
        column_identifier (str): name of column in table corresponding to the
            unique identifier of records for each point
        column_name (str): name of column in table corresponding to the name of
            records for each point
        column_response_markers (str): name of column in table corresponding to
            values for representation as color of markers for individual
            points; either categorical groups or values on a quantitative,
            continuous scale of measurement
        column_group_ellipses (str): name of column in table corresponding to
            names or categories that designate groups of observations for which
            to create confidence ellipses
        column_abscissa (str): name of column in table corresponding to values
            for representation on the abscissa horizontal axis
        column_ordinate (str): name of column in table corresponding to values
            for representation on the ordinate vertical axis
        type_response (bool): type of the response feature, either 'continuity'
            or 'category'
        title_chart (str): title of the chart
        title_response (str): title of the feature to represent as colors of
            individual points
        title_abscissa (str): title of the feature to represent on the abscissa
            horizontal axis
        title_ordinate (str): title of the feature to represent on the ordinate
            vertical axis
        identifiers_emphasis (list<str>): identifiers corresponding to a
            special selection of records for which to emphasize points on
            chart and for which to create individual text labels adjacent to
            the points on the chart
        size_marker (int): size of markers for points representing values
        size_edge_marker (float): size for line width of edge
        size_edge_ellipse (float): size for line width of edge
        factor_confidence_ellipse (float): factor product for standard
            deviations in the confidence ellipses
        colors_fill_markers (list<tuple>): definitions of color properties
        colors_fill_ellipses (list<tuple>): definitions of color properties
        color_edge_markers (tuple): definition of color properties
        color_edge_ellipses (tuple): definition of color properties
        color_emphasis (tuple): definition of color properties
        show_confidence_ellipse (bool): whether to show confidence ellipse;
            show empty ellipse edge for continuous response; show color-matched
            distinct ellipses for discrete categorical groups
        show_emphasis_marker (bool): whether to create special markers to emphasize
            a special selection of points
        show_emphasis_label (bool): whether to create text labels adjacent to
            the special selection of points for special emphasis
        show_legend_bar (bool): whether to show legend or scale bar on chart
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object from MatPlotLib

    """

    ##########
    # Organize information for plot.

    # Copy information in table.
    table = table.copy(deep=True)

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    #colors = pplot.define_color_properties()
    # Create figure.
    figure = plot_dot_forest_category_ordinate_three_series(
        table=table,
        column_feature=column_response_identifier,
        column_feature_name=column_response_name,
        column_value_primary=column_effect_primary,
        column_interval_low_primary=column_interval_low_primary,
        column_interval_high_primary=column_interval_high_primary,
        column_value_secondary=column_effect_secondary,
        column_interval_low_secondary=column_interval_low_secondary,
        column_interval_high_secondary=column_interval_high_secondary,
        column_value_tertiary=column_effect_tertiary,
        column_interval_low_tertiary=column_interval_low_tertiary,
        column_interval_high_tertiary=column_interval_high_tertiary,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        title_legend=title_legend,
        legend_series_primary=label_effect_primary,
        legend_series_secondary=label_effect_secondary,
        legend_series_tertiary=label_effect_tertiary,
        aspect="portrait",
        minimum_abscissa=minimum_abscissa,
        center_abscissa=center_abscissa,
        maximum_abscissa=maximum_abscissa,
        position_line_origin=center_abscissa,
        factor_space_series=None, # if not zero or None, overrides space
        space_between_series=0.33,
        size_title_chart="nine",
        size_title_legend="eleven",
        size_title_abscissa="ten",
        size_title_ordinate="ten",
        size_label_abscissa="eleven",
        size_label_ordinate="eleven",
        size_label_legend="fourteen",
        size_marker_primary=size_marker,
        size_marker_secondary=float(size_marker * 0.7),
        size_marker_tertiary=float(size_marker * 0.7),
        size_edge_marker=size_edge_marker,
        size_line_origin=3.5,
        size_line_interval=2, # 3.5
        colors_fill_markers=colors_fill_markers,
        color_edge_markers=color_edge_markers,
        color_edge_intervals=color_edge_intervals,
        fonts=fonts,
        show_legend=show_legend,
        report=report,
    )

    # Write product information to file.

    # Bundle information.
    pail_write_plot = dict()
    pail_write_plot[name_chart] = figure

    # Write figure object to file.
    if True:
        pplot.write_product_plots_parent_directory(
            pail_write=pail_write_plot,
            format="jpg", # jpg, png, svg
            resolution=96, # 72, 96, 300
            path_directory=path_directory_parent,
        )
    if False:
        pplot.write_product_plots_parent_directory(
            pail_write=pail_write_plot,
            format="png", # jpg, png, svg
            resolution=150,
            path_directory=path_directory_parent,
        )
    if False:
        pplot.write_product_plots_parent_directory(
            pail_write=pail_write_plot,
            format="svg", # jpg, png, svg
            resolution=150,
            path_directory=path_directory_parent,
        )

    # Return information.
    return figure




################################################################################
# Procedure


##########
# Call main procedure.


def execute_procedure(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_features_observations=None,
    column_filter_binary=None,
    column_category=None,
    column_response_identifier=None,
    column_response_name=None,
    column_effect_primary=None,
    column_interval_low_primary=None,
    column_interval_high_primary=None,
    column_effect_secondary=None,
    column_interval_low_secondary=None,
    column_interval_high_secondary=None,
    column_effect_tertiary=None,
    column_interval_low_tertiary=None,
    column_interval_high_tertiary=None,
    name_chart=None,
    title_chart=None,
    title_legend=None,
    title_abscissa=None,
    title_ordinate=None,
    label_effect_primary=None,
    label_effect_secondary=None,
    label_effect_tertiary=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    size_marker=None,
    size_edge_marker=None,
    colors_fill_markers=None,
    color_edge_markers=None,
    color_edge_intervals=None,
    show_legend=None,
    name_batch=None,
    categories_batch=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Date, revision or review: 25 February 2026

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_file_source_table_features_observations (str): path to source file
        column_identifier_observation (str): name of column in source table
        column_identifier_signal (str): name of column in source table
        ...
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parse parameters.
    pail_parameters = parse_text_parameters(
        path_directory_dock=path_directory_dock,
        path_directory_dock_pail=path_directory_dock_pail,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        column_filter_binary=column_filter_binary,
        column_category=column_category,
        column_response_identifier=column_response_identifier,
        column_response_name=column_response_name,
        column_effect_primary=column_effect_primary,
        column_interval_low_primary=column_interval_low_primary,
        column_interval_high_primary=column_interval_high_primary,
        column_effect_secondary=column_effect_secondary,
        column_interval_low_secondary=column_interval_low_secondary,
        column_interval_high_secondary=column_interval_high_secondary,
        column_effect_tertiary=column_effect_tertiary,
        column_interval_low_tertiary=column_interval_low_tertiary,
        column_interval_high_tertiary=column_interval_high_tertiary,
        name_chart=name_chart,
        title_chart=title_chart,
        title_legend=title_legend,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        label_effect_primary=label_effect_primary,
        label_effect_secondary=label_effect_secondary,
        label_effect_tertiary=label_effect_tertiary,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        size_marker=size_marker,
        size_edge_marker=size_edge_marker,
        colors_fill_markers=colors_fill_markers,
        color_edge_markers=color_edge_markers,
        color_edge_intervals=color_edge_intervals,
        show_legend=show_legend,
        name_batch=name_batch,
        categories_batch=categories_batch,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_forest_effect.py"
        )
        print(str("module: " + module))
        print("function: execute_procedure()")
        print("system: local")
        putly.print_terminal_partition(level=5)
        pass


    ##########
    # Read source information from file.
    pail_source = read_source(
        path_directory_dock=pail_parameters["path_directory_dock"],
        path_directory_dock_pail=pail_parameters["path_directory_dock_pail"],
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_file_source_table_features_observations=(
            pail_parameters["path_file_source_table_features_observations"]
        ),
        report=pail_parameters["report"],
    )

    # Organize and filter table of parameters by categories in current batch.
    pail_batch = organize_table_features_observations(
        table=pail_source["table_features_observations"],
        column_filter_binary=pail_parameters["column_filter_binary"],
        column_category=pail_parameters["column_category"],
        name_batch=pail_parameters["name_batch"],
        categories_batch=pail_parameters["categories_batch"],
        filter_rows=True,
        report=pail_parameters["report"],
    )

    # Define paths to directories.
    path_directory_chart = os.path.join(
        pail_parameters["path_directory_product"],
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_chart,
    )

    # Create plot chart and write to file.
    create_write_plot_chart_dot_forest(
        path_directory_parent=path_directory_chart,
        name_chart=pail_parameters["name_chart"],
        table=pail_batch["table"],
        column_response_identifier=(
            pail_parameters["column_response_identifier"]
        ),
        column_response_name=pail_parameters["column_response_name"],
        column_effect_primary=pail_parameters["column_effect_primary"],
        column_interval_low_primary=(
            pail_parameters["column_interval_low_primary"]
        ),
        column_interval_high_primary=(
            pail_parameters["column_interval_high_primary"]
        ),
        column_effect_secondary=pail_parameters["column_effect_secondary"],
        column_interval_low_secondary=(
            pail_parameters["column_interval_low_secondary"]
        ),
        column_interval_high_secondary=(
            pail_parameters["column_interval_high_secondary"]
        ),
        column_effect_tertiary=pail_parameters["column_effect_tertiary"],
        column_interval_low_tertiary=(
            pail_parameters["column_interval_low_tertiary"]
        ),
        column_interval_high_tertiary=(
            pail_parameters["column_interval_high_tertiary"]
        ),
        title_chart=pail_parameters["title_chart"],
        title_legend=pail_parameters["title_legend"],
        title_abscissa=pail_parameters["title_abscissa"],
        title_ordinate=pail_parameters["title_ordinate"],
        label_effect_primary=pail_parameters["label_effect_primary"],
        label_effect_secondary=pail_parameters["label_effect_secondary"],
        label_effect_tertiary=pail_parameters["label_effect_tertiary"],
        minimum_abscissa=pail_parameters["minimum_abscissa"],
        center_abscissa=0.0,
        maximum_abscissa=pail_parameters["maximum_abscissa"],
        size_marker=pail_parameters["size_marker"],
        size_edge_marker=pail_parameters["size_edge_marker"],
        colors_fill_markers=pail_parameters["colors_fill_markers"],
        color_edge_markers=pail_parameters["color_edge_markers"],
        color_edge_intervals=pail_parameters["color_edge_intervals"],
        show_legend=pail_parameters["show_legend"],
        report=pail_parameters["report"],
    )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_dock = sys.argv[1]
    path_directory_dock_pail = sys.argv[2]
    path_directory_source = sys.argv[3]
    path_directory_product = sys.argv[4]
    path_file_source_table_features_observations = sys.argv[5]
    column_filter_binary = sys.argv[6]
    column_category = sys.argv[7]
    column_response_identifier = sys.argv[8]
    column_response_name = sys.argv[9]
    column_effect_primary = sys.argv[10]
    column_interval_low_primary = sys.argv[11]
    column_interval_high_primary = sys.argv[12]
    column_effect_secondary = sys.argv[13]
    column_interval_low_secondary = sys.argv[14]
    column_interval_high_secondary = sys.argv[15]
    column_effect_tertiary = sys.argv[16]
    column_interval_low_tertiary = sys.argv[17]
    column_interval_high_tertiary = sys.argv[18]
    name_chart = sys.argv[19]
    title_chart = sys.argv[20]
    title_legend = sys.argv[21]
    title_abscissa = sys.argv[22]
    title_ordinate = sys.argv[23]
    label_effect_primary = sys.argv[24]
    label_effect_secondary = sys.argv[25]
    label_effect_tertiary = sys.argv[26]
    minimum_abscissa = sys.argv[27]
    maximum_abscissa = sys.argv[28]
    size_marker = sys.argv[29]
    size_edge_marker = sys.argv[30]
    colors_fill_markers = sys.argv[31]
    color_edge_markers = sys.argv[32]
    color_edge_intervals = sys.argv[33]
    show_legend = sys.argv[34]
    name_batch = sys.argv[35]
    categories_batch = sys.argv[36]
    report = sys.argv[37]

    # Call function for procedure.
    execute_procedure(
        path_directory_dock=path_directory_dock,
        path_directory_dock_pail=path_directory_dock_pail,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        column_filter_binary=column_filter_binary,
        column_category=column_category,
        column_response_identifier=column_response_identifier,
        column_response_name=column_response_name,
        column_effect_primary=column_effect_primary,
        column_interval_low_primary=column_interval_low_primary,
        column_interval_high_primary=column_interval_high_primary,
        column_effect_secondary=column_effect_secondary,
        column_interval_low_secondary=column_interval_low_secondary,
        column_interval_high_secondary=column_interval_high_secondary,
        column_effect_tertiary=column_effect_tertiary,
        column_interval_low_tertiary=column_interval_low_tertiary,
        column_interval_high_tertiary=column_interval_high_tertiary,
        name_chart=name_chart,
        title_chart=title_chart,
        title_legend=title_legend,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        label_effect_primary=label_effect_primary,
        label_effect_secondary=label_effect_secondary,
        label_effect_tertiary=label_effect_tertiary,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        size_marker=size_marker,
        size_edge_marker=size_edge_marker,
        colors_fill_markers=colors_fill_markers,
        color_edge_markers=color_edge_markers,
        color_edge_intervals=color_edge_intervals,
        show_legend=show_legend,
        name_batch=name_batch,
        categories_batch=categories_batch,
        report=report,
    )

    pass



#
