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
# Date, initialization: 28 October 2025
# Date, review or revision: 28 October 2025
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
#import partner.regression as preg
import partner.plot as pplot

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Colors.


def create_divergent_color_map(
    value_minimum=None,
    value_center=None,
    value_maximum=None,
    color_minimum=None,
    color_center=None,
    color_maximum=None,
):
    """
    Create a custom divergent color map.

    References:
    1. https://matplotlib.org/stable/users/explain/colors/
       colormap-manipulation.html

https://matplotlib.org/stable/users/explain/colors/colormap-manipulation.html

    Review or revision: TCW; 13 November 2025

    arguments:
        ...
    raises:

    returns:
        (object): MatPlotLib color map

    """

    # Organize anchors and colors for color map.
    anchors = [0.0, 0.5, 1.0,]
    colors = [color_minimum, color_center, color_maximum,]
    # Create color map.
    map = matplotlib.colors.LinearSegmentedColormap.from_list(
        "map_tcw_custom",
        list(zip(anchors, colors,)),
        N=256,
        gamma=1.0,
    )
    # Create normalization scale to fit the values to the color map.
    scale = matplotlib.colors.TwoSlopeNorm(
        vmin=value_minimum,
        vcenter=value_center,
        vmax=value_maximum,
    )
    # Bundle information.
    pail = dict()
    pail["map"] = map
    pail["scale"] = scale
    # Return object.
    return pail




##########
# Heatmap.
# 1.
# function: plot_heatmap_signal_features_observations_labels()
# signal: individual signals or summaries (such as mean or median) for small
# counts of features across observations or groups of observations
# design: explicit labels for individual features on vertical ordinate axis and
# for individual observations or groups on horizontal abscissa axis
# 2.
# function: plot_heatmap_signal_features_sets_observations_labels()
# signal: individual signals or summaries (such as mean or median) for small or
# large counts of features which belong to sets across observations or groups
# of observations
# design: graphical representation of allocation of each feature to overlapping
# sets and explicit labels for individual observations or groups
# 3.
# function: plot_heatmap_signal_features_sets_observations_groups()
# signal: individual signals for large counts of features which belong to sets
# across observations which belong to groups
# design: graphical representation of allocation of each feature to overlapping
# sets and graphical representation of observations in mutually exclusive
# groups


def extract_prepare_table_signals_categories_for_heatmap(
    table=None,
    format_table=None,
    index_columns=None,
    index_rows=None,
    column_group=None,
    transpose_table=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_center=None,
    value_maximum=None,
    labels_ordinate_categories=None,
    labels_abscissa_categories=None,
    report=None,
):
    """
    Extract from table and prepare signals and categorical labels for heatmap.

    Format of source table

    Option 1.
    Format of source table is in wide format with floating-point values of
    signal intensities or a single, specific type of descriptive statistics
    (usually either mean or median) corresponding to features across rows and
    observations or groups of observations across columns.
    ----------
    observation_or_group group_1 group_2 group_3 group_4
    feature
    feature_1            0.01    0.001   0.001   0.015
    feature_2            0.01    0.001   0.001   0.015
    feature_3            -0.01   0.001   0.001   0.015
    feature_4            -0.01   0.001   0.001   0.015
    feature_5            -0.01   0.001   0.001   0.015
    ----------

    Option 2.
    Format of source table is in wide format with floating-point values of
    signal intensities corresponding to features across columns and distinct
    observations across rows. A special column gives identifiers corresponding
    to each observation across rows. Another special column provides names
    of categorical groups of observations.
    ----------
    observation     group   feature_1 feature_2 feature_3 feature_4 feature_5
    observation_1   group_1 0.001     0.001     0.001     0.001     0.001
    observation_2   group_1 0.001     0.001     0.001     0.001     0.001
    observation_3   group_2 0.001     0.001     0.001     0.001     0.001
    observation_4   group_2 0.001     0.001     0.001     0.001     0.001
    observation_5   group_3 0.001     0.001     0.001     0.001     0.001
    ----------

    For versatility, the source table does not have explicitly defined indices
    across columns or rows.

    This function preserves the original sequence of features. This function
    also preserves the original sequence of groups and observations within
    groups.

    Review: 27 December 2024

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            for features across sample observations or groups of sample
            observations
        format_table (int): value 1 for features across rows and observations
            or groups of observations across columns, value 2 for features
            across columns and observations across rows with potential special
            column for groups of observations
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        column_group (str): name of column in table to use for groups
        transpose_table (bool): whether to transpose matrix from table
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_missing_fill (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        value_minimum (float): minimal value for threshold constraint on
            signals and for anchor on scale of projection to visual
            representation in color
        value_center (float): central value for anchor on scale of projection
            to visual representation in color
        value_maximum (float): maximal value for threshold constraint on
            signals and for anchor on scale of projection to visual
            representation in color
        labels_ordinate_categories (list<str>): optional, explicit labels for
            ordinate or vertical axis
        labels_abscissa_categories (list<str>): optional, explicit labels for
            abscissa or horizontal axis
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information for figure

    """

    # Copy information in table.
    table_signal = table.copy(deep=True)
    # Organize indices in table.
    table_signal.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_signal.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    if (format_table == 1):
        table_signal.set_index(
            index_rows,
            append=False,
            drop=True,
            inplace=True
        )
    elif ((format_table == 2) and (column_group is None)):
        table_signal.set_index(
            index_rows,
            append=False,
            drop=True,
            inplace=True
        )
    elif ((format_table == 2) and (column_group is not None)):
        table_signal.set_index(
            [index_rows, column_group,],
            append=False,
            drop=True,
            inplace=True
        )
        pass
    # Extract categorical names from original indices in table.
    #labels_columns = copy.deepcopy(table_signal.columns.to_list())
    #labels_rows = table.index.to_list()
    #labels_rows = table_signal["group_primary"].to_list()
    labels_index_columns = copy.deepcopy(
        table_signal.columns.get_level_values(
            index_columns
        ).to_list()
    )
    labels_index_rows = copy.deepcopy(
        table_signal.index.get_level_values(
            index_rows
        ).to_list()
    )
    if ((format_table == 2) and (column_group is not None)):
        labels_group_rows = copy.deepcopy(
            table_signal.index.get_level_values(
                column_group
            ).to_list()
        )
        labels_group_unique = putly.collect_unique_elements(
            elements=labels_group_rows,
        )
    else:
        labels_group_rows = None
        labels_group_unique = None
        pass
    # Extract counts of individual observations in each group.
    if ((format_table == 2) and (column_group is not None)):
        groups_counts = table_signal.groupby(
            level=column_group,
        ).size().to_dict()
        pass
    else:
        groups_counts = None
        pass

    # Transpose table.
    if (transpose_table):
        table_signal = table_signal.transpose(copy=True)
        pass
    # Extract values.
    #matrix_signal = numpy.transpose(numpy.copy(table_signal.to_numpy()))
    matrix_signal = numpy.copy(table_signal.to_numpy())
    # Extract minimal and maximal values.
    if (
        (not constrain_signal_values) and
        (value_minimum is None) and
        (value_maximum is None)
    ):
        round_offset = (numpy.nanmin(matrix_signal) * 0.10)
        value_minimum = round((numpy.nanmin(matrix_signal) - round_offset), 3)
        value_maximum = round((numpy.nanmax(matrix_signal) + round_offset), 3)
    # Fill missing values.
    if fill_missing:
        matrix_signal = numpy.nan_to_num(
            matrix_signal,
            copy=True,
            nan=value_missing_fill,
            posinf=value_maximum, # or + 1.0 for correlations
            neginf=value_minimum, # or - 1.0 for correlations
        )
        pass
    # Constrain values.
    if constrain_signal_values:
        matrix_signal[matrix_signal < value_minimum] = value_minimum
        matrix_signal[matrix_signal > value_maximum] = value_maximum
        pass
    # Determine labels for axes.
    if (transpose_table):
        if (
            (format_table == 2) and
            (column_group is not None) and
            (
                (labels_abscissa_categories is None) or
                (len(labels_abscissa_categories) < 2)
            )
        ):
            labels_abscissa_categories = labels_index_rows # vertical_axis
            #labels_abscissa_categories = labels_group_rows # vertical_axis
            labels_abscissa_groups = labels_group_unique # vertical_axis

        elif (
            (
                (labels_abscissa_categories is None) or
                (len(labels_abscissa_categories) < 2)
            )
        ):
            labels_abscissa_categories = labels_index_rows # vertical axis
            labels_abscissa_groups = None
        if (
            (labels_ordinate_categories is None) or
            (len(labels_ordinate_categories) < 2)
        ):
            labels_ordinate_categories = labels_index_columns # vertical axis
            labels_ordinate_groups = None
    else:
        if (
            (labels_abscissa_categories is None) or
            (len(labels_abscissa_categories) < 2)
        ):
            labels_abscissa_categories = labels_index_columns # horizontal axis
            labels_abscissa_groups = None,
        if (
            (format_table == 2) and
            (column_group is not None) and
            (
                (labels_ordinate_categories is None) or
                (len(labels_ordinate_categories) < 2)
            )
        ):
            labels_ordinate_categories = labels_index_rows # vertical_axis
            #labels_ordinate_categories = labels_group_rows # vertical_axis
            labels_ordinate_groups = labels_group_unique # vertical_axis
        elif (
            (
                (labels_ordinate_categories is None) or
                (len(labels_ordinate_categories) < 2)
            )
        ):
            labels_ordinate_categories = labels_index_rows # vertical axis
            labels_ordinate_groups = None
        pass

    # Define discrete numerical representation of categorical groups.
    if ((format_table == 2) and (column_group is not None)):
        groups_indices = dict()
        indices_groups = dict()
        index = 0
        for name in labels_group_unique:
            groups_indices[name] = index
            indices_groups[index] = name
            index += 1
            pass
        integers_group_rows = list(map(
            lambda name: groups_indices[name], labels_group_rows
        ))
        #groups_representation = dict()
        #groups_representation["names"] = labels_group_rows
        #groups_representation["integers"] = integers_group_rows
        #table_groups_representation = pandas.DataFrame(
        #    data=groups_representation,
        #)
        # Organize the integer representations of discrete categorical groups
        # as a matrix.
        matrix_group_integers = numpy.array(integers_group_rows).reshape(
            1, len(integers_group_rows)
        )
    else:
        groups_indices = None
        indices_groups = None
        integers_group_rows = None
        matrix_group_integers = None
        pass

    ##########
    # 9. Collect information.
    pail = dict()
    pail["matrix_signal"] = matrix_signal
    pail["matrix_group_integers"] = matrix_group_integers
    pail["value_minimum"] = value_minimum
    pail["value_center"] = value_center
    pail["value_maximum"] = value_maximum
    pail["labels_index_columns"] = labels_index_columns
    pail["labels_index_rows"] = labels_index_rows
    pail["labels_group_rows"] = labels_group_rows
    pail["integers_group_rows"] = integers_group_rows
    pail["labels_group_unique"] = labels_group_unique
    pail["groups_indices"] = groups_indices
    pail["indices_groups"] = indices_groups
    pail["groups_counts"] = groups_counts
    pail["labels_ordinate_categories"] = labels_ordinate_categories
    pail["labels_ordinate_groups"] = labels_ordinate_groups
    pail["labels_abscissa_categories"] = labels_abscissa_categories
    pail["labels_abscissa_groups"] = labels_abscissa_groups

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.plot.py")
        function = str(
            "extract_prepare_table_signals_categories_for_heatmap" +
            "()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=4)
        print("source table of signals:")
        print("rows in table: " + str(table_signal.shape[0]))
        print("columns in table: " + str(table_signal.shape[1]))
        putly.print_terminal_partition(level=4)
        print("matrix of signals:")
        count_rows = copy.deepcopy(matrix_signal.shape[0])
        count_columns = copy.deepcopy(matrix_signal.shape[1])
        print("matrix rows (dimension 0): " + str(count_rows))
        print("matrix columns (dimension 1): " + str(count_columns))
        putly.print_terminal_partition(level=4)
        print("abscissa, horizontal axis labels:")
        #print(labels_abscissa_categories)
        print("ordinate, vertical axis labels:")
        #print(labels_ordinate_categories)
        putly.print_terminal_partition(level=4)
        print("labels_group_rows:")
        print(labels_group_rows)
        print("integers_group_rows:")
        print(integers_group_rows)
        putly.print_terminal_partition(level=4)

    # Return information.
    return pail


def determine_size_axis_labels_categories(
    count_labels=None,
):
    """
    Determine appropriate size of font for text labels on an axis to represent
    explicit names of features or categories.

    The maximal count of labels is 200.

    These sizes of font correspond to definitions in the function below.
    package: partner
    module: plot.py
    function: define_font_properties()

    Review: 24 October 2024

    arguments:
        count_labels (int): count of labels for axis

    raises:

    returns:
        (str): text designation of size for font in labels

    """

    # Determine appropriate size of font for text labels of explicit names of
    # features or categories along axis.
    if (150 <= count_labels and count_labels < 201):
        size_label = "twenty-one"
    elif (120 <= count_labels and count_labels < 150):
        size_label = "twenty"
    elif (100 <= count_labels and count_labels < 120):
        size_label = "nineteen"
    elif (80 <= count_labels and count_labels < 100):
        size_label = "eighteen"
    elif (60 <= count_labels and count_labels < 80):
        size_label = "seventeen"
    elif (40 <= count_labels and count_labels < 60):
        size_label = "sixteen"
    elif (30 <= count_labels and count_labels < 40):
        size_label = "fifteen"
    elif (20 <= count_labels and count_labels < 30):
        size_label = "fourteen"
    elif (10 <= count_labels and count_labels < 20):
        size_label = "thirteen"
    elif (1 <= count_labels and count_labels < 10):
        size_label = "eleven"
    else:
        size_label = None
        pass
    # Return information.
    return size_label


def plot_heatmap_signal_features_observations_labels(
    table=None,
    format_table=None,
    index_columns=None,
    index_rows=None,
    transpose_table=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_center=None,
    value_maximum=None,
    title_ordinate=None,
    title_abscissa=None,
    title_bar=None,
    labels_ordinate_categories=None,
    labels_abscissa_categories=None,
    size_title_ordinate=None,
    size_title_abscissa=None,
    size_title_bar=None,
    size_label_ordinate=None,
    size_label_abscissa=None,
    size_label_bar=None,
    show_labels_ordinate=None,
    show_labels_abscissa=None,
    show_scale_bar=None,
    aspect=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Heatmap.

    features of this chart design...
    labels of categorical groups on both axes: True
    labels of significance on individual cells: False
    clustering: False

    Format of source table

    Format of source table is in wide format with floating-point values of
    signal intensities or a single, specific type of descriptive statistics
    (usually either mean or median) corresponding to features across rows and
    observations or groups of observations across columns.
    ----------
    observation_or_group group_1 group_2 group_3 group_4
    feature
    feature_1            0.01    0.001   0.001   0.015
    feature_2            0.01    0.001   0.001   0.015
    feature_3            -0.01   0.001   0.001   0.015
    feature_4            -0.01   0.001   0.001   0.015
    feature_5            -0.01   0.001   0.001   0.015
    ----------

    For versatility, the source table does not have explicitly defined indices
    across columns or rows.

    This function preserves the original sequence of features. This function
    also preserves the original sequence of groups and observations within
    groups.

    MatPlotLib color maps.
    https://matplotlib.org/stable/tutorials/colors/colormaps.html


    Review or revision: 13 November 2025
    Review or revision: 30 December 2024

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            for features in rows across sample observations or groups of
            sample observations in columns
        format_table (int): value 1 for features across rows and observations
            or groups of observations across columns, value 2 for features
            across columns and observations across rows with potential special
            column for groups of observations
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        transpose_table (bool): whether to transpose matrix from table
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_missing_fill (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        value_minimum (float): minimal value for threshold constraint on
            signals and for anchor on scale of projection to visual
            representation in color
        value_center (float): central value for anchor on scale of projection
            to visual representation in color
        value_maximum (float): maximal value for threshold constraint on
            signals and for anchor on scale of projection to visual
            representation in color
        title_ordinate (str): title for ordinate vertical axis
        title_abscissa (str): title for abscissa horizontal axis
        title_bar (str): title for scale bar
        labels_ordinate_categories (list<str>): optional, explicit labels for
            ordinate or vertical axis
        labels_abscissa_categories (list<str>): optional, explicit labels for
            abscissa or horizontal axis
        size_title_ordinate (str): font size
        size_title_abscissa (str): font size
        size_title_bar (str): font size
        size_label_ordinate (str): font size
        size_label_abscissa (str): font size
        size_label_bar (str): font size
        show_labels_ordinate (bool): whether to show on vertical ordinate axis
            of plot chart explicit text labels for individual categories
        show_labels_abscissa (bool): whether to show on horizontal abscissa
            axis of plot chart explicit text labels for individual categories
        show_scale_bar (bool): whether to create scale bar
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): MatPlotLib figure object

    """

    ##########
    # Prepare information for figure.
    pail = extract_prepare_table_signals_categories_for_heatmap(
        table=table,
        format_table=format_table, # 1: features in rows, observations in columns
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=None,
        transpose_table=transpose_table,
        fill_missing=fill_missing,
        value_missing_fill=value_missing_fill,
        constrain_signal_values=constrain_signal_values,
        value_minimum=value_minimum,
        value_center=value_center,
        value_maximum=value_maximum,
        labels_ordinate_categories=labels_ordinate_categories,
        labels_abscissa_categories=labels_abscissa_categories,
        report=report,
    )

    ##########
    # Create figure.
    figure = pplot.initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Create axes.
    #axes_main = matplotlib.pyplot.axes()
    axes_main = figure.add_subplot(111)
    # Keep axes, ticks, and labels, but remove border.
    for position in ['right', 'top', 'bottom', 'left']:
        matplotlib.pyplot.gca().spines[position].set_visible(False)
    # Adjust margins.
    figure.subplots_adjust(
        left=0.02,
        right=0.99,
        top=0.98,
        bottom=0.15,
    )

    # Plot values as a grid of color on continuous scale.
    # This function represents values acros matrix dimension 0 as vertical
    # rows.
    # This function represents values across matrix dimension 1 as horizontal
    # columns.
    # Diverging color maps: "PRGn", "PRGn_r", "PiYG", "PiYG_r",
    # Diverging color maps: "PuOr", "PuOr_r",
    # Diverging color maps: "PuOr", "PuOr_r", "RdBu", "RdBu_r", "BrBG",
    # Sequential color maps: "Reds", "Reds_r", "Oranges", "Oranges_r",
    # site: https://montoliu.naukas.com/2021/11/18/color-blindness-purple-and-
    #     orange-are-the-solution/
    #image_main = axes_main.imshow(
    #    pail["matrix_signal"],
    #    cmap=matplotlib.colormaps["PuOr"], # binary, Reds, RdBu_r, PuOr, PuOr_r
    #    vmin=pail["value_minimum"],
    #    vmax=pail["value_maximum"],
    #    aspect="auto", # "auto", "equal",
    #    origin="lower",
    #    # Extent: (left, right, bottom, top)
    #    #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    #)
    pail_color_map = create_divergent_color_map(
        value_minimum=pail["value_minimum"],
        value_center=pail["value_center"],
        value_maximum=pail["value_maximum"],
        color_minimum=(0.05,0.05,0.6,1.0,), # "navy blue";(red: 13; green: 13; blue: 153) # TCW; 13 November 2025
        color_center=(1.0,1.0,1.0,1.0,), # "white"
        color_maximum=(0.9,0.5,0.05,1.0,), # "orange"; (red: 230, green: 128, blue: 13) # TCW; 13 November 2025
    )
    image_main = axes_main.imshow(
        pail["matrix_signal"],
        cmap=pail_color_map["map"],
        norm=pail_color_map["scale"],
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set titles for axes.
    if (len(title_ordinate) > 0):
        axes_main.set_ylabel(
            ylabel=title_ordinate,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_ordinate]
        )
        pass
    if (len(title_abscissa) > 0):
        axes_main.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
        pass
    # Set tick parameters for axes.
    axes_main.tick_params(
        axis="both", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=2.5, # 5.0
        width=1.5, # 3.5
        pad=5, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=False,
        bottom=False,
        left=False,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=False,
    )
    # Manage labels for categories along the axes.
    # Determine feasibility and appropriate font size for representing labels
    # of individual categorical features or observations along the axes.
    count_labels_ordinate = len(pail["labels_ordinate_categories"])
    count_labels_abscissa = len(pail["labels_abscissa_categories"])
    if (size_label_ordinate is None):
        size_label_ordinate = determine_size_axis_labels_categories(
            count_labels=count_labels_ordinate,
        )
        pass
    if (size_label_abscissa is None):
        size_label_abscissa = determine_size_axis_labels_categories(
            count_labels=count_labels_abscissa,
        )
        pass
    # Determine whether to show labels for features or categories along the
    # vertical ordinate axis.
    if (
        (show_labels_ordinate) and
        (size_label_ordinate is not None) and
        (pail["labels_ordinate_categories"] is not None) and
        (count_labels_ordinate > 1)
    ):
        # Set tick positions and labels on vertical ordinate axis.
        axes_main.set_yticks(
            numpy.arange(pail["matrix_signal"].shape[0]),
        )
        axes_main.set_yticklabels(
            pail["labels_ordinate_categories"],
            #minor=False,
            ha="right", # horizontal alignment
            va="center", # vertical alignment
            alpha=1.0,
            rotation=0.0,
            rotation_mode="anchor",
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_label_ordinate]
        )
        # Set tick parameters for vertical ordinate axis.
        axes_main.tick_params(
            axis="y", # "y", "x", or "both"
            which="both", # "major", "minor", or "both"
            length=2.5, # 5.0
            width=1.5, # 3.5
            pad=5, # 7.5
            direction="out",
            color=colors["black"],
            labelcolor=colors["black"],
            left=True,
            right=False,
            labelleft=True,
            labelright=False,
        )
        pass
    # Determine whether to show labels for features or categories along the
    # horizontal abscissa axis.
    if (
        (show_labels_abscissa) and
        (size_label_abscissa is not None) and
        (pail["labels_abscissa_categories"] is not None) and
        (count_labels_abscissa > 1)
    ):
        # Set tick positions and labels on horizontal abscissa axis.
        axes_main.set_xticks(
            numpy.arange(pail["matrix_signal"].shape[1]),
        )
        axes_main.set_xticklabels(
            pail["labels_abscissa_categories"],
            #minor=False,
            ha="left", # horizontal alignment
            va="top", # vertical alignment
            alpha=1.0,
            rotation=-60,
            rotation_mode="anchor",
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_label_abscissa]
        )
        # Set tick parameters for horizontal abscissa axis.
        axes_main.tick_params(
            axis="x", # "y", "x", or "both"
            which="both", # "major", "minor", or "both"
            length=2.5, # 5.0
            width=1.5, # 3.5
            pad=7, # 7.5 - 17.5
            direction="out",
            color=colors["black"],
            labelcolor=colors["black"],
            top=False,
            bottom=True,
            labeltop=False,
            labelbottom=True,
        )
        pass

    # Create legend for scale of color grid.
    if show_scale_bar:
        bar = axes_main.figure.colorbar(
            image_main,
            orientation="vertical",
            ax=axes_main,
            location="right",
            shrink=0.5, # 0.7; factor for dimensions of the Scale Bar.
        )
        if (len(title_bar) > 0):
            bar.ax.set_ylabel(
                title_bar,
                rotation=-90,
                va="bottom",
                labelpad=5, # 5
                alpha=1.0,
                backgroundcolor=colors["white"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_title_bar],
            )
        bar.ax.tick_params(
            axis="both",
            which="both", # major, minor, or both
            direction="out",
            length=5, # 5.0, 7.5
            width=2.5, # 2.5, 5.0
            color=colors["black"],
            pad=5, # 5, 7
            labelsize=fonts["values"][size_label_bar]["size"],
            labelcolor=colors["black"],
        )

    # Return figure.
    return figure


def plot_heatmap_signal_features_sets_observations_labels(
    table_signal=None,
    table_feature_sets=None,
    format_table_signal=None,
    index_columns=None,
    index_rows=None,
    transpose_table=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_maximum=None,
    title_ordinate=None,
    title_abscissa=None,
    title_bar=None,
    labels_abscissa_categories=None,
    size_title_ordinate=None,
    size_title_abscissa=None,
    size_title_bar=None,
    size_label_feature_set=None,
    size_label_abscissa=None,
    size_label_bar=None,
    show_labels_abscissa=None,
    show_scale_bar=None,
    aspect=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Heat map.

    Notice that this chart design does not show explicit categorical labels for
    features on the vertical ordinate axis or for observations or groups of
    observations on the horizontal abscissa axis.

    Format of source table of signals
    Format of source table is in wide format with floating-point values of
    signal intensities or a single, specific type of descriptive statistics
    (usually either mean or median) corresponding to features across rows and
    observations or groups of observations across columns.
    ----------
    observation_or_group group_1 group_2 group_3 group_4
    feature
    feature_1            0.01    0.001   0.001   0.015
    feature_2            0.01    0.001   0.001   0.015
    feature_3            -0.01   0.001   0.001   0.015
    feature_4            -0.01   0.001   0.001   0.015
    feature_5            -0.01   0.001   0.001   0.015
    ----------

    Format of source table of feature allocations to sets
    ----------
    feature     set_1 set_2 set_3 set_4 set_5
    feature_1   1     0     0     0     1
    feature_2   1     1     0     0     0
    feature_3   0     1     1     0     0
    feature_4   0     0     1     1     0
    feature_5   0     0     0     1     1
    ----------

    For versatility, the source tables do not have explicitly defined indices
    across columns or rows.

    This function preserves the original sequence of features. This function
    also preserves the original sequence of groups and observations within
    groups.

    This function assumes that the table has many more observations than
    features, and for this reason, the design orients features across the
    vertical axis and observations across the horizontal axis. With a landscape
    aspect ratio, this design allows more space for the horizontal axis.

    ----------

    Reference:

    Review: 30 December 2024

    arguments:
        table_signal (object): Pandas data-frame table of values of signal
            intensity for features in columns across sample observations or
            groups of sample observations in rows
        table_feature_sets (object): Pandas data-frame table of indications of
            allocation of features to sets in a sort sequence that matches the
            sequence of features across columns in table of signals
        format_table_signal (int): value 1 for features across rows and
            observations or groups of observations across columns, value 2 for
            features across columns and observations across rows with potential
            special column for groups of observations
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        transpose_table (bool): whether to transpose matrix from table
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_missing_fill (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
        title_ordinate (str): title for ordinate vertical axis
        title_abscissa (str): title for abscissa horizontal axis
        title_bar (str): title for scale bar
        labels_abscissa_categories (list<str>): optional, explicit labels for
            abscissa or horizontal axis
        size_title_ordinate (str): font size
        size_title_abscissa (str): font size
        size_title_bar (str): font size
        size_label_feature_set (str): font size
        size_label_abscissa (str): font size
        size_label_bar (str): font size
        show_labels_abscissa (bool): whether to show on horizontal abscissa
            axis of plot chart explicit text labels for individual categories
        show_scale_bar (bool): whether to create scale bar
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): MatPlotLib figure object

    """

    ##########
    # Prepare information for figure.
    pail = extract_prepare_table_signals_categories_for_heatmap(
        table=table_signal,
        format_table=format_table_signal, # 1: features in rows, observations in columns
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=None,
        transpose_table=transpose_table,
        fill_missing=fill_missing,
        value_missing_fill=value_missing_fill,
        constrain_signal_values=constrain_signal_values,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        labels_ordinate_categories=None,
        labels_abscissa_categories=labels_abscissa_categories,
        report=report,
    )
    # Copy information in table.
    table_feature_sets = table_feature_sets.copy(deep=True)
    # Organize indices in table.
    table_feature_sets.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_feature_sets.set_index(
        index_rows,
        append=False,
        drop=True,
        inplace=True
    )
    table_feature_sets.columns.rename(
        index_columns,
        inplace=True,
    ) # single-dimensional index
    # Extract categorical names from original indices in table.
    labels_feature_sets_index_columns = copy.deepcopy(
        table_feature_sets.columns.get_level_values(
            index_columns
        ).to_list()
    )
    labels_feature_sets_index_rows = copy.deepcopy(
        table_feature_sets.index.get_level_values(
            index_rows
        ).to_list()
    )
    # Extract values.
    #matrix_signal = numpy.transpose(numpy.copy(table_signal.to_numpy()))
    matrix_feature_sets = numpy.copy(table_feature_sets.to_numpy())

    ##########
    # Create figure.

    ##########
    # Initialize figure.
    figure = pplot.initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Initialize grid within figure.
    #matplotlib.gridspec.GridSpec()
    #figure.add_gridspec()
    #sharey=True, # sharey="row"
    #sharex=True, # sharex="col",
    grid = matplotlib.gridspec.GridSpec(
        nrows=1,
        ncols=4,
        wspace=0.005, # horizontal width space between grid blocks for subplots
        hspace=0.005, # vertical height space between grid blocks for subplots
        width_ratios=(40,55,1,4),
        height_ratios=(100,),
    )
    grid.update(
        wspace=0.005, # horizontal width space between grid blocks for subplots
        hspace=0.005, # vertical height space between grid blocks for subplots
    )
    # Initialize axes within grid within figure.
    axes_set = figure.add_subplot(grid[0,0]) # first row, first column
    axes_main = figure.add_subplot(grid[0,1]) # first row, second column
    axes_space = figure.add_subplot(grid[0,2]) # first row, third column
    axes_bar = figure.add_subplot(grid[0,3]) # first row, fourt column
    axes_set.clear()
    axes_main.clear()
    axes_space.clear()
    axes_bar.clear()
    # Set axes to empty as a space holder.
    axes_space.axis("off")
    # Keep axes, ticks, and labels, but remove border.
    for position in ['right', 'top', 'bottom', 'left']:
        matplotlib.pyplot.gca().spines[position].set_visible(False)
    # Adjust margins.
    # Method "tight_layout()" does not function properly with object
    # "GridSpec()".
    #grid.tight_layout(
    #    figure,
    #    #pad=1.0,
    #    #h_pad=1.0,
    #    #w_pad=1.0,
    #    rect=[0,0.05,1.0,1.0], # left, bottom, right, top
    #)
    #grid.update(
    #    wspace=0.005, # horizontal width space between grid blocks for subplots
    #    hspace=0.005, # vertical height space between grid blocks for subplots
    #)
    figure.subplots_adjust(
        left=0.01,
        right=0.90,
        top=0.99,
        bottom=0.20,
    )

    ##########
    # axes: main
    # Plot values as a grid of color on continuous scale.
    # This function represents values acros matrix dimension 0 as vertical
    # rows.
    # This function represents values across matrix dimension 1 as horizontal
    # columns.
    # Diverging color maps: "PRGn", "PRGn_r", "PiYG", "PiYG_r",
    # Diverging color maps: "PuOr", "PuOr_r",
    # Diverging color maps: "PuOr", "PuOr_r", "RdBu", "RdBu_r", "BrBG",
    # Sequential color maps: "Reds", "Reds_r", "Oranges", "Oranges_r",
    # site: https://montoliu.naukas.com/2021/11/18/color-blindness-purple-and-
    #     orange-are-the-solution/
    image_main = axes_main.imshow(
        pail["matrix_signal"],
        cmap=matplotlib.colormaps["PuOr"], # binary, Reds, RdBu_r, PuOr, PuOr_r
        vmin=pail["value_minimum"],
        vmax=pail["value_maximum"],
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set titles for axes.
    if (len(title_ordinate) > 0):
        axes_main.set_ylabel(
            ylabel=title_ordinate,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_ordinate]
        )
        pass
    if (len(title_abscissa) > 0):
        axes_main.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
        pass
    # Set tick parameters for axes.
    axes_main.tick_params(
        axis="both", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=5.0, # 5.0
        width=3.5, # 3.5
        pad=10, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=False,
        bottom=False,
        left=False,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=False,
    )
    # Manage labels for categories along the axes.
    # Determine feasibility and appropriate font size for representing labels
    # of individual categorical features or observations along the axes.
    count_labels_abscissa = len(pail["labels_abscissa_categories"])
    if (size_label_abscissa is None):
        size_label_abscissa = determine_size_axis_labels_categories(
            count_labels=count_labels_abscissa,
        )
        pass
    # Determine whether to show labels for features or categories along the
    # horizontal abscissa axis.
    if (
        (show_labels_abscissa) and
        (size_label_abscissa is not None) and
        (pail["labels_abscissa_categories"] is not None) and
        (count_labels_abscissa > 1)
    ):
        # Set tick positions and labels on horizontal abscissa axis.
        axes_main.set_xticks(
            numpy.arange(pail["matrix_signal"].shape[1]),
        )
        axes_main.set_xticklabels(
            pail["labels_abscissa_categories"],
            #minor=False,
            ha="left", # horizontal alignment
            va="top", # vertical alignment
            alpha=1.0,
            rotation=-60,
            rotation_mode="anchor",
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_label_abscissa]
        )
        # Set tick parameters for horizontal abscissa axis.
        axes_main.tick_params(
            axis="x", # "y", "x", or "both"
            which="both", # "major", "minor", or "both"
            length=5.0, # 5.0
            width=3.5, # 3.5
            pad=10, # 7.5
            direction="out",
            color=colors["black"],
            labelcolor=colors["black"],
            top=False,
            bottom=True,
            labeltop=False,
            labelbottom=True,
        )
        pass

    ##########
    # axes: set
    # Define color map for discrete, binary integer representation of
    # allocation to sets.
    # https://matplotlib.org/3.1.1/gallery/color/named_colors.html
    # "white", "black"
    # "white", "dimgray"
    color_map_set = matplotlib.colors.ListedColormap([
        "white", "dimgray"
    ])
    # Plot values as a grid of color on discrete scale.
    image = axes_set.imshow(
        matrix_feature_sets,
        cmap=color_map_set,
        vmin=0,
        vmax=1,
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set tick positions and labels on axes.
    axes_set.set_xticks(
        numpy.arange(matrix_feature_sets.shape[1]),
    )
    axes_set.set_xticklabels(
        labels_feature_sets_index_columns,
        #minor=False,
        ha="center", # horizontal alignment
        va="top", # vertical alignment
        alpha=1.0,
        rotation=90, # negative: clockwise; positive: count-clockwise
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"][size_label_feature_set],
    )
    # Set tick parameters for axes.
    axes_set.tick_params(
        axis="y", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        left=False,
        right=False,
        labelleft=False,
        labelright=False,
    )
    axes_set.tick_params(
        axis="x", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=5.0, # 5.0
        width=3.5, # 3.5
        pad=10, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=False,
        bottom=True,
        labeltop=False,
        labelbottom=True,
    )

    ##########
    # axes: bar
    # Create legend for scale of color grid.
    # Notice that use of the "cax" argument causes to ignore the "shrink"
    # argument.
    if show_scale_bar:
        # Create scale bar.
        bar_main = axes_bar.figure.colorbar(
            image_main,
            orientation="vertical",
            cax=axes_bar,
            location="right",
            #shrink=0.5, # 0.7; factor for dimensions of the Scale Bar.
        )
        if (len(title_bar) > 0):
            bar_main.ax.set_ylabel(
                title_bar,
                rotation=-90,
                loc="center",
                va="center",
                labelpad=15, # 5
                alpha=1.0,
                backgroundcolor=colors["white"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_title_bar],
            )
        bar_main.ax.tick_params(
            axis="both",
            which="both", # major, minor, or both
            direction="out",
            length=7.5, # 5.0, 7.5
            width=3, # 2.5, 5.0
            color=colors["black"],
            pad=3, # 5, 7
            labelsize=fonts["values"][size_label_bar]["size"],
            labelcolor=colors["black"],
        )
        pass

    ##########
    # Return figure.
    return figure


def plot_heatmap_signal_features_labels_observations_groups(
    table_signal=None,
    format_table_signal=None,
    index_columns=None,
    index_rows=None,
    column_group=None,
    transpose_table=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_maximum=None,
    title_ordinate=None,
    title_abscissa=None,
    title_bar=None,
    labels_ordinate_categories=None,
    size_title_ordinate=None,
    size_title_abscissa=None,
    size_title_bar=None,
    size_label_ordinate=None,
    size_label_legend_observation_group=None,
    size_label_bar=None,
    show_labels_ordinate=None,
    show_scale_bar=None,
    aspect=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Heat map.

    Notice that this chart design shows explicit categorical labels for
    features on the vertical ordinate axis but not for observations or groups
    of observations on the horizontal abscissa axis.

    Format of source table of signals
    Format of source table is in wide format with floating-point values of
    signal intensities corresponding to features across columns and distinct
    observations across rows. A special column gives identifiers corresponding
    to each observation across rows. Another special column provides names
    of categorical groups of observations.
    ----------
    observation     group   feature_1 feature_2 feature_3 feature_4 feature_5
    observation_1   group_1 0.001     0.001     0.001     0.001     0.001
    observation_2   group_1 0.001     0.001     0.001     0.001     0.001
    observation_3   group_2 0.001     0.001     0.001     0.001     0.001
    observation_4   group_2 0.001     0.001     0.001     0.001     0.001
    observation_5   group_3 0.001     0.001     0.001     0.001     0.001
    ----------

    For versatility, the source tables do not have explicitly defined indices
    across columns or rows.

    This function preserves the original sequence of features. This function
    also preserves the original sequence of groups and observations within
    groups.

    This function assumes that the table has many more observations than
    features, and for this reason, the design orients features across the
    vertical axis and observations across the horizontal axis. With a landscape
    aspect ratio, this design allows more space for the horizontal axis.

    ----------

    Reference:
    https://www.kaggle.com/code/sgalella/correlation-heatmaps-with-
        hierarchical-clustering
    https://www.kaggle.com/code/chinoysen/beginners-guide-to-heatmaps-cluster-
        heatmap
    https://seaborn.pydata.org/generated/seaborn.clustermap.html
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.heatmap.html
    https://gist.github.com/peterk87/5505691
    https://matplotlib.org/stable/tutorials/colors/colormaps.html

    Maybe useful for problem-solving:
    https://matplotlib.org/stable/gallery/subplots_axes_and_figures/
        subplots_demo.html
        - subplots and gridspec
    http://www.acgeospatial.co.uk/colour-bar-for-discrete-rasters-with-matplotlib/
        - color maps for discrete values
    https://matplotlib.org/stable/users/explain/colors/colorbar_only.html
        #colorbar-with-arbitrary-colors
        - color bar with custom dimensions of discrete, categorical intervals

    https://matplotlib.org/stable/gallery/color/colormap_reference.html
        - color maps in MatPlotLib

    https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
        - categorical color bars
        - helpful reference for figuring out how to set the thresholds between
          the discrete levels of the color bar to correspond properly to the
          columns of the main heatmap...
    https://stackoverflow.com/questions/9707676/defining-a-discrete-colormap-for-imshow
    https://stackoverflow.com/questions/7229971/2d-grid-data-visualization-in-python
    https://matplotlib.org/stable/users/explain/colors/colormapnorms.html
        - normalization of information for color maps???

    Review: 3 October 2025
    Review: 30 December 2024

    arguments:
        table_signal (object): Pandas data-frame table of values of signal
            intensity for features in columns across sample observations or
            groups of sample observations in rows
        format_table_signal (int): value 1 for features across rows and
            observations or groups of observations across columns, value 2 for
            features across columns and observations across rows with potential
            special column for groups of observations
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        column_group (str): name of column in table to use for groups
        transpose_table (bool): whether to transpose matrix from table
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_missing_fill (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
        title_ordinate (str): title for ordinate vertical axis
        title_abscissa (str): title for abscissa horizontal axis
        title_bar (str): title for scale bar
        labels_ordinate_categories (list<str>): optional, explicit labels for
            ordinate or vertical axis
        size_title_ordinate (str): font size
        size_title_abscissa (str): font size
        size_title_bar (str): font size
        size_label_ordinate (str): font size
        size_label_legend_observation_group (str): font size
        size_label_bar (str): font size
        show_labels_ordinate (bool): whether to show on vertical ordinate axis
            of plot chart explicit text labels for individual categories
        show_scale_bar (bool): whether to create scale bar
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): MatPlotLib figure object

    """

    ##########
    # Prepare information for figure.
    pail = extract_prepare_table_signals_categories_for_heatmap(
        table=table_signal,
        format_table=format_table_signal, # 1: features in rows, observations in columns
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=column_group,
        transpose_table=transpose_table,
        fill_missing=fill_missing,
        value_missing_fill=value_missing_fill,
        constrain_signal_values=constrain_signal_values,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        report=report,
    )

    ##########
    # Create figure.

    ##########
    # Initialize figure.
    figure = pplot.initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Initialize grid within figure.
    #matplotlib.gridspec.GridSpec()
    #figure.add_gridspec()
    #sharey=True, # sharey="row"
    #sharex=True, # sharex="col",
    grid = matplotlib.gridspec.GridSpec(
        nrows=4,
        ncols=2,
        wspace=0.005, # horizontal width space between grid blocks for subplots
        hspace=0.005, # vertical height space between grid blocks for subplots
        width_ratios=(20,80),
        height_ratios=(90,3,5,2),
    )
    grid.update(
        wspace=0.005, # horizontal width space between grid blocks for subplots
        hspace=0.005, # vertical height space between grid blocks for subplots
    )
    # Initialize axes within grid within figure.
    axes_space_1 = figure.add_subplot(grid[0,0]) # first row, first column
    axes_main = figure.add_subplot(grid[0,1]) # first row, second column
    axes_space_2 = figure.add_subplot(grid[1,0]) # second row, first column
    axes_group = figure.add_subplot(grid[1,1]) # second row, second column
    axes_space_3 = figure.add_subplot(grid[2,0]) # third row, first column
    axes_space_4 = figure.add_subplot(grid[2,1]) # third row, second column
    axes_space_5 = figure.add_subplot(grid[3,0]) # fourth row, first column
    axes_bar = figure.add_subplot(grid[3,1]) # fourth row, second column
    axes_main.clear()
    axes_group.clear()
    axes_bar.clear()
    axes_space_1.clear()
    axes_space_2.clear()
    axes_space_3.clear()
    axes_space_4.clear()
    axes_space_5.clear()
    # Set axes to empty as a space holder.
    axes_space_1.axis("off")
    axes_space_2.axis("off")
    axes_space_3.axis("off")
    axes_space_4.axis("off")
    axes_space_5.axis("off")
    #axes_set.axis("off")
    #axes_group.axis("off")
    #axes_bar.axis("off")
    #axes_main.axis("off")
    # Keep axes, ticks, and labels, but remove border.
    for position in ['right', 'top', 'bottom', 'left']:
        matplotlib.pyplot.gca().spines[position].set_visible(False)
    # Adjust margins.
    # Method "tight_layout()" does not function properly with object
    # "GridSpec()".
    #grid.tight_layout(
    #    figure,
    #    #pad=1.0,
    #    #h_pad=1.0,
    #    #w_pad=1.0,
    #    rect=[0,0.05,1.0,1.0], # left, bottom, right, top
    #)
    #grid.update(
    #    wspace=0.005, # horizontal width space between grid blocks for subplots
    #    hspace=0.005, # vertical height space between grid blocks for subplots
    #)
    figure.subplots_adjust(
        left=0.02,
        right=0.98,
        top=0.98,
        bottom=0.05,
    )

    ##########
    # axes: main
    # Plot values as a grid of color on continuous scale.
    # This function represents values acros matrix dimension 0 as vertical
    # rows.
    # This function represents values across matrix dimension 1 as horizontal
    # columns.
    # Diverging color maps: "PRGn", "PRGn_r", "PiYG", "PiYG_r",
    # Diverging color maps: "PuOr", "PuOr_r",
    # Diverging color maps: "PuOr", "PuOr_r", "RdBu", "RdBu_r", "BrBG",
    # Sequential color maps: "Reds", "Reds_r", "Oranges", "Oranges_r",
    # site: https://montoliu.naukas.com/2021/11/18/color-blindness-purple-and-
    #     orange-are-the-solution/
    image_main = axes_main.imshow(
        pail["matrix_signal"],
        cmap=matplotlib.colormaps["PuOr"], # binary, Reds, RdBu_r, PuOr, PuOr_r
        vmin=pail["value_minimum"],
        vmax=pail["value_maximum"],
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set titles for axes.
    if (len(title_ordinate) > 0):
        axes_main.set_ylabel(
            ylabel=title_ordinate,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_ordinate]
        )
        pass
    if (len(title_abscissa) > 0):
        axes_main.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
        pass
    # Set tick parameters for axes.
    axes_main.tick_params(
        axis="both", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=5.0, # 5.0
        width=3.5, # 3.5
        pad=10, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=False,
        bottom=False,
        left=False,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=False,
    )
    # Manage labels for categories along the axes.
    # Determine feasibility and appropriate font size for representing labels
    # of individual categorical features or observations along the axes.
    count_labels_ordinate = len(pail["labels_ordinate_categories"])
    if (size_label_ordinate is None):
        size_label_ordinate = determine_size_axis_labels_categories(
            count_labels=count_labels_ordinate,
        )
        pass
    # Determine whether to show labels for features or categories along the
    # vertical ordinate axis.
    if (
        (show_labels_ordinate) and
        (size_label_ordinate is not None) and
        (pail["labels_ordinate_categories"] is not None) and
        (count_labels_ordinate > 1)
    ):
        # Set tick positions and labels on vertical ordinate axis.
        axes_main.set_yticks(
            numpy.arange(pail["matrix_signal"].shape[0]),
        )
        axes_main.set_yticklabels(
            pail["labels_ordinate_categories"],
            #minor=False,
            ha="right", # horizontal alignment
            va="center", # vertical alignment
            alpha=1.0,
            rotation=0.0,
            rotation_mode="anchor",
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_label_ordinate]
        )
        # Set tick parameters for vertical ordinate axis.
        axes_main.tick_params(
            axis="y", # "y", "x", or "both"
            which="both", # "major", "minor", or "both"
            length=5.0, # 5.0
            width=3.5, # 3.5
            pad=10, # 7.5
            direction="out",
            color=colors["black"],
            labelcolor=colors["black"],
            left=True,
            right=False,
            labelleft=True,
            labelright=False,
        )
        pass

    ##########
    # axes: group
    # Define color map for discrete, integer representation of categorical
    # groups.
    #discrete_minimum = 0
    #discrete_maximum = len(pail["labels_group_unique"])
    discrete_minimum = numpy.nanmin(pail["matrix_group_integers"])
    discrete_maximum = numpy.nanmax(pail["matrix_group_integers"])
    color_map_group = matplotlib.pyplot.get_cmap(
        "tab10", # "Set1", "Set2", "Dark2", "tab10",
        ((discrete_maximum - discrete_minimum) + 1)
    )
    # Plot values as a grid of color on discrete scale.
    image_group = axes_group.imshow(
        pail["matrix_group_integers"],
        cmap=color_map_group,
        vmin=discrete_minimum,
        vmax=discrete_maximum,
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(
        #    0.0,
        #    (pail["matrix_group_integers"].shape[1]),
        #    0.0,
        #    (3*(pail["matrix_group_integers"].shape[0]))
        #),
    )
    axes_group.tick_params(
        top=False,
        bottom=False,
        left=False,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=False,
    )
    # Create legend.
    # See nested function "create_legend_elements()" within main function
    # "plot_scatter_factor_groups()".
    handles = [matplotlib.patches.Patch(
        color=color_map_group(i),
        label=pail["indices_groups"][i]
    ) for i in pail["indices_groups"].keys()]
    axes_group.legend(
        handles=handles,
        prop=fonts["properties"][size_label_legend_observation_group],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=4,
    )

    ##########
    # axes: bar
    # Create legend for scale of color grid.
    # Notice that use of the "cax" argument causes to ignore the "shrink"
    # argument.
    if show_scale_bar:
        # Set definitive scale to avoid conflicts with global objects.
        #color_map_main = matplotlib.pyplot.cm.ScalarMappable(
        #    cmap="PuOr",
        #    norm=matplotlib.pyplot.Normalize(
        #        vmin=pail["value_minimum"],
        #        vmax=pail["value_maximum"],
        #    ),
        #)
        #color_map_main.set_array(pail["matrix_signal"])
        # Create scale bar.
        bar_main = axes_bar.figure.colorbar(
            image_main, # image_main, color_map_main
            orientation="horizontal",
            cax=axes_bar,
            location="bottom",
            #shrink=0.5, # 0.7; factor for dimensions of the Scale Bar.
        )
        if (len(title_bar) > 0):
            bar_main.ax.set_xlabel(
                title_bar,
                rotation=0,
                loc="center",
                va="bottom",
                labelpad=20, # 5
                alpha=1.0,
                backgroundcolor=colors["white"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_title_bar],
            )
        bar_main.ax.tick_params(
            axis="both",
            which="both", # major, minor, or both
            direction="out",
            length=7.5, # 5.0, 7.5
            width=3, # 2.5, 5.0
            color=colors["black"],
            pad=3, # 5, 7
            labelsize=fonts["values"][size_label_bar]["size"],
            labelcolor=colors["black"],
        )
        pass

    ##########
    # Return figure.
    return figure


def plot_heatmap_signal_features_sets_observations_groups(
    table_signal=None,
    table_feature_sets=None,
    format_table_signal=None,
    index_columns=None,
    index_rows=None,
    column_group=None,
    transpose_table=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_maximum=None,
    title_ordinate=None,
    title_abscissa=None,
    title_bar=None,
    size_title_ordinate=None,
    size_title_abscissa=None,
    size_title_bar=None,
    size_label_feature_set=None,
    size_label_legend_observation_group=None,
    size_label_bar=None,
    show_scale_bar=None,
    aspect=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Heat map.

    Notice that this chart design does not show explicit categorical labels for
    features on the vertical ordinate axis or for observations or groups of
    observations on the horizontal abscissa axis.

    Format of source table of signals
    Format of source table is in wide format with floating-point values of
    signal intensities corresponding to features across columns and distinct
    observations across rows. A special column gives identifiers corresponding
    to each observation across rows. Another special column provides names
    of categorical groups of observations.
    ----------
    observation     group   feature_1 feature_2 feature_3 feature_4 feature_5
    observation_1   group_1 0.001     0.001     0.001     0.001     0.001
    observation_2   group_1 0.001     0.001     0.001     0.001     0.001
    observation_3   group_2 0.001     0.001     0.001     0.001     0.001
    observation_4   group_2 0.001     0.001     0.001     0.001     0.001
    observation_5   group_3 0.001     0.001     0.001     0.001     0.001
    ----------

    Format of source table of feature allocations to sets
    ----------
    feature     set_1 set_2 set_3 set_4 set_5
    feature_1   1     0     0     0     1
    feature_2   1     1     0     0     0
    feature_3   0     1     1     0     0
    feature_4   0     0     1     1     0
    feature_5   0     0     0     1     1
    ----------

    For versatility, the source tables do not have explicitly defined indices
    across columns or rows.

    This function preserves the original sequence of features. This function
    also preserves the original sequence of groups and observations within
    groups.

    This function assumes that the table has many more observations than
    features, and for this reason, the design orients features across the
    vertical axis and observations across the horizontal axis. With a landscape
    aspect ratio, this design allows more space for the horizontal axis.

    ----------

    Reference:
    https://www.kaggle.com/code/sgalella/correlation-heatmaps-with-
        hierarchical-clustering
    https://www.kaggle.com/code/chinoysen/beginners-guide-to-heatmaps-cluster-
        heatmap
    https://seaborn.pydata.org/generated/seaborn.clustermap.html
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.heatmap.html
    https://gist.github.com/peterk87/5505691
    https://matplotlib.org/stable/tutorials/colors/colormaps.html

    Maybe useful for problem-solving:
    https://matplotlib.org/stable/gallery/subplots_axes_and_figures/
        subplots_demo.html
        - subplots and gridspec
    http://www.acgeospatial.co.uk/colour-bar-for-discrete-rasters-with-matplotlib/
        - color maps for discrete values
    https://matplotlib.org/stable/users/explain/colors/colorbar_only.html
        #colorbar-with-arbitrary-colors
        - color bar with custom dimensions of discrete, categorical intervals

    https://matplotlib.org/stable/gallery/color/colormap_reference.html
        - color maps in MatPlotLib

    https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
        - categorical color bars
        - helpful reference for figuring out how to set the thresholds between
          the discrete levels of the color bar to correspond properly to the
          columns of the main heatmap...
    https://stackoverflow.com/questions/9707676/defining-a-discrete-colormap-for-imshow
    https://stackoverflow.com/questions/7229971/2d-grid-data-visualization-in-python
    https://matplotlib.org/stable/users/explain/colors/colormapnorms.html
        - normalization of information for color maps???

    Review: 30 December 2024

    arguments:
        table_signal (object): Pandas data-frame table of values of signal
            intensity for features in columns across sample observations or
            groups of sample observations in rows
        table_feature_sets (object): Pandas data-frame table of indications of
            allocation of features to sets in a sort sequence that matches the
            sequence of features across columns in table of signals
        format_table_signal (int): value 1 for features across rows and
            observations or groups of observations across columns, value 2 for
            features across columns and observations across rows with potential
            special column for groups of observations
        index_columns (str): name to define an index corresponding to
            information across columns in source table
        index_rows (str): name of a column in source table which defines an
            index corresponding to information across rows
        column_group (str): name of column in table to use for groups
        transpose_table (bool): whether to transpose matrix from table
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_missing_fill (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
        title_ordinate (str): title for ordinate vertical axis
        title_abscissa (str): title for abscissa horizontal axis
        title_bar (str): title for scale bar
        size_title_ordinate (str): font size
        size_title_abscissa (str): font size
        size_title_bar (str): font size
        size_label_feature_set (str): font size
        size_label_legend_observation_group (str): font size
        size_label_bar (str): font size
        show_scale_bar (bool): whether to create scale bar
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): MatPlotLib figure object

    """

    ##########
    # Prepare information for figure.
    pail = extract_prepare_table_signals_categories_for_heatmap(
        table=table_signal,
        format_table=format_table_signal, # 1: features in rows, observations in columns
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=column_group,
        transpose_table=transpose_table,
        fill_missing=fill_missing,
        value_missing_fill=value_missing_fill,
        constrain_signal_values=constrain_signal_values,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        report=report,
    )
    # Copy information in table.
    table_feature_sets = table_feature_sets.copy(deep=True)
    # Organize indices in table.
    table_feature_sets.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_feature_sets.set_index(
        index_columns,
        append=False,
        drop=True,
        inplace=True
    )
    table_feature_sets.columns.rename(
        "feature_sets",
        inplace=True,
    ) # single-dimensional index
    # Extract categorical names from original indices in table.
    labels_feature_sets_index_columns = copy.deepcopy(
        table_feature_sets.columns.get_level_values(
            "feature_sets"
        ).to_list()
    )
    labels_feature_sets_index_rows = copy.deepcopy(
        table_feature_sets.index.get_level_values(
            index_columns
        ).to_list()
    )
    # Extract values.
    #matrix_signal = numpy.transpose(numpy.copy(table_signal.to_numpy()))
    matrix_feature_sets = numpy.copy(table_feature_sets.to_numpy())

    ##########
    # Create figure.

    ##########
    # Initialize figure.
    figure = pplot.initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Initialize grid within figure.
    #matplotlib.gridspec.GridSpec()
    #figure.add_gridspec()
    #sharey=True, # sharey="row"
    #sharex=True, # sharex="col",
    grid = matplotlib.gridspec.GridSpec(
        nrows=4,
        ncols=2,
        wspace=0.005, # horizontal width space between grid blocks for subplots
        hspace=0.005, # vertical height space between grid blocks for subplots
        width_ratios=(15,85), # first column 1/10th width of second column
        height_ratios=(90,3,5,2), # first row 30 times the height of second row
    )
    grid.update(
        wspace=0.005, # horizontal width space between grid blocks for subplots
        hspace=0.005, # vertical height space between grid blocks for subplots
    )
    # Initialize axes within grid within figure.

    # Initialize axes within grid within figure.
    axes_set = figure.add_subplot(grid[0,0]) # first row, first column
    axes_main = figure.add_subplot(grid[0,1]) # first row, second column
    axes_space_1 = figure.add_subplot(grid[1,0]) # second row, first column
    axes_group = figure.add_subplot(grid[1,1]) # second row, second column
    axes_space_2 = figure.add_subplot(grid[2,0]) # third row, first column
    axes_space_3 = figure.add_subplot(grid[2,1]) # third row, second column
    axes_space_4 = figure.add_subplot(grid[3,0]) # fourth row, first column
    axes_bar = figure.add_subplot(grid[3,1]) # fourth row, second column
    axes_set.clear()
    axes_main.clear()
    axes_group.clear()
    axes_bar.clear()
    axes_space_1.clear()
    axes_space_2.clear()
    axes_space_3.clear()
    axes_space_4.clear()
    # Set axes to empty as a space holder.
    axes_space_1.axis("off")
    axes_space_2.axis("off")
    axes_space_3.axis("off")
    axes_space_4.axis("off")
    #axes_set.axis("off")
    #axes_group.axis("off")
    #axes_bar.axis("off")
    #axes_main.axis("off")
    # Keep axes, ticks, and labels, but remove border.
    for position in ['right', 'top', 'bottom', 'left']:
        matplotlib.pyplot.gca().spines[position].set_visible(False)
    # Adjust margins.
    # Method "tight_layout()" does not function properly with object
    # "GridSpec()".
    #grid.tight_layout(
    #    figure,
    #    #pad=1.0,
    #    #h_pad=1.0,
    #    #w_pad=1.0,
    #    rect=[0,0.05,1.0,1.0], # left, bottom, right, top
    #)
    #grid.update(
    #    wspace=0.005, # horizontal width space between grid blocks for subplots
    #    hspace=0.005, # vertical height space between grid blocks for subplots
    #)
    figure.subplots_adjust(
        left=0.02,
        right=0.98,
        top=0.98,
        bottom=0.05,
    )

    ##########
    # axes: main
    # Plot values as a grid of color on continuous scale.
    # This function represents values acros matrix dimension 0 as vertical
    # rows.
    # This function represents values across matrix dimension 1 as horizontal
    # columns.
    # Diverging color maps: "PRGn", "PRGn_r", "PiYG", "PiYG_r",
    # Diverging color maps: "PuOr", "PuOr_r",
    # Diverging color maps: "PuOr", "PuOr_r", "RdBu", "RdBu_r", "BrBG",
    # Sequential color maps: "Reds", "Reds_r", "Oranges", "Oranges_r",
    # site: https://montoliu.naukas.com/2021/11/18/color-blindness-purple-and-
    #     orange-are-the-solution/
    image_main = axes_main.imshow(
        pail["matrix_signal"],
        cmap=matplotlib.colormaps["PuOr"], # binary, Reds, RdBu_r, PuOr, PuOr_r
        vmin=pail["value_minimum"],
        vmax=pail["value_maximum"],
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set titles for axes.
    if (len(title_ordinate) > 0):
        axes_main.set_ylabel(
            ylabel=title_ordinate,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_ordinate]
        )
        pass
    if (len(title_abscissa) > 0):
        axes_main.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
        pass
    # Manage tick gradations, labels, and titles on axes.
    axes_main.tick_params(
        axis="both", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        left=False,
        right=False,
        labelleft=False,
        labelright=False,
    )

    ##########
    # axes: set
    # Define color map for discrete, binary integer representation of
    # allocation to sets.
    # https://matplotlib.org/3.1.1/gallery/color/named_colors.html
    # "white", "black"
    # "white", "dimgray"
    color_map_set = matplotlib.colors.ListedColormap([
        "white", "dimgray"
    ])
    # Plot values as a grid of color on discrete scale.
    image_set = axes_set.imshow(
        matrix_feature_sets,
        cmap=color_map_set,
        vmin=0,
        vmax=1,
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )
    # Set tick positions and labels on axes.
    axes_set.set_xticks(
        numpy.arange(matrix_feature_sets.shape[1]),
    )
    axes_set.set_xticklabels(
        labels_feature_sets_index_columns,
        #minor=False,
        ha="center", # horizontal alignment
        va="top", # vertical alignment
        alpha=1.0,
        rotation=90, # negative: clockwise; positive: count-clockwise
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"][size_label_feature_set],
    )
    # Set tick parameters for axes.
    axes_set.tick_params(
        axis="y", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        left=False,
        right=False,
        labelleft=False,
        labelright=False,
    )
    axes_set.tick_params(
        axis="x", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=5.0, # 5.0
        width=3.5, # 3.5
        pad=10, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=False,
        bottom=True,
        labeltop=False,
        labelbottom=True,
    )

    ##########
    # axes: group
    # Define color map for discrete, integer representation of categorical
    # groups.
    #discrete_minimum = 0
    #discrete_maximum = len(pail["labels_group_unique"])
    discrete_minimum = numpy.nanmin(pail["matrix_group_integers"])
    discrete_maximum = numpy.nanmax(pail["matrix_group_integers"])
    color_map_group = matplotlib.colormaps.get_cmap(
        "tab10", # "Set1", "Set2", "Dark2", "tab10",
        ((discrete_maximum - discrete_minimum) + 1)
    )
    # Plot values as a grid of color on discrete scale.
    image_group = axes_group.imshow(
        pail["matrix_group_integers"],
        cmap=color_map_group,
        vmin=discrete_minimum,
        vmax=discrete_maximum,
        aspect="auto", # "auto", "equal",
        origin="lower",
        # Extent: (left, right, bottom, top)
        #extent=(
        #    0.0,
        #    (pail["matrix_group_integers"].shape[1]),
        #    0.0,
        #    (3*(pail["matrix_group_integers"].shape[0]))
        #),
    )
    axes_group.tick_params(
        top=False,
        bottom=False,
        left=False,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=False,
    )
    # Create legend.
    # See nested function "create_legend_elements()" within main function
    # "plot_scatter_factor_groups()".
    handles = [matplotlib.patches.Patch(
        color=color_map_group(i),
        label=pail["indices_groups"][i]
    ) for i in pail["indices_groups"].keys()]
    axes_group.legend(
        handles=handles,
        prop=fonts["properties"][size_label_legend_observation_group],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=4,
    )

    ##########
    # axes: bar
    # Create legend for scale of color grid.
    # Notice that use of the "cax" argument causes to ignore the "shrink"
    # argument.
    if show_scale_bar:
        # Set definitive scale to avoid conflicts with global objects.
        #color_map_main = matplotlib.pyplot.cm.ScalarMappable(
        #    cmap="PuOr",
        #    norm=matplotlib.pyplot.Normalize(
        #        vmin=pail["value_minimum"],
        #        vmax=pail["value_maximum"],
        #    ),
        #)
        #color_map_main.set_array(pail["matrix_signal"])
        # Create scale bar.
        bar_main = axes_bar.figure.colorbar(
            image_main, # image_main, color_map_main
            orientation="horizontal",
            cax=axes_bar,
            location="bottom",
            #shrink=0.5, # 0.7; factor for dimensions of the Scale Bar.
        )
        if (len(title_bar) > 0):
            bar_main.ax.set_xlabel(
                title_bar,
                rotation=0,
                loc="center",
                va="bottom",
                labelpad=20, # 5
                alpha=1.0,
                backgroundcolor=colors["white"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_title_bar],
            )
        bar_main.ax.tick_params(
            axis="both",
            which="both", # major, minor, or both
            direction="out",
            length=7.5, # 5.0, 7.5
            width=3, # 2.5, 5.0
            color=colors["black"],
            pad=3, # 5, 7
            labelsize=fonts["values"][size_label_bar]["size"],
            labelcolor=colors["black"],
        )
        pass

    ##########
    # Return figure.
    return figure


##########
# Scatter points.


def create_confidence_ellipse(
    axes=None,
    values_abscissa=None,
    values_ordinate=None,
    factor_confidence_ellipse=None,
    size_edge=None,
    color_edge=None,
    color_fill=None,
    **kwargs,
):
    """
    Create and plot an ellipse as a representation of the covariance confidence
    corresponding to values of features represented on the abscissa (x) and
    ordinate (y) axes.

    References:
    1. https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html
    2. https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.
        Ellipse.html#matplotlib.patches.Ellipse
    3. https://stats.stackexchange.com/questions/391706/drawing-95-ellipse-
        over-scatter-plot

    Review or revision: TCW; 4 November 2025
    Review or revision: TCW; 31 October 2025

    arguments:
        axes (object): MatPlotLib axes object
        values_abscissa (object): Numpy array of values for representation on
            the horizontal abscissa 'x' axis
        values_ordinate (object): Numpy array of values for representation on
            the vertical ordinate 'y' axis
        factor_confidence_ellipse (float): factor by which to scale the
            dimensions of the confidence ellipse as multiples of the standard
            deviation
        size_edge (float): size for line width of edge
        color_edge (object): color for edge of the ellipse
        color_fill (object): color to fill the ellipse
        kwargs (object): additional key-word arguments to forward to function
            or class 'matplotlib.patches.Ellipse()'
    raises:
        warning about dimesions of the the two NumPy arrays of values

    returns:
        (object): chart figure object 'matplotlib.patches.Ellipse'

    """

    # Warning.
    if (values_abscissa.size != values_ordinate.size):
        raise ValueError("x and y must be the same size")
        pass

    # Calculate covariance.
    cov = numpy.cov(values_abscissa, values_ordinate)
    pearson = cov[0, 1]/numpy.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = numpy.sqrt(1 + pearson)
    ell_radius_y = numpy.sqrt(1 - pearson)
    ellipse = matplotlib.patches.Ellipse(
        (0, 0),
        width=(ell_radius_x * 2),
        height=(ell_radius_y * 2),
        linestyle="--", # "--", ":",
        linewidth=size_edge,
        edgecolor=color_edge,
        facecolor=color_fill,
        alpha=0.35,
        zorder=0,
        **kwargs,
    )

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = numpy.sqrt(cov[0, 0]) * factor_confidence_ellipse
    mean_x = numpy.mean(values_abscissa)

    # calculating the standard deviation of y ...
    scale_y = numpy.sqrt(cov[1, 1]) * factor_confidence_ellipse
    mean_y = numpy.mean(values_ordinate)

    # Transform coordinates.
    transf = matplotlib.transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)
    ellipse.set_transform(transf + axes.transData)
    # Return object.
    return axes.add_patch(ellipse)


def plot_scatter_point_color_response_discrete_or_continuous(
    table=None,
    column_identifier=None,
    column_name=None,
    column_response_markers=None,
    column_group_ellipses=None,
    column_abscissa=None,
    column_ordinate=None,
    type_response=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    title_chart=None,
    title_response=None,
    title_abscissa=None,
    title_ordinate=None,
    identifiers_emphasis=None,
    size_title_chart=None,
    size_title_abscissa=None,
    size_title_ordinate=None,
    size_title_legend_bar=None,
    size_label_abscissa=None,
    size_label_ordinate=None,
    size_label_emphasis=None,
    size_label_legend_bar=None,
    size_marker=None,
    size_edge_marker=None,
    size_edge_ellipse=None,
    factor_confidence_ellipse=None,
    aspect=None,
    fonts=None,
    colors_fill_markers=None,
    colors_fill_ellipses=None,
    color_edge_markers=None,
    color_edge_ellipses=None,
    color_emphasis=None,
    set_axis_limits=None,
    show_confidence_ellipse=None,
    show_lines_origin=None,
    show_line_diagonal=None,
    show_emphasis_marker=None,
    show_emphasis_label=None,
    show_legend_bar=None,
    report=None,
):
    """
    Create a plot chart of type scatter point with color representation of a
    third feature variable on either a discrete on continuous measurement
    scale.

    Scatter Plot with color representation of a third feature for points
    chart type: scatter
    response
       - representation: color of points
    abscissa
       - representation: position along horizontal abscissa x axis
    ordinate
       - representation: position along vertical ordinate y axis

    Review or revision: TCW; 4 November 2025
    Review or revision: TCW; 31 October 2025
    Review or revision: TCW; 1 October 2025
    Review or revision: TCW; 2 December 2024

    arguments:
        table (object): Pandas data-frame table of features across columns and
            values for observations across rows
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
        minimum_abscissa (float): value for minimal limit to represent on the
            abscissa horizontal axis
        maximum_abscissa (float): value for maximal limit to represent on the
            abscissa horizontal axis
        minimum_ordinate (float): value for minimal limit to represent on the
            ordinate horizontal axis
        maximum_ordinate (float): value for maximal limit to represent on the
            ordinate horizontal axis
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
        size_title_chart (str): font size for title of chart
        size_title_abscissa (str): font size for title on abscissa horizontal
            axis
        size_title_ordinate (str): font size for title on ordinate vertical
            axis
        size_title_legend_bar (str): font size for title on the legend for
            discrete categorical values or on the scale bar for continuous
            values
        size_label_abscissa (str): font size for labels on abscissa horizontal
            axis
        size_label_ordinate (str): font size for labels on ordinate vertical
            axis
        size_label_emphasis (str): font size for labels adjacent to points for
            special emphasis
        size_label_legend_bar (str): font size for labels on the legend or on
            the scale bar
        size_marker (int): size of markers for points representing values
        size_edge_marker (float): size for line width of edge
        size_edge_ellipse (float): size for line width of edge
        factor_confidence_ellipse (float): factor product for standard
            deviations in the confidence ellipses
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): definitions of font properties
        colors_fill_markers (list<tuple>): definitions of color properties
        colors_fill_ellipses (list<tuple>): definitions of color properties
        color_edge_markers (tuple): definition of color properties
        color_edge_ellipses (tuple): definition of color properties
        color_emphasis (tuple): definition of color properties
        set_axis_limits (bool): whether to set explicity limits on axes
        show_confidence_ellipse (bool): whether to show confidence ellipse;
            show empty ellipse edge for continuous response; show color-matched
            distinct ellipses for discrete categorical groups
        show_lines_origin (bool): whether to draw vertical and horizontal lines to
            represent origin (zero) on abscissa and ordinate axes, respectively
        show_line_diagonal (bool): whether to draw diagonal line for equality
            between abscissa and ordinate
        show_emphasis_marker (bool): whether to create special markers to emphasize
            a special selection of points
        show_emphasis_label (bool): whether to create text labels adjacent to
            the special selection of points for special emphasis
        show_legend_bar (bool): whether to show legend or scale bar on chart
        report (bool): whether to print reports

    raises:

    returns:
        (object): chart figure object

    """

    ##########
    # Organize information for chart.

    # Copy information.
    table = table.copy(deep=True)
    # Filter columns in table.
    #table = table.loc[
    #    :, table.columns.isin(columns_sequence)
    #]
    table = table.filter(
        items=[
            column_identifier,
            column_name,
            column_response_markers,
            column_group_ellipses,
            column_abscissa,
            column_ordinate,
        ],
        axis="columns",
    )
    # Organize information in table.
    table.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    limits = [
        {
            "type": "low",
            "value": minimum_abscissa,
            "column": column_abscissa,
        },
        {
            "type": "high",
            "value": maximum_abscissa,
            "column": column_abscissa,
        },
        {
            "type": "low",
            "value": minimum_ordinate,
            "column": column_ordinate,
        },
        {
            "type": "high",
            "value": maximum_ordinate,
            "column": column_ordinate,
        },
    ]
    for limit in limits:
        if (limit["value"] is not None):
            if (limit["type"] == "low"):
                table = table.loc[
                    (table[limit["column"]] >= limit["value"]), :
                ].copy(deep=True)
            elif (limit["type"] == "high"):
                table = table.loc[
                    (table[limit["column"]] <= limit["value"]), :
                ].copy(deep=True)
                pass
            pass
        pass
    # Extract information about ranges of first and second features for
    # representation as position.
    values_abscissa_raw = table[column_abscissa].to_numpy(
        dtype="float64",
        na_value=numpy.nan,
        copy=True,
    )
    values_ordinate_raw = table[column_ordinate].to_numpy(
        dtype="float64",
        na_value=numpy.nan,
        copy=True,
    )
    if (minimum_abscissa is None):
        minimum_abscissa = (numpy.nanmin(values_abscissa_raw))
    if (maximum_abscissa is None):
        maximum_abscissa = (numpy.nanmax(values_abscissa_raw))
    if (minimum_ordinate is None):
        minimum_ordinate = (numpy.nanmin(values_ordinate_raw))
    if (maximum_ordinate is None):
        maximum_ordinate = (numpy.nanmax(values_ordinate_raw))
        pass
    center_abscissa = ((maximum_abscissa + minimum_abscissa)/2)
    center_ordinate = ((maximum_ordinate + minimum_ordinate)/2)

    # Extract information about selection of records for special emphasis.
    if (
        (identifiers_emphasis is not None) and
        (len(identifiers_emphasis) > 0)
    ):
        table_standard = table.loc[
            ~table[column_identifier].isin(identifiers_emphasis), :
        ].copy(deep=True)
        table_special = table.loc[
            table[column_identifier].isin(identifiers_emphasis), :
        ].copy(deep=True)
        pass
    else:
        table_standard = table.copy(deep=True)
        table_special = pandas.DataFrame()
        pass

    # Extract information about range of third response feature for
    # representation as color of markers for points.
    if (
        (column_response_markers is not None) and
        (column_response_markers != "") and
        (type_response is not None) and
        (type_response == "continuity")
    ):
        values_response_raw = table[column_response_markers].to_numpy(
            dtype="float64",
            na_value=numpy.nan,
            copy=True,
        )
        minimum_response = numpy.nanmin(values_response_raw)
        maximum_response = numpy.nanmax(values_response_raw)
        pass

    # Copy information.
    table_response = table.copy(deep=True)
    table_ellipse = table.copy(deep=True)

    # Stratify groups of observations by categorical factor response feature.
    if (
        (column_response_markers is not None) and
        (column_response_markers != "") and
        (type_response is not None) and
        (type_response == "category")
    ):
        # Count unique categorical values.
        #categories = copy.deepcopy(
        #    table_response[column_category].unique().tolist()
        #)
        #count_category_response = len(categories)
        count_categories_response = int(
            table_response[column_response_markers].nunique(dropna=True)
        )
        # Organize indices in table.
        table_response.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table_response.set_index(
            [column_identifier, column_response_markers],
            append=False,
            drop=True,
            inplace=True
        )
        # Split rows within table by factor columns.
        groups_response = table_response.groupby(
            level=column_response_markers,
        )
        pass

    # Stratify groups of observations by categorical factor ellipse feature.
    if (
        (column_group_ellipses is not None) and
        (column_group_ellipses != "") and
        (show_confidence_ellipse)
    ):
        # Count unique categorical values.
        #categories = copy.deepcopy(
        #    table_ellipse[column_category].unique().tolist()
        #)
        #count_category_ellipse = len(categories)
        count_categories_ellipse = (
            table_ellipse[column_group_ellipses].nunique(dropna=True)
        )
        # Organize indices in table.
        table_ellipse.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table_ellipse.set_index(
            [column_identifier, column_group_ellipses],
            append=False,
            drop=True,
            inplace=True
        )
        # Split rows within table by factor columns.
        groups_ellipse = table_ellipse.groupby(
            level=column_group_ellipses,
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_special.py"
        )
        print(str("module: " + module))
        function = "plot_scatter_point_color_response_discrete_or_continuous()"
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Organize parameters for visual representation.

    # Colors.
    # Determine whether the parameters define adequate custom colors.
    # Markers.
    if (
        (column_response_markers is not None) and
        (column_response_markers != "") and
        (type_response is not None) and
        (type_response == "category") and
        (count_categories_response > 0)
    ):
        if (
            (colors_fill_markers is None) or
            (len(colors_fill_markers) < count_categories_response)
        ):
            # Define colors.
            # Create discrete color map for categorical values.
            # color_map = matplotlib.pyplot.get_cmap()
            color_map_fill_markers = matplotlib.colormaps.get_cmap("tab10") # "Set1", "Set2", "Dark2", "tab10",
        else:
            color_map_fill_markers = matplotlib.colors.ListedColormap(
                colors_fill_markers
            )
        pass
    if (color_edge_markers is None):
        color_set_edge_markers = matplotlib.colors.to_rgba("black", 1.0)
    elif (not isinstance(color_edge_markers, tuple)):
        color_set_edge_markers = matplotlib.colors.to_rgba(
            color_edge_markers, 1.0
        )
    else:
        color_set_edge_markers = color_edge_markers
        pass

    # Ellipses.
    if (
        (column_group_ellipses is not None) and
        (column_group_ellipses != "") and
        (show_confidence_ellipse) and
        (count_categories_ellipse > 0)
    ):
        if (
            (colors_fill_ellipses is None) or
            (len(colors_fill_ellipses) < count_categories_ellipse)
        ):
            # Define colors.
            # Create discrete color map for categorical values.
            # color_map = matplotlib.pyplot.get_cmap()
            color_map_fill_ellipses = matplotlib.colormaps.get_cmap(
                "tab10",
                count_categories_ellipse,
            ) # "Set1", "Set2", "Dark2", "tab10",
        else:
            color_map_fill_ellipses = matplotlib.colors.ListedColormap(
                colors_fill_ellipses
            )
        pass
    if (
        (color_edge_ellipses is None)
    ):
        color_set_edge_ellipses = matplotlib.colors.to_rgba("black", 1.0)
    elif (not isinstance(color_edge_ellipses, tuple)):
        color_set_edge_ellipses = matplotlib.colors.to_rgba(
            color_edge_ellipses, 1.0
        )
    else:
        color_set_edge_ellipses = color_edge_ellipses
        pass

    # Emphasis.
    if (
        (color_emphasis is None)
    ):
        color_set_emphasis = matplotlib.colors.to_rgba("orange", 1.0)
    elif (not isinstance(color_emphasis, tuple)):
        color_set_emphasis = matplotlib.colors.to_rgba(
            color_emphasis, 1.0
        )
    else:
        color_set_emphasis = color_emphasis
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
    #axes.margins(
    #    x=1,
    #    y=1,
    #    tight=True,
    #)
    # Define limits for axes.
    if (set_axis_limits):
        if (minimum_abscissa is not None):
            axes.set_xlim(xmin=minimum_abscissa)
        if (maximum_abscissa is not None):
            axes.set_xlim(xmax=maximum_abscissa)
        if (minimum_ordinate is not None):
            axes.set_ylim(ymin=minimum_ordinate)
        if (maximum_ordinate is not None):
            axes.set_ylim(ymax=maximum_ordinate)
            pass
        pass
    # Include title label on chart.
    if len(title_chart) > 0:
        axes.set_title(
            title_chart,
            fontproperties=fonts["properties"][size_title_chart],
            loc="center",
            horizontalalignment="center",
            verticalalignment="top",
            pad=5,
        )
        pass
    # Set titles for axes.
    if (len(title_abscissa) > 0):
        axes.set_xlabel(
            xlabel=title_abscissa,
            labelpad=20,
            alpha=1.0,
            #loc="left",
            backgroundcolor=matplotlib.colors.to_rgba("white", 1.0), # colors["white"]
            color=matplotlib.colors.to_rgba("black", 1.0), # colors["black"]
            fontproperties=fonts["properties"][size_title_abscissa]
        )
    if (len(title_ordinate) > 0):
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=20,
            alpha=1.0,
            backgroundcolor=matplotlib.colors.to_rgba("white", 1.0), # colors["white"]
            color=matplotlib.colors.to_rgba("black", 1.0), # colors["black"]
            fontproperties=fonts["properties"][size_title_ordinate]
        )
    # Define parameters for tick labels on axes.
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

    # Create lines to represent origins (zero) on axes.
    if (show_lines_origin):
        axes.axhline(
            y=0.0,
            #xmin=minimum_abscissa,
            #xmax=maximum_abscissa,
            alpha=1.0,
            color=matplotlib.colors.to_rgba("black", 1.0),
            linestyle="--",
            linewidth=2.5,
        )
        axes.axvline(
            x=0.0,
            #ymin=minimum_ordinate,
            #ymax=maximum_ordinate,
            alpha=1.0,
            color=matplotlib.colors.to_rgba("black", 1.0),
            linestyle="--",
            linewidth=2.5,
        )
        pass

    # Create diagonal line to represent equality between abscissa and ordinate.
    # Notice that the current definition does not actually correspond to a 1:1
    # relationship between abscissa and ordinate scales. Instead, it depends on
    # the relative ranges of the abscissa and ordinate axes.
    if (show_line_diagonal):
        axes.plot(
            [0, 1,],
            [0, 1,],
            transform=axes.transAxes,
            alpha=1.0,
            color=matplotlib.colors.to_rgba("black", 1.0),
            linestyle="--",
            linewidth=2,
        )
        axes.text(
            0.02,
            0.025,
            str(
                "diagonal line for visual reference only; " +
                "not proportional to respective ranges of axes"
            ),
            transform=axes.transAxes, # positions in coordinates relative to scope of axes
            transform_rotates_text=True, # match rotation to scope of axes
            rotation_mode="anchor",
            rotation=45,
            backgroundcolor=matplotlib.colors.to_rgba("white", 0.75), # colors["white_faint"],
            color=matplotlib.colors.to_rgba("black", 1.0),
            fontproperties=fonts["properties"]["seventeen"],
        )
        pass

    ##########
    # Represent information on the chart figure object.

    # Markers for points.

    # For table of standard selection of records, determine whether to set
    # color of points to represent a third response feature.
    if (
        (column_response_markers is None) or
        (column_response_markers == "") or
        (type_response is None) or
        (type_response == "")
    ):
        # Plot points to represent main observations.
        handle_standard = axes.plot(
            table[column_abscissa].values,
            table[column_ordinate].values,
            linestyle="",
            linewidth=size_edge_marker,
            marker="o",
            markersize=size_marker,
            markeredgecolor=color_set_edge_markers,
            markerfacecolor=color_map_fill_markers(0),
        )

    elif (
        (column_response_markers is not None) and
        (column_response_markers != "") and
        (type_response is not None) and
        (type_response == "continuity")
    ):
        # Plot points to represent main observations.
        handle_standard = axes.scatter(
            table[column_abscissa].values,
            table[column_ordinate].values,
            c=table[column_response_markers].values,
            s=(size_marker**2), # scale with area (dimension squared)
            norm="linear",
            cmap="binary", # 'binary', 'plasma', 'viridis', 'civids'
            vmin=minimum_response,
            vmax=maximum_response,
            alpha=1,
            marker="o",
            linewidths=size_edge_marker,
            edgecolors=color_set_edge_markers,
            #linestyle="",
        )

    elif (
        (column_response_markers is not None) and
        (column_response_markers != "") and
        (type_response is not None) and
        (type_response == "category") and
        (count_categories_response > 0)
    ):
        # An alternative would be to map categorical values to discrete
        # integers that respectively map to colors. Then create discrete color
        # bar as a reference in place of a legend.
        # Collect information.
        labels_groups_response = []
        counter = 0
        # Iterate on groups of records from table.
        for name_group, table_group in groups_response:
            # Copy information in table.
            table_group = table_group.copy(deep=True)
            # Collect information.
            label_group = str(name_group)
            labels_groups_response.append(label_group)
            # Plot points to represent main observations.
            handle_group = axes.plot(
                table_group[column_abscissa].values,
                table_group[column_ordinate].values,
                linestyle="",
                marker="o",
                markersize=size_marker,
                markeredgewidth=size_edge_marker,
                markeredgecolor=color_set_edge_markers, # colors_groups[counter]
                markerfacecolor=color_map_fill_markers(counter),
            )
            # Update index counter.
            counter += 1
            pass
        pass

    # Confidence ellipses.
    # Create confidence ellipse for groups of observations.
    if (
        (show_confidence_ellipse) and
        (count_categories_ellipse > 0)
    ):
        # Collect information.
        labels_groups_ellipse = []
        counter = 0
        # Iterate on groups of records from table.
        for name_group, table_group in groups_ellipse:
            # Copy information in table.
            table_group = table_group.copy(deep=True)
            # Collect information.
            label_group = str(name_group)
            labels_groups_ellipse.append(label_group)
            # Create ellipse to represent confidence in range of observations
            # in  group.
            create_confidence_ellipse(
                axes=axes,
                values_abscissa=table_group[column_abscissa].values,
                values_ordinate=table_group[column_ordinate].values,
                factor_confidence_ellipse=factor_confidence_ellipse,
                size_edge=size_edge_ellipse,
                color_edge=color_set_edge_ellipses,
                color_fill=color_map_fill_ellipses(counter),
            )
            # Update index counter.
            counter += 1
            pass
        pass

    # For table of special selection of records, create text labels adjacent
    # to points.
    if (
        (identifiers_emphasis is not None) and
        (len(identifiers_emphasis) > 0) and
        (not table_special.empty) and
        (show_emphasis_label)
    ):
        for index, row in table_special.iterrows():
            # Determine position coordinates of label.
            abscissa_label_raw = row[column_abscissa]
            if (abscissa_label_raw >= center_abscissa):
                alignment_horizontal = "right"
                abscissa_label = (
                    abscissa_label_raw - (
                        0.025 * (maximum_abscissa - minimum_abscissa)
                    )
                )
            if (abscissa_label_raw < center_abscissa):
                alignment_horizontal = "left"
                abscissa_label = (
                    abscissa_label_raw + (
                        0.025 * (maximum_abscissa - minimum_abscissa)
                    )
                )
            ordinate_label = row[column_ordinate]
            # Create label on chart.
            axes.text(
                abscissa_label,
                ordinate_label,
                str(row[column_name]),
                horizontalalignment=alignment_horizontal,
                verticalalignment="center",
                backgroundcolor=matplotlib.colors.to_rgba("white", 0.75),
                color=matplotlib.colors.to_rgba("black", 1.0),
                fontproperties=fonts["properties"][size_label_emphasis],
            )
            pass
        pass

    # For table of special selection of records, plot points with a special
    # color for emphasis.
    if (
        (identifiers_emphasis is not None) and
        (len(identifiers_emphasis) > 0) and
        (not table_special.empty) and
        (show_emphasis_marker)
    ):
        handle_special = axes.plot(
            table_special[column_abscissa].values,
            table_special[column_ordinate].values,
            linestyle="",
            marker="o",
            markersize=(size_marker*2),
            markeredgecolor=color_set_emphasis,
            markeredgewidth=3.0,
            markerfacecolor="None"
        )

    # Create legend or reference bar for color map.
    if (
        (show_legend_bar) and
        (column_response_markers is not None) and
        (column_response_markers != "") and
        (type_response is not None) and
        (type_response == "continuity")
    ):
        bar = axes.figure.colorbar(
            handle_standard,
            orientation="vertical",
            ax=axes,
            location="right",
            shrink=0.9, # 0.7; factor for dimensions of the Scale Bar.
        )
        if (len(title_response) > 0):
            bar.ax.set_ylabel(
                title_response,
                rotation=-90,
                va="bottom",
                labelpad=5, # 5
                alpha=1.0,
                backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
                color=matplotlib.colors.to_rgba("black", 1.0),
                fontproperties=fonts["properties"][size_title_legend_bar],
            )
        else:
            bar.ax.set_ylabel(
                column_response_markers,
                rotation=-90,
                va="bottom",
                labelpad=5, # 5
                alpha=1.0,
                backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
                color=matplotlib.colors.to_rgba("black", 1.0),
                fontproperties=fonts["properties"][size_title_legend_bar],
            )
            pass
        bar.ax.tick_params(
            axis="both",
            which="both", # major, minor, or both
            direction="out",
            length=7.5, # 5.0, 7.5
            width=3, # 2.5, 5.0
            color=matplotlib.colors.to_rgba("black", 1.0),
            pad=5, # 5, 7
            labelsize=fonts["values"][size_label_legend_bar]["size"],
            labelcolor=matplotlib.colors.to_rgba("black", 1.0),
        )
    elif (
        (show_legend_bar) and
        (column_response_markers is not None) and
        (column_response_markers != "") and
        (type_response is not None) and
        (type_response == "category")
    ):
        # Create legend.
        handles_legend = [matplotlib.patches.Patch(
            color=color_map_fill_markers(i), # colors_groups[i]
            label=labels_groups_response[i]#,
            #size?
        ) for i in range(count_categories_response)]
        axes.legend(
            handles=handles_legend,
            loc="upper right", # "upper right", "upper center", "lower right"
            #bbox_to_anchor=(0.5, -0.05),
            #ncol=4,
            prop=fonts["properties"][size_label_legend_bar],
            title=title_response,
            title_fontsize=fonts["values"][size_label_legend_bar]["size"],
        )
        pass

    ##########
    # Return figure.
    return figure


def create_write_plot_chart_scatter_point_response(
    path_directory_parent=None,
    name_chart=None,
    table=None,
    column_identifier=None,
    column_name=None,
    column_response_markers=None,
    column_group_ellipses=None,
    column_abscissa=None,
    column_ordinate=None,
    type_response=None,
    title_chart=None,
    title_response=None,
    title_abscissa=None,
    title_ordinate=None,
    identifiers_emphasis=None,
    size_marker=None,
    size_edge_marker=None,
    size_edge_ellipse=None,
    factor_confidence_ellipse=None,
    colors_fill_markers=None,
    colors_fill_ellipses=None,
    color_edge_markers=None,
    color_edge_ellipses=None,
    color_emphasis=None,
    show_confidence_ellipse=None,
    show_emphasis_marker=None,
    show_emphasis_label=None,
    show_legend_bar=None,
    report=None,
):
    """
    Create and plot a chart of the scatter point type.

    Review or revision: TCW; 4 November 2025
    Review or revision: TCW; 31 October 2025
    Review or revision: TCW; 4 November 2025

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
    figure = plot_scatter_point_color_response_discrete_or_continuous(
        table=table,
        column_identifier=column_identifier,
        column_name=column_name,
        column_response_markers=column_response_markers,
        column_group_ellipses=column_group_ellipses,
        column_abscissa=column_abscissa,
        column_ordinate=column_ordinate,
        type_response=type_response,
        minimum_abscissa=None,
        maximum_abscissa=None,
        minimum_ordinate=None,
        maximum_ordinate=None,
        title_chart=title_chart,
        title_response=title_response,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        identifiers_emphasis=identifiers_emphasis,
        size_title_chart="nine",
        size_title_abscissa="nine",
        size_title_ordinate="nine",
        size_title_legend_bar="thirteen",
        size_label_abscissa="thirteen",
        size_label_ordinate="thirteen",
        size_label_emphasis="sixteen",
        size_label_legend_bar="thirteen",
        size_marker=size_marker,
        size_edge_marker=size_edge_marker,
        size_edge_ellipse=size_edge_ellipse,
        factor_confidence_ellipse=factor_confidence_ellipse, # 2.0 * standard deviation (95% of normal distribution)
        aspect="landscape",
        fonts=fonts,
        colors_fill_markers=colors_fill_markers,
        colors_fill_ellipses=colors_fill_ellipses,
        color_edge_markers=color_edge_markers,
        color_edge_ellipses=color_edge_ellipses,
        color_emphasis=color_emphasis,
        set_axis_limits=False,
        show_confidence_ellipse=show_confidence_ellipse,
        show_lines_origin=True,
        show_line_diagonal=False, # diagonal is not proportional to respective ranges of axes
        show_emphasis_marker=show_emphasis_marker,
        show_emphasis_label=show_emphasis_label,
        show_legend_bar=show_legend_bar,
        report=None,
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
# Currently, this module is not directly executable.

################################################################################
# End
