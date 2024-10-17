"""
Supply basic plotting functionality.

This module 'plot_cluster' is part of the 'partner' package.

This module is not directly executable.

This subpackage 'partner' provides executable functionality under the
management of a higher level package. Importation paths require this hierarchy.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This module file is part of the project package directory 'partner'
    (https://github.com/tcameronwaller/partner/).

    Project 'partner' supports data analysis in multiple other projects.
    Copyright (C) 2024 Thomas Cameron Waller

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

###############################################################################
# Notes



###############################################################################
# Installation and importation

# Standard
import os
import pickle
import itertools
import statistics
import copy
import math

# Relevant
import pandas
import numpy
import scipy
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
import partner.organization as porg
import partner.plot as pplot

#dir()
#importlib.reload()


###############################################################################
# Functionality



##########
# Heatmap plots

# Label a small set of features.
# Cluster observations within categorical groups.
# Label categorical groups either on axis or by color key in legend.

# Arguments:
# table
# column_observation
# column_group
# show_legend_groups
# label_axis_groups



# Determine ColorBrewer colors for categorical groups.
# Create color bar for categorical groups.
# This color bar will have the same count of discrete values as there are
# observations across rows in the original table for representation in cells of
# the main heat map chart.
# Option 1
# Create the color bar as a heat map (matplotlib imshow) with discrete mappings
# of values to colors. This option guarantees that the colors have the same
# dimensions as the corresponding groups of cells in the main heat map. Create
# a categorical legend to interpret the colors and groups.
# Option 2
# Create the color bar as a color bar with custom bin thresholds to align with
# the corresponding groups of cells in the main heat map.


# TODO: TCW; 17 October 2024
# Determine if the "tight_layout" is a problem with gridspec...

# TODO: TCW; 17 October 2024
# Create legend for discrete representation of groups
# Need to map from the integers back to the categories

def plot_heatmap_signal_label_features_groups_of_observations(
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
    value_maximum=None,
    labels_ordinate_categories=None,
    labels_abscissa_categories=None,
    show_scale_bar=None,
    title_ordinate=None,
    title_abscissa=None,
    title_bar=None,
    size_title_ordinate=None,
    size_title_abscissa=None,
    size_label_ordinate=None,
    size_label_abscissa=None,
    size_title_bar=None,
    size_label_bar=None,
    aspect=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Heat map.

    Format of source table

    Format of source table is in wide format with floating-point values of
    signal intensities corresponding to features across columns and distinct
    observations across rows. A special column gives identifiers corresponding
    to each observation across rows. Another special column provides names
    of categorical groups of observations. For convenience in separating the
    values from their labels or identifiers, the source table must have
    explicitly defined indices across columns and rows.
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



    https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
        - categorical color bars
        - helpful reference for figuring out how to set the thresholds between
          the discrete levels of the color bar to correspond properly to the
          columns of the main heatmap...
    https://stackoverflow.com/questions/9707676/defining-a-discrete-colormap-for-imshow
    https://stackoverflow.com/questions/7229971/2d-grid-data-visualization-in-python
    https://matplotlib.org/stable/users/explain/colors/colormapnorms.html
        - normalization of information for color maps???

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            for features in columns across sample observations or groups of
            sample observations in rows
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
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
        labels_ordinate_categories (list<str>): optional, explicit labels for
            ordinate or vertical axis
        labels_abscissa_categories (list<str>): optional, explicit labels for
            abscissa or horizontal axis

        ...
        ...
        ...

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
    pail = pplot.extract_prepare_table_signals_categories_for_heatmap(
        table=table,
        format_table=format_table, # 1: features in rows, observations in columns
        index_columns=index_columns,
        index_rows=index_rows,
        column_group=column_group,
        transpose_table=transpose_table,
        fill_missing=fill_missing,
        value_missing_fill=value_missing_fill,
        constrain_signal_values=constrain_signal_values,
        value_minimum=value_minimum,
        value_maximum=value_maximum,
        labels_ordinate_categories=labels_ordinate_categories,
        labels_abscissa_categories=labels_abscissa_categories,
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
        nrows=2,
        ncols=1,
        wspace=0.0, # horizontal width space between grid blocks for subplots
        hspace=0.0, # vertical height space between grid blocks for subplots
        width_ratios=None, # all columns have same width
        height_ratios=(1,15), # first row 1/10th height of second row
    )
    # Initialize axes within grid within figure.
    axes_group = figure.add_subplot(grid[0,0]) # first row, first column
    axes_main = figure.add_subplot(grid[1,0]) # second row, first column

    # Keep axes, ticks, and labels, but remove border.
    for position in ['right', 'top', 'bottom', 'left']:
        matplotlib.pyplot.gca().spines[position].set_visible(False)

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
    image = axes_main.imshow(
        pail["matrix_signal"],
        cmap=matplotlib.colormaps["PuOr"], # RdBu_r, PuOr_r
        vmin=pail["value_minimum"],
        vmax=pail["value_maximum"],
        aspect="auto", # "auto", "equal",
        origin="upper",
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
    if (len(title_abscissa) > 0):
        axes_main.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
    # Set tick parameters for axes.
    axes_main.tick_params(
        axis="y", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=7.5, # 5.0
        width=5, # 3.5
        pad=10, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=False,
        bottom=False,
        left=True,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=True,
        labelright=False,
    )
    # Set tick positions and labels on axes.
    axes_main.set_yticks(
        numpy.arange(pail["matrix_signal"].shape[0]),
    )
    axes_main.set_yticklabels(
        pail["labels_ordinate_categories"],
        #minor=False,
        ha="right", # horizontal alignment
        va="center", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"][size_label_ordinate]
    )
    axes_main.tick_params(
        top=False,
        labeltop=False,
        bottom=False,
        labelbottom=False,
    )

    ##########
    # axes: group

    # Define color map for discrete, integer representation of categorical
    # groups.
    #discrete_minimum = 0
    #discrete_maximum = len(pail["labels_group_unique"])
    discrete_minimum = numpy.nanmin(pail["matrix_group_integers"])
    discrete_maximum = numpy.nanmax(pail["matrix_group_integers"])
    color_map = matplotlib.pyplot.get_cmap(
        "tab10",
        ((discrete_maximum - discrete_minimum) + 1)
    )
    # Plot values as a grid of color on discrete scale.
    image = axes_group.imshow(
        pail["matrix_group_integers"],
        cmap=color_map,
        vmin=discrete_minimum,
        vmax=discrete_maximum,
        aspect="auto", # "auto", "equal",
        origin="upper",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
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
    ##########
    # Plot data on figure.





    # https://gist.github.com/Prukutu/7bfd668d792f76b8f99f7a5c9fc877a4
    # - maybe helpful example of creating a color bar in a grid spec subplot
    # axis


    # TODO: TCW; 11 October 2024; I need to figure out how to set the
    # thresholds between the discrete levels of the color bar to align properly
    # with the columns of the main signal heatmap...
    # http://www.acgeospatial.co.uk/colour-bar-for-discrete-rasters-with-matplotlib/

    ##########
    # Return figure.
    return figure




###############################################################################
# Procedure
# This module is not executable.
