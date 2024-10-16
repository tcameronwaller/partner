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




def prepare_data_label_features_groups_cluster_observations(
    table=None,
    column_observation=None,
    column_group=None,
    report=None,
):
    """
    Prepare source information for plotting a heat map.

    Format of source table is in wide format with floating-point values of
    signal intensities corresponding to features across columns and distinct
    observations across rows. A special column provides names or values of
    categorical groups of observations. The table does not need to have an
    explicit index across rows.

    observation   group   feature_1 feature_2 feature_3 feature_4 feature_5
    observation_1 group_1 0.001     0.001     0.001     0.001     0.001
    observation_2 group_1 0.001     0.001     0.001     0.001     0.001
    observation_3 group_2 0.001     0.001     0.001     0.001     0.001
    observation_4 group_2 0.001     0.001     0.001     0.001     0.001
    observation_5 group_3 0.001     0.001     0.001     0.001     0.001

    In the final chart, this function preserves the original sequence of
    features across columns. This function also preserves the original sequence
    of groups of observations; however, this function clusters and changes the
    sequence of observations within groups.

    MatPlotLib color maps.
    https://matplotlib.org/stable/tutorials/colors/colormaps.html

    arguments:
        table (object): Pandas data-frame table of values of signal intensity
            for features in rows across sample observations or groups of
            sample observations in columns
        column_observation (str): name of column in table for identifiers that
            correspond to observations in rows with values of signal intensity
            across features in other columns
        column_group (str): name of column in table for categorical names that
            correspond to groups of observations in rows with values of signal
            intensity across features in other columns
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information for plot

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Extract labels for features from names of columns in table.
    columns_all = copy.deepcopy(table.columns.to_list())
    features = list(filter(
        lambda column: (column not in [column_observation, column_group]),
        columns_all
    ))
    # Extract labels for categorical groups of observations from values in the
    # special column for these categorical groups.
    groups = copy.deepcopy(table[column_group].unique().tolist())
    # Cluster rows in table within groups.
    table_cluster = porg.cluster_table_rows_by_group(
        table=table,
        group=column_group,
        index=column_observation,
    )
    # Collect information.
    pail = dict()
    pail["table_cluster"] = table_cluster
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.organization.py")
        function = "prepare_data_label_features_groups_cluster_observations"
        print(str("function: " + function + "()"))
        print("version check: 7")
        putly.print_terminal_partition(level=5)
        print(pail["table_cluster"])
        putly.print_terminal_partition(level=4)
    # Return information.
    return pail

# TODO: don't require the plotting function to perform the clustering
# run the clustering BEFORE passing the table to the plot function




def plot_heat_map_label_features_groups_cluster_observations(
    table=None,
    transpose_table=None,
    index_group_columns=None,
    index_group_rows=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_maximum=None,
    show_scale_bar=None,
    labels_ordinate_categories=None,
    labels_abscissa_categories=None,
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

    Format of source table is in wide format with floating-point values of
    signal intensities corresponding to features across columns and distinct
    observations across rows. A special column provides names or values of
    categorical groups of observations. The table does not need to have an
    explicit index across rows.

    observation   group   feature_1 feature_2 feature_3 feature_4 feature_5
    observation_1 group_1 0.001     0.001     0.001     0.001     0.001
    observation_2 group_1 0.001     0.001     0.001     0.001     0.001
    observation_3 group_2 0.001     0.001     0.001     0.001     0.001
    observation_4 group_2 0.001     0.001     0.001     0.001     0.001
    observation_5 group_3 0.001     0.001     0.001     0.001     0.001

    In the final chart, this function preserves the original sequence of
    features across columns. This function also preserves the original sequence
    of groups of observations; however, this function clusters and changes the
    sequence of observations within groups.

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
    https://matplotlib.org/stable/gallery/subplots_axes_and_figures/subplots_demo.html
       - subplots and gridspec
    https://stackoverflow.com/questions/42562440/multiple-heatmaps-with-fixed-grid-size
    http://www.acgeospatial.co.uk/colour-bar-for-discrete-rasters-with-matplotlib/
        - discrete color maps
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

        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): MatPlotLib figure object

    """

    ##########
    # Prepare data.


    ##########
    # Create figure.


    ##########
    # Initialize figure.
    figure = initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Initialize grid within figure.
    #matplotlib.gridspec.GridSpec()
    grid = figure.add_gridspec(
        nrows=2,
        ncols=1,
        hspace=0.01,
        wspace=0,
        height_ratios=(1,10),
        width_ratios=None, # all columns same width
        sharex=True, # sharex="col",

    )
    # Initialize axes within grid within figure.
    axes_group = figure.add_subplot(grid[0,0]) # first row, first column
    axes_main = figure.add_subplot(grid[1,0]) # second row, first column

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
    #return figure
    pass







###############################################################################
# Procedure
# This module is not executable.
