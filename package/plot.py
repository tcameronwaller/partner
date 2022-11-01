"""
Supply basic plotting functionality.

This module is not directly executable.

This subpackage 'promiscuity' provides executable functionality under the
management of a higher level package. Importation paths represent this
hierarchy.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of Promiscuity
    (https://github.com/tcameronwaller/promiscuity/).

    Promiscuity supports data analysis in multiple other projects.
    Copyright (C) 2022 Thomas Cameron Waller

    Promiscuity is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    Promiscuity is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with Promiscuity. If not, see <http://www.gnu.org/licenses/>.
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
import promiscuity.utility as utility

#dir()
#importlib.reload()


###############################################################################
# Functionality


def define_font_properties():
    """
    Defines font properties.

    arguments:

    raises:

    returns:
        (dict<object>): references to definitions of font properties

    """

    # Define font values.
    values_1 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 1000, # '1000' is the maximal permissable font stretch
        "weight": 1000, # '1000' is the maximal permissable font weight
        "size": 70,
    }
    values_2 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 1000, # '1000' is the maximal permissable font stretch
        "weight": 1000, # '1000' is the maximal permissable font weight
        "size": 60,
    }
    values_3 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 1000, # '1000' is the maximal permissable font stretch
        "weight": 1000, # '1000' is the maximal permissable font weight
        "size": 50,
    }
    values_4 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 1000, # '1000' is the maximal permissable font stretch
        "weight": 1000, # '1000' is the maximal permissable font weight
        "size": 45,
    }
    values_5 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 950, # '1000' is the maximal permissable font stretch
        "weight": 950, # '1000' is the maximal permissable font weight
        "size": 40,
    }
    values_6 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 900,
        "weight": 900,
        "size": 30,
    }
    values_7 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 750,
        "weight": 750,
        "size": 20,
    }
    values_8 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 500,
        "weight": 500,
        "size": 17,
    }
    values_9 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 400,
        "weight": 400,
        "size": 15,
    }
    values_10 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 300,
        "weight": 300,
        "size": 13,
    }
    values_11 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 300,
        "weight": 300,
        "size": 10,
    }
    values_12 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 200,
        "weight": 200,
        "size": 7,
    }
    values_13 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 150,
        "weight": 150,
        "size": 5,
    }
    values_14 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 100,
        "weight": 100,
        "size": 3,
    }
    # Define font properties.
    properties_1 = matplotlib.font_manager.FontProperties(
        family=values_1["family"],
        style=values_1["style"],
        variant=values_1["variant"],
        stretch=values_1["stretch"],
        weight=values_1["weight"],
        size=values_1["size"]
    )
    properties_2 = matplotlib.font_manager.FontProperties(
        family=values_2["family"],
        style=values_2["style"],
        variant=values_2["variant"],
        stretch=values_2["stretch"],
        weight=values_2["weight"],
        size=values_2["size"]
    )
    properties_3 = matplotlib.font_manager.FontProperties(
        family=values_3["family"],
        style=values_3["style"],
        variant=values_3["variant"],
        stretch=values_3["stretch"],
        weight=values_3["weight"],
        size=values_3["size"]
    )
    properties_4 = matplotlib.font_manager.FontProperties(
        family=values_4["family"],
        style=values_4["style"],
        variant=values_4["variant"],
        stretch=values_4["stretch"],
        weight=values_4["weight"],
        size=values_4["size"]
    )
    properties_5 = matplotlib.font_manager.FontProperties(
        family=values_5["family"],
        style=values_5["style"],
        variant=values_5["variant"],
        stretch=values_5["stretch"],
        weight=values_5["weight"],
        size=values_5["size"]
    )
    properties_6 = matplotlib.font_manager.FontProperties(
        family=values_6["family"],
        style=values_6["style"],
        variant=values_6["variant"],
        stretch=values_6["stretch"],
        weight=values_6["weight"],
        size=values_6["size"]
    )
    properties_7 = matplotlib.font_manager.FontProperties(
        family=values_7["family"],
        style=values_7["style"],
        variant=values_7["variant"],
        stretch=values_7["stretch"],
        weight=values_7["weight"],
        size=values_7["size"]
    )
    properties_8 = matplotlib.font_manager.FontProperties(
        family=values_8["family"],
        style=values_8["style"],
        variant=values_8["variant"],
        stretch=values_8["stretch"],
        weight=values_8["weight"],
        size=values_8["size"]
    )
    properties_9 = matplotlib.font_manager.FontProperties(
        family=values_9["family"],
        style=values_9["style"],
        variant=values_9["variant"],
        stretch=values_9["stretch"],
        weight=values_9["weight"],
        size=values_9["size"]
    )
    properties_10 = matplotlib.font_manager.FontProperties(
        family=values_10["family"],
        style=values_10["style"],
        variant=values_10["variant"],
        stretch=values_10["stretch"],
        weight=values_10["weight"],
        size=values_10["size"]
    )
    properties_11 = matplotlib.font_manager.FontProperties(
        family=values_11["family"],
        style=values_11["style"],
        variant=values_11["variant"],
        stretch=values_11["stretch"],
        weight=values_11["weight"],
        size=values_11["size"]
    )
    properties_12 = matplotlib.font_manager.FontProperties(
        family=values_12["family"],
        style=values_12["style"],
        variant=values_12["variant"],
        stretch=values_12["stretch"],
        weight=values_12["weight"],
        size=values_12["size"]
    )
    properties_13 = matplotlib.font_manager.FontProperties(
        family=values_13["family"],
        style=values_13["style"],
        variant=values_13["variant"],
        stretch=values_13["stretch"],
        weight=values_13["weight"],
        size=values_13["size"]
    )
    properties_14 = matplotlib.font_manager.FontProperties(
        family=values_14["family"],
        style=values_14["style"],
        variant=values_14["variant"],
        stretch=values_14["stretch"],
        weight=values_14["weight"],
        size=values_14["size"]
    )
    # Compile and return references.
    return {
        "values": {
            "one": values_1,
            "two": values_2,
            "three": values_3,
            "four": values_4,
            "five": values_5,
            "six": values_6,
            "seven": values_7,
            "eight": values_8,
            "nine": values_9,
            "ten": values_10,
            "eleven": values_11,
            "twelve": values_12,
            "thirteen": values_13,
            "fourteen": values_14,
        },
        "properties": {
            "one": properties_1,
            "two": properties_2,
            "three": properties_3,
            "four": properties_4,
            "five": properties_5,
            "six": properties_6,
            "seven": properties_7,
            "eight": properties_8,
            "nine": properties_9,
            "ten": properties_10,
            "eleven": properties_11,
            "twelve": properties_12,
            "thirteen": properties_13,
            "fourteen": properties_14,
        }
    }


def define_color_properties():
    """
    Defines color properties.

    https://www.w3schools.com/colors/colors_rgb.asp

    arguments:

    raises:

    returns:
        (dict<tuple>): references to definitions of color properties

    """
    # (0.1 * 255) = 25.5
    # (0.2 * 255) = 51.0
    # (0.3 * 255) = 76.5
    # (0.4 * 255) = 102.0
    # (0.5 * 255) = 127.5
    # (0.6 * 255) = 153.0
    # (0.7 * 255) = 178.5
    # (0.8 * 255) = 204.0
    # (0.9 * 255) = 229.5

    # Black.
    black = (0.0, 0.0, 0.0, 1.0)
    # Gray.
    gray = (0.5, 0.5, 0.5, 1.0)
    gray_dark = (0.3, 0.3, 0.3, 1.0) # (red: 75; green: 75; blue: 75)
    gray_light = (0.7, 0.7, 0.7, 1.0) # (red: 180; green: 180; blue: 180)
    # White.
    white = (1.0, 1.0, 1.0, 1.0)
    white_faint = (1.0, 1.0, 1.0, 0.75)
    # Clear.
    clear = (1.0, 1.0, 1.0, 0.0)
    clear_faint = (1.0, 1.0, 1.0, 0.25)

    blue_navy = (0.039, 0.196, 0.588, 1.0) # (red: 10; green: 50; blue: 150)
    blue_navy_faint = (0.039, 0.196, 0.588, 0.75)
    blue_navy_light = (0.118, 0.314, 0.706, 1.0) # (r: 30; g: 80; b: 180)
    blue_sky = (0.196, 0.588, 1.0, 1.0) # (red: 50; green: 150; blue: 255)
    blue_navy_bright = (0.784, 0.824, 1.0, 1.0) # (r: 200; g: 210; b: 255)

    purple = (0.510, 0.039, 0.510, 1.0) # (red: 130; green: 10; blue: 130)
    purple_light = (0.588, 0.196, 0.588, 1.0) # (red: 150; green: 50; blue: 150)
    magenta = (0.784, 0.275, 0.784, 1.0) # (red: 200; green: 70; blue: 200)

    red_brick = (0.667, 0.196, 0.039, 1.0) # (red: 170; green: 50; blue: 10)

    orange = (1.0, 0.588, 0.039, 1.0) # (red: 255; green: 150; blue: 10)
    orange_faint = (1.0, 0.588, 0.039, 0.75)

    # Green (red: 64; green: 191; blue: 64).
    green = (0.25, 0.75, 0.25, 1.0) # a brighter green, similar to "Malachite"
    # Yellow.
    yellow = (1.0, 0.8, 0.0, 1.0)
    # Compile and return references.
    return {
        "black": black,
        "gray": gray,
        "gray_dark": gray_dark,
        "gray_light": gray_light,
        "white": white,
        "white_faint": white_faint,
        "clear": clear,
        "clear_faint": clear_faint,
        "blue_navy": blue_navy,
        "blue_navy_faint": blue_navy_faint,
        "blue_navy_light": blue_navy_light,
        "blue_sky": blue_sky,
        "purple": purple,
        "magenta": magenta,
        "green": green,
        "orange": orange,
        "orange_faint": orange_faint,
        "yellow": yellow,
    }


# TODO: pass variable for label on scale bar...

def plot_heatmap_symmetric_diverge(
    data=None,
    minimum=None,
    maximum=None,
    label_rows=None,
    label_columns=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type heatmap.

    arguments:
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap
        minimum (float): minimal value
        maximum (float): maximal value
        label_rows (bool): whether to include explicit labels on heatmap's rows
        label_columns (bool): whether to include explicit labels on heatmap's
            columns
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    # Map data's columns to heatmap's rows.
    # Map data's rows to heatmap's columns.
    labels_rows = data.columns.to_list()
    labels_columns = data.index.to_list()
    matrix = numpy.transpose(data.to_numpy())

    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    axes = matplotlib.pyplot.axes()

    # Represent data as a color grid.
    image = axes.imshow(
        matrix,
        cmap="PuOr", # "PuOr", "RdBu" diverging color map
        vmin=minimum,
        vmax=maximum,
        aspect="equal",
        origin="upper",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )

    # Create legend for color map.
    label_bar = "Spearman correlation (FDR <= 0.05)"
    bar = axes.figure.colorbar(
        image,
        orientation="vertical",
        ax=axes,
    )
    bar.ax.set_ylabel(
        label_bar,
        rotation=-90,
        va="bottom",
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["three"],
    )
    bar.ax.tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=2.0,
        color=colors["black"],
        pad=2,
        labelsize=fonts["values"]["five"]["size"],
        labelcolor=colors["black"],
    )

    # Create ticks and labels for each grid.
    # Let the horizontal axes labeling appear on top.
    axes.tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=1.0,
        color=colors["black"],
        pad=5,
        labelcolor=colors["black"],
        top=True,
        bottom=False, # False
        left=True,
        right=False,
        labeltop=True,
        labelbottom=False,
        labelleft=True,
        labelright=False,
    )

    # Create ticks and labels.
    if (label_columns and data.shape[1] <= 100):
        if (100 >= data.shape[1] and data.shape[1] >= 90):
            size_column = "ten"
        elif (90 > data.shape[1] and data.shape[1] >= 75):
            size_column = "nine"
        elif (75 > data.shape[1] and data.shape[1] >= 50):
            size_column = "eight"
        elif (50 > data.shape[1] and data.shape[1] >= 25):
            size_column = "seven"
        elif (25 > data.shape[1]):
            size_column = "six"
        axes.set_xticks(numpy.arange(matrix.shape[1]))
        axes.set_xticklabels(
            labels_columns,
            #minor=False,
            rotation=-60,
            rotation_mode="anchor",
            ha="right", # horizontal alignment
            va="bottom", # vertical alignment
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_column]
        )
    if (label_rows and data.shape[0] <= 100):
        if (100 >= data.shape[0] and data.shape[0] >= 90):
            size_row = "nine"
        elif (90 > data.shape[0] and data.shape[0] >= 75):
            size_row = "eight"
        elif (75 > data.shape[0] and data.shape[0] >= 50):
            size_row = "six"
        elif (50 > data.shape[0] and data.shape[0] >= 25):
            size_row = "five"
        elif (25 > data.shape[0]):
            size_row = "three"
        axes.set_yticks(numpy.arange(matrix.shape[0]))
        axes.set_yticklabels(
            labels_rows,
            #minor=False,
            ha="right", # horizontal alignment
            va="center", # vertical alignment
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_row]
        )
    if False:
        # Rotate the tick labels and set their alignment.
        matplotlib.pyplot.setp(
            axes.get_xticklabels(),
            rotation=-30,
            ha="right",
            rotation_mode="anchor"
        )
        pass
    # Return figure.
    return figure


# TODO: accommodate different variable types... categorical, continuous divergent, continuous, ordinal, etc...


def organize_heatmap_asymmetric_data(
    data=None,
):
    """
    Organizes data for chart.

    arguments:
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap

    raises:

    returns:
        (dict): collection of information for chart

    """

    # Copy data.
    data = data.copy(deep=True)
    # Map data's columns to heatmap's rows.
    # Map data's rows to heatmap's columns.
    labels_rows = data.columns.to_list()
    labels_columns = data.index.to_list()
    matrix = numpy.transpose(data.to_numpy())
    # Determine maximal and minimal values.
    array = matrix.flatten()
    minimum = round(numpy.nanmin(array), 2)
    maximum = round(numpy.nanmax(array), 2)
    # Collect information.
    bin = dict()
    bin["data"] = data
    bin["matrix"] = matrix
    bin["minimum"] = minimum
    bin["maximum"] = maximum
    bin["labels_rows"] = labels_rows
    bin["labels_columns"] = labels_columns
    # Return information.
    return bin


def plot_heatmap_asymmetric(
    data=None,
    title=None,
    label_scale=None,
    type=None,
    label_rows=None,
    label_columns=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type heatmap.

    arguments:
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap
        title (str): title for top left of figure
        label_scale (str): label for heatmap scale bar
        type (str): type of property, category, binary, ordinal, continuous, or
            continuous_divergent
        label_rows (bool): whether to include explicit labels on heatmap's rows
        label_columns (bool): whether to include explicit labels on heatmap's
            columns
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    bucket = organize_heatmap_asymmetric_data(
        data=data,
    )
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811), # 40 cm, 30 cm
        tight_layout=True,
    )
    figure.suptitle(
        title,
        x=0.01,
        y=0.99,
        ha="left",
        va="top",
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"]
    )
    axes = matplotlib.pyplot.axes()
     # Main chart in bottom panel.
    organize_heatmap_asymmetric_master_main_bottom(
        label_scale=label_scale,
        type=type,
        matrix=bucket["matrix"],
        minimum=bucket["minimum"],
        maximum=bucket["maximum"],
        labels_rows=bucket["labels_rows"],
        labels_columns=bucket["labels_columns"],
        fonts=fonts,
        colors=colors,
        axis_main=axes,
        axis_scale=axes,
        figure=figure,
     )

    return figure


# Asymmetric heatmap with master and main panels
# Sort by master variable and cluster by similarities
# TODO: would this set of functions be more appropriate for the utility module?

def organize_data_master_main_sort_cluster(
    type_master=None,
    sequence=None,
    data=None,
    index=None,
):
    """
    Organize sort and cluster of data.

    arguments:
        type_master (str): type of master variable, category, binary, ordinal,
            or continuous
        data (object): Pandas data frame of main and master feature variables
            across instances
        index (str): name of index that is common between data_master and
            data_main
        sequence (str): method to determine sequence of instances, sort by
            master variable's values or cluster by similarities across features


    raises:

    returns:
        (object): Pandas data frame of pan-tissue signals across genes and
            persons with property

    """

    if (sequence == "sort"):
        data.sort_values(
            by=["master"],
            axis="index",
            ascending=True,
            inplace=True,
        )
        data_cluster_columns = utility.cluster_data_columns(
            data=data,
        )
        data_cluster_columns.reset_index(
            level=None,
            inplace=True
        )
        data_cluster_columns.set_index(
            [index],
            append=False,
            drop=True,
            inplace=True
        )
        if (type_master in ["category", "binary", "ordinal"]):
            data_cluster_columns.reset_index(
                level=None,
                inplace=True
            )
            data_cluster_columns.set_index(
                [index, "master"],
                append=False,
                drop=True,
                inplace=True
            )
            data_cluster_rows = utility.cluster_data_rows_by_group(
                group="master",
                index=index,
                data=data_cluster_columns,
            )
            # Organize data.
            data_cluster_rows.reset_index(
                level=None,
                inplace=True
            )
            data_cluster_rows.set_index(
                [index],
                append=False,
                drop=True,
                inplace=True
            )
            data_sequence = data_cluster_rows.copy(deep=True)
        else:
            data_sequence = data_cluster_columns.copy(deep=True)
    elif (sequence == "cluster"):
        data_cluster_columns = utility.cluster_data_columns(
            data=data,
        )
        data_cluster_rows = utility.cluster_data_rows(
            data=data_cluster_columns,
        )
        data_sequence = data_cluster_rows.copy(deep=True)
    # Return information.
    return data_sequence


def organize_data_master_main(
    data_master=None,
    master=None,
    type_master=None,
    data_main=None,
    type_main=None,
    scale_unit_main=None,
    columns_main_scale_unit=None,
    fill_missing=None,
    index=None,
    sequence=None,
):
    """
    Organize information to sort multiple main data variables by a single
    master variable.

    Data have features across columns and instances across rows.

    Sequence of features depends on hierarchical cluster by their
    similarities across instances.
    Sequence of instances depends either on sort by values of master variable
    or on hierarchical cluster by their similarities across features.

    arguments:
        data_master (object): Pandas data frame including master variable
            across instances that match those of data_main
        master (str): name of feature in data_master that is master variable
        type_master (str): type of master variable, category, binary, ordinal,
            or continuous
        data_main (object): Pandas data frame of feature variables across
            instances that match those of data_master
        type_main (str): type of main variables, category, binary, ordinal,
            continuous, or continuous_divergent
        scale_unit_main (bool): whether to scale columns in data_main
        columns_main_scale_unit (list<str>): names of columns in data_main to
            scale to unit distribution between 0 and 1
        fill_missing (bool): whether to fill missing values with zero
        index (str): name of index that is common between data_master and
            data_main
        sequence (str): method to determine sequence of instances, sort by
            master variable's values or cluster by similarities across features

    raises:

    returns:
        (dict): information for chart

    """

    # Copy data.
    data_master = data_master.copy(deep=True)
    data_main = data_main.copy(deep=True)

    # Organize main variables.
    # Drop any rows or columns in main data with only missing values.
    data_main.dropna(
        axis="index",
        how="all",
        inplace=True,
    )
    data_main.dropna(
        axis="columns",
        how="all",
        inplace=True,
    )
    # Scale main variable values.
    if scale_unit_main:
        data_main_scale = data_main.apply(
            lambda column:
                sklearn.preprocessing.minmax_scale(
                    column.to_numpy(),
                    feature_range=(0, 1),
                    axis=0,
                    copy=True,
                ) if (column.name in columns_main_scale_unit) else column,
            axis="index",
        )
    else:
        data_main_scale = data_main
    # Replace missing values in main data with zero.
    if fill_missing:
        # Set any infinite values to missing.
        data_main_scale[data_main_scale == numpy.inf] = numpy.nan
        data_main_scale.fillna(
            value=0.0,
            #axis="columns",
            inplace=True,
        )
    else:
        data_main_scale.fillna(
            value=-1.0,
            #axis="columns",
            inplace=True,
        )
    # Organize master variable.
    if (type_master in ["category", "binary", "ordinal"]):
        data_master["master"], labels_categories_master = pandas.factorize(
            data_master[master],
            sort=True
        )
        labels_categories_master = list(labels_categories_master)
    elif type_master == "continuous":
        data_master["master"] = data_master[master]
        labels_categories_master = list()
    data_master = data_master.loc[
        :, data_master.columns.isin(["master", master])
    ]
    # Join master and main data.
    data_hybrid = data_main_scale.join(
        data_master,
        how="left",
        on=index
    )
    data_hybrid.drop(
        labels=[master],
        axis="columns",
        inplace=True
    )
    # Drop any rows in hybrid data with missing values in master.
    data_hybrid.dropna(
        axis="index",
        how="all",
        subset=["master"],
        inplace=True,
    )
    # Sort and cluster data.
    if sequence != "none":
        data_hybrid_sequence = organize_data_master_main_sort_cluster(
            type_master=type_master,
            sequence=sequence,
            data=data_hybrid,
            index=index,
        )
    else:
        data_hybrid_sequence = data_hybrid
    # Compile information.
    bin = dict()
    bin["labels_categories_master"] = labels_categories_master
    bin["data"] = data_hybrid_sequence
    # Return information
    return bin


def organize_heatmap_asymmetric_master_main_data(
    data=None,
):
    """
    Organizes data for chart.

    arguments:
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap

    raises:

    returns:
        (dict): collection of information for chart

    """

    # Copy data.
    data_master = data.copy(deep=True)
    data_main = data.copy(deep=True)
    # Groups.
    values_master = data_master["master"].to_list()
    values_master_unique = utility.collect_unique_elements(
        elements_original=values_master,
    )
    # Map data's columns to heatmap's rows.
    # Map data's rows to heatmap's columns.
    data_main.drop(
        labels=["master"],
        axis="columns",
        inplace=True
    )
    labels_rows_main = data_main.columns.to_list()
    labels_columns_main = data_main.index.to_list()
    matrix_main = numpy.transpose(data_main.to_numpy())
    # Determine maximal and minimal values.
    minimum_master = round(min(values_master), 2)
    maximum_master = round(max(values_master), 2)
    array_main = matrix_main.flatten()
    minimum_main = round(numpy.nanmin(array_main), 2)
    maximum_main = round(numpy.nanmax(array_main), 2)

    # Collect information.
    bin = dict()
    bin["data_master"] = data_master
    bin["values_master"] = values_master
    bin["values_master_unique"] = values_master_unique
    bin["minimum_master"] = minimum_master
    bin["maximum_master"] = maximum_master
    bin["data_main"] = data_main
    bin["matrix_main"] = matrix_main
    bin["minimum_main"] = minimum_main
    bin["maximum_main"] = maximum_main
    bin["labels_rows_main"] = labels_rows_main
    bin["labels_columns_main"] = labels_columns_main
    # Return information.
    return bin


def initialize_heatmap_asymmetric_master_main_figure_axes(
    title=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type heatmap.

    arguments:
        title (str): title for chart
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object, object): instances of figure and axes objects

    """

    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811), # 40 cm, 30 cm
        #tight_layout=True,
    )
    figure.suptitle(
        title,
        x=0.01,
        y=0.99,
        ha="left",
        va="top",
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"]
    )
    axes = figure.subplots(
        nrows=2,
        ncols=2,
        sharex=False,
        sharey=False,
        #squeeze=True,
        gridspec_kw=dict(
            hspace=0.05,
            wspace=0.05,
            height_ratios=[1, 10],
            width_ratios=[50, 1],
            left=0.11,
            right=0.94,
            top=0.94,
            bottom=0.05,
        ),
    )
    # Return information.
    return figure, axes


def organize_heatmap_asymmetric_master_main_top(
    title_subordinate=None,
    label=None,
    type=None,
    values=None,
    values_unique=None,
    labels_categories=None,
    minimum=None,
    maximum=None,
    fonts=None,
    colors=None,
    axes=None,
    figure=None,
):
    """
    Organizes top panel of figure.

    arguments:
        title_subordinate (str): subordinate title for top right of figure
        label (str): label for master heatmap
        type (str): type of property, category or continuity
        values (list): values of property
        values_unique (list): unique values of property
        labels_categories (list<str>): labels for categorical ticks
        minimum (float): minimal value of property
        maximum (float): maximal value of property
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        axes (object): instance axis object
        figure (object): instance figure object

    raises:

    returns:
        (object): figure object

    """

    # Define color map and labels.
    if ((type == "category") or (type == "binary") or (type == "ordinal")):
        #color_map = matplotlib.colors.ListedColormap(
        #    [colors["blue_navy"], colors["orange"]], 2
        #)
        color_map = matplotlib.pyplot.get_cmap("GnBu", len(values_unique))
        ticks = values_unique
        labels_ticks = labels_categories
    elif (type == "continuous"):
        color_map = "GnBu"
        ticks = [minimum, maximum]
        labels_ticks = [minimum, maximum]
    # Initialize chart.
    image = axes[0, 0].imshow(
        [values],
        cmap=color_map,
        vmin=minimum,
        vmax=maximum,
        aspect="auto",
        origin="upper",
        # Extent: (left, right, bottom, top)
        extent=(-0.5, (len(values) - 0.5), (1 + 0.5), -0.5),
    )
    axes[0, 0].set_yticks(numpy.arange(1))
    axes[0, 0].set_yticklabels(
        [str(label)],
        #minor=False,
        ha="right", # horizontal alignment
        va="center", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["four"]
    )
    axes[0, 0].set_xlabel(
        str(title_subordinate),
        rotation=0,
        ha="right",
        va="bottom",
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["three"],
    )
    #bar_group.ax.xaxis.set_label_position("top")
    axes[0, 0].xaxis.set_label_coords(1.0, 1.2)
    axes[0, 0].tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=3.0,
        pad=7,
        top=False,
        bottom=False,
        left=True,
        right=False,
        labeltop=False,
        labelbottom=False,
        labelleft=True,
        labelright=False,
        #color=colors["black"],
        #labelsize=fonts["values"]["four"]["size"],
        #labelcolor=colors["black"],
    )
    # Create legend for color map.
    bar = figure.colorbar(
        image,
        cax=axes[0, 1],
        ticks=ticks,
        orientation="vertical",
        use_gridspec=True,
    )
    bar.ax.set_yticklabels(
        labels_ticks,
        #minor=False,
        ha="left", # horizontal alignment
        va="bottom", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["four"]
    )
    bar.ax.tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=3.0,
        pad=7,
        top=False,
        bottom=False,
        left=False,
        right=True,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=True,
        #color=colors["black"],
        #labelsize=fonts["values"]["four"]["size"],
        #labelcolor=colors["black"],
    )
    pass


def organize_heatmap_asymmetric_master_main_bottom(
    label_scale=None,
    type=None,
    matrix=None,
    minimum=None,
    maximum=None,
    labels_rows=None,
    labels_columns=None,
    fonts=None,
    colors=None,
    axis_main=None,
    axis_scale=None,
    figure=None,
):
    """
    Organizes top panel of figure.

    arguments:
        label_scale (str): label for heatmap scale bar
        type (str): type of property, category, binary, ordinal, continuous, or
            continuous_divergent
        matrix (array<array>): values of properties
        minimum (float): minimal value of property
        maximum (float): maximal value of property
        labels_rows (list<str>): labels for matrix rows
        labels_columns (list<str>): labels for matrix columns
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        axis_main (object): instance of axis for main plot
        axis_scale (object): instance of axis for scale bar plot
        figure (object): instance figure object

    raises:

    returns:
        (object): figure object

    """

    # Define color map and labels.
    if (
        (type == "category") or
        (type == "binary") or
        (type == "ordinal") or
        (type == "continuous")
    ):
        color_map = "RdPu"
        scale = matplotlib.colors.Normalize(
            vmin=minimum,
            vmax=maximum,
        )
        ticks = [minimum, maximum]
        labels_ticks = [minimum, maximum]
    elif (type == "continuous_divergent"):
        color_map = "PuOr"
        scale = matplotlib.colors.TwoSlopeNorm(
            vmin=minimum,
            vcenter=0.0,
            vmax=maximum,
        )
        ticks = [minimum, 0.0, maximum]
        labels_ticks = [minimum, 0.0, maximum]

    # Initialize chart.
    image = axis_main.imshow(
        matrix,
        cmap=color_map, # sequential: "RdPu", diverging: "PuOr"
        aspect="auto", # "auto", "equal"
        origin="upper", # "lower" or "upper"
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
        alpha=1.0,
        norm=scale,
        filternorm=True,
        resample=True,
    )
    #axes[2, 0].set_xlim(-10, (matrix.shape[1] + 10))
    #axes[2, 0].set_ylim(-3, (matrix.shape[0] + 3))

    # Create ticks and labels for each grid.
    # Let the horizontal axes labeling appear on top.
    if matrix.shape[0] <= 100:
        if (100 >= matrix.shape[0] and matrix.shape[0] >= 90):
            size_count = "ten"
        elif (90 > matrix.shape[0] and matrix.shape[0] >= 75):
            size_count = "nine"
        elif (75 > matrix.shape[0] and matrix.shape[0] >= 50):
            size_count = "eight"
        elif (50 > matrix.shape[0] and matrix.shape[0] >= 25):
            size_count = "seven"
        elif (25 > matrix.shape[0]):
            size_count = "three"
        axis_main.tick_params(
            axis="both",
            which="both",
            direction="out",
            length=5.0,
            width=3.0,
            pad=7,
            top=False,
            bottom=True,
            left=True,
            right=False,
            labeltop=False,
            labelbottom=True,
            labelleft=True,
            labelright=False,
            color=colors["black"],
            labelsize=fonts["values"][size_count]["size"],
            labelcolor=colors["black"],
        )
        axis_main.set_yticks(numpy.arange(matrix.shape[0]))
        axis_main.set_yticklabels(
            labels_rows,
            #minor=False,
            ha="right", # horizontal alignment
            va="center", # vertical alignment
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_count]
        )
        pass
    # Create legend for color map.
    if axis_scale is axis_main:
        bar = axis_main.figure.colorbar(
            image,
            orientation="vertical",
            ax=axis_main,
            ticks=ticks,
        )
    else:
        bar = figure.colorbar(
            image,
            cax=axis_scale,
            ticks=ticks,
            orientation="vertical",
            use_gridspec=True,
        )
        pass
    bar.ax.set_ylabel(
        label_scale,
        rotation=-90,
        ha="center",
        va="bottom",
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["three"],
    )
    bar.ax.yaxis.set_label_coords(2.5, 0.5)
    bar.ax.set_yticklabels(
        labels_ticks,
        #minor=False,
        ha="left", # horizontal alignment
        va="center", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["four"]
    )
    bar.ax.tick_params(
        axis="both",
        which="both", # major, minor, or both
        direction="out",
        length=5.0,
        width=3.0,
        pad=7,
        top=False,
        bottom=False,
        left=False,
        right=True,
        labeltop=False,
        labelbottom=False,
        labelleft=False,
        labelright=True,
        #color=colors["black"],
        #labelsize=fonts["values"]["four"]["size"],
        #labelcolor=colors["black"],
    )
    pass


def plot_heatmap_asymmetric_master_main_top_bottom(
    title=None,
    title_subordinate=None,
    label_master=None,
    labels_categories_master=None,
    label_main=None,
    type_master=None,
    type_main=None,
    data=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type heatmap.

    Data must have observations across rows.
    Data's observations must already be in sort order.
    Data must have features across columns.
    Data must also have a single column of name "master".
    Values in column "master" must be integers.

    The chart will organize features across rows and observations across
    columns.
    The chart will represent integer values of the group of each observation in
    a separate chart across columns.


    arguments:
        title (str): title for top left of figure
        title_subordinate (str): subordinate title for top right of figure
        label_master (str): label for left of master heatmap row
        labels_categories_master (list<str>): labels for scale ticks of
            categorical master variable
        label_main (str): label for main heatmap scale
        type_master (str): type of master variable, category, binary, ordinal,
            or continuous
        type_main (str): type of main variables, category, binary, ordinal,
            continuous, or continuous_divergent
        data (object): Pandas data frame of quantitative values with mappings
            to columns and rows that will be transposed in heatmap
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    bin_data = organize_heatmap_asymmetric_master_main_data(
        data=data,
    )
    # Initialize figure.
    figure, axes = initialize_heatmap_asymmetric_master_main_figure_axes(
        title=title,
        fonts=fonts,
        colors=colors,
    )

    # Master chart in top panel.
    organize_heatmap_asymmetric_master_main_top(
        title_subordinate=title_subordinate,
        label=label_master,
        type=type_master,
        values=bin_data["values_master"],
        values_unique=bin_data["values_master_unique"],
        labels_categories=labels_categories_master,
        minimum=bin_data["minimum_master"],
        maximum=bin_data["maximum_master"],
        fonts=fonts,
        colors=colors,
        axes=axes,
        figure=figure,
     )

     # Main chart in bottom panel.
    organize_heatmap_asymmetric_master_main_bottom(
        label_scale=label_main,
        type=type_main,
        matrix=bin_data["matrix_main"],
        minimum=bin_data["minimum_main"],
        maximum=bin_data["maximum_main"],
        labels_rows=bin_data["labels_rows_main"],
        labels_columns=bin_data["labels_columns_main"],
        fonts=fonts,
        colors=colors,
        axis_main=axes[1, 0],
        axis_scale=axes[1, 1],
        figure=figure,
     )

    # Return figure.
    return figure


def plot_distribution_histogram(
    array=None,
    title=None,
    bin_method=None,
    bin_count=None,
    bar_width=None,
    label_bins=None,
    label_counts=None,
    fonts=None,
    colors=None,
    line=None,
    line_position=None,
    label_title=None,
    label_report=None,
):
    """
    Creates a figure of a chart of type histogram to represent the frequency
    distribution of a single series of values.

    Consider utility of numpy.digitize(), numpy.bincount(), numpy.unique(),
    and numpy.histogram().

    arguments:
        array (object): NumPy array of values
        title (str): title for figure
        bin_method (str): method to define bins, "count" or "auto"
        bin_count (int): count of bins to define and populate
        bar_width (float): proportional of width of bar relative to bin
        label_bins (str): label for bins
        label_counts (str): label for counts
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        line (bool): whether to draw a vertical line
        line_position (float): position for vertical line
        label_title (str): text label title to include on figure
        label_report (bool): whether to include label for report of count and
            mean of values

    raises:

    returns:
        (object): figure object

    """

    # Determine count, mean, and median of values in array.
    count_values = int(array.size)
    mean_values = round(numpy.nanmean(array), 3)
    median_values = round(numpy.nanmedian(array), 3)

    # TODO: TCW; 11 October 2022
    # TODO: I do not remember what is the rationale to call
    # TODO: numpy.histogram_bin_edges() directly. ?

    # Define and populate bins.
    # Bin method "auto" is useful.
    #values, edges = numpy.histogram(series, bins=count_bins)
    if bin_method == "count":
        bin_edges = numpy.histogram_bin_edges(array, bins=bin_count)
    else:
        bin_edges = numpy.histogram_bin_edges(array, bins=bin_method)

    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    axes = matplotlib.pyplot.axes()
    values, bins, patches = axes.hist(
        array,
        bins=bin_edges,
        histtype="bar",
        align="left",
        orientation="vertical",
        rwidth=bar_width, # 0.35, 0.50
        log=False,
        color=colors["blue_navy"],
        label=title,
        stacked=False
    )
    if False:
        axes.set_title(
            title,
            #fontdict=fonts["properties"]["one"],
            loc="center",
        )
        axes.legend(
            loc="upper right",
            markerscale=2.5,
            markerfirst=True,
            prop=fonts["properties"]["one"],
            edgecolor=colors["black"]
        )
    axes.set_xlabel(
        xlabel=label_bins,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"]
    )
    axes.set_ylabel(
        ylabel=label_counts,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["two"]["size"],
        labelcolor=colors["black"]
    )
    if line:
        axes.axvline(
            x=line_position,
            ymin=0,
            ymax=1,
            alpha=1.0,
            color=colors["orange"],
            linestyle="--",
            linewidth=7.5, # 3.0, 5.0, 7.5
        )
    if len(label_title) > 0:
        matplotlib.pyplot.text(
            0.99,
            0.99,
            label_title,
            horizontalalignment="right",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )
    if label_report:
        matplotlib.pyplot.text(
            0.99,
            0.90,
            str("(count: " + str(count_values) + ")"),
            horizontalalignment="right",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )
        matplotlib.pyplot.text(
            0.99,
            0.85,
            str("(mean: " + str(mean_values) + ")"),
            horizontalalignment="right",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )
        matplotlib.pyplot.text(
            0.99,
            0.80,
            str("(median: " + str(median_values) + ")"),
            horizontalalignment="right",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )
    return figure


def plot_bar_stack(
    data=None,
    label_vertical=None,
    label_horizontal=None,
    fonts=None,
    colors=None,
    color_count=None,
    rotation=None,
    legend=None,
):
    """
    Creates a figure of a chart of type histogram to represent the frequency
    distribution of a single series of values.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        label_vertical (str): label for vertical axis
        label_horizontal (str): label for horizontal axis
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        color_count (int): count of colors for subcategories (stacks)
        rotation (str): value for rotation of labels
        legend (bool): whether to include a legend for series on the chart

    raises:

    returns:
        (object): figure object

    """

    # Define groups and series.
    groups = list(data.loc[:, "groups"])
    data.set_index(
        ["groups"],
        append=False,
        drop=True,
        inplace=True
    )
    series_names = list(data.columns)

    if color_count == 1:
        colors_series = [colors["blue_navy"]]
    elif color_count == 2:
        colors_series = [colors["blue_navy"], colors["orange"]]
    elif color_count > 2:
        colors_series = list(
            seaborn.color_palette("hls", n_colors=color_count)
        )
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()

    # Initialize bases.
    bases = list(itertools.repeat(0, len(groups)))
    # Initialize bars.
    bars = []
    # Iterate on series.
    index = 0
    for series_name in series_names:
        # Access series values.
        values = list(data.loc[:, series_name])
        print(series_name)
        print(values)
        # Create charts on axes.
        series_bars = axes.bar(
            range(len(groups)),
            values,
            width=0.8,
            bottom=bases,
            align="center",
            color=colors_series[index],
        )
        bars.append(series_bars[0])
        # Restore bases.
        bases = list(map(sum, zip(bases, values)))
        print(bases)
        index += 1
        pass

    axes.set_ylabel(
        ylabel=label_vertical,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"]
    )
    axes.set_xlabel(
        xlabel=label_horizontal,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["two"],
        rotation="horizontal",
    )
    axes.set_xticks(range(len(groups)), minor=False)
    axes.set_xticklabels(
        groups,
        rotation=-30,
        rotation_mode="anchor",
        ha="left", # horizontal alignment
        va="top", # vertical alignment
        minor=False,
        #rotation=rotation,
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["two"]["size"],
        labelcolor=colors["black"]
    )
    if legend:
        axes.legend(
            bars,
            series_names,
            loc="upper left",
            markerscale=2.5,
            markerfirst=True,
            prop=fonts["properties"]["one"],
            edgecolor=colors["black"]
        )
    return figure


# https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.violinplot.html#matplotlib.pyplot.violinplot


def plot_boxes_groups(
    values_groups=None,
    title_ordinate=None,
    title_abscissa=None,
    titles_abscissa_groups=None,
    colors_groups=None,
    label_top_center=None,
    label_top_left=None,
    label_top_right=None,
    aspect=None,
    orientation_box=None,
    axis_linear_minimum=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type box plot to represent the distribution
    of multiple series of values.

    https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.boxplot.html
    https://towardsdatascience.com/understanding-boxplots-5e2df7bcbd51

    arguments:
        values_groups (list<array>): NumPy arrays of non-missing values for each
            group
        title_ordinate (str): title for ordinate or vertical axis
        title_abscissa (str): title for abscissa or horizontal axis
        titles_abscissa_groups (list<str>): titles for groups on abscissa or
            horizontal axis in same sort order as 'cohorts_groups'
        colors_groups (list<tuple>): color parameters for boxes of groups in
            same sort order as 'values_groups'
        label_top_center (str): label for top center of plot area
        label_top_left (str): label for top left of plot area
        label_top_right (str): label for top right of plot area
        aspect (str): orientation and aspect ratio of figure: either 'portrait',
            'portrait_half_width', 'landscape', or 'landscape_half_height'
        orientation_box (str): whether the orientation of boxes is 'horizontal'
            or 'vertical'
        axis_linear_minimum (float): minimal value for range of linear axis
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    #colors_groups = list(seaborn.color_palette("hls", n_colors=color_count))

    # Create figure.
    if aspect == "portrait":
        figure = matplotlib.pyplot.figure(
            figsize=(11.811, 15.748), # aspect 3 X 4; 15.748 inches = 39.999 cm
            tight_layout=True
        )
    elif aspect == "portrait_half_width":
        figure = matplotlib.pyplot.figure(
            figsize=(5.906, 15.748), # aspect 1.5 X 4; 5.906 inches = 15.001 cm
            tight_layout=True
        )
    elif aspect == "landscape":
        figure = matplotlib.pyplot.figure(
            figsize=(15.748, 11.811), # aspect 4 X 3; 11.811 inches = 29.999 cm
            tight_layout=True
        )
    elif aspect == "landscape_half_height":
        figure = matplotlib.pyplot.figure(
            figsize=(15.748, 5.906), # aspect 4 X 1.5; 5.906 inches = 15.001 cm
            tight_layout=True
        )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    # Create boxes.
    if orientation_box == "horizontal":
        boxes_vertical = False
    elif orientation_box == "vertical":
        boxes_vertical = True
    handle = axes.boxplot(
        values_groups,
        notch=False, # whether to draw notch at center of box
        showfliers=False, # whether to show flier (outlier) points
        vert=boxes_vertical, # whether box groups across horizontal axis
        widths=0.7,
        patch_artist=True,
        labels=titles_abscissa_groups,
        manage_ticks=True,
        showmeans=True, # whether to show marker (point or line) for mean
        meanline=False, # whether to show line for mean
        boxprops={
            "linestyle": "solid",
            "linewidth": 1.0,
            "color": colors["black"],
        },
        medianprops={
            "linestyle": "solid",
            "linewidth": 5.0, # 1.0, 2.5, 5.0
            "color": colors["black"],
        },
        meanprops={
            "marker": "D", # diamond
            "markersize": 15.0, # 10.0, 20.0
            "markeredgecolor": colors["orange"],
            "markerfacecolor": colors["orange"],
        },
        whiskerprops={
            "linestyle": "solid",
            "linewidth": 3.5,
            "color": colors["black"],
        },
    )

    # Set minimum value of linear axis.
    if orientation_box == "horizontal":
        axes.set_xlim(
            xmin=axis_linear_minimum,
        )
    elif orientation_box == "vertical":
        axes.set_ylim(
            ymin=axis_linear_minimum,
        )

    # Fill boxes with colors.
    for box_patch, color_box in zip(handle["boxes"], colors_groups):
        box_patch.set_facecolor(color_box)
        pass
    # Label axes.
    if len(title_abscissa) > 0:
        axes.set_xlabel(
            xlabel=title_abscissa,
            labelpad=20,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"],
            rotation="horizontal",
        )
    if len(title_ordinate) > 0:
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=20,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )

    # Set tick parameters for axes.
    if orientation_box == "horizontal":
        # No label rotation necessary for horizontal orientation.
        axes.tick_params(
            axis="y", # "y", "x", or "both"
            which="both",
            direction="out",
            length=10.0, # 5.0, 10.0, 15.0
            width=7.5, # 3.0, 7.5, 11.0
            color=colors["black"],
            pad=10,
            labelsize=fonts["values"]["five"]["size"],
            labelcolor=colors["black"],
        )
        axes.tick_params(
            axis="x", # "y", "x", or "both"
            which="both",
            direction="out",
            length=10.0, # 5.0, 10.0, 15.0
            width=7.5, # 3.0, 7.5, 11.0
            color=colors["black"],
            pad=10,
            labelsize=fonts["values"]["five"]["size"],
            labelcolor=colors["black"],
        )
    elif orientation_box == "vertical":
        # Label rotation necessary for vertical orientation.
        axes.tick_params(
            axis="y", # "y", "x", or "both"
            which="both",
            direction="out",
            length=10.0, # 5.0, 10.0, 15.0
            width=7.5, # 3.0, 7.5, 11.0
            color=colors["black"],
            pad=10,
            labelsize=fonts["values"]["five"]["size"],
            labelcolor=colors["black"],
        )
        axes.tick_params(
            axis="x", # "y", "x", or "both"
            which="both",
            direction="out",
            length=10.0, # 5.0, 10.0, 15.0
            width=7.5, # 3.0, 7.5, 11.0
            color=colors["black"],
            pad=10,
            labelsize=fonts["values"]["five"]["size"],
            labelcolor=colors["black"],
            #labelrotation=45.0, # 45.0, 60,0
            #rotation_mode="anchor",
            #horizontalalignment="right", # not supported in current version
            #position=(-0.25, 0.0), # (-0.5, 0.0)
        )
        #axes.set_xticklabels(axes.get_xticklabels(), ha="right")
        #axes.set_xticklabels(axes.get_xticklabels(), rotation=45.0, ha="right")
        #axes.axis["bottom"].major_ticklabels.set_ha("right")

        # Access and modify the tick labels as text objects.
        for text_tick in axes.get_xticklabels():
            # set_horizontalalignment()
            # set_rotation()
            # set_rotation_mode()
            text_tick.set(
                rotation=45.0,
                rotation_mode="anchor",
                horizontalalignment="right",
                verticalalignment="top",
            )

    # Include label or labels on plot area.
    if len(label_top_center) > 0:
        matplotlib.pyplot.text(
            0.5,
            0.9,
            label_top_center,
            horizontalalignment="center",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["eight"]
        )
    if len(label_top_left) > 0:
        matplotlib.pyplot.text(
            0.1,
            0.9,
            label_top_left,
            horizontalalignment="left",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["eight"]
        )
    if len(label_top_right) > 0:
        matplotlib.pyplot.text(
            0.9,
            0.9,
            label_top_right,
            horizontalalignment="right",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["eight"]
        )
    return figure


def plot_scatter_factor_groups(
    data=None,
    abscissa=None,
    ordinate=None,
    label_horizontal=None,
    label_vertical=None,
    factor=None,
    label_factor=None,
    labels_factors=None,
    fonts=None,
    colors=None,
    point_size=None,
    plot_factor_labels=None,
    legend=None,
):
    """
    Creates a figure of a chart of type scatter to represent the association
    of two variables.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with variable for horizontal (x)
            axis
        ordinate (str): name of data column with variable for vertical (y) axis
        label_horizontal (str): label for horizontal axis
        label_vertical (str): label for vertical axis
        factor (str): name of data column with categorical variable to
            distinguish groups of instances
        label_factor (str): label to describe the factor
        labels_factors (list<str>): labels to describe the factor's values
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        point_size (float): size for scatter points
        plot_factor_labels (bool): whether to plot factor labels on chart
        legend (bool): whether to include a legend for series on the chart

    raises:

    returns:
        (object): figure object

    """

    ##########
    # Organize data.
    data = data.copy(deep=True)
    data.reset_index(
        level=None,
        inplace=True
    )
    data.set_index(
        factor,
        append=False,
        drop=True,
        inplace=True
    )
    groups = data.groupby(level=[factor])
    colors_series = list(seaborn.color_palette("hls", n_colors=len(groups)))

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    axes.set_xlabel(
        xlabel=label_horizontal,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel=label_vertical,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["one"]["size"],
        labelcolor=colors["black"]
    )
    # Plot points for values from each group.
    index = 0
    for name, group in groups:
        values_x = group[abscissa].to_list()
        values_y = group[ordinate].to_list()
        handle = axes.plot(
            values_x,
            values_y,
            linestyle="",
            marker="o",
            markersize=point_size,
            markeredgecolor=colors_series[index],
            markerfacecolor=colors_series[index]
        )
        index += 1
        pass
    # Plot labels for each group.
    labels = []
    index = 0
    for name, group in groups:
        if plot_factor_labels:
            values_x = group[abscissa].to_list()
            mean_x = statistics.mean(values_x)
            values_y = group[ordinate].to_list()
            mean_y = statistics.mean(values_y)
            axes.text(
                mean_x,
                mean_y,
                str(index+1),
                backgroundcolor=colors["white_faint"],
                color=colors["black"],
                fontproperties=fonts["properties"]["three"],
                horizontalalignment="center",
                verticalalignment="center"
            )
        label = str(str(index+1) + ": " + str(labels_factors[index]))
        labels.append(label)
        index += 1
        pass
    # Create legend.
    # Create custome elements for the legend.
    if legend:
        elements = create_legend_elements(
            colors=colors_series,
            labels=labels
        )
        axes.legend(
            handles=elements,
            loc="upper right",
            prop=fonts["properties"]["four"],
            title=label_factor,
            title_fontsize=fonts["values"]["three"]["size"]
        )
    return figure


def plot_scatter_points_discrete_abscissa_ordinate_error_bars(
    table=None,
    abscissa=None,
    ordinate=None,
    ordinate_error_low=None,
    ordinate_error_high=None,
    title_abscissa=None,
    title_ordinate=None,
    fonts=None,
    colors=None,
    size_marker=None,
    label_title=None,
):
    """
    Creates a figure of a chart of type scatter.

    This plot's abscissa (horizontal, x axis) fits best for a variable with
    discrete values.
    This plot's ordinate (vertical, y axis) fits best for a variable with
    continuous values.

    arguments:
        table (object): Pandas data frame of feature variables across columns
            and observation records across rows
        abscissa (str): name of table's column with variable for horizontal (x)
            axis
        ordinate (str): name of table's column with variable for vertical (y)
            axis
        ordinate_error_low (str): name of table's column for the extent of the
            bottom error bar below the dot
        ordinate_error_high (str): name of table's column for the extent of the
            top error bar above the dot
        title_abscissa (str): title for abscissa on horizontal axis
        title_ordinate (str): title for ordinate on vertical axis
        factor (str): name of data column with groups or factors of samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        size_marker (int): size of marker
        label_title (str): text label title to include on figure

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    table = table.copy(deep=True)
    columns = [abscissa, ordinate, ordinate_error_low, ordinate_error_high]
    table_selection = table.loc[
        :, table.columns.isin(columns)
    ]
    table_selection.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    values_abscissa = table_selection[abscissa].to_numpy()
    values_ordinate = table_selection[ordinate].to_numpy()
    errors_ordinate_low = table_selection[ordinate_error_low].to_list()
    errors_ordinate_high = table_selection[ordinate_error_high].to_list()
    # Shape (n, 2)
    #errors_ordinate = numpy.array(list(zip(
    #    errors_ordinate_low, errors_ordinate_high
    #)))
    # Shape (2, n)
    errors_ordinate = numpy.array([errors_ordinate_low, errors_ordinate_high])

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    #axes.set_ylim(ymin=0)
    #axes.set_xlim(
    #    xmin=abscissa_minimum,
    #    xmax=abscissa_maximum,
    #)
    # Set titles for axes.
    if (len(title_abscissa) > 0):
        axes.set_xlabel(
            xlabel=title_abscissa,
            labelpad=20,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )
    if (len(title_ordinate) > 0):
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=20,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )
    # Set tick parameters for axes.
    axes.tick_params(
        axis="y", # "y", "x", or "both"
        which="both",
        direction="out",
        length=10.0, # 5.0
        width=7.0, # 3.0
        color=colors["black"],
        pad=10,
        labelsize=fonts["values"]["three"]["size"],
        labelcolor=colors["black"]
    )
    axes.tick_params(
        axis="x", # "y", "x", or "both"
        which="both",
        direction="out",
        length=10.0, # 5.0
        width=7.0, # 3.0
        color=colors["black"],
        pad=10,
        labelsize=fonts["values"]["three"]["size"],
        labelcolor=colors["black"]
    )
    # Plot points for values from each group.
    handle = axes.errorbar(
        values_abscissa,
        values_ordinate,
        yerr=errors_ordinate,
        xerr=None,
        ecolor=colors["gray_dark"],
        elinewidth=7.0, # 10.0 is a little too thick with 30 points
        barsabove=False, # whether to print error bars in layer above points
        linestyle="",
        marker="o",
        markersize=size_marker, # 15, 25
        markeredgecolor=colors["blue_navy"],
        markerfacecolor=colors["blue_navy"],
    )
    # Include title label on plot.
    if (len(label_title) > 0):
        matplotlib.pyplot.text(
            0.99,
            0.99,
            label_title,
            horizontalalignment="right",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["eight"]
        )
    # Return figure.
    return figure


def plot_scatter(
    data=None,
    abscissa=None,
    ordinate=None,
    title_abscissa=None,
    title_ordinate=None,
    fonts=None,
    colors=None,
    size=None,
):
    """
    Creates a figure of a chart of type scatter.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with variable for horizontal (x)
            axis
        ordinate (str): name of data column with variable for vertical (y) axis
        title_abscissa (str): title for abscissa on horizontal axis
        title_ordinate (str): title for ordinate on vertical axis
        factor (str): name of data column with groups or factors of samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        size (int): size of marker

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    data_copy = data.copy(deep=True)
    data_selection = data_copy.loc[:, [abscissa, ordinate]]
    data_selection.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    values_abscissa = data_selection[abscissa].to_numpy()
    values_ordinate = data_selection[ordinate].to_nympy()

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    axes.set_xlabel(
        xlabel=title_abscissa,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel=title_ordinate,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["one"]["size"],
        labelcolor=colors["black"]
    )
    # Plot points for values from each group.
    handle = axes.plot(
        values_abscissa,
        values_ordinate,
        linestyle="",
        marker="o",
        markersize=size, # 5, 15
        markeredgecolor=colors["blue_navy"],
        markerfacecolor=colors["blue_navy"]
    )

    # Return figure.
    return figure


def plot_scatter_points_forest_category_ordinate_two_groups(
    table=None,
    column_group=None,
    column_ordinate_label=None,
    column_ordinate_sort=None,
    column_abscissa_value=None,
    column_abscissa_interval_below=None,
    column_abscissa_interval_above=None,
    group_one=None,
    group_two=None,
    abscissa_minimum=None,
    abscissa_maximum=None,
    ordinate_title=None,
    abscissa_title=None,
    label_chart=None,
    label_size_ordinate_categories=None,
    label_size_abscissa_values=None,
    size_marker_one=None,
    size_marker_two=None,
    color_marker_one=None,
    color_marker_two=None,
    space_groups=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type scatter with discrete, categorical,
    textual labels on the vertical ordinate axis and continuous points with
    error bars on the horizontal abscissa axis.

    Common name of this chart design is "Forest Plot".

    This function accommodates exactly two groups of values (series).

    This plot's ordinate (vertical, y axis) fits best for a variable with
    discrete categorical values that serve as names.

    This plot's abscissa (horizontal, x axis) fits best for a variable with
    continuous, ratio-scale values that center on zero.

    MatPlotLib accepts intervals, not ranges for error bars. The function does
    the arithmetic to calculate ranges below and above the central value.

    arguments:
        table (object): Pandas data frame of feature variables across columns
            and observation records across rows
        column_group (str): name of table's column with names of groups or
            series of values
        column_ordinate_label (str): name of table's column with textual
            categorical values for representation as labels on the vertical
            ordinate (y) axis
        column_ordinate_sort (str): name of table's column with integer values
            for sort order on categorical values for ordinate
        column_abscissa_value (str): name of table's column with floating point
            continuous scale coefficient values for representation as points on
            the horizontal abscissa (x) axis
        column_abscissa_interval_below (str): name of table's column for the
            extent (interval) of the error bar below the point value
        column_abscissa_interval_above (str): name of table's column for the
            extent (interval) of the error bar above the point value
        group_one (str): textual categorical name of first group of values
            (series) within table's column with name 'column_group'
        group_two (str): textual categorical name of second group of values
            (series) within table's column with name 'group'
        abscissa_minimum (float): minimal value for range of abscissa axis
        abscissa_maximum (float): maximal value for range of abscissa axis
        ordinate_title (str): title for ordinate or vertical axis
        abscissa_title (str): title for abscissa or horizontal axis
        label_chart (str): text label title of chart to include on figure
        label_size_ordinate_categories (str): label size for categories on
            vertical axis
        label_size_abscissa_values (str): label size for values on
            horizontal axis
        size_marker_one (int): size of markers for group one
        size_marker_two (int): size of markers for group two
        color_marker_one (str): name of color for markers of group one
        color_marker_two (str): name of color for markers of group two
        space_groups (float): vertical spacing between markers for groups
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize information from table.
    # It is important that the sorts orders of labels and values are identical
    # for series one and series two.
    table = table.copy(deep=True)
    columns = [
        column_group, column_ordinate_label, column_ordinate_sort,
        column_abscissa_value,
        column_abscissa_interval_below, column_abscissa_interval_above,
    ]
    table_columns = table.loc[
        :, table.columns.isin(columns)
    ]
    table_columns.sort_values(
        by=[column_group, column_ordinate_sort],
        axis="index",
        ascending=False, # Non-ascending sort order is more intuitive here
        inplace=True,
    )
    # Accommodate missing values.
    table_columns[column_abscissa_value].fillna(
        value=-3,
        inplace=True,
    )
    table_columns[column_abscissa_interval_below].fillna(
        value=0,
        inplace=True,
    )
    table_columns[column_abscissa_interval_above].fillna(
        value=0,
        inplace=True,
    )
    #print("plot function")
    #print(table_columns)
    # Stratify.
    table_group_one = table_columns.loc[
        (
            (table_columns[column_group] == group_one)
        ), :
    ]
    table_group_two = table_columns.loc[
        (
            (table_columns[column_group] == group_two)
        ), :
    ]

    # Organize information for categorical labels on ordinate vertical axis.
    # Assign positions for Group One to be above center point.
    # Assign positions for Group Two to be below center point.
    ordinate_labels = table_group_one[column_ordinate_label].to_list()
    ordinate_positions_center = range(len(ordinate_labels))
    ordinate_positions_one = list(map(
        lambda position: (position + space_groups),
        ordinate_positions_center
    ))
    ordinate_positions_two = list(map(
        lambda position: (position - space_groups),
        ordinate_positions_center
    ))
    # Extract information for labels, values, and error bars.
    abscissa_positions_one = table_group_one[column_abscissa_value].to_numpy()
    abscissa_intervals_one_below = (
        table_group_one[column_abscissa_interval_below].to_numpy()
    )
    abscissa_intervals_one_above = (
        table_group_one[column_abscissa_interval_above].to_numpy()
    )
    abscissa_positions_two = table_group_two[column_abscissa_value].to_numpy()
    abscissa_intervals_two_below = (
        table_group_two[column_abscissa_interval_below].to_numpy()
    )
    abscissa_intervals_two_above = (
        table_group_two[column_abscissa_interval_above].to_numpy()
    )
    # Shape (n, 2)
    #errors_ordinate = numpy.array(list(zip(
    #    errors_ordinate_low, errors_ordinate_high
    #)))
    # Shape (2, n)
    abscissa_intervals_one = numpy.array(
        [abscissa_intervals_one_below, abscissa_intervals_one_above]
    )
    abscissa_intervals_two = numpy.array(
        [abscissa_intervals_two_below, abscissa_intervals_two_above]
    )

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(11.811, 15.748),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    axes.set_xlim(
        xmin=abscissa_minimum,
        xmax=abscissa_maximum,
    )
    # Set titles for axes.
    if (len(abscissa_title) > 0):
        axes.set_xlabel(
            xlabel=abscissa_title,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["one"]
        )
    if (len(ordinate_title) > 0):
        axes.set_ylabel(
            ylabel=ordinate_title,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["one"]
        )
    # Set tick parameters for axes.
    axes.tick_params(
        axis="y", # "y", "x", or "both"
        which="both",
        direction="out",
        length=15.0, # 5.0
        width=11.0, # 3.0
        color=colors["black"],
        pad=25,
        labelsize=fonts["values"][label_size_ordinate_categories]["size"],
        labelcolor=colors["black"]
    )
    axes.tick_params(
        axis="x", # "y", "x", or "both"
        which="both",
        direction="out",
        length=15.0, # 5.0
        width=11.0, # 3.0
        color=colors["black"],
        pad=25,
        labelsize=fonts["values"][label_size_abscissa_values]["size"],
        labelcolor=colors["black"]
    )
    # Set explicit tick positions and labels on vertical ordinate axis.
    # (https://matplotlib.org/3.5.1/api/_as_gen/
    # matplotlib.axes.Axes.set_yticks.html)
    axes.set_xticks(
        [round(abscissa_minimum, 1), 0.0, round(abscissa_maximum, 1)],
        labels=None,
        minor=False,
    )
    axes.set_yticks(
        ordinate_positions_center, # center positions with even spacing
        labels=ordinate_labels, # place labels at center positions
        minor=False,
    )
    # Plot dashed line at zero.
    # Consider the "dashes" argument for fine control of line dash pattern.
    axes.axvline(
        x=0,
        ymin=0,
        ymax=1,
        alpha=1.0,
        color=colors["black"],
        linestyle="--",
        linewidth=11.0, # 7.5
    )
    # Plot points and error bars for values and ranges from each group.
    # First plot markers for group two so that these are below.
    # Second plot markers for group one so that these are above.
    # (https://matplotlib.org/3.5.1/api/_as_gen/
    # matplotlib.axes.Axes.errorbar.html)
    handle_two = axes.errorbar(
        abscissa_positions_two,
        ordinate_positions_two,
        yerr=None,
        xerr=abscissa_intervals_two,
        ecolor=colors["gray_light"],
        elinewidth=11.0, # 7.5
        barsabove=False, # whether to print error bars in layer above points,
        linestyle="",
        marker="D", # "^" marker shape: up triangle
        markersize=size_marker_two, # 5, 15
        markeredgecolor=colors[color_marker_two], # colors["green"],
        markerfacecolor=colors[color_marker_two], # colors["green"],
    )
    handle_one = axes.errorbar(
        abscissa_positions_one,
        ordinate_positions_one,
        yerr=None,
        xerr=abscissa_intervals_one,
        ecolor=colors["gray_dark"],
        elinewidth=13.0, # 7.5
        barsabove=False, # whether to print error bars in layer above points
        linestyle="",
        marker="o", # marker shape: circle
        markersize=size_marker_one, # 5, 15, 50, 70
        markeredgecolor=colors[color_marker_one], # colors["purple"],
        markerfacecolor=colors[color_marker_one], # colors["purple"],
    )
    # Include title label on plot.
    if len(label_chart) > 0:
        matplotlib.pyplot.text(
            0.99,
            0.99,
            label_chart,
            horizontalalignment="right",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["eight"]
        )

    # Return figure.
    return figure


def plot_scatter_points_forest_category_ordinate(
    table=None,
    abscissa=None,
    ordinate=None,
    abscissa_interval_low=None,
    abscissa_interval_high=None,
    title_ordinate=None,
    title_abscissa=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    label_size_ordinate_categories=None,
    fonts=None,
    colors=None,
    size=None,
    label_title=None,
):
    """
    Creates a figure of a chart of type scatter.

    Common name is "Forest Plot".

    This plot's ordinate (vertical, y axis) fits best for a variable with
    discrete categorical values that serve as names.

    This plot's abscissa (horizontal, x axis) fits best for a variable with
    continuous, ratio-scale values that center on zero.

    arguments:
        table (object): Pandas data frame of feature variables across columns
            and observation records across rows
        abscissa (str): name of table's column with variable for horizontal (x)
            axis
        ordinate (str): name of table's column with variable for vertical (y)
            axis
        abscissa_interval_low (str): name of table's column for the extent of
            the error bar below the dot
        abscissa_interval_high (str): name of table's column for the extent of
            the error bar above the dot
        title_ordinate (str): title for ordinate on vertical axis
        title_abscissa (str): title for abscissa on horizontal axis
        minimum_abscissa (float): minimal value for abscissa axis' range
        maximum_abscissa (float): maximal value for abscissa axis' range
        label_size_ordinate_categories (str): label size for categories on
            vertical axis
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        size (int): size of marker
        label_title (str): text label title to include on figure

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    table = table.copy(deep=True)
    columns = [
        abscissa, ordinate, abscissa_interval_low, abscissa_interval_high,
    ]
    table_selection = table.loc[
        :, table.columns.isin(columns)
    ]
    #table_selection.dropna(
    #    axis="index",
    #    how="any",
    #    inplace=True,
    #)
    values_abscissa = table_selection[abscissa].to_numpy()
    values_ordinate = table_selection[ordinate].to_list()
    intervals_abscissa_low = table_selection[abscissa_interval_low].to_numpy()
    intervals_abscissa_high = table_selection[abscissa_interval_high].to_numpy()
    # Shape (n, 2)
    #errors_ordinate = numpy.array(list(zip(
    #    errors_ordinate_low, errors_ordinate_high
    #)))
    # Shape (2, n)
    intervals_abscissa = numpy.array(
        [intervals_abscissa_low, intervals_abscissa_high]
    )

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(11.811, 15.748),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    axes.set_xlim(
        xmin=minimum_abscissa,
        xmax=maximum_abscissa,
    )
    axes.set_xlabel(
        xlabel=title_abscissa,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel=title_ordinate,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.tick_params(
        axis="y",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"][label_size_ordinate_categories]["size"],
        labelcolor=colors["black"]
    )
    axes.tick_params(
        axis="x",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["one"]["size"],
        labelcolor=colors["black"]
    )
    # Plot dashed line at zero.
    axes.axvline(
        x=0,
        ymin=0,
        ymax=1,
        alpha=1.0,
        color=colors["orange"],
        linestyle="--",
        linewidth=7.5,
    )
    # Plot points for values from each group.
    # https://matplotlib.org/3.5.1/api/_as_gen/matplotlib.axes.Axes.errorbar.html
    handle = axes.errorbar(
        values_abscissa,
        values_ordinate,
        yerr=None,
        xerr=intervals_abscissa,
        elinewidth=7.5,
        barsabove=True,
        linestyle="",
        marker="o",
        markersize=size, # 5, 15
        markeredgecolor=colors["blue_navy"],
        markerfacecolor=colors["blue_navy"],
    )
    # Include title label on plot.
    if len(label_title) > 0:
        matplotlib.pyplot.text(
            0.99,
            0.99,
            label_title,
            horizontalalignment="right",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )

    # Return figure.
    return figure


def plot_scatter_points_dot_category_abscissa(
    table=None,
    abscissa=None,
    ordinate=None,
    ordinate_interval_low=None,
    ordinate_interval_high=None,
    title_ordinate=None,
    title_abscissa=None,
    label_size_abscissa_categories=None,
    fonts=None,
    colors=None,
    size=None,
    label_title=None,
):
    """
    Creates a figure of a chart of type scatter.

    Common name is "Dot Plot" or "Point Plot".

    This plot's abscissa (horizontal, x axis) fits best for a variable with
    discrete categorical values that serve as names.

    This plot's ordinate (vertical, y axis) fits best for a variable with
    continuous, ratio-scale values that base on zero.

    arguments:
        table (object): Pandas data frame of feature variables across columns
            and observation records across rows
        abscissa (str): name of table's column with variable for horizontal (x)
            axis
        ordinate (str): name of table's column with variable for vertical (y)
            axis
        ordinate_interval_low (str): name of table's column for the extent of the
            left error bar below the dot
        ordinate_interval_high (str): name of table's column for the extent of the
            right error bar above the dot
        title_ordinate (str): title for ordinate on vertical axis
        title_abscissa (str): title for abscissa on horizontal axis
        minimum_abscissa (float): minimal value for abscissa axis' range
        maximum_abscissa (float): minimal value for abscissa axis' range
        label_size_abscissa_categories (str): label size for categories on
            horizontal axis
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        size (int): size of marker
        label_title (str): text label title to include on figure

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    table = table.copy(deep=True)
    columns = [
        abscissa, ordinate, ordinate_interval_low, ordinate_interval_high,
    ]
    table_selection = table.loc[
        :, table.columns.isin(columns)
    ]
    #table_selection.dropna(
    #    axis="index",
    #    how="any",
    #    inplace=True,
    #)
    values_abscissa = table_selection[abscissa].to_list()
    values_ordinate = table_selection[ordinate].to_numpy()
    intervals_ordinate_low = table_selection[ordinate_interval_low].to_numpy()
    intervals_ordinate_high = table_selection[ordinate_interval_high].to_numpy()
    # Shape (n, 2)
    #errors_ordinate = numpy.array(list(zip(
    #    errors_ordinate_low, errors_ordinate_high
    #)))
    # Shape (2, n)
    intervals_ordinate = numpy.array(
        [intervals_ordinate_low, intervals_ordinate_high]
    )

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    #axes.set_xlim(
    #    xmin=minimum_abscissa,
    #    xmax=maximum_abscissa,
    #)
    axes.set_xlabel(
        xlabel=title_abscissa,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel=title_ordinate,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.tick_params(
        axis="y",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["one"]["size"],
        labelcolor=colors["black"]
    )
    axes.tick_params(
        axis="x",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"][label_size_abscissa_categories]["size"],
        labelcolor=colors["black"],
        labelrotation=45,
    )
    # Plot points for values from each group.
    handle = axes.errorbar(
        values_abscissa,
        values_ordinate,
        yerr=intervals_ordinate,
        xerr=None,
        elinewidth=7.5,
        barsabove=True,
        linestyle="",
        marker="o",
        markersize=size, # 5, 15
        markeredgecolor=colors["blue_navy"],
        markerfacecolor=colors["blue_navy"],
    )
    # Include title label on plot.
    if len(label_title) > 0:
        matplotlib.pyplot.text(
            0.3, # horizontal position (0: left, 1: right)
            0.99, # vertical position
            label_title,
            horizontalalignment="left",
            verticalalignment="top",
            transform=axes.transAxes,
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["four"]
        )

    # Return figure.
    return figure



# TODO: probably obsolete?
def plot_scatter_threshold(
    data=None,
    abscissa=None,
    ordinate=None,
    threshold_abscissa=None,
    selection_abscissa=None,
    threshold_ordinate=None,
    selection_ordinate=None,
    title_abscissa=None,
    title_ordinate=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type scatter with thresholds on each
        dimension.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with independent variable
        ordinate (str): name of data column with dependent variable
        threshold_abscissa (float): threshold for abscissa
        selection_abscissa (str): selection criterion for abscissa's values
            against threshold
        threshold_ordinate (float): threshold for ordinate
        selection_ordinate (str): selection criterion for ordinate's values
            against threshold
        title_abscissa (str): title for abscissa on horizontal axis
        title_ordinate (str): title for ordinate on vertical axis
        factor (str): name of data column with groups or factors of samples
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    data = data.copy(deep=True)
    data = data.loc[:, [abscissa, ordinate]]
    data.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    # Divide values by whether they pass thresholds on both dimensions.
    collection = utility.segregate_data_two_thresholds(
        data=data,
        abscissa=abscissa,
        ordinate=ordinate,
        threshold_abscissa=threshold_abscissa,
        selection_abscissa=selection_abscissa,
        threshold_ordinate=threshold_ordinate,
        selection_ordinate=selection_ordinate,
    )

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    axes.set_xlabel(
        xlabel=title_abscissa,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel=title_ordinate,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["one"]["size"],
        labelcolor=colors["black"]
    )
    # Plot points for values from each group.
    handle = axes.plot(
        collection["fail"][abscissa].values,
        collection["fail"][ordinate].values,
        linestyle="",
        marker="o",
        markersize=2.5,
        markeredgecolor=colors["gray"],
        markerfacecolor=colors["gray"]
    )
    handle = axes.plot(
        collection["pass"][abscissa].values,
        collection["pass"][ordinate].values,
        linestyle="",
        marker="o",
        markersize=5,
        markeredgecolor=colors["blue_navy"],
        markerfacecolor=colors["blue_navy"]
    )

    # Plot lines for each threshold value...
    # Create lines for thresholds.
    axes.axvline(
        x=threshold_abscissa,
        ymin=0,
        ymax=1,
        alpha=1.0,
        color=colors["orange"],
        linestyle="--",
        linewidth=3.0,
    )
    axes.axhline(
        y=threshold_ordinate,
        xmin=0,
        xmax=1,
        alpha=1.0,
        color=colors["orange"],
        linestyle="--",
        linewidth=3.0,
    )

    # Return figure.
    return figure


def plot_scatter_label_emphasis_points(
    emphasis_keys=None,
    label_keys=None,
    column_key=None,
    column_label=None,
    line_abscissa=None,
    line_ordinate=None,
    line_ordinate_origin=None,
    data=None,
    abscissa=None,
    ordinate=None,
    title_abscissa=None,
    title_ordinate=None,
    fonts=None,
    colors=None,
):
    """
    Creates a figure of a chart of type scatter with thresholds on each
        dimension.

    arguments:
        emphasis_keys (list<str>): keys for rows of points for emphasis
        label_keys (list<str>): keys for rows of points for labels
        column_key (str): name of column in data with keys
        column_label (str): name of column in data with labels
        line_abscissa (float): point on abscissa for horizontal line
        line_ordinate (float): point on ordinate for vertical line
        line_ordinate_origin (bool): whether to draw vertical origin line
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with independent variable
        ordinate (str): name of data column with dependent variable
        title_abscissa (str): title for abscissa on horizontal axis
        title_ordinate (str): title for ordinate on vertical axis
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize data.
    data = data.copy(deep=True)
    ordinate_minimum = min(data[ordinate].to_list())
    ordinate_maximum = max(data[ordinate].to_list())
    data_bore = data.loc[~data[column_key].isin(emphasis_keys), :]
    data_emphasis = data.loc[data[column_key].isin(emphasis_keys), :]
    data_label = data.loc[data[column_key].isin(label_keys), :]
    data_bore.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    data_emphasis.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    data_label.dropna(
        axis="index",
        how="any",
        inplace=True,
    )

    ##########
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    axes.set_xlabel(
        xlabel=title_abscissa,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.set_ylabel(
        ylabel=title_ordinate,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["one"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["one"]["size"],
        labelcolor=colors["black"]
    )
    # Plot points for values from each group.
    handle = axes.plot(
        data_bore[abscissa].to_numpy(),
        data_bore[ordinate].to_numpy(),
        linestyle="",
        marker="o",
        markersize=5,
        markeredgecolor=colors["gray"],
        markerfacecolor=colors["gray"]
    )
    handle = axes.plot(
        data_emphasis[abscissa].to_numpy(),
        data_emphasis[ordinate].to_numpy(),
        linestyle="",
        marker="o",
        markersize=7,
        markeredgecolor=colors["blue_navy"],
        markerfacecolor=colors["blue_navy"]
    )
    handle = axes.plot(
        data_label[abscissa].to_numpy(),
        data_label[ordinate].to_numpy(),
        linestyle="",
        marker="o",
        markersize=9,
        markeredgecolor=colors["orange"],
        markerfacecolor=colors["orange"]
    )

    # Plot lines for each threshold value...
    # Create lines for thresholds.
    if line_ordinate_origin:
        axes.axvline(
            x=0.0,
            ymin=0,
            ymax=1,
            alpha=1.0,
            color=colors["black"],
            linestyle="--",
            linewidth=3.0,
        )
    if line_abscissa is not None:
        axes.axvline(
            x=line_abscissa,
            ymin=0,
            ymax=1,
            alpha=1.0,
            color=colors["orange"],
            linestyle="--",
            linewidth=5.0,
        )
    if line_ordinate is not None:
        axes.axhline(
            y=line_ordinate,
            xmin=0,
            xmax=1,
            alpha=1.0,
            color=colors["orange"],
            linestyle="--",
            linewidth=5.0,
        )

    # Plot labels.
    # Place bottom of label above point by 5% of maximal y value.
    for label_key in label_keys:
        data_label = data.loc[data[column_key].isin([label_key]), :]
        if (data_label.shape[0] > 0):
            for index_point in data_label.index.to_list():
                axes.text(
                    (data_label.at[index_point, abscissa]),
                    (
                        data_label.at[index_point, ordinate] +
                        (0.02 * ordinate_maximum)
                    ),
                    data_label[column_label].to_list()[0],
                    backgroundcolor=colors["white_faint"],
                    color=colors["black"],
                    fontproperties=fonts["properties"]["three"],
                    horizontalalignment="center",
                    verticalalignment="bottom"
                )
                pass
            pass
        pass
    # Return figure.
    return figure


def create_legend_elements(
    colors=None,
    labels=None,
):
    """
    Creates custom elements for legend.

    arguments:
        colors (list<dict>): colors
        labels (str): name of data column with independent variable

    raises:

    returns:
        (list<object>): elements for legend

    """

    elements = []
    for index in range(len(labels)):
        element = matplotlib.lines.Line2D(
            [0],
            [0],
            marker="o",
            color=colors[index],
            label=labels[index],
            markerfacecolor=colors[index],
            markersize=15,
        )
        elements.append(element)
    return elements


def plot_overlap_sets(
    sets=None,
    fonts=None,
    colors=None,
):
    """
    Creates a Venn Diagram to represent overlap between multiple sets.

    arguments:
        sets (dict<list<str>>): values in sets
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Organize information.
    names = list(sets.keys())
    # Create figure.
    figure = matplotlib.pyplot.figure(
        figsize=(15.748, 11.811),
        tight_layout=True
    )
    if len(names) == 2:
        venn = matplotlib_venn.venn2(
            subsets=[
                set(sets[names[0]]),
                set(sets[names[1]]),
            ],
            set_labels=names,
        )
    elif len(names) == 3:
        venn = matplotlib_venn.venn3(
            subsets=[
                set(sets[names[0]]),
                set(sets[names[1]]),
                set(sets[names[2]]),
            ],
            set_labels=names,
        )
    return figure


def write_figure(
    figure=None,
    format=None,
    resolution=None,
    path=None,
):
    """
    Writes figure to file.

    arguments:
        figure (object): figure object
        format (str): format suffix, "svg", "png"
        resolution (int): dots per inch (dpi) density for bitmap image, set to
            "None" with format "svg", otherwise "300" or "600"
        path (str): path to directory and file

    raises:

    returns:

    """

    # Write information to file.
    figure.savefig(
        path,
        format=format,
        dpi=resolution,
        facecolor="w",
        edgecolor="w",
        transparent=False
    )
    pass


def write_product_plot_figure(
    name=None,
    figure=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        name (str): base name for file
        figure (object): figure object to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_file = os.path.join(
        path_parent, str(name + ".png")
    )
    # Write information to file.
    write_figure(
        figure=figure,
        format="png",
        resolution=150, # dots per inch: 150, 300, 600
        path=path_file,
    )
    pass


def write_product_plots_child_directories(
    pail_write=None,
    path_parent=None,
):
    """
    Writes product information to file.

    First dictionary tier names the child directory.
    Second dictionary tier names the file.

    arguments:
        pail_write (dict<dict<object>>): information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Structure of "plots" collection is "dict<dict<object>>".
    # First level of plots dictionary tree gives names for child directories.
    # Second level of plots dictionary tree gives names of plots.
    # Iterate across child directories.
    for name_directory in pail_write.keys():
        # Define paths to directories.
        path_child = os.path.join(
            path_parent, name_directory
        )
        # Initialize directories.
        utility.create_directories(
            path=path_child
        )
        # Iterate across charts.
        for name_file in pail_write[name_directory].keys():
            # Write chart object to file in child directory.
            write_product_plot_figure(
                name=name_file,
                figure=pail_write[name_directory][name_file],
                path_parent=path_child,
            )
            pass
        pass
    pass


def write_product_plots_child_child_directories(
    pail_write=None,
    path_parent=None,
):
    """
    Writes product information to file.

    First dictionary tier names the first child directory.
    Second dictionary tier names the second child directory.
    Third dictionary tier names the file.

    arguments:
        pail_write (dict<dict<dict<object>>>): information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Structure of "plots" collection is "dict<dict<dict<object>>>".
    # First level of plots dictionary tree gives names for first child
    # directories.
    # Second level of plots dictionary tree gives names for second child
    # directories.
    # Third level of plots dictionary tree gives names of plots.
    # Iterate across child directories.
    for name_directory in pail_write.keys():
        # Define paths to directories.
        path_child = os.path.join(
            path_parent, name_directory
        )
        # Initialize directories.
        utility.create_directories(
            path=path_child
        )
        # Parse the second level and write files.
        write_product_plots_child_directories(
            pail_write=pail_write[name_directory],
            path_parent=path_child,
        )
        pass
    pass


###############################################################################
# Drivers

# size_marker_one
# size_marker_two

def drive_iterate_plot_forest_two_groups(
    pail_tables=None,
    column_group=None,
    column_ordinate_label=None,
    column_ordinate_sort=None,
    column_abscissa_value=None,
    column_abscissa_interval_below=None,
    column_abscissa_interval_above=None,
    group_one=None,
    group_two=None,
    abscissa_minimum=None,
    abscissa_maximum=None,
    ordinate_title=None,
    abscissa_title=None,
    print_label_on_chart=None,
    label_chart_prefix=None,
    label_size_ordinate_categories=None,
    label_size_abscissa_values=None,
    size_marker_one=None,
    size_marker_two=None,
    color_marker_one=None,
    color_marker_two=None,
    space_groups=None,
):
    """
    Creates a figure of a chart of type scatter with discrete, categorical,
    textual labels on the vertical ordinate axis and continuous points with
    error bars on the horizontal abscissa axis.

    Common name of this chart design is "Forest Plot".

    This function accommodates exactly two groups of values (series).

    This plot's ordinate (vertical, y axis) fits best for a variable with
    discrete categorical values that serve as names.

    This plot's abscissa (horizontal, x axis) fits best for a variable with
    continuous, ratio-scale values that center on zero.

    MatPlotLib accepts intervals, not ranges for error bars. The function does
    the arithmetic to calculate ranges below and above the central value.

    arguments:
        pail_tables (dict<object>): collection of Pandas data-frame tables with
            feature variables across columns and observation records across rows
        column_group (str): name of table's column with names of groups or
            series of values
        column_ordinate_label (str): name of table's column with textual
            categorical values for representation as labels on the vertical
            ordinate (y) axis
        column_ordinate_sort (str): name of table's column with integer values
            for sort order on categorical values for ordinate
        column_abscissa_value (str): name of table's column with floating point
            continuous scale coefficient values for representation as points on
            the horizontal abscissa (x) axis
        column_abscissa_interval_below (str): name of table's column for the
            extent (interval) of the error bar below the point value
        column_abscissa_interval_above (str): name of table's column for the
            extent (interval) of the error bar above the point value
        group_one (str): textual categorical name of first group of values
            (series) within table's column with name 'column_group'
        group_two (str): textual categorical name of second group of values
            (series) within table's column with name 'group'
        abscissa_minimum (float): minimal value for range of abscissa axis
        abscissa_maximum (float): maximal value for range of abscissa axis
        ordinate_title (str): title for ordinate or vertical axis
        abscissa_title (str): title for abscissa or horizontal axis
        print_label_on_chart (bool): whether to print label on chart
        label_chart_prefix (str): text label title of chart to include on figure
        label_size_ordinate_categories (str): label size for categories on
            vertical axis
        label_size_abscissa_values (str): label size for values on
            horizontal axis
        size_marker_one (int): size of markers for group one
        size_marker_two (int): size of markers for group two
        color_marker_one (str): name of color for markers of group one
        color_marker_two (str): name of color for markers of group two
        space_groups (float): vertical spacing between markers for groups
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Collect figures.
    pail_figures = dict()
    # Iterate on tables for figures.
    for name_table in pail_tables.keys():
        # Organize table for figure.
        table = pail_tables[name_table]
        print(label_chart_prefix)
        print(name_table)
        # Chart label.
        if (print_label_on_chart):
            label_chart = str(label_chart_prefix + "_" + name_table)
        else:
            label_chart = ""
        figure = plot_scatter_points_forest_category_ordinate_two_groups(
            table=table,
            column_group=column_group,
            column_ordinate_label=column_ordinate_label,
            column_ordinate_sort=column_ordinate_sort,
            column_abscissa_value=column_abscissa_value,
            column_abscissa_interval_below=column_abscissa_interval_below,
            column_abscissa_interval_above=column_abscissa_interval_above,
            group_one=group_one, # markers for group one are above center
            group_two=group_two, # markers for group two are below center
            abscissa_minimum=abscissa_minimum,
            abscissa_maximum=abscissa_maximum,
            ordinate_title=ordinate_title,
            abscissa_title=abscissa_title,
            label_chart=label_chart,
            label_size_ordinate_categories=label_size_ordinate_categories,
            label_size_abscissa_values=label_size_abscissa_values,
            size_marker_one=size_marker_one,
            size_marker_two=size_marker_two,
            color_marker_one=color_marker_one,
            color_marker_two=color_marker_two,
            space_groups=space_groups, # vertical space between groups' markers
            fonts=fonts,
            colors=colors,
        )
        # Collect figures.
        pail_figures[name_table] = figure
        pass
    # Return information.
    return pail_figures


# TODO: TCW, 06 April 2022... OBSOLETE?
def split_mean_interval_table_plot_dot_category(
    column_group=None,
    column_category=None,
    column_mean=None,
    column_interval_low=None,
    column_interval_high=None,
    table=None,
    title_ordinate=None,
    title_abscissa=None,
    label_size_abscissa_categories=None,
):
    """
    Splits a table by groups and plots correlation Forest Charts for each group.

    arguments:
        column_group (str): name of table's column by which to group and split
            records for each plot
        column_category (str): name of table's column for labels for discrete
            categories
        column_mean (str): name of table's column for means
        column_interval_low (str): name of table's column for lesser interval
        column_interval_high (str): name of table's column for greater interval
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        title_ordinate (str): title for ordinate on vertical axis
        title_abscissa (str): title for abscissa on horizontal axis
        label_size_abscissa_categories (str): label size for categories on
            horizontal axis

    raises:

    returns:
        (dict): collection of figure objects from MatPlotLib

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Organize information for plot.
    # Copy information.
    table = table.copy(deep=True)
    # Select relevant information.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    columns = [
        column_group, column_category, column_mean,
        column_interval_low, column_interval_high,
    ]
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    # Split table by groups.
    table.set_index(
        [column_group],
        append=False,
        drop=True,
        inplace=True
    )
    groups = table.groupby(level=[column_group], axis="index",)
    # Collect tables for groups.
    pail = dict()
    for name, group in groups:
        table_group = group.reset_index(
            level=None,
            inplace=False,
            drop=False,
        )
        print(table_group)
        # Create figure.
        pail[name] = plot_scatter_points_dot_category_abscissa(
            table=table_group,
            abscissa=column_category,
            ordinate=column_mean,
            ordinate_interval_low=column_interval_low,
            ordinate_interval_high=column_interval_high,
            title_ordinate=title_ordinate,
            title_abscissa=title_abscissa,
            label_size_abscissa_categories=label_size_abscissa_categories,
            fonts=fonts,
            colors=colors,
            size=25,
            label_title=str("pheno: " + name),
        )
        pass
    # Return.
    return pail


# TODO: TCW, 06 April 2022... OBSOLETE?
def split_correlation_table_plot_forest_category(
    column_group=None,
    column_category=None,
    column_coefficient=None,
    column_interval_low=None,
    column_interval_high=None,
    table=None,
    title_ordinate=None,
    title_abscissa=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    label_size_ordinate_categories=None,
):
    """
    Splits a table by groups and plots correlation Forest Charts for each group.

    arguments:
        column_group (str): name of table's column by which to group and split
            records for each plot
        column_category (str): name of table's column for labels for discrete
            categories
        column_coefficient (str): name of table's column for coefficient
        column_interval_low (str): name of table's column for lesser interval
        column_interval_high (str): name of table's column for greater interval
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        title_ordinate (str): title for ordinate on vertical axis
        title_abscissa (str): title for abscissa on horizontal axis
        minimum_abscissa (float): minimal value for abscissa axis' range
        maximum_abscissa (float): minimal value for abscissa axis' range
        label_size_ordinate_categories (str): label size for categories on
            vertical axis

    raises:

    returns:
        (dict): collection of figure objects from MatPlotLib

    """

    # Define fonts.
    fonts = define_font_properties()
    # Define colors.
    colors = define_color_properties()

    # Organize information for plot.
    # Copy information.
    table = table.copy(deep=True)
    # Select relevant information.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    columns = [
        column_group, column_category, column_coefficient,
        column_interval_low, column_interval_high,
    ]
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    # Split table by groups.
    table.set_index(
        [column_group],
        append=False,
        drop=True,
        inplace=True
    )
    groups = table.groupby(level=[column_group], axis="index",)
    # Collect tables for groups.
    pail = dict()
    for name, group in groups:
        table_group = group.reset_index(
            level=None,
            inplace=False,
            drop=False,
        )
        print(table_group)
        # Create figure.
        pail[name] = plot_scatter_points_forest_category_ordinate(
            table=table_group,
            abscissa=column_coefficient,
            ordinate=column_category,
            abscissa_interval_low=column_interval_low,
            abscissa_interval_high=column_interval_high,
            title_ordinate=title_ordinate,
            title_abscissa=title_abscissa,
            minimum_abscissa=minimum_abscissa,
            maximum_abscissa=maximum_abscissa,
            label_size_ordinate_categories=label_size_ordinate_categories,
            fonts=fonts,
            colors=colors,
            size=25,
            label_title=str("pheno: " + name),
        )
        pass
    # Return.
    return pail


###############################################################################
# Procedure
# This module is not executable.
