"""
Supply basic plotting functionality.

This module 'plot' is part of the 'partner' package.

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

# TODO: TCW; 1 March 2024
# Plan to learn and implement more advanced designs that include groups of
# variables (especially on boxplots, dotplots, and forestplots) with labels
# for clarity. Also, Forest Plots look especially good when integrated with
# tables that align. PyPi forestplot has an implementation for Forest Plots
# aligned within a table of values.



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
# Utility general purpose


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
        "stretch": 950, # '1000' is the maximal permissable font stretch
        "weight": 950, # '1000' is the maximal permissable font weight
        "size": 50,
    }
    values_4 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 900, # '1000' is the maximal permissable font stretch
        "weight": 900, # '1000' is the maximal permissable font weight
        "size": 45,
    }
    values_5 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 850, # '1000' is the maximal permissable font stretch
        "weight": 850, # '1000' is the maximal permissable font weight
        "size": 40,
    }
    values_6 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 800,
        "weight": 800,
        "size": 35,
    }
    values_7 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 750,
        "weight": 750,
        "size": 30,
    }
    values_8 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 700,
        "weight": 700,
        "size": 27,
    }
    values_9 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 600,
        "weight": 600,
        "size": 25,
    }
    values_10 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 550,
        "weight": 550,
        "size": 23,
    }
    values_11 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 500,
        "weight": 500,
        "size": 20,
    }
    values_12 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 450,
        "weight": 450,
        "size": 17,
    }
    values_13 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 400,
        "weight": 400,
        "size": 15,
    }
    values_14 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 300,
        "weight": 300,
        "size": 13,
    }
    values_15 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 250,
        "weight": 250,
        "size": 10,
    }
    values_16 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 200,
        "weight": 200,
        "size": 7,
    }
    values_17 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 175,
        "weight": 175,
        "size": 6,
    }
    values_18 = {
        "family": "sans-serif",
        "style": "normal",
        "variant": "normal",
        "stretch": 150,
        "weight": 150,
        "size": 5,
    }
    values_19 = {
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
    properties_15 = matplotlib.font_manager.FontProperties(
        family=values_15["family"],
        style=values_15["style"],
        variant=values_15["variant"],
        stretch=values_15["stretch"],
        weight=values_15["weight"],
        size=values_15["size"]
    )
    properties_16 = matplotlib.font_manager.FontProperties(
        family=values_16["family"],
        style=values_16["style"],
        variant=values_16["variant"],
        stretch=values_16["stretch"],
        weight=values_16["weight"],
        size=values_16["size"]
    )
    properties_17 = matplotlib.font_manager.FontProperties(
        family=values_17["family"],
        style=values_17["style"],
        variant=values_17["variant"],
        stretch=values_17["stretch"],
        weight=values_17["weight"],
        size=values_17["size"]
    )
    properties_18 = matplotlib.font_manager.FontProperties(
        family=values_18["family"],
        style=values_18["style"],
        variant=values_18["variant"],
        stretch=values_18["stretch"],
        weight=values_18["weight"],
        size=values_18["size"]
    )
    properties_19 = matplotlib.font_manager.FontProperties(
        family=values_19["family"],
        style=values_19["style"],
        variant=values_19["variant"],
        stretch=values_19["stretch"],
        weight=values_19["weight"],
        size=values_19["size"]
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
            "fifteen": values_15,
            "sixteen": values_16,
            "seventeen": values_17,
            "eighteen": values_18,
            "nineteen": values_19,
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
            "fifteen": properties_15,
            "sixteen": properties_16,
            "seventeen": properties_17,
            "eighteen": properties_18,
            "nineteen": properties_19,
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
    # (0.15 * 255) = 38.25
    # (0.2 * 255) = 51.0
    # (0.3 * 255) = 76.5
    # (0.4 * 255) = 102.0
    # (0.5 * 255) = 127.5
    # 0.55 * 255) = 140.25
    # (0.6 * 255) = 153.0
    # (0.7 * 255) = 178.5
    # (0.8 * 255) = 204.0
    # (0.9 * 255) = 229.5

    # Collect information.
    pail = dict()

    # Black.
    pail["black"] = (0.0, 0.0, 0.0, 1.0)
    # Gray.
    pail["gray"] = (0.5, 0.5, 0.5, 1.0)
    pail["gray_dark"] = (0.3, 0.3, 0.3, 1.0) # (red: 75; green: 75; blue: 75)
    pail["gray_light"] = (0.7, 0.7, 0.7, 1.0) # (red: 180; green: 180; blue: 180)
    pail["gray_faint"] = (0.5, 0.5, 0.5, 0.5) # (red: 180; green: 180; blue: 180)
    # White.
    pail["white"] = (1.0, 1.0, 1.0, 1.0)
    pail["white_faint"] = (1.0, 1.0, 1.0, 0.75)
    # Clear.
    pail["clear"] = (1.0, 1.0, 1.0, 0.0)
    pail["clear_faint"] = (1.0, 1.0, 1.0, 0.25)

    pail["red_brick"] = (0.667, 0.196, 0.039, 1.0) # (red: 170; green: 50; blue: 10)
    pail["red_burgundy"] = (0.510, 0.039, 0.118, 1.0) # (red: 130; green: 10; blue: 30)
    pail["red_crimson"] = (0.784, 0.118, 0.196, 1.0) # (red: 200; green: 30; blue: 50)

    pail["orange"] = (1.0, 0.667, 0.039, 1.0) # (red: 255; green: 170; blue: 10)
    pail["orange_burnt"] = (0.902, 0.510, 0.039, 1.0) # (red: 230; green: 130; blue: 10)
    pail["orange_faint"] = (1.0, 0.588, 0.039, 0.75)

    pail["yellow"] = (1.0, 0.8, 0.0, 1.0) # (red: 255; green: 200; blue: 0)
    pail["yellow_sunshine"] = (1.0, 1.0, 0.196, 1.0) # (red: 255; green: 255; blue: 50)
    pail["yellow_sunflower"] = (1.0, 0.784, 0.0, 1.0) # (red: 255; green: 200; blue: 0)
    pail["yellow_butter"] = (1.0, 1.0, 0.392, 1.0) # (red: 255; green: 255; blue: 100)
    pail["yellow_cream"] = (1.0, 1.0, 0.784, 1.0) # (red: 255; green: 255; blue: 200)

    pail["green"] = (0.25, 0.75, 0.25, 1.0) # (red: 64; green: 191; blue: 64)
    pail["green_forest"] = (0.15, 0.55, 0.15, 1.0) # (red: 38; green: 140; blue: 38)
    pail["green_kelly"] = (0.3, 0.7, 0.1, 1.0) # (red: 77; green: 179; blue: 26)

    pail["blue_navy"] = (0.039, 0.196, 0.588, 1.0) # (red: 10; green: 50; blue: 150)
    pail["blue_navy_faint"] = (0.039, 0.196, 0.588, 0.75)
    pail["blue_navy_light"] = (0.118, 0.314, 0.706, 1.0) # (r: 30; g: 80; b: 180)
    pail["blue_sky"] = (0.196, 0.588, 1.0, 1.0) # (red: 50; green: 150; blue: 255)
    pail["blue_navy_bright"] = (0.784, 0.824, 1.0, 1.0) # (r: 200; g: 210; b: 255)
    pail["blue_periwinkle"] = (0.667, 0.667, 1.0, 1.0) # (red: 170; green: 170; blue: 255)
    pail["blue_indigo"] = (0.196, 0.0, 0.588, 1.0) # (red: 50; green: 0; blue: 150)

    pail["purple"] = (0.510, 0.039, 0.510, 1.0) # (red: 130; green: 10; blue: 130)
    pail["purple_light"] = (0.588, 0.196, 0.588, 1.0) # (red: 150; green: 50; blue: 150)
    pail["magenta"] = (0.784, 0.275, 0.784, 1.0) # (red: 200; green: 70; blue: 200)
    pail["purple_lavender"] = (0.784, 0.784, 1.0, 1.0) # (red: 200; green: 200; blue: 255)
    pail["purple_lilac"] = (0.784, 0.667, 0.784, 1.0) # (red: 200; green: 170; blue: 200)

    # Return information.
    return pail


def initialize_matplotlib_figure_aspect(
    aspect=None,
):
    """
    Initialize a MatPlotLib figure with a specific aspect ratio.

    MatPlotLib Pyplot is an interface that can create artifacts relating to
    global memory environment. For example, on 30 December 2024, TCW had to
    avoid Pyplot to resolve an error of it drawing a rogue dendogram on an
    otherwise empty figure.

    arguments:
        aspect (str): aspect ratio for MatPlotLib figure

    raises:

    returns:
        (object): MatPlotLib figure object

    """

    ##########
    # Create figure.
    # matplotlib.pyplot.figure()
    # matplotlib.figure.Figure()
    if aspect == "portrait":
        figure = matplotlib.figure.Figure(
            figsize=(11.811, 15.748), # aspect 3 X 4; 15.748 inches = 39.999 cm
            tight_layout=True
        )
    elif aspect == "portrait_half_width":
        figure = matplotlib.figure.Figure(
            figsize=(5.906, 15.748), # aspect 1.5 X 4; 5.906 inches = 15.001 cm
            tight_layout=True
        )
        # 80 cm = 31.496 inches
        # 30 cm = 11.811 inches
        #figure = matplotlib.pyplot.figure(
        #    figsize=(11.811, 31.496), # aspect 1.5 X 4; 5.906 inches = 15.001 cm
        #    tight_layout=True
        #)
    elif aspect == "landscape":
        figure = matplotlib.figure.Figure(
            figsize=(15.748, 11.811), # aspect 4 X 3; 11.811 inches = 29.999 cm
            tight_layout=True
        )
    elif aspect == "landscape_half_height":
        figure = matplotlib.figure.Figure(
            figsize=(15.748, 5.906), # aspect 4 X 1.5; 5.906 inches = 15.001 cm
            tight_layout=True
        )
    elif aspect == "square":
        figure = matplotlib.figure.Figure(
            figsize=(15.748, 15.748), # aspect 1 X 1; 15.748 inches = 39.999 cm
            tight_layout=True
        )
    # Initialize.
    #matplotlib.pyplot.clf()
    figure.clear()
    #matplotlib.pyplot.close(figure)

    # Return information.
    return figure



##########
# Heatmap plots


# Current applications:
# 1. genetic correlations
# 2. regression coefficients against polygenic risk scores
# review: TCW; 30 December 2024
# This design has advantages for representations of correlations or regression
# coefficients along with their significance.
def plot_heat_map_few_signal_significance_labels(
    table=None,
    transpose_table=None,
    index_group_columns=None,
    index_group_rows=None,
    fill_missing=None,
    value_missing_fill=None,
    constrain_signal_values=None,
    value_minimum=None,
    value_maximum=None,
    significance_p=None,
    significance_q=None,
    thresholds_p=None,
    thresholds_q=None,
    show_legend=None,
    show_scale_bar=None,
    label_legend=None,
    labels_ordinate_categories=None,
    labels_abscissa_categories=None,
    title_chart_top_right=None,
    title_ordinate=None,
    title_abscissa=None,
    title_bar=None,
    size_label_sig_p=None,
    size_label_sig_q=None,
    size_title_ordinate=None,
    size_title_abscissa=None,
    size_label_ordinate=None,
    size_label_abscissa=None,
    size_label_legend=None,
    size_title_bar=None,
    size_label_bar=None,
    aspect=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Heat map.

    features of this chart design...
    labels of categorical groups on both axes: True
    labels of significance on individual cells: True
    clustering: False

    Format of source table in long format with floating-point values of signals
    and p-values organized with categorical indices across columns and rows
    that will serve as labels:

    group secondary                group_2_a   group_2_b   group_2_c
    group_primary       type_value
    group_1_a           signal     -0.15       -0.2        -0.25
    group_1_a           p_value    0.001       0.001       0.001
    group_1_a           q_value    0.01        0.01        0.01
    group_1_b           signal     0.15        0.2         0.25
    group_1_b           p_value    0.001       0.001       0.001
    group_1_b           q_value    0.01        0.01        0.01

    This function does not support the transposed shape below.

    group secondary     group_1_a                    ...
    type_value          signal   p_value   q_value   ...
    group_primary                                    ...
    group_2_a           -0.15    0.001     0.01      ...
    group_2_b           -0.2     0.001     0.01      ...
    group_2_c           -0.25    0.001     0.01      ...

    This function creates labels on cells of the heatmap to indicate
    significance on the basis of the either the q-value or the p-value.

             dagger: q < 0.05
      double-dagger: q < 0.01
           asterisk: p < 0.05
    double-asterisk: p < 0.01

    str("$\u2020$") # dagger U+2020
    str("$\u2021$") # double dagger U+2021
    str("$\u002A$") # asterisk U+002A
    str("$\u2051$") # double asterisk U+2051

    The advantage of using the symbols dagger, double-dagger, asterisk, and
    double-asterisk is that they all occupy a single character space, which
    helps keep the horizontal alignment clean.

    MatPlotLib color maps.
    https://matplotlib.org/stable/tutorials/colors/colormaps.html

    arguments:
        table (object): Pandas data-frame table in long format with
            floating-point values of signal and p-values.
        transpose_table (bool): whether to transpose matrix from table
        index_group_columns (str): name of index across table's columns that
            provides labels for groups across horizontal axis (abscissa) after
            any transpose on the table
        index_group_rows (str): name of index across table's rows that provides
            labels for groups across vertical axis (ordinate) after any
            transpose on the table
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_missing_fill (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
        significance_p (bool): whether to consider p-values to indicate
            significance in labels on the heatmap
        significance_q (bool): whether to consider q-values to indicate
            significance in labels on the heatmap
        thresholds_p (list<float>): values of thresholds on p-values for
            determination of significance, or actually whether to draw the
            labels on the heatmap
        thresholds_q (list<float>): values of thresholds on q-values for
            determination of significance, or actually whether to draw the
            labels on the heatmap
        show_legend (bool): whether to create and show legend for significnce
            labels
        show_scale_bar (bool): whether to create scale bar
        label_legend (str): text character string for legend label
        labels_ordinate_categories (list<str>): optional, explicit labels for
            ordinate or vertical axis
        labels_abscissa_categories (list<str>): optional, explicit labels for
            abscissa or horizontal axis
        title_chart_top_right (str):
        title_ordinate (str): title for ordinate vertical axis
        title_abscissa (str): title for abscissa horizontal axis
        title_bar (str): title for scale bar
        size_label_sig_p (str): font size
        size_label_sig_q (str): font size
        size_title_ordinate (str): font size
        size_title_abscissa (str): font size
        size_label_ordinate (str): font size
        size_label_abscissa (str): font size
        size_label_legend (str): font size
        size_title_bar (str): font size
        size_label_bar (str): font size
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): MatPlotLib figure object

    """

    ##########
    # Organize information for figure.
    # Copy information in table.
    table = table.copy(deep=True)
    # Separate relevant information.
    if True:
        # Extract relevant information.
        # Extract values.
        # For this extraction, the "type_value" index must be oriented across
        # the table's rows.
        #table_signal = table.loc[
        #    (table["type_value"] == "signal"), :
        #].copy(deep=True)
        table_signal = table[
            table.index.get_level_values("type_value").isin(["signal"])
        ].copy(deep=True)
        table_p = table[
            table.index.get_level_values("type_value").isin(["p_value"])
        ].copy(deep=True)
        table_q = table[
            table.index.get_level_values("type_value").isin(["q_value"])
        ].copy(deep=True)
    else:
        # For this extraction, the "type_value" index must be oriented across
        # the table's columns.
        #table_signal = table.loc[
        #    (table["type_value"] == "signal"), :
        #].copy(deep=True)
        table_signal = table[
            table.columns.get_level_values("type_value").isin(["signal"])
        ].copy(deep=True)
        table_p = table[
            table.columns.get_level_values("type_value").isin(["p_value"])
        ].copy(deep=True)
        table_q = table[
            table.columns.get_level_values("type_value").isin(["q_value"])
        ].copy(deep=True)
        pass
    # Transpose table.
    if (transpose_table):
        table_signal = table_signal.transpose(copy=True)
        table_p = table_p.transpose(copy=True)
        table_q = table_q.transpose(copy=True)
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Rows in table of signals: " + str(table_signal.shape[0]))
        print("Columns in table of signals: " + str(table_signal.shape[1]))
        putly.print_terminal_partition(level=4)
        print("Rows in table of p-values or q-values: " + str(table_p.shape[0]))
        print(
            "Columns in table of p-values or q-values: "
            + str(table_p.shape[1])
        )
        putly.print_terminal_partition(level=4)
    # Extract values.
    #matrix_signal = numpy.transpose(numpy.copy(table_signal.to_numpy()))
    #matrix_p = numpy.transpose(numpy.copy(table_p.to_numpy()))
    #matrix_q = numpy.transpose(numpy.copy(table_q.to_numpy()))
    matrix_signal = numpy.copy(table_signal.to_numpy())
    matrix_p = numpy.copy(table_p.to_numpy())
    matrix_q = numpy.copy(table_q.to_numpy())

    # Organize signals in matrix.
    # Replace missing values.
    if fill_missing:
        matrix_signal = numpy.nan_to_num(
            matrix_signal,
            copy=True,
            nan=value_missing_fill,
            posinf=value_maximum, # or + 1.0 for correlations
            neginf=value_minimum, # or - 1.0 for correlations
        )
    # Constrain values.
    if constrain_signal_values:
        matrix_signal[matrix_signal < value_minimum] = value_minimum
        matrix_signal[matrix_signal > value_maximum] = value_maximum

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Matrix of signals:")
        count_rows = copy.deepcopy(matrix_signal.shape[0])
        count_columns = copy.deepcopy(matrix_signal.shape[1])
        print("Matrix rows (dimension 0): " + str(count_rows))
        print("Matrix columns (dimension 1): " + str(count_columns))
        putly.print_terminal_partition(level=4)

    # Extract categorical names for labels.
    if (
        (
            (labels_ordinate_categories is None) or
            (len(labels_ordinate_categories) < 2)
        ) or
        (
            (labels_abscissa_categories is None) or
            (len(labels_abscissa_categories) < 2)
        )
    ):
        #labels_columns = copy.deepcopy(table_signal.columns.to_list())
        #table_signal.reset_index(
        #    level=None,
        #    inplace=True,
        #    drop=False, # remove index; do not move to regular columns
        #)
        #labels_rows = table.index.to_list()
        #labels_rows = table_signal["group_primary"].to_list()
        labels_columns = copy.deepcopy(
            table_signal.columns.get_level_values(
                index_group_columns
            ).unique().to_list()
        )
        labels_rows = copy.deepcopy(
            table_signal.index.get_level_values(
                index_group_rows
            ).unique().to_list()
        )
        labels_ordinate_categories = labels_rows # vertical axis
        labels_abscissa_categories = labels_columns # horizontal axis
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Column labels:")
        print(labels_ordinate_categories)
        print("Row labels:")
        print(labels_abscissa_categories)
        putly.print_terminal_partition(level=4)

    ##########
    # Create figure.
    figure = initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
    #figure.subplots_adjust(bottom=0.2) # Does this work?
    #axes.margins(
    #    x=1,
    #    y=1,
    #    tight=True,
    #)

    # Create labels or titles on figure.
    # https://stackoverflow.com/questions/34937048/make-part-of-a-title-bold-and-a-different-color
    # ... how to make part of the title string bold-face
    if show_legend:
        axes.set_title(
            str(label_legend),
            #loc="right",
            x=1.050, # 1.0 - 1.2; 1 unit is horizontal dimension of 1 cell?
            y=0.080, # 1.0 - 1.3; 1 unit is vertical dimension of 1 cell? 0.985
            ha="left",
            va="top", # top, center, bottom
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_label_legend],
            bbox=dict(
                facecolor="none",
                edgecolor="black",
                boxstyle="round",
            )
        )
    # This scrap is obsolete and for reference only.
    if False:
        text = axes.text(
            ((matrix_signal.shape[1]) + 0.75),
            -3,
            label_figure_top_right_1,
            horizontalalignment="left",
            verticalalignment="top",
            fontproperties=fonts["properties"][size_label_legend],
            backgroundcolor=colors["white"],
            color=colors["black"],
        )
        text = axes.text(
            ((matrix_signal.shape[1]) + 0.75),
            -2.25,
            label_figure_top_right_2,
            horizontalalignment="left",
            verticalalignment="top",
            fontproperties=fonts["properties"][size_label_legend],
            backgroundcolor=colors["white"],
            color=colors["black"],
        )
        text = axes.text(
            ((matrix_signal.shape[1]) + 0.75),
            -1.5,
            label_figure_top_right_3,
            horizontalalignment="left",
            verticalalignment="top",
            fontproperties=fonts["properties"][size_label_legend],
            backgroundcolor=colors["white"],
            color=colors["black"],
        )

    # Plot values as a color grid.
    # This function represents values acros matrix dimension 0 as vertical rows.
    # This function represents values across matrix dimension 1 as horizontal
    # columns.
    # Diverging color maps: "PRGn", "PRGn_r", "PiYG", "PiYG_r",
    # Diverging color maps: "PuOr", "PuOr_r",
    # Diverging color maps: "PuOr", "PuOr_r", "RdBu", "RdBu_r", "BrBG",
    # Sequential color maps: "Reds", "Reds_r", "Oranges", "Oranges_r",
    # site: https://montoliu.naukas.com/2021/11/18/color-blindness-purple-and-orange-are-the-solution/
    image = axes.imshow(
        matrix_signal,
        cmap=matplotlib.colormaps["PuOr"], # RdBu_r, PuOr_r
        vmin=value_minimum,
        vmax=value_maximum,
        aspect="auto", # "auto", "equal",
        origin="upper",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )

    # Set titles for axes.
    if (len(title_ordinate) > 0):
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_ordinate]
        )
    if (len(title_abscissa) > 0):
        axes.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
    # Set tick parameters for axes.
    axes.tick_params(
        axis="both", # "y", "x", or "both"
        which="both", # "major", "minor", or "both"
        length=7.5, # 5.0
        width=5, # 3.5
        pad=10, # 7.5
        direction="out",
        color=colors["black"],
        labelcolor=colors["black"],
        top=True,
        bottom=False,
        left=True,
        right=False,
        labeltop=True,
        labelbottom=False,
        labelleft=True,
        labelright=False,
    )
    # Set tick positions and labels on axes.
    axes.set_yticks(
        numpy.arange(matrix_signal.shape[0]),
    )
    axes.set_yticklabels(
        labels_ordinate_categories,
        #minor=False,
        ha="right", # horizontal alignment
        va="center", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"][size_label_ordinate]
    )
    axes.set_xticks(
        numpy.arange(matrix_signal.shape[1]),
    )
    axes.set_xticklabels(
        labels_abscissa_categories,
        #minor=False,
        rotation=-60,
        rotation_mode="anchor",
        ha="left", # horizontal alignment
        va="top", # vertical alignment
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"][size_label_abscissa]
    )
    axes.tick_params(
        top=False,
        labeltop=False,
        bottom=True,
        labelbottom=True,
    )
    # Keep axes, ticks, and labels, but remove border.
    for position in ['right', 'top', 'bottom', 'left']:
        matplotlib.pyplot.gca().spines[position].set_visible(False)

    # Create value labels on individual cells on the chart.
    # Iterate across values in matrices.
    for index_row in range(matrix_signal.shape[0]):
        for index_column in range(matrix_signal.shape[1]):
            # Extract signal, p-value, and q-value from respective matrices.
            signal_value = matrix_signal[index_row, index_column]
            signal_text = "{:.2f}".format(round(signal_value, 2))
            p_value = matrix_p[index_row, index_column]
            q_value = matrix_q[index_row, index_column]
            if ((signal_value > 0.5) or (signal_value < -0.5)):
                #color="white"
                color="gray_light"
                #color="gray"
            else:
                #color="black"
                #color="gray_dark"
                color="gray"
                pass

            # Create labels on cells to represent p-values and q-values.
            #str("$\u2020$") # dagger U+2020
            #str("$\u2021$") # double dagger U+2021
            #str("$\u002A$") # asterisk U+002A
            #str("$\u2051$") # double asterisk U+2051

            # Determine whether to create labels for p-values on cells.
            if (
                (significance_p) and
                (not numpy.isnan(p_value)) and
                (p_value < thresholds_p[1])
            ):
                #label_cell = str(str(signal_text) + "**")
                label_cell = str("$\u2021$") # double dagger U+2021
                text = axes.text(
                    index_column,
                    index_row,
                    label_cell,
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontproperties=fonts["properties"][size_label_sig_p],
                    fontweight="extra bold",
                    backgroundcolor=colors["clear"],
                    color=colors[color],
                )
            elif (
                (significance_p) and
                (not numpy.isnan(p_value)) and
                (p_value < thresholds_p[0])
            ):
                label_cell = str("$\u2020$") # dagger U+2020
                text = axes.text(
                    index_column,
                    index_row,
                    label_cell,
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontproperties=fonts["properties"][size_label_sig_p],
                    fontweight="extra bold",
                    backgroundcolor=colors["clear"],
                    color=colors[color],
                )
                pass

            # Determine whether to create labels for q-values on cells.
            # Use the 'annotate' method for more versatility in position.
            if (
                False and
                (significance_q) and
                (not numpy.isnan(q_value)) and
                (q_value < thresholds_q[1])
            ):
                #label_cell = str(str(signal_text) + "**")
                label_cell = str("$\u2051$") # double asterisk U+2051
                text = axes.annotate(
                    label_cell,
                    xy=(index_column,index_row),
                    xycoords="data", # use coordinate system of object
                    xytext=(10,15), # coordinates for offset of text from main
                    textcoords="offset points", # coordinates for 'xytext'
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontproperties=fonts["properties"][size_label_sig_q],
                    backgroundcolor=colors["clear"],
                    color=colors[color],
                )
            elif (
                (significance_q) and
                (not numpy.isnan(q_value)) and
                (q_value < thresholds_q[0])
            ):
                label_cell = str("$\u002A$") # asterisk U+002A
                text = axes.annotate(
                    label_cell,
                    xy=(index_column,index_row),
                    xycoords="data", # use coordinate system of object
                    xytext=(7,10), # coordinates for offset of text from main
                    textcoords="offset points", # coordinates for 'xytext'
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontproperties=fonts["properties"][size_label_sig_q],
                    backgroundcolor=colors["clear"],
                    color=colors[color],
                )
                pass
            pass
        pass

    # Create legend for scale of color grid.
    if show_scale_bar:
        bar = axes.figure.colorbar(
            image,
            orientation="vertical",
            ax=axes,
            location="right",
            shrink=0.9, # 0.7; factor for dimensions of the Scale Bar.
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
            length=7.5, # 5.0, 7.5
            width=3, # 2.5, 5.0
            color=colors["black"],
            pad=5, # 5, 7
            labelsize=fonts["values"][size_label_bar]["size"],
            labelcolor=colors["black"],
        )
    # Create labels for chart.
    # TODO: TCW; 18 January 2024
    # TODO: I need to troubleshoot this... I want the label to be on the figure
    # to the top right of the chart area... maybe need to create text label
    # on "figure" not "axes".
    if False:
        if len(label_top_right) > 0:
            matplotlib.pyplot.text(
                0.99,
                0.99,
                label_top_right,
                horizontalalignment="right",
                verticalalignment="top",
                transform=axes.transAxes,
                backgroundcolor=colors["white_faint"],
                color=colors["black"],
                fontproperties=fonts["properties"]["four"]
            )
    # Return figure.
    return figure


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
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
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
        size_label = "nineteen"
    elif (120 <= count_labels and count_labels < 150):
        size_label = "eighteen"
    elif (100 <= count_labels and count_labels < 120):
        size_label = "seventeen"
    elif (80 <= count_labels and count_labels < 100):
        size_label = "sixteen"
    elif (60 <= count_labels and count_labels < 80):
        size_label = "fifteen"
    elif (40 <= count_labels and count_labels < 60):
        size_label = "fourteen"
    elif (30 <= count_labels and count_labels < 40):
        size_label = "thirteen"
    elif (20 <= count_labels and count_labels < 30):
        size_label = "twelve"
    elif (10 <= count_labels and count_labels < 20):
        size_label = "eleven"
    elif (1 <= count_labels and count_labels < 10):
        size_label = "ten"
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

    Review: 30 December 2024

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
        value_minimum (float): minimal value for constraint on signals and
            scale
        value_maximum (float): maximal value for constraint on signals and
            scale
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
        value_maximum=value_maximum,
        labels_ordinate_categories=labels_ordinate_categories,
        labels_abscissa_categories=labels_abscissa_categories,
        report=report,
    )

    ##########
    # Create figure.
    figure = initialize_matplotlib_figure_aspect(
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
        right=0.98,
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

    # Create legend for scale of color grid.
    if show_scale_bar:
        bar = axes_main.figure.colorbar(
            image_main,
            orientation="vertical",
            ax=axes_main,
            location="right",
            shrink=0.7, # 0.7; factor for dimensions of the Scale Bar.
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
            length=7.5, # 5.0, 7.5
            width=3, # 2.5, 5.0
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
    figure = initialize_matplotlib_figure_aspect(
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
        width_ratios=(15,80,2,3),
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
        left=0.02,
        right=0.93,
        top=0.98,
        bottom=0.15,
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
                labelpad=5, # 5
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
    figure = initialize_matplotlib_figure_aspect(
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
    axes_set = figure.add_subplot(grid[0,0]) # first row, first column
    axes_main = figure.add_subplot(grid[0,1]) # first row, second column
    axes_space_first = figure.add_subplot(grid[1,0]) # second row, first column
    axes_group = figure.add_subplot(grid[1,1]) # second row, second column
    axes_space_second = figure.add_subplot(grid[2,1]) # third row, second column
    axes_bar = figure.add_subplot(grid[3,1]) # fourth row, second column
    axes_set.clear()
    axes_main.clear()
    axes_space_first.clear()
    axes_group.clear()
    axes_space_second.clear()
    axes_bar.clear()
    # Set axes to empty as a space holder.
    axes_space_first.axis("off")
    axes_space_second.axis("off")
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




#######################################
########### Move this older scrap code to the "bimodality" legacy package

##########
# Heatmap plots with many cells... old
# These heatmap plots support options such as clustering.

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

# TODO: TCW; 11 October 2024
# TODO: currently unused or obsolete, but useful for scrap reference?
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
        data_cluster_columns = porg.cluster_data_columns(
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
            data_cluster_rows = porg.cluster_data_rows_by_group(
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
        data_cluster_columns = porg.cluster_data_columns(
            data=data,
        )
        data_cluster_rows = porg.cluster_data_rows(
            data=data_cluster_columns,
        )
        data_sequence = data_cluster_rows.copy(deep=True)
    # Return information.
    return data_sequence

# TODO: TCW; 11 October 2024
# TODO: currently unused or obsolete, but useful for scrap reference?
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
    values_master_unique = putly.collect_unique_elements(
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

########################3
#############################################################



##########
# Histogram


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
        bar_width (float): proportional width of bar relative to bin
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

    ##########
    # Create figure.
    figure = initialize_matplotlib_figure_aspect(
        aspect="landscape",
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
            prop=fonts["properties"]["five"], # one before
            edgecolor=colors["black"]
        )
    axes.set_xlabel(
        xlabel=label_bins,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["five"] # two before
    )
    axes.set_ylabel(
        ylabel=label_counts,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["five"] # two before
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["five"]["size"], # two before
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
            fontproperties=fonts["properties"]["six"] # four before
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
            fontproperties=fonts["properties"]["seven"] # four before
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
            fontproperties=fonts["properties"]["seven"] # four before
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
            fontproperties=fonts["properties"]["seven"] # four before
        )
    return figure


##########
# Bar Chart


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
    Creates a figure of a chart of type bar chart to represent the values on
    a continuous scale.

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


# Dot plot
# set minima and maxima definitively... refer to the Forest plots...



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
    report=None,
):
    """
    Creates a figure of a chart of type box plot to represent the distribution
    of multiple series of values.

         Q1-1.5IQR   Q1   median  Q3   Q3+1.5IQR
                      |-----:-----|
      o      |--------|     :     |--------|    o  o
                      |-----:-----|
    flier             <----------->            fliers
                           IQR

    lines in boxes: median
    diamond on boxes: mean
    box length: quartile 1 (Q1) to quartile 3 (Q3)
    whiskers: 1.5 * interquartile range (Q3 - Q1)

    https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.boxplot.html
    https://towardsdatascience.com/understanding-boxplots-5e2df7bcbd51

    Review: 17 October 2024

    arguments:
        values_groups (list<array>): NumPy arrays of non-missing values for
            each group
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
        report (bool): whether to print reports

    raises:

    returns:
        (object): figure object

    """

    #colors_groups = list(seaborn.color_palette("hls", n_colors=color_count))

    ##########
    # Create figure.
    figure = initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Create axes.
    #axes = matplotlib.pyplot.axes()
    axes = figure.add_subplot(111)
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
            "linewidth": 2.5, # 1.0, 2.5, 5.0
            "color": colors["black"],
        },
        meanprops={
            "marker": "D", # diamond
            "markersize": 10.0, # 10.0, 20.0
            "markeredgecolor": colors["black"], # orange_burnt
            "markerfacecolor": colors["black"],
        },
        whiskerprops={
            "linestyle": "solid",
            "linewidth": 2.5,
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
    if colors_groups is None:
        colors_groups = (
            matplotlib.colormaps["tab10"].colors[0:len(values_groups)]
        )
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
            fontproperties=fonts["properties"]["seven"],
            rotation="horizontal",
        )
    if len(title_ordinate) > 0:
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=20,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"]["seven"]
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
            labelsize=fonts["values"]["nine"]["size"],
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
            labelsize=fonts["values"]["nine"]["size"],
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
            labelsize=fonts["values"]["seven"]["size"],
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
            labelsize=fonts["values"]["seven"]["size"],
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


##########
# Chart type: scatter


def plot_scatter_fold_change_volcano(
    table=None,
    column_identifier=None,
    column_name=None,
    column_plot_fold=None,
    column_plot_p=None,
    column_fold_change=None,
    column_significance=None,
    threshold_fold_change=None,
    threshold_significance=None,
    identifiers_emphasis=None,
    lines_thresholds=None,
    emphasis_label=None,
    count_label=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    title_abscissa=None,
    title_ordinate=None,
    size_title_abscissa=None,
    size_title_ordinate=None,
    size_label_abscissa=None,
    size_label_ordinate=None,
    size_label_emphasis=None,
    size_label_count=None,
    aspect=None,
    fonts=None,
    colors=None,
    report=None,
):
    """
    Creates a chart figure of basic type scatter that represents fold changes
    with thresholds on two dimensions.

    Volcano Plot
    chart type: scatter
    abscissa
       - axis: horizontal, x
       - representation: fold change on base-two logarithmic scale
    ordinate
       - axis: vertical, y
       - representation: p-value or q-value for estimate of fold change on
          negative base-ten logarithmic scale

    Review: TCW; 3 October 2024

    arguments:
        table (object): Pandas data-frame table of features across columns and
            values for observations across rows
        column_identifier (str): name of column in table for the unique
            identifier corresponding to the fold change
        column_name (str): name of column in table for the name corresponding
            to the fold change
        column_plot_fold (str): name of column in table for values of fold
            change for representation on the plot
        column_plot_p (str): name of column in table for values of p-value or
            representation of significance for representation on the plot
        column_fold_change (str): name of column in table on which to apply the
            threshold for the fold change
        column_significance (str): name of column in table on which to apply
            the threshold for significance, corresponding to the p-value or
            q-value corresponding to the estimate of fold change
        threshold_fold_change (float): value for threshold on fold change
            (fold change > threshold) that is on the same scale, such as
            base-two logarithm, as the actual values themselves
        threshold_significance (float): value for threshold on p-values or
            q-values (p-value or q-value < threshold) that is not on a scale of
            the negative logarithm
        identifiers_emphasis (list<str>): identifiers corresponding to a
            special selection of fold changes for which to emphasize points on
            chart and for which to create text labels adjacent to the points
            on the chart
        lines_thresholds (bool): whether to draw lines to represent thresholds
        emphasis_label (bool): whether to create text labels adjacent to
            the special selection of points for special emphasis
        count_label (bool): whether to create text labels on chart to report
            counts of fold changes that pass thresholds
        minimum_abscissa (float): value for minimal limit to represent on the
            abscissa horizontal axis
        maximum_abscissa (float): value for maximal limit to represent on the
            abscissa horizontal axis
        minimum_ordinate (float): value for minimal limit to represent on the
            ordinate horizontal axis
        maximum_ordinate (float): value for maximal limit to represent on the
            ordinate horizontal axis
        title_abscissa (str): title for abscissa horizontal axis
        title_ordinate (str): title for ordinate vertical axis
        size_title_abscissa (str): font size for title on abscissa horizontal
            axis
        size_title_ordinate (str): font size for title on ordinate vertical
            axis
        size_label_abscissa (str): font size for labels on abscissa horizontal
            axis
        size_label_ordinate (str): font size for labels on ordinate vertical
            axis
        size_label_emphasis (str): font size for labels adjacent to points for
            special emphasis
        size_label_count (str): font size for labels to report counts of
            fold changes that pass thresholds
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): definitions of font properties
        colors (dict<tuple>): definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): chart figure object

    """

    ##########
    # Organize information for chart.

    # Copy information in table.
    table = table.copy(deep=True)
    # Filter columns in table.
    #table = table.loc[
    #    :, table.columns.isin(columns_sequence)
    #]
    table = table.filter(
        items=[
            column_identifier,
            column_name,
            column_plot_fold,
            column_plot_p,
            column_fold_change,
            column_significance,
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
            "column": column_plot_fold,
        },
        {
            "type": "high",
            "value": maximum_abscissa,
            "column": column_plot_fold,
        },
        {
            "type": "low",
            "value": minimum_ordinate,
            "column": column_plot_p,
        },
        {
            "type": "high",
            "value": maximum_ordinate,
            "column": column_plot_p,
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
    # Filter rows in table for segregation of values by thresholds.
    pail_threshold = porg.segregate_fold_change_values_by_thresholds(
        table=table,
        column_fold_change=column_fold_change,
        column_significance=column_significance,
        threshold_fold_change=threshold_fold_change,
        threshold_significance=threshold_significance,
        report=False,
    )
    # Extract information about selection of fold changes for special emphasis.
    #table_bore = table.loc[
    #    ~table[column_identifier].isin(identifiers_emphasis), :
    #].copy(deep=True)
    table_source = pail_threshold["table_pass_any"]
    table_emphasis = table_source.loc[
        table_source[column_identifier].isin(identifiers_emphasis), :
    ].copy(deep=True)
    # Extract counts.
    count_pass = int(pail_threshold["table_pass_any"].shape[0])
    count_pass_up = int(pail_threshold["table_pass_up"].shape[0])
    count_pass_down = int(pail_threshold["table_pass_down"].shape[0])
    count_fail = int(pail_threshold["table_fail"].shape[0])

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: partner.plot.py")
        print("function: plot_scatter_fold_change_volcano()")
        putly.print_terminal_partition(level=5)
        print(
            "count of fold changes that pass segregation thresholds: " +
            str(count_pass)
        )
        print("count ups: " + str(count_pass_up))
        print("count downs: " + str(count_pass_down))
        print(
            "count of fold changes that fail segregation thresholds: " +
            str(count_fail)
        )
        putly.print_terminal_partition(level=5)
        print("table of fold changes for special emphasis:")
        print(table_emphasis)
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Create and initialize figure chart object.

    # Create figure.
    figure = initialize_matplotlib_figure_aspect(
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
    if (minimum_abscissa is not None):
        axes.set_xlim(xmin=minimum_abscissa)
    if (maximum_abscissa is not None):
        axes.set_xlim(xmax=maximum_abscissa)
    if (minimum_ordinate is not None):
        axes.set_ylim(ymin=minimum_ordinate)
    if (maximum_ordinate is not None):
        axes.set_ylim(ymax=maximum_ordinate)

    # Set titles for axes.
    if (len(title_abscissa) > 0):
        axes.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
    if (len(title_ordinate) > 0):
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
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

    # Create lines to represent threshold values.
    if (lines_thresholds and (count_pass > 0)):
        # Determine minimal value of significance for line.
        threshold_significance_scale = numpy.nanmin(
            pail_threshold["table_pass_any"][column_plot_p].to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            )
        )
        axes.axhline(
            y=threshold_significance_scale,
            xmin=minimum_abscissa,
            xmax=maximum_abscissa,
            alpha=1.0,
            color=colors["black"],
            linestyle="--",
            linewidth=2.5,
        )
        axes.axvline(
            x=threshold_fold_change,
            ymin=minimum_ordinate,
            ymax=maximum_ordinate,
            alpha=1.0,
            color=colors["black"],
            linestyle="--",
            linewidth=2.5,
        )
        if (threshold_fold_change != 0):
            axes.axvline(
                x=(-1*threshold_fold_change),
                ymin=minimum_ordinate,
                ymax=maximum_ordinate,
                alpha=1.0,
                color=colors["black"],
                linestyle="--",
                linewidth=2.5,
            )
            pass
        pass

    ##########
    # Represent information on the chart figure object.

    # Plot points for values from each group.
    handle = axes.plot(
        pail_threshold["table_fail"][column_plot_fold].values,
        pail_threshold["table_fail"][column_plot_p].values,
        linestyle="",
        marker="o",
        markersize=2.5,
        markeredgecolor=colors["gray"],
        markerfacecolor=colors["gray"]
    )
    handle = axes.plot(
        pail_threshold["table_pass_any"][column_plot_fold].values,
        pail_threshold["table_pass_any"][column_plot_p].values,
        linestyle="",
        marker="o",
        markersize=5,
        markeredgecolor=colors["blue_navy"],
        markerfacecolor=colors["blue_navy"]
    )
    handle = axes.plot(
        table_emphasis[column_plot_fold].values,
        table_emphasis[column_plot_p].values,
        linestyle="",
        marker="o",
        markersize=5,
        markeredgecolor=colors["orange"],
        markerfacecolor=colors["orange"]
    )

    # Plot labels to report counts of fold changes that pass thresholds.
    if count_label:
        # Determine position coordinates of labels.
        center_abscissa_down = (
            minimum_abscissa + (
                abs(minimum_abscissa - (-1*threshold_fold_change))/2.5
        ))
        center_abscissa_up = (
            maximum_abscissa - ((maximum_abscissa - threshold_fold_change)/2.5)
        )
        center_ordinate = ((maximum_ordinate - 0)/2)
        #count_pass_up
        #count_pass_down
        # Create labels on chart.
        axes.text(
            center_abscissa_up,
            center_ordinate,
            str("count: " + str(count_pass_up)),
            horizontalalignment="right",
            verticalalignment="center",
            backgroundcolor=colors["white_faint"],
            color=colors["orange"],
            fontproperties=fonts["properties"][size_label_count],
        )
        axes.text(
            center_abscissa_down,
            center_ordinate,
            str("count: " + str(count_pass_down)),
            horizontalalignment="left",
            verticalalignment="center",
            backgroundcolor=colors["white_faint"],
            color=colors["orange"],
            fontproperties=fonts["properties"][size_label_count],
        )
        pass

    # Plot labels adjacent to the selection of points for special emphasis.
    # Align bottom left of label with an offset of 1% of the chart coordinate
    # dimensions.
    if emphasis_label:
        for index, row in table_emphasis.iterrows():
            # Determine position coordinates of label.
            abscissa_raw = row[column_plot_fold]
            if (abscissa_raw > 0):
                alignment_horizontal = "right"
                abscissa = (
                    abscissa_raw - (
                        0.01 * (maximum_abscissa - minimum_abscissa)
                    )
                )
            if (abscissa_raw < 0):
                alignment_horizontal = "left"
                abscissa = (
                    abscissa_raw + (
                        0.01 * (maximum_abscissa - minimum_abscissa)
                    )
                )
            ordinate = row[column_plot_p]
            # Create label on chart.
            axes.text(
                abscissa,
                ordinate,
                str(row[column_name]),
                horizontalalignment=alignment_horizontal,
                verticalalignment="center",
                backgroundcolor=colors["white_faint"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_label_emphasis],
            )
            pass
        pass

    ##########
    # Return figure.
    return figure


def plot_scatter_point_color_response_discrete_or_continuous(
    table=None,
    column_identifier=None,
    column_name=None,
    column_response=None,
    column_abscissa=None,
    column_ordinate=None,
    type_response=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    set_axis_limits=None,
    title_response=None,
    title_abscissa=None,
    title_ordinate=None,
    identifiers_emphasis=None,
    emphasis_label=None,
    line_diagonal=None,
    size_title_abscissa=None,
    size_title_ordinate=None,
    size_title_legend_bar=None,
    size_label_abscissa=None,
    size_label_ordinate=None,
    size_label_emphasis=None,
    size_label_legend_bar=None,
    aspect=None,
    fonts=None,
    colors=None,
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

    Review: TCW; 2 December 2024

    arguments:
        table (object): Pandas data-frame table of features across columns and
            values for observations across rows
        column_identifier (str): name of column in table corresponding to the
            unique identifier of records for each point
        column_name (str): name of column in table corresponding to the name of
            records for each point
        column_response (str): name of column in table corresponding to values
            for representation as color of individual points
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
        set_axis_limits (bool): whether to set explicity limits on axes
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
        emphasis_label (bool): whether to create text labels adjacent to
            the special selection of points for special emphasis
        line_diagonal (bool): whether to draw diagonal line for equality
            between abscissa and ordinate
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
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): definitions of font properties
        colors (dict<tuple>): definitions of color properties
        report (bool): whether to print reports

    raises:

    returns:
        (object): chart figure object

    """

    ##########
    # Organize information for chart.

    # Copy information in table.
    table = table.copy(deep=True)
    # Filter columns in table.
    #table = table.loc[
    #    :, table.columns.isin(columns_sequence)
    #]
    table = table.filter(
        items=[
            column_identifier,
            column_name,
            column_response,
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
        minimum_abscissa = numpy.nanmin(values_abscissa_raw)
    if (maximum_abscissa is None):
        maximum_abscissa = numpy.nanmax(values_abscissa_raw)
    if (minimum_ordinate is None):
        minimum_ordinate = numpy.nanmin(values_ordinate_raw)
    if (maximum_ordinate is None):
        maximum_ordinate = numpy.nanmax(values_ordinate_raw)
        pass
    center_abscissa = ((maximum_abscissa - minimum_abscissa)/2)
    center_ordinate = ((maximum_ordinate - minimum_ordinate)/2)
    # Extract information about range of third response feature for
    # representation as color.
    if (
        (column_response is not None) and
        (column_response != "") and
        (type_response is not None) and
        (type_response == "continuity")
    ):
        values_response_raw = table[column_response].to_numpy(
            dtype="float64",
            na_value=numpy.nan,
            copy=True,
        )
        minimum_response = numpy.nanmin(values_response_raw)
        maximum_response = numpy.nanmax(values_response_raw)
        pass

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

    # Stratify groups of records by categorical factor response feature.
    if (
        (column_response is not None) and
        (column_response != "") and
        (type_response is not None) and
        (type_response == "category")
    ):
        # Copy information in table.
        table_standard_group = table_standard.copy(deep=True)
        # Count unique categorical values.
        count_category_response = (
            table_standard[column_response].nunique(dropna=True)
        )
        # Organize indices in table.
        table_standard_group.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table_standard_group.set_index(
            [column_identifier, column_response],
            append=False,
            drop=True,
            inplace=True
        )
        # Split rows within table by factor columns.
        groups_response = table_standard_group.groupby(
            level=column_response,
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: plot.py")
        function = "plot_scatter_point_color_response_discrete_or_continuous()"
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Create and initialize figure chart object.

    # Create figure.
    figure = initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Create axes.
    axes = matplotlib.pyplot.axes()
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

    # Set titles for axes.
    if (len(title_abscissa) > 0):
        axes.set_xlabel(
            xlabel=title_abscissa,
            labelpad=30,
            alpha=1.0,
            loc="left",
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
    if (len(title_ordinate) > 0):
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=30,
            alpha=1.0,
            backgroundcolor=colors["white"],
            color=colors["black"],
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

    # Create diagonal line to represent equality between abscissa and ordinate.
    # Notice that the current definition does not actually correspond to a 1:1
    # relationship between abscissa and ordinate scales. Instead, it depends on
    # the relative ranges of the abscissa and ordinate axes.
    if (line_diagonal):
        axes.plot(
            [0, 1,],
            [0, 1,],
            transform=axes.transAxes,
            alpha=1.0,
            color=colors["black"],
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
            backgroundcolor=colors["white_faint"],
            color=colors["black"],
            fontproperties=fonts["properties"]["seventeen"],
        )
        pass

    ##########
    # Represent information on the chart figure object.
    # For table of standard selection of records, determine whether to set
    # color of points to represent a third response feature.
    if (
        (column_response is None) or
        (column_response == "") or
        (type_response is None)
    ):
        handle_standard = axes.plot(
            table_standard[column_abscissa].values,
            table_standard[column_ordinate].values,
            linestyle="",
            marker="o",
            markersize=10,
            markeredgecolor=colors["blue_navy"],
            markerfacecolor=colors["blue_navy"],
        )
    elif (
        (column_response is not None) and
        (column_response != "") and
        (type_response is not None) and
        (type_response == "continuity")
    ):
        handle_standard = axes.scatter(
            table_standard[column_abscissa].values,
            table_standard[column_ordinate].values,
            c=table_standard[column_response].values,
            s=300, # scale with area (dimension squared)
            norm="linear",
            cmap="binary", # 'binary', 'plasma', 'viridis', 'civids'
            vmin=minimum_response,
            vmax=maximum_response,
            alpha=1,
            marker="o",
            linewidths=1,
            edgecolors=colors["black"],
            #linestyle="",
        )
    elif (
        (column_response is not None) and
        (column_response != "") and
        (type_response is not None) and
        (type_response == "category")
    ):
        # An alternative would be to map categorical values to discrete
        # integers that respectively map to colors. Then create discrete color
        # bar as a reference in place of a legend.
        # Create discrete color map for categorical values.
        color_map = matplotlib.cm.get_cmap("tab10")
        colors_groups = [color_map(i) for i in range(count_category_response)]
        # Collect information.
        labels_groups = []
        counter = 0
        # Iterate on groups of records from table.
        for name_group, table_group in groups_response:
            # Copy information in table.
            table_group = table_group.copy(deep=True)
            # Collect information.
            label_group = str(name_group)
            labels_groups.append(label_group)
            # Create points on plot chart.
            handle_group = axes.plot(
                table_group[column_abscissa].values,
                table_group[column_ordinate].values,
                linestyle="",
                marker="o",
                markersize=10,
                markeredgecolor=colors_groups[counter],
                markerfacecolor=colors_groups[counter],
            )
            # Update index counter.
            counter += 1
            pass
        pass

    # For table of special selection of records, plot points with a special
    # color for emphasis.
    if (
        (identifiers_emphasis is not None) and
        (len(identifiers_emphasis) > 0) and
        (not table_special.empty)
    ):
        handle_special = axes.plot(
            table_special[column_abscissa].values,
            table_special[column_ordinate].values,
            linestyle="",
            marker="o",
            markersize=5,
            markeredgecolor=colors["red_crimson"],
            markerfacecolor=colors["red_crimson"]
        )
    # For table of special selection of records, create text labels adjacent
    # to points.
    if emphasis_label:
        for index, row in table_special.iterrows():
            # Determine position coordinates of label.
            abscissa_label_raw = row[column_abscissa]
            if (abscissa_label_raw >= center_abscissa):
                alignment_horizontal = "right"
                abscissa_label = (
                    abscissa_label_raw - (
                        0.01 * (maximum_abscissa - minimum_abscissa)
                    )
                )
            if (abscissa_label_raw < center_abscissa):
                alignment_horizontal = "left"
                abscissa_label = (
                    abscissa_label_raw + (
                        0.01 * (maximum_abscissa - minimum_abscissa)
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
                backgroundcolor=colors["white_faint"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_label_emphasis],
            )
            pass
        pass

    # Create legend or reference bar for color map.
    if (
        (column_response is not None) and
        (column_response != "") and
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
                backgroundcolor=colors["white"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_title_legend_bar],
            )
        else:
            bar.ax.set_ylabel(
                column_response,
                rotation=-90,
                va="bottom",
                labelpad=5, # 5
                alpha=1.0,
                backgroundcolor=colors["white"],
                color=colors["black"],
                fontproperties=fonts["properties"][size_title_legend_bar],
            )
            pass
        bar.ax.tick_params(
            axis="both",
            which="both", # major, minor, or both
            direction="out",
            length=7.5, # 5.0, 7.5
            width=3, # 2.5, 5.0
            color=colors["black"],
            pad=5, # 5, 7
            labelsize=fonts["values"][size_label_legend_bar]["size"],
            labelcolor=colors["black"],
        )
    elif (
        (column_response is not None) and
        (column_response != "") and
        (type_response is not None) and
        (type_response == "category")
    ):
        # Create legend.
        handles_legend = [matplotlib.patches.Patch(
            color=colors_groups[i],
            label=labels_groups[i]
        ) for i in range(count_category_response)]
        figure.legend(
            handles=handles_legend,
            loc="lower right",
            #bbox_to_anchor=(0.5, -0.05),
            #ncol=4,
            prop=fonts["properties"][size_label_legend_bar],
        )
        pass

    ##########
    # Return figure.
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
    values_ordinate = data_selection[ordinate].to_numpy()

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


def plot_scatter_qq_gwas(
    probabilities=None,
    title=None,
    label=None,
    fonts=None,
    colors=None,
):
    """
    Create a figure for a QQ scatter plot.

    The QQ scatter plots observed probabilities from a GWAS against expected
    probabilities.

    References:
    "https://github.com/ShujiaHuang/qmplot"
    "https://www.broadinstitute.org/diabetes-genetics-initiative/plotting-
    genome-wide-association-results"
    "https://www.broadinstitute.org/files/shared/diabetes/scandinavs/qqplot.R

    arguments:
        probabilities (object): NumPy array of probability values (p-values)
            from summary statistics for a GWAS already on negative, base-ten
            logarithm (-log10) scale
        title (str): title for figure
        label (str): label or title for plot area on figure
        fonts (dict<object>): references to definitions of font properties
        colors (dict<tuple>): references to definitions of color properties

    raises:

    returns:
        (object): figure object

    """

    # Suppress error warnings from NumPy.
    numpy.seterr(divide='ignore')
    numpy.seterr(invalid='ignore')

    # Prepare probability values.
    #probabilities_log = numpy.log10(probabilities_sort)
    #probabilities_neglog = numpy.multiply(probabilities_log, -1)
    probabilities = numpy.array(probabilities, dtype=numpy.float32)
    probabilities = numpy.copy(probabilities)
    probabilities = probabilities[~numpy.isnan(probabilities)]
    #probabilities = probabilities[(probabilities > 1.0e-37)]
    probabilities_sort_up = numpy.sort(
        probabilities,
        axis=0,
        kind="stable",
    )
    probabilities_sort_down = probabilities_sort_up[::-1]
    # Prepare expectation values.
    #count = int(len(probabilities_neglog))
    count = int(probabilities_sort_down.size)
    offset = 0.5 # value in (0, 1); normally 0.5; 0.375 for small sample size
    #expectations_range = range(1, (count + 1), 1)
    expectations_range = numpy.arange(1, (count + 1), 1) # range: 1 to count + 1
    expectations_offset = numpy.subtract(expectations_range, offset)
    expectations_ratio = numpy.divide(
        expectations_offset, ((count + 1) - (2*offset))
    )
    expectations_log = numpy.log10(expectations_ratio)
    expectations_neglog = numpy.multiply(expectations_log, -1)
    # Organize information for figure.
    values_abscissa = expectations_neglog # abscissa or X axis
    values_ordinate = probabilities_sort_down # ordinate or Y axis
    title_abscissa = "-log10(expected p-values)"
    title_ordinate = "-log10(observed p-values)"

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
        fontproperties=fonts["properties"]["three"]
    )
    axes.set_ylabel(
        ylabel=title_ordinate,
        labelpad=20,
        alpha=1.0,
        backgroundcolor=colors["white"],
        color=colors["black"],
        fontproperties=fonts["properties"]["three"]
    )
    axes.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5.0,
        width=3.0,
        color=colors["black"],
        pad=5,
        labelsize=fonts["values"]["three"]["size"],
        labelcolor=colors["black"]
    )
    # Plot points for values from each group.
    handle = axes.plot(
        values_abscissa,
        values_ordinate,
        linestyle="",
        marker="o",
        markersize=5, # 5, 15
        markeredgecolor=colors["blue_navy"],
        markerfacecolor=colors["blue_navy"]
    )

    # Return figure.
    return figure


# Dot, Point, Forest plots


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
    dimension. <-- outdated description? (TCW; 26 November 2024)

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


# Venn


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


# Write to file


def write_figure(
    figure=None,
    format=None,
    resolution=None,
    path_file=None,
):
    """
    Writes figure to file.

    arguments:
        figure (object): figure object
        format (str): format suffix, 'jpg', 'png', 'pdf', 'svg', etc
        resolution (int): dots per inch (dpi) density for raster image; set to
            '300' or '600' for raster or 'None' for format 'svg' or 'pdf'
        path_file (str): path to directory and file

    raises:

    returns:

    """

    # Write information to file.
    figure.savefig(
        path_file,
        format=format,
        dpi=resolution,
        facecolor="w",
        edgecolor="w",
        transparent=False
    )
    pass


def write_product_plot_figure(
    figure=None,
    format=None,
    resolution=None,
    name_file=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        figure (object): figure object to write to file
        format (str): format suffix, 'svg', 'pdf', 'png'
        resolution (int): dots per inch (dpi) density for raster image; set to
            '300' or '600' for raster or 'None' for format 'svg' or 'pdf'
        name_file (str): base name for file
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_file = os.path.join(
        path_directory, str(name_file + "." + str(format))
    )
    # Write information to file.
    write_figure(
        figure=figure,
        format=format,
        resolution=resolution, # dots per inch: 100, 150, 300, 600
        #transparent=True,
        path_file=path_file,
    )
    pass


def write_product_plots_parent_directory(
    pail_write=None,
    format=None,
    resolution=None,
    path_directory=None,
):
    """
    Writes product information to file.

    Keys of dictionary entries specify the name of the figure from which to
    derive the name of the file.
    Values of dictionary entries contain the figure objects.

    arguments:
        pail_write (dict<object>): information to write to file
        format (str): format suffix, 'svg', 'pdf', 'png'
        resolution (int): dots per inch (dpi) density for raster image; set to
            '300' or '600' for raster or 'None' for format 'svg' or 'pdf'
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    # Iterate across figures.
    for name_file in pail_write.keys():
        # Write figure object to file within parent directory.
        write_product_plot_figure(
            figure=pail_write[name_file],
            format=format,
            resolution=resolution,
            name_file=name_file,
            path_directory=path_directory,
        )
        pass
    pass


def write_product_plots_child_directories(
    pail_write=None,
    format=None,
    resolution=None,
    path_parent=None,
):
    """
    Writes product information to file.

    First dictionary tier names the child directory.
    Second dictionary tier names the file.

    arguments:
        pail_write (dict<dict<object>>): information to write to file
        format (str): format suffix, 'svg', 'pdf', 'png'
        resolution (int): dots per inch (dpi) density for raster image; set to
            '300' or '600' for raster or 'None' for format 'svg' or 'pdf'
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
        putly.create_directories(
            path=path_child
        )
        # Iterate across figure objects.
        # Write chart objects to file in child directory.
        write_product_plots_parent_directory(
            pail_write=pail_write[name_directory],
            format=format,
            resolution=resolution,
            path_directory=path_child,
        )
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
        putly.create_directories(
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
