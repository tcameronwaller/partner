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
# Date, initialization: 2 July 2025
# Date, review or revision: 5 March 2026
# Date, review or revision: 16 July 2025
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


def parse_text_parameters(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_effects=None,
    path_file_source_list_emphasis=None,
    identifiers_emphasis=None,
    column_effect_identifier=None,
    column_effect_estimate=None,
    column_effect_error=None,
    column_effect_p=None,
    column_effect_q=None,
    column_feature_identifier=None,
    column_feature_name=None,
    column_plot_abscissa=None,
    column_plot_ordinate=None,
    column_significance=None,
    name_chart=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    threshold_abscissa=None,
    threshold_ordinate=None,
    threshold_significance=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    size_marker=None,
    size_marker_emphasis=None,
    size_edge_marker=None,
    size_font_title_abscissa=None,
    size_font_title_ordinate=None,
    size_font_label_abscissa=None,
    size_font_label_ordinate=None,
    size_font_label_emphasis=None,
    size_font_label_count=None,
    color_fill_markers_subtle=None,
    color_fill_markers_significant=None,
    color_fill_markers_emphasis=None,
    color_edge_markers=None,
    color_edge_markers_emphasis=None,
    line_threshold_abscissa=None,
    line_threshold_ordinate=None,
    label_emphasis=None,
    label_count=None,
    report=None,
):
    """
    Parse parameters from text.

    Date, revision or review: 5 March 2026

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
    pail["path_file_source_table_effects"] = str(
        path_file_source_table_effects
    ).strip()
    pail["path_file_source_list_emphasis"] = str(
        path_file_source_list_emphasis
    ).strip()

    # Lists, simple text.
    # Iterate on individual parameters for names and categories.
    lists = {
        "identifiers_emphasis": identifiers_emphasis,
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

    # Names and categories.
    # It is problematic to pass any white space in parameters from a script in
    # Bash. Designate the hash symbol "#" as a substitute for white space.
    # It is also problematic to pass an empty string in parameters from a
    # script in Bash. Designate the word "none" as a substitute for missing or
    # empty.
    # Iterate on individual names that could be empty or missing.
    names_categories = {
        "column_effect_identifier": column_effect_identifier,
        "column_effect_estimate": column_effect_estimate,
        "column_effect_error": column_effect_error,
        "column_effect_p": column_effect_p,
        "column_effect_q": column_effect_q,
        "column_feature_identifier": column_feature_identifier,
        "column_feature_name": column_feature_name,
        "column_plot_abscissa": column_plot_abscissa,
        "column_plot_ordinate": column_plot_ordinate,
        "column_significance": column_significance,
        "name_chart": name_chart,
        "title_chart": title_chart,
        "title_abscissa": title_abscissa,
        "title_ordinate": title_ordinate,
        "size_font_title_abscissa": size_font_title_abscissa,
        "size_font_title_ordinate": size_font_title_ordinate,
        "size_font_label_abscissa": size_font_label_abscissa,
        "size_font_label_ordinate": size_font_label_ordinate,
        "size_font_label_emphasis": size_font_label_emphasis,
        "size_font_label_count": size_font_label_count,
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
    pail["threshold_abscissa"] = float(str(threshold_abscissa).strip())
    pail["threshold_ordinate"] = float(str(threshold_ordinate).strip())
    pail["threshold_significance"] = float(str(threshold_significance).strip())
    pail["minimum_abscissa"] = float(str(minimum_abscissa).strip())
    pail["maximum_abscissa"] = float(str(maximum_abscissa).strip())
    pail["minimum_ordinate"] = float(str(minimum_ordinate).strip())
    pail["maximum_ordinate"] = float(str(maximum_ordinate).strip())
    pail["size_marker"] = float(str(size_marker).strip())
    pail["size_marker_emphasis"] = float(str(size_marker_emphasis).strip())
    pail["size_edge_marker"] = float(str(size_edge_marker).strip())

    # Iterate on individual colors.
    colors = {
        "color_fill_markers_subtle": color_fill_markers_subtle,
        "color_fill_markers_significant": color_fill_markers_significant,
        "color_fill_markers_emphasis": color_fill_markers_emphasis,
        "color_edge_markers": color_edge_markers,
        "color_edge_markers_emphasis": color_edge_markers_emphasis,
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
        "line_threshold_abscissa": line_threshold_abscissa,
        "line_threshold_ordinate": line_threshold_ordinate,
        "label_emphasis": label_emphasis,
        "label_count": label_count,
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
        print("module: plot_chart_volcano_effect.py")
        print("function: parse_text_parameters()")
        putly.print_terminal_partition(level=5)
        print("parameters:")
        print(pail)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_effects=None,
    path_file_source_list_emphasis=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Date, revision or review: 5 March 2026

    arguments:
        path_file_source_table_effects (str): path to source file
        path_file_source_list_emphasis (str): path to source file
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Determine whether paths point to a directory or file that exist.
    #existence_directory = os.path.exists(path_directory)
    existence_effects = os.path.exists(path_file_source_table_effects)
    existence_features = os.path.exists(path_file_source_list_emphasis)

    # Bundle information.
    pail = dict()

    # Read information from file.
    if (existence_effects):
        pail["table_effects"] = pandas.read_csv(
            path_file_source_table_effects,
            sep="\t",
            header=0,
            #dtype=types_columns,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        pail["table_effects"] = None
        pass
    if (existence_features):
        # Read information from file.
        pail["features_emphasis"] = putly.read_file_text_list(
            path_file=path_file_source_list_emphasis,
            delimiter="\n",
            unique=True,
        )
    else:
        pail["features_emphasis"] = list()
        pass

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: plot_chart_volcano_effect.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("table of effects:")
        print(pail["table_effects"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def plot_scatter_fold_change_volcano(
    table=None,
    features_emphasis=None,
    column_effect_identifier=None,
    column_effect_estimate=None,
    column_effect_error=None,
    column_effect_p=None,
    column_effect_q=None,
    column_feature_identifier=None,
    column_feature_name=None,
    column_plot_abscissa=None,
    column_plot_ordinate=None,
    column_significance=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    threshold_abscissa=None,
    threshold_ordinate=None,
    threshold_significance=None,
    aspect=None,
    minimum_abscissa=None,
    center_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    size_marker=None,
    size_marker_emphasis=None,
    size_edge_marker=None,
    size_line_threshold=None,
    size_font_title_abscissa=None,
    size_font_title_ordinate=None,
    size_font_label_abscissa=None,
    size_font_label_ordinate=None,
    size_font_label_emphasis=None,
    size_font_label_count=None,
    color_fill_markers_subtle=None,
    color_fill_markers_significant=None,
    color_fill_markers_emphasis=None,
    color_edge_markers=None,
    color_edge_markers_emphasis=None,
    fonts=None,
    line_threshold_abscissa=None,
    line_threshold_ordinate=None,
    label_emphasis=None,
    label_count=None,
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

    Date, revision or review: 5 March 2026
    Date, revision or review: 3 October 2024

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
        title_abscissa (str): title for abscissa horizontal axis
        title_ordinate (str): title for ordinate vertical axis
        lines_thresholds (bool): whether to draw lines to represent thresholds
        label_emphasis (bool): whether to create text labels adjacent to
            the special selection of points for special emphasis
        label_count (bool): whether to create text labels on chart to report
            counts of fold changes that pass thresholds
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
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: plot_chart_forest_effect.py")
        print("function: plot_scatter_fold_change_volcano()")
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Organize information for chart.

    # Copy information.
    table = table.copy(deep=True)
    columns_available = copy.deepcopy(table.columns.to_list())

    ##########
    # Organize information for chart.

    # Filter and sort columns in table.
    columns_sequence = [
        column_feature_identifier,
        column_feature_name,
        column_effect_identifier,
        column_effect_estimate,
        column_effect_error,
        column_effect_p,
        column_effect_q,
        column_plot_abscissa,
        column_plot_ordinate,
        column_significance,
    ]
    columns_sequence = putly.collect_unique_items(
        items=columns_sequence,
    )
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=columns_sequence,
        report=report,
    )
    # Filter rows in table.
    table.dropna(
        axis="index",
        how="any",
        inplace=True,
    )
    limits = [
        {
            "type": "low",
            "value": minimum_abscissa,
            "column": column_plot_abscissa,
        },
        {
            "type": "high",
            "value": maximum_abscissa,
            "column": column_plot_abscissa,
        },
        {
            "type": "low",
            "value": minimum_ordinate,
            "column": column_plot_ordinate,
        },
        {
            "type": "high",
            "value": maximum_ordinate,
            "column": column_plot_ordinate,
        },
    ]
    for limit in limits:
        if (limit["value"] is not None):
            if (limit["type"] == "low"):
                if True:
                    # Clip values in table's column.
                    table[str(limit["column"])] = (
                        table[str(limit["column"])].clip(
                            lower=limit["value"],
                            #inplace=True,
                    ))
                else:
                    # Filter rows in table.
                    table = table.loc[
                        (table[str(limit["column"])] >= limit["value"])#, :
                    ].copy(deep=True)
                    pass
            elif (limit["type"] == "high"):
                if True:
                    # Clip values in table's column.
                    table[str(limit["column"])] = (
                        table[str(limit["column"])].clip(
                            upper=limit["value"],
                            #inplace=True,
                    ))
                else:
                    # Filter rows in table.
                    table = table.loc[
                        (table[str(limit["column"])] < limit["value"])#, :
                    ].copy(deep=True)
                    pass
                pass
            pass
        pass
    pass
    # Filter rows in table for segregation of values by thresholds.
    pail_threshold = porg.segregate_effects_by_thresholds(
        table=table,
        column_effect=str(column_effect_estimate),
        column_significance=str(column_significance),
        threshold_effect=threshold_abscissa,
        threshold_significance=threshold_significance,
        report=False,
    )

    # Extract information about selection of fold changes for special emphasis.
    #table_bore = table.loc[
    #    ~table[column_identifier].isin(identifiers_emphasis), :
    #].copy(deep=True)
    if True:
        # Extract emphasis from whole table.
        table_emphasis = table.loc[
            table[column_feature_identifier].isin(features_emphasis), :
        ].copy(deep=True)
    else:
        # Extract emphasis from table of effects that pass thresholds.
        table_source = pail_threshold["table_pass_any"]
        table_emphasis = table_source.loc[
            table_source[column_feature_identifier].isin(features_emphasis), :
        ].copy(deep=True)
        pass
    # Extract counts.
    count_pass = int(pail_threshold["table_pass_any"].shape[0])
    count_pass_up = int(pail_threshold["table_pass_up"].shape[0])
    count_pass_down = int(pail_threshold["table_pass_down"].shape[0])
    count_fail = int(pail_threshold["table_fail"].shape[0])

    # Report.
    if report:
        putly.print_terminal_partition(level=5)
        print(
            "count of effects that pass segregation thresholds: " +
            str(count_pass)
        )
        print("count ups: " + str(count_pass_up))
        print("count downs: " + str(count_pass_down))
        print(
            "count of effects that fail segregation thresholds: " +
            str(count_fail)
        )
        putly.print_terminal_partition(level=5)
        print("table of fold changes for special emphasis:")
        print(table_emphasis)
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Organize parameters for visual representation.
    # Colors.
    # Individual colors.
    # Bundle information.
    pail_colors = dict()
    # Iterate on individual parameters for red, blue, green, and alpha channels
    # in colors.
    colors = {
        "fill_subtle": color_fill_markers_subtle,
        "fill_significant": color_fill_markers_significant,
        "fill_emphasis": color_fill_markers_emphasis,
        "edge": color_edge_markers,
        "edge_emphasis": color_edge_markers_emphasis,
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
    #pail_colors["fill_subtle"]
    #pail_colors["fill_significant"]
    #pail_colors["fill_emphasis"]
    #pail_colors["edge"]
    #pail_colors["edge_emphasis"]

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
    if (minimum_abscissa is not None):
        axes.set_xlim(xmin=minimum_abscissa)
    if (maximum_abscissa is not None):
        axes.set_xlim(xmax=maximum_abscissa)
    if (minimum_ordinate is not None):
        axes.set_ylim(ymin=minimum_ordinate)
    if (maximum_ordinate is not None):
        axes.set_ylim(ymax=maximum_ordinate)

    # Include title label on chart.
    if len(title_chart) > 0:
        axes.set_title(
            title_chart,
            fontproperties=fonts["properties"][size_font_title_chart],
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
            fontproperties=fonts["properties"][size_font_title_abscissa]
        )
    if (len(title_ordinate) > 0):
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=15,
            alpha=1.0,
            backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
            color=matplotlib.colors.to_rgba("black", 1.0),
            fontproperties=fonts["properties"][size_font_title_ordinate]
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
        length=5.0, # 5.0
        width=2.5, # 3.0, 5.0
        pad=10.0, # 5.0, 7.5
        color=matplotlib.colors.to_rgba("black", 1.0),
        labelcolor=matplotlib.colors.to_rgba("black", 1.0),
    )
    axes.tick_params(
        axis="x",
        which="both",
        labelsize=fonts["values"][size_font_label_abscissa]["size"],
    )
    axes.tick_params(
        axis="y",
        which="both",
        labelsize=fonts["values"][size_font_label_ordinate]["size"],
    )
    # Keep axes, ticks, and labels, but remove border.
    # ["left", "top", "right", "bottom",]
    for position in ["top", "right",]:
        matplotlib.pyplot.gca().spines[position].set_visible(False)

    # Create lines to represent threshold values.
    if (line_threshold_abscissa and (count_pass > 0)):
        axes.axvline(
            x=threshold_abscissa,
            ymin=minimum_ordinate,
            ymax=maximum_ordinate,
            alpha=1.0,
            color=matplotlib.colors.to_rgba("black", 1.0),
            linestyle="--",
            linewidth=size_line_threshold,
        )
        if (threshold_abscissa != 0):
            axes.axvline(
                x=(-1*threshold_abscissa),
                ymin=minimum_ordinate,
                ymax=maximum_ordinate,
                alpha=1.0,
                color=matplotlib.colors.to_rgba("black", 1.0),
                linestyle="--",
                linewidth=size_line_threshold,
            )
            pass
    if (line_threshold_ordinate and (count_pass > 0)):
        # Determine minimal value of significance for line.
        threshold_significance_scale = numpy.nanmin(
            pail_threshold["table_pass_any"][column_plot_ordinate].to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            )
        )
        axes.axhline(
            y=threshold_significance_scale, # or threshold_ordinate
            xmin=minimum_abscissa,
            xmax=maximum_abscissa,
            alpha=1.0,
            color=matplotlib.colors.to_rgba("black", 1.0),
            linestyle="--",
            linewidth=size_line_threshold,
        )
        pass

    ##########
    # Represent information on the chart figure object.
    # Plot points for values from each group.
    handle = axes.plot(
        pail_threshold["table_fail"][column_plot_abscissa].values,
        pail_threshold["table_fail"][column_plot_ordinate].values,
        linestyle="",
        linewidth=size_edge_marker,
        marker="o",
        markersize=(size_marker * 0.75),
        markeredgecolor=pail_colors["edge"],
        markerfacecolor=pail_colors["fill_subtle"],
    )
    handle = axes.plot(
        pail_threshold["table_pass_any"][column_plot_abscissa].values,
        pail_threshold["table_pass_any"][column_plot_ordinate].values,
        linestyle="",
        linewidth=size_edge_marker,
        marker="o",
        markersize=size_marker, # 5, 7.5, 10
        markeredgecolor=pail_colors["edge"],
        markerfacecolor=pail_colors["fill_significant"],
    )
    handle = axes.plot(
        table_emphasis[column_plot_abscissa].values,
        table_emphasis[column_plot_ordinate].values,
        linestyle="",
        linewidth=size_edge_marker,
        marker="o",
        markersize=size_marker_emphasis, # 5, 7.5, 10
        markeredgecolor=pail_colors["edge_emphasis"],
        markerfacecolor=pail_colors["fill_emphasis"],
    )

    # Plot labels to report counts of fold changes that pass thresholds.
    if label_count:
        # Determine position coordinates of labels.
        center_abscissa_down = (
            minimum_abscissa + (
                abs(minimum_abscissa - (-1*threshold_abscissa))/2.5
        ))
        center_abscissa_up = (
            maximum_abscissa - ((maximum_abscissa - threshold_abscissa)/2.5)
        )
        center_ordinate = ((maximum_ordinate - 0)/1.5)
        #count_pass_up
        #count_pass_down
        # Create labels on chart.
        axes.text(
            center_abscissa_up,
            center_ordinate,
            str("count: " + str(count_pass_up)),
            horizontalalignment="right",
            verticalalignment="center",
            backgroundcolor=matplotlib.colors.to_rgba("white", 0.33),
            color=matplotlib.colors.to_rgba("gray", 1.0),
            fontproperties=fonts["properties"][size_font_label_count],
        )
        axes.text(
            center_abscissa_down,
            center_ordinate,
            str("count: " + str(count_pass_down)),
            horizontalalignment="left",
            verticalalignment="center",
            backgroundcolor=matplotlib.colors.to_rgba("white", 0.33),
            color=matplotlib.colors.to_rgba("gray", 1.0),
            fontproperties=fonts["properties"][size_font_label_count],
        )
        pass

    # Plot labels adjacent to the selection of points for special emphasis.
    # Align bottom left of label with an offset of 1% of the chart coordinate
    # dimensions.
    if label_emphasis:
        for index, row in table_emphasis.iterrows():
            # Determine position coordinates of label.
            abscissa_raw = row[column_plot_abscissa]
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
            ordinate = row[column_plot_ordinate]
            # Create label on chart.
            axes.text(
                abscissa,
                ordinate,
                str(row[column_feature_name]),
                horizontalalignment=alignment_horizontal,
                verticalalignment="center",
                backgroundcolor=matplotlib.colors.to_rgba("white", 0.33),
                color=matplotlib.colors.to_rgba("black", 1.0),
                fontproperties=fonts["properties"][size_font_label_emphasis],
            )
            pass
        pass

    ##########
    # Return figure.
    return figure


def create_write_plot_chart_volcano(
    path_directory_parent=None,
    name_chart=None,
    table=None,
    features_emphasis=None,
    column_effect_identifier=None,
    column_effect_estimate=None,
    column_effect_error=None,
    column_effect_p=None,
    column_effect_q=None,
    column_feature_identifier=None,
    column_feature_name=None,
    column_plot_abscissa=None,
    column_plot_ordinate=None,
    column_significance=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    threshold_abscissa=None,
    threshold_ordinate=None,
    threshold_significance=None,
    minimum_abscissa=None,
    center_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    size_marker=None,
    size_marker_emphasis=None,
    size_edge_marker=None,
    size_font_title_abscissa=None,
    size_font_title_ordinate=None,
    size_font_label_abscissa=None,
    size_font_label_ordinate=None,
    size_font_label_emphasis=None,
    size_font_label_count=None,
    color_fill_markers_subtle=None,
    color_fill_markers_significant=None,
    color_fill_markers_emphasis=None,
    color_edge_markers=None,
    color_edge_markers_emphasis=None,
    line_threshold_abscissa=None,
    line_threshold_ordinate=None,
    label_emphasis=None,
    label_count=None,
    report=None,
):
    """
    Create and write to file a plot chart of type volcano.

    Date, revision or review: 5 March 2026
    Date, revision or review: 16 July 2025

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
        name_chart (str): name of product file
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
    #colors = pplot.define_color_properties()

    # Create plot chart figure.
    figure = plot_scatter_fold_change_volcano(
        table=table,
        features_emphasis=features_emphasis,
        column_effect_identifier=column_effect_identifier,
        column_effect_estimate=column_effect_estimate,
        column_effect_error=column_effect_error,
        column_effect_p=column_effect_p,
        column_effect_q=column_effect_q,
        column_feature_identifier=column_feature_identifier,
        column_feature_name=column_feature_name,
        column_plot_abscissa=column_plot_abscissa,
        column_plot_ordinate=column_plot_ordinate,
        column_significance=column_significance,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        threshold_abscissa=threshold_abscissa,
        threshold_ordinate=threshold_ordinate,
        threshold_significance=threshold_significance,
        aspect="portrait",
        minimum_abscissa=minimum_abscissa,
        center_abscissa=center_abscissa,
        maximum_abscissa=maximum_abscissa,
        minimum_ordinate=minimum_ordinate,
        maximum_ordinate=maximum_ordinate,
        size_marker=size_marker,
        size_marker_emphasis=size_marker_emphasis,
        size_edge_marker=size_edge_marker,
        size_line_threshold=3.5,
        size_font_title_abscissa=size_font_title_abscissa,
        size_font_title_ordinate=size_font_title_ordinate,
        size_font_label_abscissa=size_font_label_abscissa,
        size_font_label_ordinate=size_font_label_ordinate,
        size_font_label_emphasis=size_font_label_emphasis,
        size_font_label_count=size_font_label_count,
        color_fill_markers_subtle=color_fill_markers_subtle,
        color_fill_markers_significant=color_fill_markers_significant,
        color_fill_markers_emphasis=color_fill_markers_emphasis,
        color_edge_markers=color_edge_markers,
        color_edge_markers_emphasis=color_edge_markers_emphasis,
        fonts=fonts,
        line_threshold_abscissa=line_threshold_abscissa,
        line_threshold_ordinate=line_threshold_ordinate,
        label_emphasis=label_emphasis,
        label_count=label_count,
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
    path_file_source_table_effects=None,
    path_file_source_list_emphasis=None,
    identifiers_emphasis=None,
    column_effect_identifier=None,
    column_effect_estimate=None,
    column_effect_error=None,
    column_effect_p=None,
    column_effect_q=None,
    column_feature_identifier=None,
    column_feature_name=None,
    column_plot_abscissa=None,
    column_plot_ordinate=None,
    column_significance=None,
    name_chart=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    threshold_abscissa=None,
    threshold_ordinate=None,
    threshold_significance=None,
    minimum_abscissa=None,
    maximum_abscissa=None,
    minimum_ordinate=None,
    maximum_ordinate=None,
    size_marker=None,
    size_marker_emphasis=None,
    size_edge_marker=None,
    size_font_title_abscissa=None,
    size_font_title_ordinate=None,
    size_font_label_abscissa=None,
    size_font_label_ordinate=None,
    size_font_label_emphasis=None,
    size_font_label_count=None,
    color_fill_markers_subtle=None,
    color_fill_markers_significant=None,
    color_fill_markers_emphasis=None,
    color_edge_markers=None,
    color_edge_markers_emphasis=None,
    line_threshold_abscissa=None,
    line_threshold_ordinate=None,
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
        name_chart (str): name of product file

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
        path_directory_dock=path_directory_dock,
        path_directory_dock_pail=path_directory_dock_pail,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_file_source_table_effects=path_file_source_table_effects,
        path_file_source_list_emphasis=path_file_source_list_emphasis,
        identifiers_emphasis=identifiers_emphasis,
        column_effect_identifier=column_effect_identifier,
        column_effect_estimate=column_effect_estimate,
        column_effect_error=column_effect_error,
        column_effect_p=column_effect_p,
        column_effect_q=column_effect_q,
        column_feature_identifier=column_feature_identifier,
        column_feature_name=column_feature_name,
        column_plot_abscissa=column_plot_abscissa,
        column_plot_ordinate=column_plot_ordinate,
        column_significance=column_significance,
        name_chart=name_chart,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        threshold_abscissa=threshold_abscissa,
        threshold_ordinate=threshold_ordinate,
        threshold_significance=threshold_significance,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        minimum_ordinate=minimum_ordinate,
        maximum_ordinate=maximum_ordinate,
        size_marker=size_marker,
        size_marker_emphasis=size_marker_emphasis,
        size_edge_marker=size_edge_marker,
        size_font_title_abscissa=size_font_title_abscissa,
        size_font_title_ordinate=size_font_title_ordinate,
        size_font_label_abscissa=size_font_label_abscissa,
        size_font_label_ordinate=size_font_label_ordinate,
        size_font_label_emphasis=size_font_label_emphasis,
        size_font_label_count=size_font_label_count,
        color_fill_markers_subtle=color_fill_markers_subtle,
        color_fill_markers_significant=color_fill_markers_significant,
        color_fill_markers_emphasis=color_fill_markers_emphasis,
        color_edge_markers=color_edge_markers,
        color_edge_markers_emphasis=color_edge_markers_emphasis,
        line_threshold_abscissa=line_threshold_abscissa,
        line_threshold_ordinate=line_threshold_ordinate,
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
        path_directory_dock=pail_parameters["path_directory_dock"],
        path_directory_dock_pail=pail_parameters["path_directory_dock_pail"],
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_file_source_table_effects=(
            pail_parameters["path_file_source_table_effects"]
        ),
        path_file_source_list_emphasis=(
            pail_parameters["path_file_source_list_emphasis"]
        ),
        report=pail_parameters["report"],
    )
    #pail_source["table_effects"]
    #pail_source["features_emphasis"]

    # Combine identifiers of features for emphasis.
    features_emphasis = copy.deepcopy(pail_source["features_emphasis"])
    features_emphasis.extend(pail_parameters["identifiers_emphasis"])

    # Define paths to directories.
    path_directory_chart = os.path.join(
        pail_parameters["path_directory_product"],
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_chart,
    )

    ##########
    # Create and write to file plot chart volcano.
    create_write_plot_chart_volcano(
        path_directory_parent=path_directory_chart,
        name_chart=pail_parameters["name_chart"],
        table=pail_source["table_effects"],
        features_emphasis=features_emphasis,
        column_effect_identifier=pail_parameters["column_effect_identifier"],
        column_effect_estimate=pail_parameters["column_effect_estimate"],
        column_effect_error=pail_parameters["column_effect_error"],
        column_effect_p=pail_parameters["column_effect_p"],
        column_effect_q=pail_parameters["column_effect_q"],
        column_feature_identifier=pail_parameters["column_feature_identifier"],
        column_feature_name=pail_parameters["column_feature_name"],
        column_plot_abscissa=pail_parameters["column_plot_abscissa"],
        column_plot_ordinate=pail_parameters["column_plot_ordinate"],
        column_significance=pail_parameters["column_significance"],
        title_chart=pail_parameters["title_chart"],
        title_abscissa=pail_parameters["title_abscissa"],
        title_ordinate=pail_parameters["title_ordinate"],
        threshold_abscissa=pail_parameters["threshold_abscissa"],
        threshold_ordinate=pail_parameters["threshold_ordinate"],
        threshold_significance=pail_parameters["threshold_significance"],
        minimum_abscissa=pail_parameters["minimum_abscissa"],
        center_abscissa=0.0,
        maximum_abscissa=pail_parameters["maximum_abscissa"],
        minimum_ordinate=pail_parameters["minimum_ordinate"],
        maximum_ordinate=pail_parameters["maximum_ordinate"],
        size_marker=pail_parameters["size_marker"],
        size_marker_emphasis=pail_parameters["size_marker_emphasis"],
        size_edge_marker=pail_parameters["size_edge_marker"],
        size_font_title_abscissa=pail_parameters["size_font_title_abscissa"],
        size_font_title_ordinate=pail_parameters["size_font_title_ordinate"],
        size_font_label_abscissa=pail_parameters["size_font_label_abscissa"],
        size_font_label_ordinate=pail_parameters["size_font_label_ordinate"],
        size_font_label_emphasis=pail_parameters["size_font_label_emphasis"],
        size_font_label_count=pail_parameters["size_font_label_count"],
        color_fill_markers_subtle=pail_parameters["color_fill_markers_subtle"],
        color_fill_markers_significant=(
            pail_parameters["color_fill_markers_significant"]
        ),
        color_fill_markers_emphasis=(
            pail_parameters["color_fill_markers_emphasis"]
        ),
        color_edge_markers=pail_parameters["color_edge_markers"],
        color_edge_markers_emphasis=(
            pail_parameters["color_edge_markers_emphasis"]
        ),
        line_threshold_abscissa=pail_parameters["line_threshold_abscissa"],
        line_threshold_ordinate=pail_parameters["line_threshold_ordinate"],
        label_emphasis=pail_parameters["label_emphasis"],
        label_count=pail_parameters["label_count"],
        report=pail_parameters["report"],
    )


    pass


# Execute program process in Python.


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_dock = sys.argv[1]
    path_directory_dock_pail = sys.argv[2]
    path_directory_source = sys.argv[3]
    path_directory_product = sys.argv[4]
    path_file_source_table_effects = sys.argv[5]
    path_file_source_list_emphasis = sys.argv[6]
    identifiers_emphasis = sys.argv[7]
    column_effect_identifier = sys.argv[8]
    column_effect_estimate = sys.argv[9]
    column_effect_error = sys.argv[10]
    column_effect_p = sys.argv[11]
    column_effect_q = sys.argv[12]
    column_feature_identifier = sys.argv[13]
    column_feature_name = sys.argv[14]
    column_plot_abscissa = sys.argv[15]
    column_plot_ordinate = sys.argv[16]
    column_significance = sys.argv[17]
    name_chart = sys.argv[18]
    title_chart = sys.argv[19]
    title_abscissa = sys.argv[20]
    title_ordinate = sys.argv[21]
    threshold_abscissa = sys.argv[22]
    threshold_ordinate = sys.argv[23]
    threshold_significance = sys.argv[24]
    minimum_abscissa = sys.argv[25]
    maximum_abscissa = sys.argv[26]
    minimum_ordinate = sys.argv[27]
    maximum_ordinate = sys.argv[28]
    size_marker = sys.argv[29]
    size_marker_emphasis = sys.argv[30]
    size_edge_marker = sys.argv[31]
    size_font_title_abscissa = sys.argv[32]
    size_font_title_ordinate = sys.argv[33]
    size_font_label_abscissa = sys.argv[34]
    size_font_label_ordinate = sys.argv[35]
    size_font_label_emphasis = sys.argv[36]
    size_font_label_count = sys.argv[37]
    color_fill_markers_subtle = sys.argv[38]
    color_fill_markers_significant = sys.argv[39]
    color_fill_markers_emphasis = sys.argv[40]
    color_edge_markers = sys.argv[41]
    color_edge_markers_emphasis = sys.argv[42]
    line_threshold_abscissa = sys.argv[43]
    line_threshold_ordinate = sys.argv[44]
    label_emphasis = sys.argv[45]
    label_count = sys.argv[46]
    report = sys.argv[47]

    # Call function for procedure.
    execute_procedure(
        path_directory_dock=path_directory_dock,
        path_directory_dock_pail=path_directory_dock_pail,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_file_source_table_effects=path_file_source_table_effects,
        path_file_source_list_emphasis=path_file_source_list_emphasis,
        identifiers_emphasis=identifiers_emphasis,
        column_effect_identifier=column_effect_identifier,
        column_effect_estimate=column_effect_estimate,
        column_effect_error=column_effect_error,
        column_effect_p=column_effect_p,
        column_effect_q=column_effect_q,
        column_feature_identifier=column_feature_identifier,
        column_feature_name=column_feature_name,
        column_plot_abscissa=column_plot_abscissa,
        column_plot_ordinate=column_plot_ordinate,
        column_significance=column_significance,
        name_chart=name_chart,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        threshold_abscissa=threshold_abscissa,
        threshold_ordinate=threshold_ordinate,
        threshold_significance=threshold_significance,
        minimum_abscissa=minimum_abscissa,
        maximum_abscissa=maximum_abscissa,
        minimum_ordinate=minimum_ordinate,
        maximum_ordinate=maximum_ordinate,
        size_marker=size_marker,
        size_marker_emphasis=size_marker_emphasis,
        size_edge_marker=size_edge_marker,
        size_font_title_abscissa=size_font_title_abscissa,
        size_font_title_ordinate=size_font_title_ordinate,
        size_font_label_abscissa=size_font_label_abscissa,
        size_font_label_ordinate=size_font_label_ordinate,
        size_font_label_emphasis=size_font_label_emphasis,
        size_font_label_count=size_font_label_count,
        color_fill_markers_subtle=color_fill_markers_subtle,
        color_fill_markers_significant=color_fill_markers_significant,
        color_fill_markers_emphasis=color_fill_markers_emphasis,
        color_edge_markers=color_edge_markers,
        color_edge_markers_emphasis=color_edge_markers_emphasis,
        line_threshold_abscissa=line_threshold_abscissa,
        line_threshold_ordinate=line_threshold_ordinate,
        label_emphasis=label_emphasis,
        label_count=label_count,
        report=report,
    )

    pass



#
