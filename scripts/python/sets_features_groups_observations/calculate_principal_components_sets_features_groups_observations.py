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
# Date, first execution: 30 September 2025
# Date, last execution or modification: 30 September 2025
# Review: TCW; 30 September 2025
################################################################################
# Note

# The specialty of this Python script is to manage a general and versatile
# design of principal components analysis.



# TODO: TCW; 1 October 2025
# Implement heatmap (filter to < 100 features and PCs) for PC loadings




##########
# Note:

# Note: TCW; 30 September 2025
# Features of interest may optionally exist within the "table_observations" or
# within the "table_signals", with or without the need for a transposition
# before merging with "table_observations".


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
#import matplotlib_venn
#import seaborn
#import sklearn

# Custom
import partner.utility as putly
#import partner.parallelization as prall
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
import partner.decomposition as pdcmp
import partner.plot as pplot

import utility_special as sutly

#dir()
#importlib.reload()

###############################################################################
# Functionality


# Organize raw parameters.


def parse_text_parameters(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_features=None,
    path_file_source_table_observations=None,
    path_file_source_table_signals=None,
    path_file_source_table_sets_features=None,
    path_file_source_table_groups_observations=None,
    column_identifier_feature=None,
    column_name_feature=None,
    column_identifier_observation=None,
    column_identifier_signal=None,
    count_components=None,
    transpose_table_signals=None,
    allow_replicate_observations=None,
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
    pail["path_directory_source"] = str(path_directory_source).strip()
    pail["path_directory_product"] = str(path_directory_product).strip()
    pail["path_directory_dock"] = str(path_directory_dock).strip()

    # Paths to files.
    pail["path_file_source_table_features"] = str(
        path_file_source_table_features
    ).strip()
    pail["path_file_source_table_observations"] = str(
        path_file_source_table_observations
    ).strip()
    pail["path_file_source_table_signals"] = str(
        path_file_source_table_signals
    ).strip()
    pail["path_file_source_table_sets_features"] = str(
        path_file_source_table_sets_features
    ).strip()
    pail["path_file_source_table_groups_observations"] = str(
        path_file_source_table_groups_observations
    ).strip()

    # Names of columns.
    pail["column_identifier_feature"] = str(
        column_identifier_feature
    ).strip()
    pail["column_name_feature"] = str(
        column_name_feature
    ).strip()
    pail["column_identifier_observation"] = str(
        column_identifier_observation
    ).strip()
    pail["column_identifier_signal"] = str(
        column_identifier_signal
    ).strip()

    # Number.
    pail["count_components"] = int(count_components)

    # Boolean, true or false.
    if (
        (transpose_table_signals is not None) and
        (str(transpose_table_signals) != "") and
        (str(transpose_table_signals) != "none") and
        (str(transpose_table_signals) == "true")
    ):
        pail["transpose_table_signals"] = True
    else:
        pail["transpose_table_signals"] = False
        pass
    if (
        (allow_replicate_observations is not None) and
        (str(allow_replicate_observations) != "") and
        (str(allow_replicate_observations) != "none") and
        (str(allow_replicate_observations) == "true")
    ):
        pail["allow_replicate_observations"] = True
    else:
        pail["allow_replicate_observations"] = False
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

    # Report.
    if pail["report"]:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
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
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_features=None,
    path_file_source_table_observations=None,
    path_file_source_table_signals=None,
    path_file_source_table_sets_features=None,
    path_file_source_table_groups_observations=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 29 September 2025

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files

        ...

        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Determine whether paths point to a directory or file that exist.
    #existence_directory = os.path.exists(path_directory)
    existence_file_signals = os.path.exists(path_file_source_table_signals)

    # Bundle information.
    pail = dict()

    # Read information from file.
    pail["table_features"] = pandas.read_csv(
        path_file_source_table_features,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_observations"] = pandas.read_csv(
        path_file_source_table_observations,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    if (existence_file_signals):
        pail["table_signals"] = pandas.read_csv(
            path_file_source_table_signals,
            sep="\t",
            header=0,
            #dtype=types_columns,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        pail["table_signals"] = None
        pass
    pail["table_sets"] = pandas.read_csv(
        path_file_source_table_sets_features,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    pail["table_groups"] = pandas.read_csv(
        path_file_source_table_groups_observations,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
        )
        print(str("module: " + module))
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def filter_combine_features_observations(
    table_main=None,
    table_supplement=None,
    column_identifier_feature=None,
    column_identifier_observation=None,
    column_identifier_signal=None,
    columns_categories=None,
    features_selection=None,
    observations_selection=None,
    filter_table_main=None,
    transpose_table_supplement=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 30 September 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table_main = table_main.copy(deep=True)
    if (table_supplement is not None):
        table_supplement = table_supplement.copy(deep=True)
        pass
    columns_categories = copy.deepcopy(columns_categories)
    features_selection = copy.deepcopy(features_selection)
    observations_selection = copy.deepcopy(observations_selection)

    # Filter columns and rows in main table for specific features and
    # observations.
    # Notice that the function below already includes the name of the column
    # for the index across rows.
    if (filter_table_main):
        # Copy information.
        categories_features_selection = copy.deepcopy(columns_categories)
        # Prepare inclusive list of columns.
        categories_features_selection.insert(0, column_identifier_observation)
        #categories_features_selection.insert(0, column_identifier_signal)
        categories_features_selection.extend(features_selection)
        # Filter columns and rows in table.
        table_main_selection = (
            porg.filter_select_table_columns_rows_by_identifiers(
                table=table_main,
                index_rows=column_identifier_signal,
                identifiers_columns=categories_features_selection,
                identifiers_rows=observations_selection,
                report=False,
        ))
    else:
        # Copy information.
        table_main_selection = table_main.copy(deep=True)
        pass

    # Determine whether there is a table of supplemental features (signals) to
    # merge into the main table of features and observations.
    if (table_supplement is not None):
        # Optional preliminary transposition.
        # Copy information.
        table_supplement_format = table_supplement.copy(deep=True)
        # Determine whether to apply optional transposition.
        if (transpose_table_supplement):
            # Organize indices in table.
            table_supplement_format = (
                porg.explicate_table_indices_columns_rows_single_level(
                    table=table_supplement_format,
                    index_columns=column_identifier_signal,
                    index_rows=column_identifier_feature,
                    explicate_indices=True,
                    report=report,
            ))
            # Transpose table.
            table_supplement_format = table_supplement_format.transpose(
                copy=True
            )
            # Organize indices in table.
            table_supplement_format.reset_index(
                level=None,
                inplace=True,
                drop=False, # remove index; do not move to regular columns
            )
            table_supplement_format.columns.rename(
                None,
                inplace=True,
            ) # single-dimensional index

        # Filter columns and rows in table for specific features and
        # observations.
        table_supplement_selection = (
            porg.filter_select_table_columns_rows_by_identifiers(
                table=table_supplement_format,
                index_rows=column_identifier_signal,
                identifiers_columns=features_selection,
                identifiers_rows=observations_selection,
                report=False,
        ))

        # Merge together features from table of supplemental signals to those
        # in table of main observations.
        table_features_observations = porg.merge_columns_two_tables(
            identifier_first=column_identifier_signal,
            identifier_second=column_identifier_signal,
            table_first=table_main_selection,
            table_second=table_supplement_selection,
            preserve_index=False,
            report=report,
        )
        pass

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
        )
        print(str("module: " + module))
        print("function: filter_combine_features_observations()")
        putly.print_terminal_partition(level=5)
        print("table of combined features and observations:")
        print(table_features_observations)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table_features_observations


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
    aspect=None,
    fonts=None,
    colors=None,
    set_axis_limits=None,
    emphasis_marker=None,
    emphasis_label=None,
    line_diagonal=None,
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

    Review: TCW; 1 October 2025
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
        aspect (str): aspect ratio for MatPlotLib chart figure
        fonts (dict<object>): definitions of font properties
        colors (dict<tuple>): definitions of color properties
        set_axis_limits (bool): whether to set explicity limits on axes
        emphasis_marker (bool): whether to create special markers to emphasize
            a special selection of points
        emphasis_label (bool): whether to create text labels adjacent to
            the special selection of points for special emphasis
        line_diagonal (bool): whether to draw diagonal line for equality
            between abscissa and ordinate
        show_legend_bar (bool): whether to show legend or scale bar on chart
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
        table_group = table.copy(deep=True)

        # Count unique categorical values.
        #categories = copy.deepcopy(
        #    table[column_category].unique().tolist()
        #)
        #count_category_response = len(categories)
        count_category_response = (
            table[column_response].nunique(dropna=True)
        )
        # Organize indices in table.
        table.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table.set_index(
            [column_identifier, column_response],
            append=False,
            drop=True,
            inplace=True
        )
        # Split rows within table by factor columns.
        groups_response = table.groupby(
            level=column_response,
        )
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
        )
        print(str("module: " + module))
        function = "plot_scatter_point_color_response_discrete_or_continuous()"
        print("function: " + function)
        putly.print_terminal_partition(level=5)
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
            labelpad=20,
            alpha=1.0,
            #loc="left",
            backgroundcolor=colors["white"],
            color=colors["black"],
            fontproperties=fonts["properties"][size_title_abscissa]
        )
    if (len(title_ordinate) > 0):
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=20,
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
            table[column_abscissa].values,
            table[column_ordinate].values,
            linestyle="",
            marker="o",
            markersize=size_marker,
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
            table[column_abscissa].values,
            table[column_ordinate].values,
            c=table_standard[column_response].values,
            s=(size_marker**2), # scale with area (dimension squared)
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
        # matplotlib.pyplot.get_cmap("tab10", count_categories)
        color_map = matplotlib.cm.get_cmap("tab10") # "Set1", "Set2", "Dark2", "tab10",
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
                markersize=size_marker,
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
        (not table_special.empty) and
        (emphasis_marker)
    ):
        handle_special = axes.plot(
            table_special[column_abscissa].values,
            table_special[column_ordinate].values,
            linestyle="",
            marker="o",
            markersize=(size_marker*1.5),
            markeredgecolor=colors["green_kelly"],
            markerfacecolor=colors["green_kelly"]
        )
    # For table of special selection of records, create text labels adjacent
    # to points.
    if (
        (identifiers_emphasis is not None) and
        (len(identifiers_emphasis) > 0) and
        (not table_special.empty) and
        (emphasis_label)
    ):
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
        (show_legend_bar) and
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
        (show_legend_bar) and
        (column_response is not None) and
        (column_response != "") and
        (type_response is not None) and
        (type_response == "category")
    ):
        # Create legend.
        handles_legend = [matplotlib.patches.Patch(
            color=colors_groups[i],
            label=labels_groups[i]#,
            #size?
        ) for i in range(count_category_response)]
        figure.legend(
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
    column_response=None,
    column_abscissa=None,
    column_ordinate=None,
    type_response=None,
    title_chart=None,
    title_response=None,
    title_abscissa=None,
    title_ordinate=None,
    identifiers_emphasis=None,
    report=None,
):
    """
    Create and plot a chart of the scatter point type.

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
        column_response (str): name of column in table corresponding to values
            for representation as color of individual points
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

        ...
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
    colors = pplot.define_color_properties()
    # Create figure.
    if False:
        figure = pplot.plot_scatter(
            data=table,
            abscissa=column_abscissa,
            ordinate=column_ordinate,
            title_abscissa=title_abscissa,
            title_ordinate=title_ordinate,
            fonts=fonts,
            colors=colors,
            size=10,
        )
    if True:
        figure = plot_scatter_point_color_response_discrete_or_continuous(
            table=table,
            column_identifier=column_identifier,
            column_name=column_name,
            column_response=column_response,
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
            size_label_emphasis="fifteen",
            size_label_legend_bar="thirteen",
            size_marker=15,
            aspect="landscape",
            fonts=fonts,
            colors=colors,
            set_axis_limits=False,
            emphasis_marker=True,
            emphasis_label=True,
            line_diagonal=False, # diagonal is not proportional to respective ranges of axes
            show_legend_bar=True,
            report=None,
        )

    # Write product information to file.

    # Define paths to directories.
    path_directory_chart = os.path.join(
        path_directory_parent, "scatter_categories",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_chart,
    )

    # Bundle information.
    pail_write_plot = dict()
    pail_write_plot[name_chart] = figure

    # Write figure object to file.
    pplot.write_product_plots_parent_directory(
        pail_write=pail_write_plot,
        format="jpg", # jpg, png, svg
        resolution=150,
        path_directory=path_directory_chart,
    )

    # Return information.
    return figure


def create_write_plot_charts_principal_component_scores(
    path_directory_parent=None,
    table_features_observations=None,
    table_scores=None,
    column_identifier_observation=None,
    column_identifier_signal=None,
    columns_categories=None,
    columns_components=None,
    name_set_features=None,
    features_response_quantitative=None,
    observations_selection=None,
    groups_observations=None,
    names_groups_observations_sequence=None,
    translations_observations=None,
    write_charts=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 1 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table_features_observations = table_features_observations.copy(deep=True)
    table_scores = table_scores.copy(deep=True)
    columns_categories = copy.deepcopy(columns_categories)
    columns_components = copy.deepcopy(columns_components)
    features_response_quantitative = copy.deepcopy(
        features_response_quantitative
    )
    observations_selection = copy.deepcopy(observations_selection)
    groups_observations = copy.deepcopy(groups_observations)
    names_groups_observations_sequence = copy.deepcopy(
        names_groups_observations_sequence
    )
    translations_observations = copy.deepcopy(translations_observations)

    # Determine translations for names of principal components.
    translations_components = dict()
    for column_component in columns_components:
        name_prefix = str(name_set_features + "_")
        translation = str(column_component).replace(name_prefix, "")
        translations_components[column_component] = translation
        pass

    # Prepare reference to sort groups of rows in table.
    sequence_groups_observations = dict()
    index = 0
    for name in names_groups_observations_sequence:
        sequence_groups_observations[name] = index
        index += 1
        pass

    # Filter columns and rows in table.
    columns_sequence = [
        column_identifier_observation,
        column_identifier_signal,
    ]
    columns_sequence.extend(columns_categories)
    columns_sequence.extend(columns_components)
    columns_sequence.extend(features_response_quantitative)

    table_selection = porg.filter_select_table_columns_rows_by_identifiers(
        table=table_features_observations,
        index_rows=column_identifier_signal,
        identifiers_columns=columns_sequence,
        identifiers_rows=observations_selection,
        report=False,
    )

    # Filter rows in table for non-missing values across relevant columns.
    table_selection.dropna(
        axis="index",
        how="any",
        subset=columns_components,
        inplace=True,
    )

    # Determine and fill groups of observations.
    # These next two functions create a new column named "group" and assign
    # categorical names corresponding to specific groups of observations.
    # Each observation can only belong to a single group.
    table_group = porg.determine_fill_table_groups_rows(
        table=table_selection,
        index_rows=column_identifier_signal,
        column_group="group",
        groups_rows=groups_observations,
        report=False,
    )
    # Sort rows in table by groups.
    table_group = porg.sort_table_rows_by_single_column_reference(
        table=table_group,
        index_rows=column_identifier_signal,
        column_reference="group",
        column_sort_temporary="sort_temporary",
        reference_sort=sequence_groups_observations,
    )

    # Define comparisons between principal components.
    comparisons = [
        ["1", "2",],
        ["1", "3",],
        ["1", "4",],
        ["1", "5",],
        ["2", "3",],
        ["2", "4",],
        ["2", "5",],
        ["3", "4",],
        ["3", "5",],
    ]
    # Iterate on comparisons.
    for comparison in comparisons:
        # Determine name for chart.
        name_chart = str(
            "chart_" +
            name_set_features +
            "_pc_" + comparison[0] + "_" + comparison[1]
        )
        # Determine names of columns for values to represent on abscissa and
        # ordinate axes.
        column_abscissa = str(name_set_features + "_pc_" + comparison[0])
        column_ordinate = str(name_set_features + "_pc_" + comparison[1])
        title_abscissa = str("PC-" + comparison[0])
        title_ordinate = str("PC-" + comparison[1])

        # Create plot chart and write to file.
        create_write_plot_chart_scatter_point_response(
            path_directory_parent=path_directory_parent,
            name_chart=name_chart,
            table=table_group,
            column_identifier=column_identifier_signal,
            column_name=column_identifier_observation,
            column_response="group",
            column_abscissa=column_abscissa,
            column_ordinate=column_ordinate,
            type_response="category",
            title_chart="",
            title_response="group",
            title_abscissa=title_abscissa,
            title_ordinate=title_ordinate,
            #identifiers_emphasis=list(),
            identifiers_emphasis=[
                "QCW-7445",
                "BNI-8563",
                "CWA-8409",
                "PZA-6684",
                "SMH-9005",
                "DCH-0868",
            ],
            report=report,
        )
        pass




    ##########
    # Plot chart: scatter with colors for categorical groups of observations


    ##########
    # Plot chart: scatter with color gradient for response feature on a
    # quantitative scale.





    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
        )
        print(str("module: " + module))
        print("function: create_write_plot_charts()")
        putly.print_terminal_partition(level=5)
        pass

    pass


def manage_create_write_plot_charts(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    table_features_observations=None,
    table_loadings=None,
    table_variances=None,
    table_scores=None,
    column_identifier_observation=None,
    column_identifier_signal=None,
    columns_categories=None,
    columns_components=None,
    name_set_features=None,
    features_selection=None,
    features_response_quantitative=None,
    observations_selection=None,
    groups_observations=None,
    names_groups_observations_sequence=None,
    translations_features=None,
    translations_observations=None,
    write_charts=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 1 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table_features_observations = table_features_observations.copy(deep=True)
    table_loadings = table_loadings.copy(deep=True)
    table_variances = table_variances.copy(deep=True)
    table_scores = table_scores.copy(deep=True)
    columns_categories = copy.deepcopy(columns_categories)
    columns_components = copy.deepcopy(columns_components)
    features_selection = copy.deepcopy(features_selection)
    observations_selection = copy.deepcopy(observations_selection)
    groups_observations = copy.deepcopy(groups_observations)
    names_groups_observations_sequence = copy.deepcopy(
        names_groups_observations_sequence
    )
    translations_features = copy.deepcopy(translations_features)
    translations_observations = copy.deepcopy(translations_observations)

    ##########
    # Write product information to file.
    # Define paths to directories.
    path_directory_instance = os.path.join(
        path_directory_product, name_set_features,
    )
    path_directory_instance_charts = os.path.join(
        path_directory_instance, "charts",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_instance_charts,
    )

    # Create and write plot charts for the main table of features and
    # observations, including principal components as features.
    create_write_plot_charts_principal_component_scores(
        path_directory_parent=path_directory_instance_charts,
        table_features_observations=table_features_observations,
        table_scores=table_scores,
        column_identifier_observation=column_identifier_observation,
        column_identifier_signal=column_identifier_signal,
        columns_categories=columns_categories,
        columns_components=columns_components,
        name_set_features=name_set_features,
        features_response_quantitative=features_response_quantitative,
        observations_selection=observations_selection,
        groups_observations=groups_observations,
        names_groups_observations_sequence=names_groups_observations_sequence,
        translations_observations=translations_observations,
        write_charts=write_charts,
        report=report,
    )

    # Create and write plot charts for the table of principal component
    # loadings.
    # TODO: TCW; 1 October 2025
    # heatmap with labeled axes


    ##########
    # Plot chart: scatter with colors for categorical groups of observations


    ##########
    # Plot chart: scatter with color gradient for response feature on a
    # quantitative scale.


    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
        )
        print(str("module: " + module))
        print("function: create_write_plot_charts()")
        putly.print_terminal_partition(level=5)
        pass

    pass


def manage_components_set_features_groups_observations(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    table=None,
    column_identifier_observation=None,
    column_identifier_signal=None,
    columns_categories=None,
    name_set_features=None,
    features_selection=None,
    observations_selection=None,
    groups_observations=None,
    names_groups_observations_sequence=None,
    translations_features=None,
    translations_observations=None,
    allow_replicate_observations=None,
    write_lists=None,
    write_tables=None,
    write_charts=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 1 October 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    table = table.copy(deep=True)
    columns_categories = copy.deepcopy(columns_categories)
    features_selection = copy.deepcopy(features_selection)
    observations_selection = copy.deepcopy(observations_selection)
    groups_observations = copy.deepcopy(groups_observations)
    names_groups_observations_sequence = copy.deepcopy(
        names_groups_observations_sequence
    )
    translations_features = copy.deepcopy(translations_features)
    translations_observations = copy.deepcopy(translations_observations)

    # Define prefix for names of principal components.
    prefix_name_components = str(
        name_set_features + "_" + "pc"
    )

    # Calculate principal components and combine with table of main features
    # and observations.
    pail_components = (
        pdcmp.calculate_principal_components_table_features_observations(
            table=table,
            name_index_columns="features",
            name_index_rows=column_identifier_signal,
            columns_selection=features_selection,
            rows_selection=observations_selection,
            prefix=prefix_name_components,
            separator="_",
            #threshold_proportion=None,
            threshold_proportion=0.70, # cumulative proportion of variance; float or None to keep all
            explicate_indices=False,
            report=report,
    ))
    #pail_components["table_merge"] = table_merge
    #pail_components["table_loadings"] = table_loadings
    #pail_components["table_variances"] = table_variances
    #pail_components["table_scores"] = table_scores
    #pail_components["columns_scores"] = columns_scores

    # Simplify main table of features and observations.
    columns_all = copy.deepcopy(
        pail_components["table_merge"].columns.to_list()
    )
    columns_sequence = list(filter(
        lambda item: item not in features_selection, columns_all
    ))
    # Filter and sort columns in table.
    table_features_observations = porg.filter_sort_table_columns(
        table=pail_components["table_merge"],
        columns_sequence=columns_sequence,
        report=report,
    )

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Lists.
    pail_write_lists = dict()
    pail_write_lists["columns_scores"] = (
        pail_components["columns_scores"]
    )
    # Tables.
    pail_write_tables = dict()
    pail_write_tables["table_loadings"] = (
        pail_components["table_loadings"]
    )
    pail_write_tables["table_variances"] = (
        pail_components["table_variances"]
    )
    pail_write_tables["table_scores"] = pail_components["table_scores"]

    ##########
    # Write product information to file.

    # Define paths to directories.
    path_directory_instance = os.path.join(
        path_directory_product, name_set_features,
    )
    path_directory_instance_lists = os.path.join(
        path_directory_instance, "lists",
    )
    path_directory_instance_tables = os.path.join(
        path_directory_instance, "tables",
    )
    path_directory_instance_charts = os.path.join(
        path_directory_instance, "charts",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_instance_lists,
    )
    putly.create_directories(
        path=path_directory_instance_tables,
    )
    putly.create_directories(
        path=path_directory_instance_charts,
    )
    # Lists.
    if (write_lists):
        putly.write_lists_to_file_text(
            pail_write=pail_write_lists,
            path_directory=path_directory_instance_lists,
            delimiter="\n",
        )
        pass
    # Tables.
    if (write_tables):
        putly.write_tables_to_file(
            pail_write=pail_write_tables,
            path_directory=path_directory_instance_tables,
            reset_index_rows=False,
            write_index_rows=False,
            write_index_columns=True,
            type="text",
            delimiter="\t",
            suffix=".tsv",
        )
        pass

    # Manage the creation and write of plot charts.
    # Calculate principal components.
    manage_create_write_plot_charts(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        table_features_observations=table_features_observations,
        table_loadings=pail_components["table_loadings"],
        table_variances=pail_components["table_variances"],
        table_scores=pail_components["table_scores"],
        column_identifier_observation=column_identifier_observation,
        column_identifier_signal=column_identifier_signal,
        columns_categories=columns_categories,
        columns_components=pail_components["columns_scores"],
        name_set_features=name_set_features,
        features_selection=features_selection,
        features_response_quantitative=list(), # place holder (TCW; 1 October 2025)
        observations_selection=observations_selection,
        groups_observations=groups_observations,
        names_groups_observations_sequence=names_groups_observations_sequence,
        translations_features=translations_features,
        translations_observations=translations_observations,
        write_charts=write_charts,
        report=report,
    )


    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
        )
        print(str("module: " + module))
        print("function: manage_components_set_features_groups_observations()")
        putly.print_terminal_partition(level=5)
        pass


    # Bundle information.
    pail_return = dict()
    pail_return["table"] = table_features_observations
    # Return information.
    return pail_return



################################################################################
# Procedure




##########
# Call main procedure.


def execute_procedure(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_features=None,
    path_file_source_table_observations=None,
    path_file_source_table_signals=None,
    path_file_source_table_sets_features=None,
    path_file_source_table_groups_observations=None,
    column_identifier_feature=None,
    column_name_feature=None,
    column_identifier_observation=None,
    column_identifier_signal=None,
    count_components=None,
    transpose_table_signals=None,
    allow_replicate_observations=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review: TCW; 30 September 2025

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_file_source_table_features (str): path to source file
        path_file_source_table_observations (str): path to source file
        path_file_source_table_signals (str): path to source file
        path_file_source_table_sets_features (str): path to source file
        path_file_source_table_groups_observations (str): path to source file
        column_identifier_feature (str): name of column in source table
        column_name_feature (str): name of column in source table
        column_identifier_observation (str): name of column in source table
        column_identifier_signal (str): name of column in source table
        count_components (int): count of principal components to keep for each
            set of features
        transpose_table_signals (bool): whether to transpose the table of
            signals to match orientation in table of observations (columns:
            features; rows: observations)
        allow_replicate_observations (bool): whether to allow replicate
            observations or to require groups to be mutually exclusive, such
            that any individual observation can only belong to one group
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parse parameters.
    pail_parameters = parse_text_parameters(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_features=path_file_source_table_features,
        path_file_source_table_observations=(
            path_file_source_table_observations
        ),
        path_file_source_table_signals=path_file_source_table_signals,
        path_file_source_table_sets_features=(
            path_file_source_table_sets_features
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        column_identifier_feature=column_identifier_feature,
        column_name_feature=column_name_feature,
        column_identifier_observation=column_identifier_observation,
        column_identifier_signal=column_identifier_signal,
        count_components=count_components,
        transpose_table_signals=transpose_table_signals,
        allow_replicate_observations=allow_replicate_observations,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "calculate_principal_components_sets_features_groups_" +
            "observations.py"
        )
        print(str("module: " + module))
        print("function: execute_procedure()")
        print("system: local")
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        path_file_source_table_features=(
            pail_parameters["path_file_source_table_features"]
        ),
        path_file_source_table_observations=(
            pail_parameters["path_file_source_table_observations"]
        ),
        path_file_source_table_signals=(
            pail_parameters["path_file_source_table_signals"]
        ),
        path_file_source_table_sets_features=(
            pail_parameters["path_file_source_table_sets_features"]
        ),
        path_file_source_table_groups_observations=(
            pail_parameters["path_file_source_table_groups_observations"]
        ),
        report=pail_parameters["report"],
    )

    # Parameters.

    features_available = sutly.determine_features_available(
        table_features=pail_source["table_features"],
        table_observations=pail_source["table_observations"],
        table_signals=pail_source["table_signals"],
        column_identifier_feature=pail_parameters["column_identifier_feature"],
        column_name_feature=pail_parameters["column_name_feature"],
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_identifier_signal=pail_parameters["column_identifier_signal"],
        transpose_table_signals=pail_parameters["transpose_table_signals"],
        report=pail_parameters["report"],
    )

    pail_sets = sutly.read_organize_parameters_sets_features(
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        table=pail_source["table_sets"],
        column_name="abbreviation",
        features_available=features_available,
        report=pail_parameters["report"],
    )
    #pail_sets["table"]
    #pail_sets["names_sets_features_sequence"]
    #pail_sets["sets_features"]
    #pail_sets["features_sets_union"]
    #pail_sets["records"]

    pail_features = sutly.organize_parameters_further_sets_features(
        table_features=pail_source["table_features"],
        column_identifier_feature=(
            pail_parameters["column_identifier_feature"]
        ),
        column_name_feature=(
            pail_parameters["column_name_feature"]
        ),
        features_selection=pail_sets["features_sets_union"],
        features_sets_union=pail_sets["features_sets_union"],
        sets_features=pail_sets["sets_features"],
        names_sets_features_sequence=pail_sets["names_sets_features_sequence"],
        prefix_name_feature="",
        report=pail_parameters["report"],
    )
    #pail_features["table_features_selection"]
    #pail_features["features_selection"]
    #pail_features["features_selection_translation"]
    #pail_features["features_selection_prefix"]
    #pail_features["translations_features"]
    #pail_features["translations_features_prefix"]
    #pail_features["names_sets_features_sequence"]
    #pail_features["sets_features"]

    pail_groups = sutly.organize_parameters_groups_observations(
        table=pail_source["table_groups"],
        column_name="abbreviation",
        report=pail_parameters["report"],
    )
    #pail_groups["table"]
    #pail_groups["names_groups_observations_sequence"]
    #pail_groups["categories_groups"]
    #pail_groups["records"]

    pail_observations = sutly.organize_parameters_further_groups_observations(
        table_observations=pail_source["table_observations"],
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_identifier_signal=(
            pail_parameters["column_identifier_signal"]
        ),
        instances_groups_observations=pail_groups["records"],
        key_name="abbreviation",
        names_groups_observations_sequence=(
            pail_groups["names_groups_observations_sequence"]
        ),
        report=report,
    )
    #pail_observations["table_observations_selection"]
    #pail_observations["observations_selection"]
    #pail_observations["translations_observations"]
    #pail_observations["names_groups_observations_sequence"]
    #pail_observations["groups_observations"]

    # Filter features and observations anc combine those from "table_signals"
    # to those in "table_observations".
    table_features_observations = filter_combine_features_observations(
        table_main=pail_source["table_observations"],
        table_supplement=pail_source["table_signals"],
        column_identifier_feature=(
            pail_parameters["column_identifier_feature"]
        ),
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_identifier_signal=pail_parameters["column_identifier_signal"],
        columns_categories=pail_groups["categories_groups"],
        features_selection=pail_features["features_selection"],
        observations_selection=pail_observations["observations_selection"],
        filter_table_main=False,
        transpose_table_supplement=(
            pail_parameters["transpose_table_signals"]
        ),
        report=pail_parameters["report"],
    )

    # Iterate on sets of features.
    for name_set_features in pail_features["names_sets_features_sequence"]:

        # Use "name_set_features":
        # 1. to access features in set
        # 2. to name the PCs
        # 3. to name the files (tables, charts)

        # select relevant features
        # need to select relevant observations from table


        # Access set of features for current instance.
        features_set = pail_features["sets_features"][name_set_features]
        # Calculate principal components.
        pail_instance = manage_components_set_features_groups_observations(
            path_directory_source=pail_parameters["path_directory_source"],
            path_directory_product=pail_parameters["path_directory_product"],
            path_directory_dock=pail_parameters["path_directory_dock"],
            table=table_features_observations,
            column_identifier_observation=(
                pail_parameters["column_identifier_observation"]
            ),
            column_identifier_signal=(
                pail_parameters["column_identifier_signal"]
            ),
            columns_categories=pail_groups["categories_groups"],
            name_set_features=name_set_features,
            features_selection=features_set,
            observations_selection=pail_observations["observations_selection"],
            groups_observations=pail_observations["groups_observations"],
            names_groups_observations_sequence=(
                pail_observations["names_groups_observations_sequence"]
            ),
            translations_features=(
                pail_features["translations_features_prefix"]
            ),
            translations_observations=(
                pail_observations["translations_observations"]
            ),
            allow_replicate_observations=(
                pail_parameters["allow_replicate_observations"]
            ),
            write_lists=True,
            write_tables=True,
            write_charts=True,
            report=pail_parameters["report"],
        )
        # Update reference to table of features and observations.
        # Collect principal components for all sets of features.
        # Copy information.
        table_features_observations = pail_instance["table"].copy(deep=True)
        pass

    # Review and write the table after introduction of all relevant PCs.

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_source = sys.argv[1]
    path_directory_product = sys.argv[2]
    path_directory_dock = sys.argv[3]
    path_file_source_table_features = sys.argv[4]
    path_file_source_table_observations = sys.argv[5]
    path_file_source_table_signals = sys.argv[6]
    path_file_source_table_sets_features = sys.argv[7]
    path_file_source_table_groups_observations = sys.argv[8]
    column_identifier_feature = sys.argv[9]
    column_name_feature = sys.argv[10]
    column_identifier_observation = sys.argv[11]
    column_identifier_signal = sys.argv[12]
    count_components = sys.argv[13]
    transpose_table_signals = sys.argv[14]
    allow_replicate_observations = sys.argv[15]
    report = sys.argv[16]

    # Call function for procedure.
    execute_procedure(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_features=path_file_source_table_features,
        path_file_source_table_observations=(
            path_file_source_table_observations
        ),
        path_file_source_table_signals=path_file_source_table_signals,
        path_file_source_table_sets_features=(
            path_file_source_table_sets_features
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        column_identifier_feature=column_identifier_feature,
        column_name_feature=column_name_feature,
        column_identifier_observation=column_identifier_observation,
        column_identifier_signal=column_identifier_signal,
        count_components=count_components,
        transpose_table_signals=transpose_table_signals,
        allow_replicate_observations=allow_replicate_observations,
        report=report,
    )

    pass



#
