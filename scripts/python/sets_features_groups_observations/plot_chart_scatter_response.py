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
# Date, review or revision: 31 October 2025
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
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_features_observations=None,
    column_identifier_observation=None,
    column_name_observation=None,
    column_response_markers=None,
    column_group_ellipses=None,
    column_abscissa=None,
    column_ordinate=None,
    identifiers_emphasis=None,
    name_chart=None,
    title_chart=None,
    title_response=None,
    title_abscissa=None,
    title_ordinate=None,
    type_response=None,
    size_marker=None,
    size_edge_marker=None,
    size_edge_ellipse=None,
    factor_confidence_ellipse=None,
    colors_fill_markers=None,
    colors_fill_ellipses=None,
    color_edge_markers=None,
    color_edge_ellipses=None,
    show_confidence_ellipse=None,
    show_emphasis_marker=None,
    show_emphasis_label=None,
    show_legend_bar=None,
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
    pail["path_file_source_table_features_observations"] = str(
        path_file_source_table_features_observations
    ).strip()

    # Names and categories.
    # It is problematic to pass any white space in parameters from a script in
    # Bash. Designate the hash symbol "#" as a substitute for white space.
    pail["column_identifier_observation"] = str(
        column_identifier_observation
    ).strip().replace("#", " ")
    pail["column_name_observation"] = str(
        column_name_observation
    ).strip().replace("#", " ")
    pail["column_response_markers"] = str(
        column_response_markers
    ).strip().replace("#", " ")
    pail["column_group_ellipses"] = str(
        column_group_ellipses
    ).strip().replace("#", " ")
    pail["column_abscissa"] = str(
        column_abscissa
    ).strip().replace("#", " ")
    pail["column_ordinate"] = str(
        column_ordinate
    ).strip().replace("#", " ")
    pail["name_chart"] = str(
        name_chart
    ).strip().replace("#", " ")
    pail["title_chart"] = str(
        title_chart
    ).strip().replace("#", " ")
    pail["title_response"] = str(
        title_response
    ).strip().replace("#", " ")
    pail["title_abscissa"] = str(
        title_abscissa
    ).strip().replace("#", " ")
    pail["title_ordinate"] = str(
        title_ordinate
    ).strip().replace("#", " ")
    pail["type_response"] = str(
        type_response
    ).strip().replace("#", " ")

    # Number.
    pail["size_marker"] = float(str(size_marker).strip())
    pail["size_edge_marker"] = float(str(size_edge_marker).strip())
    pail["size_edge_ellipse"] = float(str(size_edge_ellipse).strip())
    pail["factor_confidence_ellipse"] = float(str(
        factor_confidence_ellipse
    ).strip())

    # List, simple text.
    identifiers_emphasis_parse = putly.parse_text_list_values(
        text=identifiers_emphasis,
        delimiter=",",
    )
    pail["identifiers_emphasis"] = putly.collect_unique_items(
        items=identifiers_emphasis_parse,
    )

    # List, objects.

    # Iterate on lists of colors.
    colors_lists = {
        "colors_fill_markers": colors_fill_markers,
        "colors_fill_ellipses": colors_fill_ellipses,
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
    # Iterate on individual of colors.
    colors = {
        "color_edge_markers": color_edge_markers,
        "color_edge_ellipses": color_edge_ellipses,
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
    if (
        (show_confidence_ellipse is not None) and
        (str(show_confidence_ellipse) != "") and
        (str(show_confidence_ellipse) != "none") and
        (str(show_confidence_ellipse) == "true")
    ):
        pail["show_confidence_ellipse"] = True
    else:
        pail["show_confidence_ellipse"] = False
        pass

    if (
        (show_emphasis_marker is not None) and
        (str(show_emphasis_marker) != "") and
        (str(show_emphasis_marker) != "none") and
        (str(show_emphasis_marker) == "true")
    ):
        pail["show_emphasis_marker"] = True
    else:
        pail["show_emphasis_marker"] = False
        pass

    if (
        (show_emphasis_label is not None) and
        (str(show_emphasis_label) != "") and
        (str(show_emphasis_label) != "none") and
        (str(show_emphasis_label) == "true")
    ):
        pail["show_emphasis_label"] = True
    else:
        pail["show_emphasis_label"] = False
        pass

    if (
        (show_legend_bar is not None) and
        (str(show_legend_bar) != "") and
        (str(show_legend_bar) != "none") and
        (str(show_legend_bar) == "true")
    ):
        pail["show_legend_bar"] = True
    else:
        pail["show_legend_bar"] = False
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
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_features_observations=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 4 November 2025

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


################################################################################
# Procedure


##########
# Call main procedure.


def execute_procedure(
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_features_observations=None,
    column_identifier_observation=None,
    column_name_observation=None,
    column_response_markers=None,
    column_group_ellipses=None,
    column_abscissa=None,
    column_ordinate=None,
    identifiers_emphasis=None,
    name_chart=None,
    title_chart=None,
    title_response=None,
    title_abscissa=None,
    title_ordinate=None,
    type_response=None,
    size_marker=None,
    size_edge_marker=None,
    size_edge_ellipse=None,
    factor_confidence_ellipse=None,
    colors_fill_markers=None,
    colors_fill_ellipses=None,
    color_edge_markers=None,
    color_edge_ellipses=None,
    show_confidence_ellipse=None,
    show_emphasis_marker=None,
    show_emphasis_label=None,
    show_legend_bar=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review: TCW; 4 November 2025

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
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        column_response_markers=column_response_markers,
        column_group_ellipses=column_group_ellipses,
        column_abscissa=column_abscissa,
        column_ordinate=column_ordinate,
        identifiers_emphasis=identifiers_emphasis,
        name_chart=name_chart,
        title_chart=title_chart,
        title_response=title_response,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        type_response=type_response,
        size_marker=size_marker,
        size_edge_marker=size_edge_marker,
        size_edge_ellipse=size_edge_ellipse,
        factor_confidence_ellipse=factor_confidence_ellipse,
        colors_fill_markers=colors_fill_markers,
        colors_fill_ellipses=colors_fill_ellipses,
        color_edge_markers=color_edge_markers,
        color_edge_ellipses=color_edge_ellipses,
        show_confidence_ellipse=show_confidence_ellipse,
        show_emphasis_marker=show_emphasis_marker,
        show_emphasis_label=show_emphasis_label,
        show_legend_bar=show_legend_bar,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_scatter_response.py"
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
        path_file_source_table_features_observations=(
            pail_parameters["path_file_source_table_features_observations"]
        ),
        report=pail_parameters["report"],
    )

    # Define paths to directories.
    path_directory_chart = os.path.join(
        pail_parameters["path_directory_product"], "scatter_categories",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_chart,
    )

    # Create plot chart and write to file.
    splot.create_write_plot_chart_scatter_point_response(
        path_directory_parent=path_directory_chart,
        name_chart=pail_parameters["name_chart"],
        table=pail_source["table_features_observations"],
        column_identifier=pail_parameters["column_identifier_observation"],
        column_name=pail_parameters["column_name_observation"],
        column_response_markers=pail_parameters["column_response_markers"],
        column_group_ellipses=pail_parameters["column_group_ellipses"],
        column_abscissa=pail_parameters["column_abscissa"],
        column_ordinate=pail_parameters["column_ordinate"],
        type_response=pail_parameters["type_response"],
        title_chart=pail_parameters["title_chart"],
        title_response=pail_parameters["title_response"],
        title_abscissa=pail_parameters["title_abscissa"],
        title_ordinate=pail_parameters["title_ordinate"],
        identifiers_emphasis=pail_parameters["identifiers_emphasis"],
        size_marker=pail_parameters["size_marker"],
        size_edge_marker=pail_parameters["size_edge_marker"],
        size_edge_ellipse=pail_parameters["size_edge_ellipse"],
        factor_confidence_ellipse=pail_parameters["factor_confidence_ellipse"],
        colors_fill_markers=pail_parameters["colors_fill_markers"],
        colors_fill_ellipses=pail_parameters["colors_fill_ellipses"],
        color_edge_markers=pail_parameters["color_edge_markers"],
        color_edge_ellipses=pail_parameters["color_edge_ellipses"],
        show_confidence_ellipse=pail_parameters["show_confidence_ellipse"],
        show_emphasis_marker=pail_parameters["show_emphasis_marker"],
        show_emphasis_label=pail_parameters["show_emphasis_label"],
        show_legend_bar=pail_parameters["show_legend_bar"],
        report=pail_parameters["report"],
    )


    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_source = sys.argv[1]
    path_directory_product = sys.argv[2]
    path_directory_dock = sys.argv[3]
    path_file_source_table_features_observations = sys.argv[4]
    column_identifier_observation = sys.argv[5]
    column_name_observation = sys.argv[6]
    column_response_markers = sys.argv[7]
    column_group_ellipses = sys.argv[8]
    column_abscissa = sys.argv[9]
    column_ordinate = sys.argv[10]
    identifiers_emphasis = sys.argv[11]
    name_chart = sys.argv[12]
    title_chart = sys.argv[13]
    title_response = sys.argv[14]
    title_abscissa = sys.argv[15]
    title_ordinate = sys.argv[16]
    type_response = sys.argv[17]
    size_marker = sys.argv[18]
    size_edge_marker = sys.argv[19]
    size_edge_ellipse = sys.argv[20]
    factor_confidence_ellipse = sys.argv[21]
    colors_fill_markers = sys.argv[22]
    colors_fill_ellipses = sys.argv[23]
    color_edge_markers = sys.argv[24]
    color_edge_ellipses = sys.argv[25]
    show_confidence_ellipse = sys.argv[26]
    show_emphasis_marker = sys.argv[27]
    show_emphasis_label = sys.argv[28]
    show_legend_bar = sys.argv[29]
    report = sys.argv[30]

    # Call function for procedure.
    execute_procedure(
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        column_response_markers=column_response_markers,
        column_group_ellipses=column_group_ellipses,
        column_abscissa=column_abscissa,
        column_ordinate=column_ordinate,
        identifiers_emphasis=identifiers_emphasis,
        name_chart=name_chart,
        title_chart=title_chart,
        title_response=title_response,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        type_response=type_response,
        size_marker=size_marker,
        size_edge_marker=size_edge_marker,
        size_edge_ellipse=size_edge_ellipse,
        factor_confidence_ellipse=factor_confidence_ellipse,
        colors_fill_markers=colors_fill_markers,
        colors_fill_ellipses=colors_fill_ellipses,
        color_edge_markers=color_edge_markers,
        color_edge_ellipses=color_edge_ellipses,
        show_confidence_ellipse=show_confidence_ellipse,
        show_emphasis_marker=show_emphasis_marker,
        show_emphasis_label=show_emphasis_label,
        show_legend_bar=show_legend_bar,
        report=report,
    )

    pass



#
