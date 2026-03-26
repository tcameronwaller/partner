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
# Date, initialization: 25 March 2026
# Date, review or revision: 25 March 2026
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
    path_file_source_list_sets=None,
    path_file_source_list_groups=None,
    name_batch=None,
    prefix_negative=None,
    prefix_positive=None,
    suffix_file_source=None,
    column_set=None,
    column_score=None,
    column_p=None,
    name_chart=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    title_legend=None,
    title_bar=None,
    label_legend=None,
    aspect=None,
    scale_minimum=None,
    scale_center=None,
    scale_maximum=None,
    threshold_p_first=None,
    threshold_p_second=None,
    threshold_q_first=None,
    threshold_q_second=None,
    size_title_chart=None,
    size_title_abscissa=None,
    size_title_ordinate=None,
    size_title_legend=None,
    size_title_bar=None,
    size_label_abscissa=None,
    size_label_ordinate=None,
    size_label_legend=None,
    size_label_bar=None,
    size_label_significance_p=None,
    size_label_significance_q=None,
    color_scale_maximum=None,
    color_scale_center=None,
    color_scale_minimum=None,
    show_significance_p=None,
    show_significance_q=None,
    show_legend=None,
    show_scale_bar=None,
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
    pail["path_file_source_list_sets"] = str(
        path_file_source_list_sets
    ).strip()
    pail["path_file_source_list_groups"] = str(
        path_file_source_list_groups
    ).strip()

    # Names and categories.
    # It is problematic to pass any white space in parameters from a script in
    # Bash. Designate the hash symbol "#" as a substitute for white space.
    # It is also problematic to pass an empty string in parameters from a
    # script in Bash. Designate the word "none" as a substitute for missing or
    # empty.
    # Iterate on individual names that could be empty or missing.
    names_categories = {
        "name_batch": name_batch,
        "prefix_negative": prefix_negative,
        "prefix_positive": prefix_positive,
        "suffix_file_source": suffix_file_source,
        "column_set": column_set,
        "column_score": column_score,
        "column_p": column_p,
        "name_chart": name_chart,
        "title_chart": title_chart,
        "title_abscissa": title_abscissa,
        "title_ordinate": title_ordinate,
        "title_legend": title_legend,
        "title_bar": title_bar,
        "label_legend": label_legend,
        "aspect": aspect,
        "size_title_chart": size_title_chart,
        "size_title_abscissa": size_title_abscissa,
        "size_title_ordinate": size_title_ordinate,
        "size_title_legend": size_title_legend,
        "size_title_bar": size_title_bar,
        "size_label_abscissa": size_label_abscissa,
        "size_label_ordinate": size_label_ordinate,
        "size_label_legend": size_label_legend,
        "size_label_bar": size_label_bar,
        "size_label_significance_p": size_label_significance_p,
        "size_label_significance_q": size_label_significance_q,
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

    # Numbers.
    # Iterate on individual numbers.
    numbers = {
        "scale_minimum": scale_minimum,
        "scale_center": scale_center,
        "scale_maximum": scale_maximum,
        "threshold_p_first": threshold_p_first,
        "threshold_p_second": threshold_p_second,
        "threshold_q_first": threshold_q_first,
        "threshold_q_second": threshold_q_second,
    }
    for key_number in numbers.keys():
        # Determine whether parameter has a valid value.
        if (
            (numbers[key_number] != "none")
        ):
            # Parse number.
            pail[key_number] = float(str(numbers[key_number]).strip())
        else:
            pail[key_number] = None
            pass
        pass

    # Lists, simple text.
    # List, objects.

    # Iterate on individual colors.
    colors = {
        "color_scale_maximum": color_scale_maximum,
        "color_scale_center": color_scale_center,
        "color_scale_minimum": color_scale_minimum,
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
        pass

    # Boolean, true or false.
    # Iterate on individual of Boolean designations.
    designations = {
        "show_significance_p": show_significance_p,
        "show_significance_q": show_significance_q,
        "show_legend": show_legend,
        "show_scale_bar": show_scale_bar,
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
            "plot_chart_heatmap_enrichments.py"
        )
        print(str("module: " + module))
        function = ("parse_text_parameters()")
        print(str("function: " + function))
        putly.print_terminal_partition(level=5)
        print("parameters:")
        print(pail)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Read source information.


def read_source_directory_files_tables(
    path_directory_parent=None,
    name_file_child_prefixes=None,
    name_file_child_suffix=None,
    name_file_child_not=None,
    report=None,
):
    """
    Read and organize source information.

    Date, revision or review: 25 March 2026

    arguments:
        path_directory_parent (str): path to parent directory in which to find
            child files
        name_file_child_prefixes (list<str>): prefixes in name by which to
            recognize relevant child files within parent directory
        name_file_child_suffix (str): suffix in name by which to recognize
            relevant child files within parent directory
        name_file_child_not (str): character string in names of files to exclude
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Bundle information.
    pail = dict()

    # Iterate on common prefixes for files in directory.
    for prefix in name_file_child_prefixes:
        # Extract and filter complete paths to child files within parent
        # directory.
        paths_prefix = putly.extract_filter_child_file_names_paths(
            path_directory=path_directory_parent,
            name_file_prefix=prefix,
            name_file_suffix=name_file_child_suffix,
            name_file_not=name_file_child_not,
            report=report,
        )
        # Iterate on paths to files.
        for path_file in paths_prefix:
            # Organize information for read.
            # Extract information from path to directory and file.
            path_directory = os.path.dirname(path_file)
            name_suffix_file = os.path.basename(path_file)
            #name_file = name_suffix_file.replace(str(name_file_child_suffix), "")
            name_file, suffix_file = os.path.splitext(name_suffix_file)
            # Read information from file.
            # Collect information.
            pail[name_file] = pandas.read_csv(
                path_file,
                sep="\t",
                header=0,
                #dtype=types_columns,
                na_values=[
                    "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
                ],
                encoding="utf-8",
            )
            pass
        pass

    # Report.
    if report:
        # Organize information.
        count_files = len(pail.keys())
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_enrichments.py"
        )
        print(str("module: " + module))
        function = ("read_source_directory_files_tables()")
        print(str("function: " + function))
        putly.print_terminal_partition(level=5)
        print("count of files in parent directory:")
        print(str(count_files))
        print(pail.keys())
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_list_sets=None,
    path_file_source_list_groups=None,
    prefix_negative=None,
    prefix_positive=None,
    suffix_file_source=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Date, revision or review: 25 March 2026

    arguments:
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_directory_dock_pail (str): path to dock directory for procedure's
            source and product directories and files
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_file_source_list_sets (str): path to source file
        path_file_source_list_groups (str): path to source file
        ...
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Read information from file.
    # Read tables of information about functional enrichments.
    pail_tables = read_source_directory_files_tables(
        path_directory_parent=path_directory_source,
        name_file_child_prefixes=[
            prefix_negative,
            prefix_positive,
        ],
        name_file_child_suffix=suffix_file_source,
        name_file_child_not="blargharium_291367",
        report=report,
    )

    # Determine whether paths point to a directory or file that exist.
    #existence_directory = os.path.exists(path_directory)
    existence_file_list_sets = os.path.exists(
        path_file_source_list_sets
    )
    existence_file_list_groups = os.path.exists(
        path_file_source_list_groups
    )

    # Bundle information.
    pail = dict()
    pail["tables"] = pail_tables
    pail["lists"] = dict()

    # Read information from file.
    if (existence_file_list_sets):
        # Collect information.
        pail["lists"]["sets"] = putly.read_file_text_list(
            path_file=path_file_source_list_sets,
            delimiter="\n",
            unique=True,
        )
    else:
        pail["lists"]["sets"] = None
        pass
    if (existence_file_list_groups):
        # Collect information.
        pail["lists"]["groups"] = putly.read_file_text_list(
            path_file=path_file_source_list_groups,
            delimiter="\n",
            unique=True,
        )
    else:
        pail["lists"]["groups"] = None
        pass

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_enrichments.py"
        )
        print(str("module: " + module))
        function = ("read_source()")
        print(str("function: " + function))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def organize_table_enrichments_wide(
    pail_tables=None,
    prefix_negative=None,
    prefix_positive=None,
    column_set=None,
    column_score=None,
    column_p=None,
    selection_sets=None,
    selection_groups=None,
    report=None,
):
    """
    Organize table of enrichments from all sets and groups in current batch.

    Date, revision or review: 25 March 2026

    arguments:
        pail_tables (dict<object>): collection of Pandas data-frame tables
        prefix_negative (str): prefix for names of tables
        prefix_positive (str): prefix for names of tables
        column_set (str): name of column in source table
        column_score (str): name of column in source table
        column_p (str): name of column in source table
        selection_sets (list<str>): names of sets in selection
        selection_groups (list<str>): names of groups in selection
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    # Copy information.
    pail_tables = copy.deepcopy(pail_tables)
    selection_sets = copy.deepcopy(selection_sets)
    selection_groups = copy.deepcopy(selection_groups)

    ##########
    # Collect information from all groups of enrichments.

    # Extract names of all tables in collection.
    names_tables = copy.deepcopy(list(pail_tables.keys()))
    # Extract unique base names of tables.
    names_tables_base = list(map(
        lambda name: str(name).strip().replace(prefix_negative, ""),
        names_tables
    ))
    names_tables_base = list(map(
        lambda name: str(name).strip().replace(prefix_positive, ""),
        names_tables_base
    ))
    names_tables_base = putly.collect_unique_items(
        items=names_tables_base,
    )

    # Collect information.
    table_collection = None
    # Iterate on names of groups corresponding to the collection of tables.
    for name_group in names_tables_base:
        # Derive names of tables.
        name_table_negative = str(prefix_negative + name_group)
        name_table_positive = str(prefix_positive + name_group)
        # Copy information.
        table_negative = pail_tables[name_table_negative].copy(deep=True)
        table_positive = pail_tables[name_table_positive].copy(deep=True)
        # Concatenate tables.
        table_group = pandas.concat(
            [table_negative, table_positive],
            ignore_index=False,
        )
        # Organize indices in table.
        table_group.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        # Organize new column to table for name of group.
        table_group["group"] = str(name_group)
        # Collect information.
        if (table_collection is not None):
            # Concatenate tables.
            table_collection = pandas.concat(
                [table_collection, table_group],
                ignore_index=True,
            )
        else:
            # Copy information.
            table_collection = table_group.copy(deep=True)
            pass
        pass

    ##########
    # Organize table of enrichments in partial wide format.

    # Copy information.
    table_wide_partial = table_collection.copy(deep=True)
    # Translate names of columns in table.
    translations = {
        column_set: "name_set",
        column_score: "score",
        column_p: "p_value",
    }
    table_wide_partial.rename(
        columns=translations,
        inplace=True,
    )
    # Organize information.
    table_wide_partial["score"] = pandas.to_numeric(
        table_wide_partial["score"],
        downcast="float",
        errors="coerce",
    )
    table_wide_partial["p_value"] = pandas.to_numeric(
        table_wide_partial["p_value"],
        downcast="float",
        errors="coerce",
    )
    # Filter rows in table by names of sets.
    if (
        (selection_sets is not None) and
        (len(selection_sets) > 0)
    ):
        table_wide_partial = table_wide_partial.loc[(
            (table_wide_partial["name_set"].isin(selection_sets))
        ), :].copy(deep=True)
    if (
        (selection_groups is not None) and
        (len(selection_groups) > 0)
    ):
        table_wide_partial = table_wide_partial.loc[(
            (table_wide_partial["group"].isin(selection_groups))
        ), :].copy(deep=True)
        pass
    # Organize indices in table.
    table_wide_partial.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Calculate Benjamini-Hochberg q-values for False-Discovery Rate (FDR).
    # Calculate q-values across all comparisons in table.
    # FDR 5% (q <= 0.05).
    table_wide_partial = (
        pdesc.calculate_table_false_discovery_rate_q_values(
            threshold=0.05, # alpha; family-wise error rate
            name_column_p_value="p_value",
            name_column_q_value="q_value",
            name_column_significance="q_significance",
            table=table_wide_partial,
    ))
    # Filter rows in table by q-value.
    if (False):
        table_wide_partial = table_wide_partial.loc[(
            (table_wide_partial["q_value"] < 0.05)
        ), :].copy(deep=True)
        pass
    # Filter and sort columns within table.
    columns_sequence = [
        "group",
        "name_set",
        "score",
        "p_value",
        "q_value",
    ]
    table_wide_partial = porg.filter_sort_table_columns(
        table=table_wide_partial,
        columns_sequence=columns_sequence,
        report=False,
    )

    ##########
    # Organize table of enrichments in full wide format.

    # Collect information.
    table_collection = None
    # Iterate on names of groups corresponding to the collection of tables.
    if (
        (selection_groups is not None) and
        (len(selection_groups) > 0)
    ):
        names_groups = copy.deepcopy(selection_groups)
    else:
        names_groups = copy.deepcopy(names_tables_base)
        pass
    for name_group in names_groups:
        # Filter rows in main table to select those for current group.
        table_group = table_wide_partial.loc[(
            (table_wide_partial["group"] == name_group)
        ), :].copy(deep=True)
        # Filter and sort columns within table.
        columns_sequence = [
            "name_set",
            "score",
            #"p_value",
            "q_value",
        ]
        table_group = porg.filter_sort_table_columns(
            table=table_group,
            columns_sequence=columns_sequence,
            report=False,
        )
        # Translate names of columns in table.
        translations = {
            "score": str(name_group + "_score"),
            #"p_value": str(name_group + "_p_value"),
            "q_value": str(name_group + "_q_value"),
        }
        table_group.rename(
            columns=translations,
            inplace=True,
        )
        # Collect information.
        if (table_collection is not None):
            # Merge tables.
            table_collection = porg.merge_columns_two_tables(
                identifier_first="name_set",
                identifier_second="name_set",
                table_first=table_collection,
                table_second=table_group,
                preserve_index=False,
                report=False,
            )
        else:
            # Copy information.
            table_collection = table_group.copy(deep=True)
            pass
        pass

    # Copy information.
    table_wide_full = table_collection.copy(deep=True)
    # Sort rows within table.
    table_wide_full.sort_values(
        by=[
            "name_set",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Organize indices in table.
    table_wide_full.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    # Bundle information.
    pail = dict()
    pail["table_wide_partial"] = table_wide_partial
    pail["table_wide_full"] = table_wide_full

    # Report.
    if report:
        # Organize information.
        # Print information.

        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_enrichments.py"
        )
        print(str("module: " + module))
        function = ("organize_table_enrichments_wide()")
        print(str("function: " + function))
        putly.print_terminal_partition(level=5)
        print("table of enrichments in partial wide format:")
        print(table_wide_partial)
        putly.print_terminal_partition(level=5)
        print("table of enrichments in full wide format:")
        print(table_wide_full)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def organize_table_enrichments_long(
    table_wide_partial=None,
    selection_sets=None,
    selection_groups=None,
    report=None,
):
    """
    Organize table of enrichments with a long format.

    Date, revision or review: 25 March 2026

    arguments:
        table_wide_partial (object): Pandas data-frame table in partial wide
            format
        selection_sets (list<str>): names of sets in selection
        selection_groups (list<str>): names of groups in selection
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information

    """

    # Copy information.
    table_wide_partial = table_wide_partial.copy(deep=True)
    selection_sets = copy.deepcopy(selection_sets)
    selection_groups = copy.deepcopy(selection_groups)

    ##########
    # Transform table from partial wide to full long format.
    # Organize indices in table.
    table_wide_partial.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_wide_partial.set_index(
        [
            "group",
            "name_set",
        ],
        append=False,
        drop=True,
        inplace=True,
    )
    table_wide_partial.columns.rename(
        "type_value",
        inplace=True,
    ) # single-dimensional index
    # Pandas dataframe methods "stack", "melt", and "wide_to_long", can all be
    # useful in this context.
    # Method "stack" converts to a multi-index series when the column index only
    # has a single level.
    # Method "wide_to_long" assumes that the information about multiple levels
    # in the column index is stored in delimited strings of compound column
    # names.
    # Transform the format of the table.
    if False:
        table_long = table_wide_partial.stack(
            level=name_row_index,
            #future_stack=True,
        )
    table_long_full = table_wide_partial.melt(
        id_vars=None,
        value_vars=None,
        var_name="type_value",
        value_name="value",
        ignore_index=False,
    )
    # Organize indices in table.
    table_long_full.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_long_full.set_index(
        [
            "group",
            "name_set",
            "type_value",
        ],
        append=False,
        drop=True,
        inplace=True,
    )

    ##########
    # Transform table from full long to partial long format.
    # Copy information.
    table_long_copy = table_long_full.copy(deep=True)
    # Organize indices in table.
    table_long_copy.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    # Pandas dataframe methods "unstack" and "pivot" can be useful in this
    # context.
    if False:
        table_long_copy.set_index(
            [name_column_index_1, name_row_index, name_column_index_2,],
            append=False,
            drop=True,
            inplace=True,
        )
        table_long_partial_q = table_long_copy.unstack(
            level=[name_row_index,],
        )
        pass
    table_long_partial = table_long_copy.pivot(
        columns=["group",],
        index=["name_set", "type_value",],
        values="value",
    )
    # Sort columns in table.
    # Alternatively, sort again by the explicit list that is a parameter for
    # this module.
    if True:
        table_long_partial = table_long_partial[[*selection_groups]]
    else:
        sequence_row_index = copy.deepcopy(
            table_wide_partial.index.get_level_values("group").to_list()
        )
        table_long_partial = table_long_partial[[*sequence_row_index]]
        pass

    # Bundle information.
    pail = dict()
    pail["table_long_full"] = table_long_full
    pail["table_long_partial"] = table_long_partial

    # Report.
    if report:
        # Organize information.
        # Print information.

        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "plot_chart_heatmap_enrichments.py"
        )
        print(str("module: " + module))
        function = ("organize_table_enrichments_long()")
        print(str("function: " + function))
        putly.print_terminal_partition(level=5)
        print("table of enrichments in full long format:")
        print(table_long_full)
        putly.print_terminal_partition(level=5)
        print("table of enrichments in partial long format:")
        print(table_long_partial)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# This design has advantages for representations of correlations or regression
# coefficients along with their significance.
def plot_heat_map_few_signal_significance_labels(
    table=None,
    transpose_table=None,
    name_index_columns_abscissa=None,
    name_index_rows_ordinate=None,
    name_index_rows_type=None,
    type_value_signal=None,
    type_value_p=None,
    type_value_q=None,
    fill_missing=None,
    value_fill_missing=None,
    constrain_signal_values=None,
    scale_minimum=None,
    scale_center=None,
    scale_maximum=None,
    show_significance_p=None,
    show_significance_q=None,
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
    size_label_significance_p=None,
    size_label_significance_q=None,
    size_title_ordinate=None,
    size_title_abscissa=None,
    size_label_ordinate=None,
    size_label_abscissa=None,
    size_label_legend=None,
    size_title_bar=None,
    size_label_bar=None,
    aspect=None,
    fonts=None,
    color_scale_maximum=None,
    color_scale_center=None,
    color_scale_minimum=None,
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
        name_index_columns_abscissa (str): name of index across table's columns that
            provides labels for groups across horizontal axis (abscissa) after
            any transpose on the table
        name_index_rows_ordinate (str): name of index across table's rows that provides
            labels for groups across vertical axis (ordinate) after any
            transpose on the table
        ...
        fill_missing (bool): whether to fill any missing values in every
            element of matrix
        value_fill_missing (float): value with which to fill any missing values
        constrain_signal_values (bool): whether to constrain all values in
            matrix
        scale_minimum (float): minimal value for constraint on signals and
            scale
        scale_maximum (float): maximal value for constraint on signals and
            scale
        show_significance_p (bool): whether to consider p-values to indicate
            significance in labels on the heatmap
        show_significance_q (bool): whether to consider q-values to indicate
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
        size_label_significance_p (str): font size
        size_label_significance_q (str): font size
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
    # Organize information for chart.

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
            table.index.get_level_values(
                name_index_rows_type
            ).isin([type_value_signal])
        ].copy(deep=True)
        table_p = table[
            table.index.get_level_values(
                name_index_rows_type
            ).isin([type_value_p])
        ].copy(deep=True)
        table_q = table[
            table.index.get_level_values(
                name_index_rows_type
            ).isin([type_value_q])
        ].copy(deep=True)
    else:
        # For this extraction, the "type_value" index must be oriented across
        # the table's columns.
        #table_signal = table.loc[
        #    (table["type_value"] == "signal"), :
        #].copy(deep=True)
        table_signal = table[
            table.columns.get_level_values(
                name_index_rows_type
            ).isin([type_value_signal])
        ].copy(deep=True)
        table_p = table[
            table.columns.get_level_values(
                name_index_rows_type
            ).isin([type_value_p])
        ].copy(deep=True)
        table_q = table[
            table.columns.get_level_values(
                name_index_rows_type
            ).isin([type_value_q])
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
    # Extract minimal and maximal values of signal intensity.
    if (
        (scale_minimum is None) or
        (scale_center is None) or
        (scale_maximum is None)
    ):
        round_offset = abs(numpy.nanmin(matrix_signal) * 0.01)
        value_minimum = round((numpy.nanmin(matrix_signal) - round_offset), 3)
        value_maximum = round((numpy.nanmax(matrix_signal) + round_offset), 3)
        scale_minimum = value_minimum
        scale_center = 0.0
        scale_maximum = value_maximum
        pass
    # Organize signals in matrix.
    # Replace missing values.
    if fill_missing:
        matrix_signal = numpy.nan_to_num(
            matrix_signal,
            copy=True,
            nan=value_fill_missing,
            posinf=scale_maximum, # or + 1.0 for correlations
            neginf=scale_minimum, # or - 1.0 for correlations
        )
    # Constrain values.
    if constrain_signal_values:
        matrix_signal[matrix_signal < scale_minimum] = scale_minimum
        matrix_signal[matrix_signal > scale_maximum] = scale_maximum

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
                name_index_columns_abscissa
            ).unique().to_list()
        )
        labels_rows = copy.deepcopy(
            table_signal.index.get_level_values(
                name_index_rows_ordinate
            ).unique().to_list()
        )
        labels_ordinate_categories = labels_rows # vertical axis
        labels_abscissa_categories = labels_columns # horizontal axis
        pass
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Column labels:")
        print(labels_ordinate_categories)
        print("Row labels:")
        print(labels_abscissa_categories)
        putly.print_terminal_partition(level=4)

    ##########
    # Organize parameters for visual representation.
    # Colors.
    pail_color_map = splot.create_color_map_divergent(
        value_minimum=scale_minimum,
        value_center=scale_center,
        value_maximum=scale_maximum,
        color_minimum=color_scale_minimum,
        color_center=color_scale_center,
        color_maximum=color_scale_maximum,
    )

    ##########
    # Create figure.
    figure = pplot.initialize_matplotlib_figure_aspect(
        aspect=aspect,
    )
    # Create axes.
    #axes = matplotlib.pyplot.axes()
    axes = figure.add_subplot(111)
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
            backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
            color=matplotlib.colors.to_rgba("black", 1.0),
        )
        text = axes.text(
            ((matrix_signal.shape[1]) + 0.75),
            -2.25,
            label_figure_top_right_2,
            horizontalalignment="left",
            verticalalignment="top",
            fontproperties=fonts["properties"][size_label_legend],
            backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
            color=matplotlib.colors.to_rgba("black", 1.0),
        )
        text = axes.text(
            ((matrix_signal.shape[1]) + 0.75),
            -1.5,
            label_figure_top_right_3,
            horizontalalignment="left",
            verticalalignment="top",
            fontproperties=fonts["properties"][size_label_legend],
            backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
            color=matplotlib.colors.to_rgba("black", 1.0),
        )

    # Plot values as a color grid.
    # This function represents values acros matrix dimension 0 as vertical
    # rows.
    # This function represents values across matrix dimension 1 as horizontal
    # columns.
    # Diverging color maps: "PRGn", "PRGn_r", "PiYG", "PiYG_r",
    # Diverging color maps: "PuOr", "PuOr_r",
    # Diverging color maps: "PuOr", "PuOr_r", "RdBu", "RdBu_r", "BrBG",
    # Sequential color maps: "Reds", "Reds_r", "Oranges", "Oranges_r",
    # site: https://montoliu.naukas.com/2021/11/18/color-blindness-purple-and-
    #   orange-are-the-solution/
    image = axes.imshow(
        matrix_signal,
        #cmap=matplotlib.colormaps["PuOr"], # RdBu_r, PuOr_r
        #vmin=scale_minimum,
        #vmax=scale_maximum,
        cmap=pail_color_map["map"],
        norm=pail_color_map["scale"],
        aspect="auto", # "auto", "equal",
        origin="upper",
        # Extent: (left, right, bottom, top)
        #extent=(-0.5, (matrix.shape[1] - 0.5), (matrix.shape[0] - 0.5), -0.5),
    )

    # Set titles for axes.
    if (len(title_ordinate) > 0):
        axes.set_ylabel(
            ylabel=title_ordinate,
            labelpad=15,
            alpha=1.0,
            backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
            color=matplotlib.colors.to_rgba("black", 1.0),
            fontproperties=fonts["properties"][size_title_ordinate]
        )
    if (len(title_abscissa) > 0):
        axes.set_xlabel(
            xlabel=title_abscissa,
            labelpad=15,
            alpha=1.0,
            backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
            color=matplotlib.colors.to_rgba("black", 1.0),
            fontproperties=fonts["properties"][size_title_abscissa]
        )
    # Set tick parameters for axes.
    axes.tick_params(
        axis="both", # "y", "x", or "both"
        which="major", # "major", "minor", or "both"
        length=5.0, # 5.0
        width=3.5, # 3.5
        pad=7.5, # 7.5
        direction="out",
        color=matplotlib.colors.to_rgba("black", 1.0),
        labelcolor=matplotlib.colors.to_rgba("black", 1.0),
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
    axes.set_xticks(
        numpy.arange(matrix_signal.shape[1]),
    )
    axes.set_yticks(
        numpy.arange(matrix_signal.shape[0]),
    )
    axes.set_xticklabels(
        labels_abscissa_categories,
        #minor=False,
        rotation=-60,
        rotation_mode="anchor",
        ha="left", # horizontal alignment
        va="top", # vertical alignment
        alpha=1.0,
        backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
        color=matplotlib.colors.to_rgba("black", 1.0),
        fontproperties=fonts["properties"][size_label_abscissa]
    )
    axes.set_yticklabels(
        labels_ordinate_categories,
        #minor=False,
        ha="right", # horizontal alignment
        va="center", # vertical alignment
        alpha=1.0,
        backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
        color=matplotlib.colors.to_rgba("black", 1.0),
        fontproperties=fonts["properties"][size_label_ordinate]
    )
    axes.tick_params(
        which="major",
        top=False,
        labeltop=False,
        bottom=True,
        labelbottom=True,
    )
    # Set grid for minor ticks without changing major ticks.
    # Shift x and y ticks by -0.5 to place grid between cells.
    axes.set_xticks(numpy.arange(-0.5, matrix_signal.shape[1], 1), minor=True)
    axes.set_yticks(numpy.arange(-0.5, matrix_signal.shape[0], 1), minor=True)
    axes.grid(
        True,
        which="minor",
        color="white",
        linestyle="-",
        linewidth=2
    )
    axes.tick_params(
        which="minor",
        top=False,
        labeltop=False,
        bottom=False,
        labelbottom=False,
        left=False,
        labelleft=False,
        right=False,
        labelright=False,
    )
    # Turn off grid for major ticks.
    axes.grid(which="major", linestyle="", color="none")
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
            #if ((signal_value > 0.5) or (signal_value < -0.5)):

            # Create labels on cells to represent p-values and q-values.
            #str("$\u2020$") # dagger U+2020
            #str("$\u2021$") # double dagger U+2021
            #str("$\u002A$") # asterisk U+002A
            #str("$\u2051$") # double asterisk U+2051

            # Determine whether to create labels for p-values on cells.
            if (
                (show_significance_p) and
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
                    fontproperties=(
                        fonts["properties"][size_label_significance_p]
                    ),
                    fontweight="extra bold",
                    backgroundcolor=matplotlib.colors.to_rgba("white", 0.0),
                    color=matplotlib.colors.to_rgba("black", 1.0),
                )
            elif (
                (show_significance_p) and
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
                    fontproperties=(
                        fonts["properties"][size_label_significance_p]
                    ),
                    fontweight="extra bold",
                    backgroundcolor=matplotlib.colors.to_rgba("white", 0.0),
                    color=matplotlib.colors.to_rgba("black", 1.0),
                )
                pass

            # Determine whether to create labels for q-values on cells.
            # Use the 'annotate' method for more versatility in position.
            if (
                (show_significance_q) and
                (not numpy.isnan(q_value)) and
                (q_value < thresholds_q[1])
            ):
                #label_cell = str(str(signal_text) + "**")
                label_cell = str("$\u2051$") # double asterisk U+2051
                text = axes.annotate(
                    label_cell,
                    xy=(index_column,index_row),
                    xycoords="data", # use coordinate system of object
                    #xytext=(10,15), # coordinates for offset of text from main
                    xytext=(0,0), # coordinates for offset of text from main
                    textcoords="offset points", # coordinates for 'xytext'
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontproperties=(
                        fonts["properties"][size_label_significance_q]
                    ),
                    backgroundcolor=matplotlib.colors.to_rgba("white", 0.0),
                    color=matplotlib.colors.to_rgba("black", 1.0),
                )
            elif (
                (show_significance_q) and
                (not numpy.isnan(q_value)) and
                (q_value < thresholds_q[0])
            ):
                label_cell = str("$\u002A$") # asterisk U+002A
                text = axes.annotate(
                    label_cell,
                    xy=(index_column,index_row),
                    xycoords="data", # use coordinate system of object
                    #xytext=(7,10), # coordinates for offset of text from main
                    xytext=(0,0), # coordinates for offset of text from main
                    textcoords="offset points", # coordinates for 'xytext'
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontproperties=(
                        fonts["properties"][size_label_significance_q]
                    ),
                    backgroundcolor=matplotlib.colors.to_rgba("white", 0.0),
                    color=matplotlib.colors.to_rgba("black", 1.0),
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
            shrink=0.5, # 0.7; factor for dimensions of the Scale Bar.
        )
        if (len(title_bar) > 0):
            bar.ax.set_ylabel(
                title_bar,
                rotation=-90,
                va="bottom",
                labelpad=5, # 5
                alpha=1.0,
                backgroundcolor=matplotlib.colors.to_rgba("white", 1.0),
                color=matplotlib.colors.to_rgba("black", 1.0),
                fontproperties=fonts["properties"][size_title_bar],
            )
        bar.ax.tick_params(
            axis="both",
            which="both", # major, minor, or both
            direction="out",
            length=5.0, # 5.0, 7.5
            width=3, # 2.5, 5.0
            color=matplotlib.colors.to_rgba("black", 1.0),
            pad=5, # 5, 7
            labelsize=fonts["values"][size_label_bar]["size"],
            labelcolor=matplotlib.colors.to_rgba("black", 1.0),
        )
        pass
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
                backgroundcolor=matplotlib.colors.to_rgba("white", 0.1),
                color=matplotlib.colors.to_rgba("black", 1.0),
                fontproperties=fonts["properties"]["four"]
            )
    # Return figure.
    return figure


def create_write_plot_chart_heatmap_enrichment(
    path_directory_parent=None,
    name_chart=None,
    table=None,
    name_index_columns_abscissa=None,
    name_index_rows_ordinate=None,
    name_index_rows_type=None,
    type_value_signal=None,
    type_value_p=None,
    type_value_q=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    title_legend=None,
    title_bar=None,
    label_legend=None,
    aspect=None,
    value_fill_missing=None,
    scale_minimum=None,
    scale_center=None,
    scale_maximum=None,
    thresholds_p=None,
    thresholds_q=None,
    size_title_chart=None,
    size_title_abscissa=None,
    size_title_ordinate=None,
    size_title_legend=None,
    size_title_bar=None,
    size_label_ordinate=None,
    size_label_abscissa=None,
    size_label_legend=None,
    size_label_bar=None,
    size_label_significance_p=None,
    size_label_significance_q=None,
    color_scale_maximum=None,
    color_scale_center=None,
    color_scale_minimum=None,
    constrain_signal_values=None,
    fill_missing=None,
    show_significance_p=None,
    show_significance_q=None,
    show_legend=None,
    show_scale_bar=None,
    report=None,
):
    """
    Create, plot, and write to file a chart of the type dot forest.

    Date, revision or review: 25 March 2026

    arguments:
        path_directory_parent (str): path to parent directory for procedure's
            product directories and files
        name_chart (str): name for writing figure object to file
        table (object): Pandas data-frame table of features across columns and
            observations across rows with values on quantitative, continuous
            interval or ratio scales of measurement
        name_index_columns_abscissa (str): name of index across columns in
            table for which to represent categories as labels on the horizontal
            abscissa axis
        name_index_rows_ordinate (str): name of column in table corresponding
            to an index across rows for which to represent categories as labels
            on the vertical ordinate axis
        ...
        title_chart (str): title of the chart
        title_abscissa (str): title of the feature to represent on the abscissa
            horizontal axis
        title_ordinate (str): title of the feature to represent on the ordinate
            vertical axis
        ...
        color_scale_maximum (tuple): definition of color properties
        color_scale_center (tuple): definition of color properties
        color_scale_minimum (tuple): definition of color properties
        ...
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

    # scale_center
    # title_legend
    # size_title_legend

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = plot_heat_map_few_signal_significance_labels(
        table=table,
        transpose_table=False,
        name_index_columns_abscissa=name_index_columns_abscissa,
        name_index_rows_ordinate=name_index_rows_ordinate,
        name_index_rows_type=name_index_rows_type,
        type_value_signal=type_value_signal,
        type_value_p=type_value_p,
        type_value_q=type_value_q,
        label_legend=label_legend,
        labels_ordinate_categories=None,
        labels_abscissa_categories=None,
        title_chart_top_right=title_chart,
        title_ordinate=title_ordinate,
        title_abscissa=title_abscissa,
        title_bar=title_bar,
        aspect=aspect,
        value_fill_missing=value_fill_missing,
        scale_minimum=scale_minimum,
        scale_center=scale_center,
        scale_maximum=scale_maximum,
        thresholds_p=thresholds_p,
        thresholds_q=thresholds_q,
        size_title_abscissa=size_title_abscissa,
        size_title_ordinate=size_title_ordinate,
        #size_title_legend=size_title_legend,
        size_title_bar=size_title_bar,
        size_label_abscissa=size_label_abscissa,
        size_label_ordinate=size_label_ordinate,
        size_label_legend=size_label_legend,
        size_label_bar=size_label_bar,
        size_label_significance_p=size_label_significance_p,
        size_label_significance_q=size_label_significance_q,
        color_scale_maximum=color_scale_maximum,
        color_scale_center=color_scale_center,
        color_scale_minimum=color_scale_minimum,
        fonts=fonts,
        constrain_signal_values=constrain_signal_values,
        fill_missing=fill_missing,
        show_significance_p=show_significance_p,
        show_significance_q=show_significance_q,
        show_legend=show_legend,
        show_scale_bar=show_scale_bar,
        report=report,
    )
    if False:
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
    path_file_source_list_sets=None,
    path_file_source_list_groups=None,
    name_batch=None,
    prefix_negative=None,
    prefix_positive=None,
    suffix_file_source=None,
    column_set=None,
    column_score=None,
    column_p=None,
    name_chart=None,
    title_chart=None,
    title_abscissa=None,
    title_ordinate=None,
    title_legend=None,
    title_bar=None,
    label_legend=None,
    aspect=None,
    scale_minimum=None,
    scale_center=None,
    scale_maximum=None,
    threshold_p_first=None,
    threshold_p_second=None,
    threshold_q_first=None,
    threshold_q_second=None,
    size_title_chart=None,
    size_title_abscissa=None,
    size_title_ordinate=None,
    size_title_legend=None,
    size_title_bar=None,
    size_label_abscissa=None,
    size_label_ordinate=None,
    size_label_legend=None,
    size_label_bar=None,
    size_label_significance_p=None,
    size_label_significance_q=None,
    color_scale_maximum=None,
    color_scale_center=None,
    color_scale_minimum=None,
    show_significance_p=None,
    show_significance_q=None,
    show_legend=None,
    show_scale_bar=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Date, revision or review: 25 March 2026

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
        path_file_source_list_sets=path_file_source_list_sets,
        path_file_source_list_groups=path_file_source_list_groups,
        name_batch=name_batch,
        prefix_negative=prefix_negative,
        prefix_positive=prefix_positive,
        suffix_file_source=suffix_file_source,
        column_set=column_set,
        column_score=column_score,
        column_p=column_p,
        name_chart=name_chart,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        title_legend=title_legend,
        title_bar=title_bar,
        label_legend=label_legend,
        aspect=aspect,
        scale_minimum=scale_minimum,
        scale_center=scale_center,
        scale_maximum=scale_maximum,
        threshold_p_first=threshold_p_first,
        threshold_p_second=threshold_p_second,
        threshold_q_first=threshold_q_first,
        threshold_q_second=threshold_q_second,
        size_title_chart=size_title_chart,
        size_title_abscissa=size_title_abscissa,
        size_title_ordinate=size_title_ordinate,
        size_title_legend=size_title_legend,
        size_title_bar=size_title_bar,
        size_label_abscissa=size_label_abscissa,
        size_label_ordinate=size_label_ordinate,
        size_label_legend=size_label_legend,
        size_label_bar=size_label_bar,
        size_label_significance_p=size_label_significance_p,
        size_label_significance_q=size_label_significance_q,
        color_scale_maximum=color_scale_maximum,
        color_scale_center=color_scale_center,
        color_scale_minimum=color_scale_minimum,
        show_significance_p=show_significance_p,
        show_significance_q=show_significance_q,
        show_legend=show_legend,
        show_scale_bar=show_scale_bar,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("system: local")
        print("package: partner")
        module = str(
            "plot_chart_heatmap_enrichments.py"
        )
        print(str("module: " + module))
        function = ("execute_procedure()")
        print(str("function: " + function))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_directory_dock=pail_parameters["path_directory_dock"],
        path_directory_dock_pail=pail_parameters["path_directory_dock_pail"],
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_file_source_list_sets=(
            pail_parameters["path_file_source_list_sets"]
        ),
        path_file_source_list_groups=(
            pail_parameters["path_file_source_list_groups"]
        ),
        prefix_negative=pail_parameters["prefix_negative"],
        prefix_positive=pail_parameters["prefix_positive"],
        suffix_file_source=pail_parameters["suffix_file_source"],
        report=pail_parameters["report"],
    )
    #pail_source["tables"]
    #pail_source["lists"]["sets"]
    #pail_source["lists"]["groups"]

    # Organize table of enrichments from all sets and groups in current batch.
    pail_wide = organize_table_enrichments_wide(
        pail_tables=pail_source["tables"],
        prefix_negative=pail_parameters["prefix_negative"],
        prefix_positive=pail_parameters["prefix_positive"],
        column_set=pail_parameters["column_set"],
        column_score=pail_parameters["column_score"],
        column_p=pail_parameters["column_p"],
        selection_sets=pail_source["lists"]["sets"],
        selection_groups=pail_source["lists"]["groups"],
        report=pail_parameters["report"],
    )

    ##########
    # Organize the table in the appropriate format for the function that will
    # plot the chart.
    # Organize table of enrichments in a long format.
    pail_long = organize_table_enrichments_long(
        table_wide_partial=pail_wide["table_wide_partial"],
        selection_sets=pail_source["lists"]["sets"],
        selection_groups=pail_source["lists"]["groups"],
        report=pail_parameters["report"],
    )

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Text.
    # Lists.
    # Tables.
    pail_write_tables = dict()
    pail_write_tables[str("table_enrichment_wide_partial")] = (
        pail_wide["table_wide_partial"]
    )
    pail_write_tables[str("table_enrichment_wide_full")] = (
        pail_wide["table_wide_full"]
    )
    pail_write_tables[str("table_enrichment_long_full")] = (
        pail_long["table_long_full"]
    )
    pail_write_tables[str("table_enrichment_long_partial")] = (
        pail_long["table_long_partial"]
    )
    ##########
    # Write product information to file.
    # Extract information from path to directory and file.
    #path_directory = os.path.dirname(path_file)
    #name_suffix_file = os.path.basename(path_file)
    #name_file, suffix_file = os.path.splitext(name_suffix_file)
    # Define paths to directories.
    path_directory_tables = os.path.join(
        pail_parameters["path_directory_product"], "tables",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_tables,
    )
    # Text.
    # Lists.
    # Tables.
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=path_directory_tables,
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )

    # Define paths to directories.
    path_directory_chart = os.path.join(
        pail_parameters["path_directory_product"], "charts"
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_chart,
    )

    # Create plot chart and write to file.
    create_write_plot_chart_heatmap_enrichment(
        path_directory_parent=path_directory_chart,
        name_chart=pail_parameters["name_chart"],
        table=pail_long["table_long_partial"],
        name_index_columns_abscissa="group",
        name_index_rows_ordinate="name_set",
        name_index_rows_type="type_value",
        type_value_signal="score",
        type_value_p="p_value",
        type_value_q="q_value",
        title_chart=pail_parameters["title_chart"],
        title_abscissa=pail_parameters["title_abscissa"],
        title_ordinate=pail_parameters["title_ordinate"],
        title_legend=pail_parameters["title_legend"],
        title_bar=pail_parameters["title_bar"],
        label_legend=pail_parameters["label_legend"],
        aspect=pail_parameters["aspect"],
        value_fill_missing=0.0,
        scale_minimum=pail_parameters["scale_minimum"],
        scale_center=pail_parameters["scale_center"],
        scale_maximum=pail_parameters["scale_maximum"],
        thresholds_p=[
            pail_parameters["threshold_p_first"],
            pail_parameters["threshold_p_second"],
        ],
        thresholds_q=[
            pail_parameters["threshold_q_first"],
            pail_parameters["threshold_q_second"],
        ],
        size_title_chart=pail_parameters["size_title_chart"],
        size_title_abscissa=pail_parameters["size_title_abscissa"],
        size_title_ordinate=pail_parameters["size_title_ordinate"],
        size_title_legend=pail_parameters["size_title_legend"],
        size_title_bar=pail_parameters["size_title_bar"],
        size_label_abscissa=pail_parameters["size_label_abscissa"],
        size_label_ordinate=pail_parameters["size_label_ordinate"],
        size_label_legend=pail_parameters["size_label_legend"],
        size_label_bar=pail_parameters["size_label_bar"],
        size_label_significance_p=pail_parameters["size_label_significance_p"],
        size_label_significance_q=pail_parameters["size_label_significance_q"],
        color_scale_maximum=pail_parameters["color_scale_maximum"],
        color_scale_center=pail_parameters["color_scale_center"],
        color_scale_minimum=pail_parameters["color_scale_minimum"],
        constrain_signal_values=False,
        fill_missing=True,
        show_significance_p=pail_parameters["show_significance_p"],
        show_significance_q=pail_parameters["show_significance_q"],
        show_legend=pail_parameters["show_legend"],
        show_scale_bar=pail_parameters["show_scale_bar"],
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
    path_file_source_list_sets = sys.argv[5]
    path_file_source_list_groups = sys.argv[6]
    name_batch = sys.argv[7]
    prefix_negative = sys.argv[8]
    prefix_positive = sys.argv[9]
    suffix_file_source = sys.argv[10]
    column_set = sys.argv[11]
    column_score = sys.argv[12]
    column_p = sys.argv[13]
    name_chart = sys.argv[14]
    title_chart = sys.argv[15]
    title_abscissa = sys.argv[16]
    title_ordinate = sys.argv[17]
    title_legend = sys.argv[18]
    title_bar = sys.argv[19]
    label_legend = sys.argv[20]
    aspect = sys.argv[21]
    scale_minimum = sys.argv[22]
    scale_center = sys.argv[23]
    scale_maximum = sys.argv[24]
    threshold_p_first = sys.argv[25]
    threshold_p_second = sys.argv[26]
    threshold_q_first = sys.argv[27]
    threshold_q_second = sys.argv[28]
    size_title_chart = sys.argv[29]
    size_title_abscissa = sys.argv[30]
    size_title_ordinate = sys.argv[31]
    size_title_legend = sys.argv[32]
    size_title_bar = sys.argv[33]
    size_label_abscissa = sys.argv[34]
    size_label_ordinate = sys.argv[35]
    size_label_legend = sys.argv[36]
    size_label_bar = sys.argv[37]
    size_label_significance_p = sys.argv[38]
    size_label_significance_q = sys.argv[39]
    color_scale_maximum = sys.argv[40]
    color_scale_center = sys.argv[41]
    color_scale_minimum = sys.argv[42]
    show_significance_p = sys.argv[43]
    show_significance_q = sys.argv[44]
    show_legend = sys.argv[45]
    show_scale_bar = sys.argv[46]
    report = sys.argv[47]

    # Call function for procedure.
    execute_procedure(
        path_directory_dock=path_directory_dock,
        path_directory_dock_pail=path_directory_dock_pail,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_file_source_list_sets=path_file_source_list_sets,
        path_file_source_list_groups=path_file_source_list_groups,
        name_batch=name_batch,
        prefix_negative=prefix_negative,
        prefix_positive=prefix_positive,
        suffix_file_source=suffix_file_source,
        column_set=column_set,
        column_score=column_score,
        column_p=column_p,
        name_chart=name_chart,
        title_chart=title_chart,
        title_abscissa=title_abscissa,
        title_ordinate=title_ordinate,
        title_legend=title_legend,
        title_bar=title_bar,
        label_legend=label_legend,
        aspect=aspect,
        scale_minimum=scale_minimum,
        scale_center=scale_center,
        scale_maximum=scale_maximum,
        threshold_p_first=threshold_p_first,
        threshold_p_second=threshold_p_second,
        threshold_q_first=threshold_q_first,
        threshold_q_second=threshold_q_second,
        size_title_chart=size_title_chart,
        size_title_abscissa=size_title_abscissa,
        size_title_ordinate=size_title_ordinate,
        size_title_legend=size_title_legend,
        size_title_bar=size_title_bar,
        size_label_abscissa=size_label_abscissa,
        size_label_ordinate=size_label_ordinate,
        size_label_legend=size_label_legend,
        size_label_bar=size_label_bar,
        size_label_significance_p=size_label_significance_p,
        size_label_significance_q=size_label_significance_q,
        color_scale_maximum=color_scale_maximum,
        color_scale_center=color_scale_center,
        color_scale_minimum=color_scale_minimum,
        show_significance_p=show_significance_p,
        show_significance_q=show_significance_q,
        show_legend=show_legend,
        show_scale_bar=show_scale_bar,
        report=report,
    )

    pass



#
