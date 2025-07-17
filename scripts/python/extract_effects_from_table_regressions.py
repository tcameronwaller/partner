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
# Date, first execution: 2 July 2025
# Date, last execution or modification: 16 July 2025
# Review: TCW; 16 July 2025
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
    path_file_source_table_regressions=None,
    path_file_source_table_features=None,
    name_file_product_table_effects=None,
    name_file_product_table_significance=None,
    name_file_product_ranks=None,
    name_file_product_list_change=None,
    name_file_product_list_negative=None,
    name_file_product_list_positive=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    column_effect_identifier=None,
    column_effect_estimate=None,
    column_effect_error=None,
    column_effect_p=None,
    column_feature_identifier=None,
    column_feature_name=None,
    columns_extra_keep=None,
    rate_false_discovery=None,
    report=None,
):
    """
    Parse parameters from text.

    arguments:
        path_file_source_table_regressions (str): path to source file
        path_file_source_table_features (str): path to source file
        name_file_product_table_effects (str): name of product file
        name_file_product_table_significance (str): name of product file
        name_file_product_ranks (str): name of product file
        name_file_product_list_change (str): name of product file
        name_file_product_list_negative (str): name of product file
        name_file_product_list_positive (str): name of product file
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        column_effect_identifier (str): name of column in source table
        column_effect_estimate (str): name of column in source table
        column_effect_error (str): name of column in source table
        column_effect_p (str): name of column in source table
        column_feature_identifier (str): name of column in source table
        column_feature_name (str): name of column in source table
        columns_extra_keep (list<str>): names of columns in source table
        rate_false_discovery (str): value for rate of acceptable false
            discoveries
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information
    """

    # Bundle information.
    pail = dict()
    # Parse information.
    pail["path_file_source_table_regressions"] = str(
        path_file_source_table_regressions
    ).strip()
    pail["path_file_source_table_features"] = str(
        path_file_source_table_features
    ).strip()
    pail["name_file_product_table_effects"] = str(
        name_file_product_table_effects
    ).strip()
    pail["name_file_product_table_significance"] = str(
        name_file_product_table_significance
    ).strip()
    pail["name_file_product_ranks"] = str(name_file_product_ranks).strip()
    pail["name_file_product_list_change"] = str(
        name_file_product_list_change
    ).strip()
    pail["name_file_product_list_positive"] = str(
        name_file_product_list_positive
    ).strip()
    pail["name_file_product_list_negative"] = str(
        name_file_product_list_negative
    ).strip()
    pail["path_directory_source"] = str(path_directory_source).strip()
    pail["path_directory_product"] = str(path_directory_product).strip()
    pail["path_directory_dock"] = str(path_directory_dock).strip()
    pail["column_effect_identifier"] = str(column_effect_identifier).strip()
    pail["column_effect_estimate"] = str(column_effect_estimate).strip()
    pail["column_effect_error"] = str(column_effect_error).strip()
    pail["column_effect_p"] = str(column_effect_p).strip()
    pail["column_feature_identifier"] = str(column_feature_identifier).strip()
    pail["column_feature_name"] = str(column_feature_name).strip()
    pail["columns_extra_keep"] = putly.parse_text_list_values(
        text=columns_extra_keep,
        delimiter=",",
    )
    pail["columns_extra_keep"] = putly.collect_unique_items(
        items=pail["columns_extra_keep"],
    )
    pail["rate_false_discovery"] = float(str(rate_false_discovery).strip())
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

    # Collect names of relevant columns in original source table.
    pail["columns_source_text"] = list()
    pail["columns_source_number"] = list()
    pail["columns_source_text"].append(pail["column_effect_identifier"])
    pail["columns_source_number"].append(pail["column_effect_estimate"])
    pail["columns_source_number"].append(pail["column_effect_error"])
    pail["columns_source_number"].append(pail["column_effect_p"])

    # Report.
    if pail["report"]:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: extract_effects_from_table_regressions.py")
        print("function: parse_text_parameters()")
        putly.print_terminal_partition(level=5)
        print("parameters:")
        print(pail)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def define_column_types_table_source_regression(
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
    path_file_source_table_regressions=None,
    path_file_source_table_features=None,
    path_directory_dock=None,
    columns_source_text=None,
    columns_source_number=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 2 July 2025

    arguments:
        path_file_source_table_regressions (str): path to source file
        path_file_source_table_features (str): path to source file
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        columns_source_text (list<str>): names of relevant columns in table
        columns_source_number (list<str>): names of relevant columns in table
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Define types of variables in columns of table.
    types_columns = define_column_types_table_source_regression(
        columns_source_text=columns_source_text,
        columns_source_number=columns_source_number,
    )
    # Determine paths point to files that exist.
    existence_regression = os.path.exists(path_file_source_table_regressions)
    existence_features = os.path.exists(path_file_source_table_features)
    # Read information from file.
    if (existence_regression):
        table_regression = pandas.read_csv(
            path_file_source_table_regressions,
            sep="\t",
            header=0,
            dtype=types_columns,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        table_regression = None
        pass
    if (existence_features):
        table_feature = pandas.read_csv(
            path_file_source_table_features,
            sep="\t",
            header=0,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        table_feature = None
        pass

    # Bundle information.
    pail = dict()
    pail["table_regression"] = table_regression
    pail["table_feature"] = table_feature

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: extract_effects_from_table_regressions.py")
        print("function: read_source_table_data()")
        putly.print_terminal_partition(level=5)
        print("table of results from regression:")
        print(pail["table_regression"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def organize_table_regression_effects(
    table_regression=None,
    table_feature=None,
    column_effect_identifier=None,
    column_effect_estimate=None,
    column_effect_error=None,
    column_effect_p=None,
    column_feature_identifier=None,
    column_feature_name=None,
    columns_extra_keep=None,
    rate_false_discovery=None,
    report=None,
):
    """
    Organize information from regressions, with effects across many entities or
    features.

    Review: TCW; 17 July 2025

    arguments:
        table_regression (object): Pandas data-frame table
        table_feature (object): Pandas data-frame table
        column_effect_identifier (str): name of column in source table
        column_effect_estimate (str): name of column in source table
        column_effect_error (str): name of column in source table
        column_effect_p (str): name of column in source table
        column_feature_identifier (str): name of column in source table
        column_feature_name (str): name of column in source table
        columns_extra_keep (list<str>): names of columns in source table
        rate_false_discovery (float): value for acceptable rate of false
            discoveries
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy information.
    table_regression = table_regression.copy(deep=True)
    table_feature = table_feature.copy(deep=True)
    #table_feature = table_feature.copy(deep=True)
    columns_extra_keep = copy.deepcopy(columns_extra_keep)

    # Determine whether there is a table with supplemental information about
    # the entities or features that correspond to the effects from regressions.
    if (table_feature is not None):
        # Extract names of columns.
        columns_feature_available = copy.deepcopy(
            table_feature.columns.to_list()
        )
        columns_feature_keep = list(filter(
            lambda feature: (feature in columns_extra_keep),
            columns_feature_available
        ))
        # Copy identifier before merge.
        table_regression["identifier_merge"] = (
            table_regression[column_effect_identifier]
        )
        #columns_relevant.insert(0, "identifier_merge")
        # Merge information about features to information about the effects
        # from regressions.
        # Notice that the merge operation does not preserve the column for the
        # identifier in the second table.
        table_merge = porg.merge_columns_two_tables(
            identifier_first="identifier_merge",
            identifier_second=column_feature_identifier,
            table_first=table_regression,
            table_second=table_feature,
            preserve_index=False,
            report=report,
        )
        # Remove from list of relevant features the name of identifiers for
        # features in the second table.
        #if (column_feature_identifier in columns_relevant):
        #    columns_relevant = list(filter(
        #        lambda feature: (feature != column_feature_identifier),
        #        columns_relevant
        #    ))
        #    pass
    else:
        # Copy information.
        table_merge = table_regression.copy(deep=True)
        columns_feature_keep = list()
        pass

    # Translate names of columns.
    translations = dict()
    translations[column_effect_identifier] = "effect_identifier"
    translations[column_effect_estimate] = "effect_estimate"
    translations[column_effect_error] = "effect_error"
    translations[column_effect_p] = "effect_p"
    table_merge.rename(
        columns=translations,
        inplace=True,
    )

    # Filter rows in table.
    table_merge = table_merge.loc[
        (
            (pandas.notna(table_merge["effect_estimate"])) &
            (pandas.notna(table_merge["effect_error"])) &
            (pandas.notna(table_merge["effect_p"])) &
            (table_merge["effect_p"] >= float(0))
        ), :
    ]
    # Fill values of zero for p-value.
    table_merge["effect_p_fill"] = table_merge.apply(
        lambda row:
            2.23E-308 if (
                float(row["effect_p"]) < 2.23E-308
            ) else row["effect_p"],
        axis="columns", # apply function to each row
    )

    # Calculate Benjamini-Hochberg q-values for False-Discovery Rate (FDR).
    # Calculate q-values across all comparisons in table.
    # FDR 5% (q <= 0.05).
    table_merge = pdesc.calculate_table_false_discovery_rate_q_values(
        threshold=rate_false_discovery, # alpha; family-wise error rate
        name_column_p_value="effect_p_fill",
        name_column_q_value="effect_q",
        name_column_significance="effect_q_significance",
        table=table_merge,
    )

    # Calculate the negative, base-ten logarithm of the p-value.
    #table_change["p_value_negative_log10"] = numpy.log10(
    #    table_change["p_value"]
    #)
    table_merge["effect_p_negative_log10"] = table_merge.apply(
        lambda row: (-1*math.log(row["effect_p_fill"], 10)),
        axis="columns", # apply function to each row
    )

    # Calculate a metric by which to rank features.
    # For differential expression in genes, this metric is useful for the
    # 'preranked' analysis in gene set enrichment (GSEA).
    # Use the p-value rather than the q-value since the q-value tends to have
    # ties between groups of comparisons.
    table_merge["rank_effect_p"] = table_merge.apply(
        lambda row: (
            (row["effect_estimate"])*(row["effect_p_negative_log10"])
        ),
        axis="columns", # apply function to each row
    )

    # Sort rows within table.
    table_merge.sort_values(
        by=[
            "rank_effect_p",
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )

    # Organize sequence of columns.
    columns_core = list()
    columns_core.append("sequence")
    columns_core.append("group")
    columns_core.append("name")
    columns_core.append("name_combination")
    columns_effect = list()
    columns_effect.append("effect_identifier")
    columns_effect.append("effect_estimate")
    columns_effect.append("effect_error")
    columns_effect.append("effect_p")
    columns_effect.append("effect_p_fill")
    columns_effect.append("effect_p_negative_log10")
    columns_effect.append("effect_q")
    columns_effect.append("effect_q_significance")
    columns_effect.append("rank_effect_p")
    columns_feature = list()
    columns_feature.append(column_feature_name)
    columns_feature.extend(columns_feature_keep)
    columns_relevant = list()
    columns_relevant.extend(columns_core)
    columns_relevant.extend(columns_effect)
    columns_relevant.extend(columns_feature)
    columns_relevant.extend(columns_extra_keep)
    columns_sequence_priority = list()
    columns_sequence_priority.extend(columns_core)
    columns_sequence_priority.extend(columns_effect)
    columns_sequence_priority.extend(columns_feature)
    columns_relevant = putly.collect_unique_items(
        items=columns_relevant,
    )
    columns_sequence_priority = putly.collect_unique_items(
        items=columns_sequence_priority,
    )
    # Filter and sort columns in table.
    table_merge = porg.filter_sort_table_columns(
        table=table_merge,
        columns_sequence=columns_relevant,
        report=report,
    )
    # Sort sequence of columns in table.
    table_merge = porg.sort_table_columns_explicit_other(
        table=table_merge,
        columns_sequence=columns_sequence_priority,
        report=report,
    )

    # Bundle product information.
    pail = dict()
    pail["table_merge"] = table_merge

    # Report.
    if report:
        # Organize.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: extract_effects_from_table_regressions.py")
        print("function: organize_table_regression_effects()")
        putly.print_terminal_partition(level=5)
        print("product regression table:")
        print(table_merge)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def extract_significant_effects(
    table_regression=None,
    column_effect_identifier=None,
    column_effect_estimate=None,
    column_effect_q=None,
    column_effect_q_significance=None,
    rate_false_discovery=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        table_regression (object): Pandas data-frame table
        column_effect_identifier (str): name of column in table for the unique
            identifier corresponding to the effect
        column_effect_estimate (str): name of column in table
        column_effect_q (str): name of column in table
        column_effect_q_significance (str): name of column in table
        rate_false_discovery (float): value for acceptable rate of false
            discoveries
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Copy information.
    table_regression = table_regression.copy(deep=True)

    # Filter rows within table.
    table_significance = table_regression.loc[
        (
            (table_regression[column_effect_q_significance])
        ), :
    ].copy(deep=True)

    # Extract identifiers of effects that were significant.
    pail_threshold = porg.segregate_effects_by_thresholds(
        table=table_regression,
        column_effect=column_effect_estimate,
        column_significance=column_effect_q,
        threshold_effect=0.0,
        threshold_significance=rate_false_discovery,
        report=False,
    )
    # Extract identifiers genes with differential expression beyond thresholds.
    effects_significance = copy.deepcopy(
        pail_threshold["table_pass_any"][column_effect_identifier].to_list()
    )
    effects_positive = copy.deepcopy(
        pail_threshold["table_pass_up"][column_effect_identifier].to_list()
    )
    effects_negative = copy.deepcopy(
        pail_threshold["table_pass_down"][column_effect_identifier].to_list()
    )
    effects_significance = putly.collect_unique_items(
        items=effects_significance,
    )
    effects_positive = putly.collect_unique_items(
        items=effects_positive,
    )
    effects_negative = putly.collect_unique_items(
        items=effects_negative,
    )

    # Bundle product information.
    pail = dict()
    pail["table_significance"] = table_significance
    pail["effects_significance"] = effects_significance
    pail["effects_positive"] = effects_positive
    pail["effects_negative"] = effects_negative

    # Report.
    if report:
        # Organize.
        count_significance = int(len(effects_significance))
        count_positive = int(len(effects_positive))
        count_negative = int(len(effects_negative))
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: extract_effects_from_table_regressions.py")
        print("function: extract_significant_effects()")
        putly.print_terminal_partition(level=5)
        print("count significant effects: " + str(count_significance))
        print("count positive effects: " + str(count_positive))
        print("count negative effects: " + str(count_negative))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def organize_table_rank_effects(
    table_regression=None,
    identifiers_exclusion=None,
    column_identifier=None,
    column_rank=None,
    report=None,
):
    """
    Prepare a table from which to extract a ranked list of effects for
    subsequent enrichment analysis. An example use case is the preranked
    analysis by gene set enrichment analysis (GSEA).

    Exercise caution in the decision to exclude effects from the list in rank
    format. Simple exclusion might disrupt the calculations of background
    statistics in GSEA or other methods. Instead of excluding effects from the
    list, consider nullifying p-values as an alternative.

    arguments:
        table_regression (object): Pandas data-frame table
        column_identifier (str): name of column in table for the unique
            identifier corresponding to the effect
        column_rank (str): name of column in table for the rank of the effect
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table_regression = table_regression.copy(deep=True)

    # Define relevant columns in sequence.
    columns_sequence = [
        column_identifier,
        column_rank,
    ]
    # Organize indices in table.
    table_regression.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Determine whether to exclude specific identifiers from selection.
    if (
        (identifiers_exclusion is not None) and
        (len(identifiers_exclusion) > 0)
    ):
        # Filter rows in table to exclude specific genes.
        table_regression = table_regression.loc[
            (
                ~table_regression[column_identifier].isin(
                    identifiers_exclusion
                )
            ), :
        ].copy(deep=True)
        pass

    # Filter rows in table for non-missing values across relevant columns.
    table_regression.dropna(
        axis="index",
        how="any",
        subset=columns_sequence,
        inplace=True,
    )

    # Remove redundant rows from table on basis of identifier.
    table_regression.drop_duplicates(
        subset=[
            column_identifier,
        ],
        keep="first",
        inplace=True,
        ignore_index=True,
    )

    # Sort rows within table.
    table_regression.sort_values(
        by=[
            column_rank,
        ],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )

    # Organize indices in table.
    table_regression.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    # Filter and sort columns within table.
    table_rank = porg.filter_sort_table_columns(
        table=table_regression,
        columns_sequence=columns_sequence,
        report=report,
    )

    # Collect information.
    pail = dict()
    pail["table_rank"] = table_rank
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: extract_effects_from_table_regressions.py")
        print("function: organize_table_rank_effects()")
        putly.print_terminal_partition(level=5)
        print(table_rank)
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return pail



################################################################################
# Procedure


##########
# Call main procedure.


def execute_procedure(
    path_file_source_table_regressions=None,
    path_file_source_table_features=None,
    name_file_product_table_effects=None,
    name_file_product_table_significance=None,
    name_file_product_ranks=None,
    name_file_product_list_change=None,
    name_file_product_list_negative=None,
    name_file_product_list_positive=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    column_effect_identifier=None,
    column_effect_estimate=None,
    column_effect_error=None,
    column_effect_p=None,
    column_feature_identifier=None,
    column_feature_name=None,
    columns_extra_keep=None,
    rate_false_discovery=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review: TCW; 15 July 2025

    arguments:
        path_file_source_table_regressions (str): path to source file
        path_file_source_table_features (str): path to source file
        name_file_product_table_effects (str): name of product file
        name_file_product_table_significance (str): name of product file
        name_file_product_ranks (str): name of product file
        name_file_product_list_change (str): name of product file
        name_file_product_list_negative (str): name of product file
        name_file_product_list_positive (str): name of product file
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        column_effect_identifier (str): name of column in source table
        column_effect_estimate (str): name of column in source table
        column_effect_error (str): name of column in source table
        column_effect_p (str): name of column in source table
        column_feature_identifier (str): name of column in source table
        column_feature_name (str): name of column in source table
        columns_extra_keep (list<str>): names of columns in source table
        rate_false_discovery (str): value for acceptable rate of false
            discoveries
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parse parameters.
    pail_parameters = parse_text_parameters(
        path_file_source_table_regressions=path_file_source_table_regressions,
        path_file_source_table_features=path_file_source_table_features,
        name_file_product_table_effects=name_file_product_table_effects,
        name_file_product_table_significance=(
            name_file_product_table_significance
        ),
        name_file_product_ranks=name_file_product_ranks,
        name_file_product_list_change=(
            name_file_product_list_change
        ),
        name_file_product_list_negative=name_file_product_list_negative,
        name_file_product_list_positive=name_file_product_list_positive,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        column_effect_identifier=column_effect_identifier,
        column_effect_estimate=column_effect_estimate,
        column_effect_error=column_effect_error,
        column_effect_p=column_effect_p,
        column_feature_identifier=column_feature_identifier,
        column_feature_name=column_feature_name,
        columns_extra_keep=columns_extra_keep,
        rate_false_discovery=rate_false_discovery,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: extract_effects_from_table_regressions.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print(
            "path_file_source_table_regressions: " +
            pail_parameters["path_file_source_table_regressions"]
        )
        print("path_directory_product: " + str(path_directory_product))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_file_source_table_regressions=(
            pail_parameters["path_file_source_table_regressions"]
        ),
        path_file_source_table_features=(
            pail_parameters["path_file_source_table_features"]
        ),
        path_directory_dock=path_directory_dock,
        columns_source_text=pail_parameters["columns_source_text"],
        columns_source_number=pail_parameters["columns_source_number"],
        report=pail_parameters["report"],
    )

    # Organize information in table.
    pail_effect = organize_table_regression_effects(
        table_regression=pail_source["table_regression"],
        table_feature=pail_source["table_feature"],
        column_effect_identifier=pail_parameters["column_effect_identifier"],
        column_effect_estimate=pail_parameters["column_effect_estimate"],
        column_effect_error=pail_parameters["column_effect_error"],
        column_effect_p=pail_parameters["column_effect_p"],
        column_feature_identifier=pail_parameters["column_feature_identifier"],
        column_feature_name=pail_parameters["column_feature_name"],
        columns_extra_keep=pail_parameters["columns_extra_keep"],
        rate_false_discovery=pail_parameters["rate_false_discovery"],
        report=pail_parameters["report"],
    )

    # Determine whether to extract identifiers of effects or names of the
    # corresponding entities or features.
    if (pail_source["table_feature"] is not None):
        column_extraction = column_feature_name
    else:
        column_extraction = "effect_identifier"
        pass

    # Extract identifiers of significant positive and negative effects.
    pail_change = extract_significant_effects(
        table_regression=pail_effect["table_merge"],
        column_effect_identifier=column_extraction,
        column_effect_estimate="effect_estimate",
        column_effect_q="effect_q",
        column_effect_q_significance="effect_q_significance",
        rate_false_discovery=pail_parameters["rate_false_discovery"],
        report=pail_parameters["report"],
    )

    # Prepare information for export and analysis using the preranked algorithm
    # in GSEA.
    pail_rank = organize_table_rank_effects(
        table_regression=pail_effect["table_merge"],
        identifiers_exclusion=list(),
        column_identifier=column_extraction,
        column_rank="rank_effect_p", # (effect * p_value_negative_log10)
        report=report,
    )

    # TODO: TCW; 2 July 2025
    # 1. organize the rank file export for GSEA
    # 2. create the volcano plot chart

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Lists.
    pail_write_lists = dict()
    pail_write_lists[name_file_product_list_change] = (
        pail_change["effects_significance"]
    )
    pail_write_lists[name_file_product_list_negative] = (
        pail_change["effects_negative"]
    )
    pail_write_lists[name_file_product_list_positive] = (
        pail_change["effects_positive"]
    )
    # Tables.
    pail_write_tables = dict()
    pail_write_tables[name_file_product_table_effects] = (
        pail_effect["table_merge"]
    )
    pail_write_tables[name_file_product_table_significance] = (
        pail_change["table_significance"]
    )
    # Ranks.
    pail_write_ranks = dict()
    pail_write_ranks[name_file_product_ranks] = (
        pail_rank["table_rank"]
    )

    ##########
    # Write product information to file.
    # Lists.
    putly.write_lists_to_file_text(
        pail_write=pail_write_lists,
        path_directory=path_directory_product,
        delimiter="\n",
    )
    # Tables.
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=path_directory_product,
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )
    putly.write_tables_to_file(
        pail_write=pail_write_ranks,
        path_directory=path_directory_product,
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=False,
        type="text",
        delimiter="\t",
        suffix=".rnk",
    )
    pass


# Execute program process in Python.


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_source_table_regressions = sys.argv[1]
    path_file_source_table_features = sys.argv[2]
    name_file_product_table_effects = sys.argv[3]
    name_file_product_table_significance = sys.argv[4]
    name_file_product_ranks = sys.argv[5]
    name_file_product_list_change = sys.argv[6]
    name_file_product_list_negative = sys.argv[7]
    name_file_product_list_positive = sys.argv[8]
    path_directory_source = sys.argv[9]
    path_directory_product = sys.argv[10]
    path_directory_dock = sys.argv[11]
    column_effect_identifier = sys.argv[12]
    column_effect_estimate = sys.argv[13]
    column_effect_error = sys.argv[14]
    column_effect_p = sys.argv[15]
    column_feature_identifier = sys.argv[16]
    column_feature_name = sys.argv[17]
    columns_extra_keep = sys.argv[18]
    rate_false_discovery = sys.argv[19]
    report = sys.argv[20]

    # Call function for procedure.
    execute_procedure(
        path_file_source_table_regressions=path_file_source_table_regressions,
        path_file_source_table_features=path_file_source_table_features,
        name_file_product_table_effects=name_file_product_table_effects,
        name_file_product_table_significance=(
            name_file_product_table_significance
        ),
        name_file_product_ranks=name_file_product_ranks,
        name_file_product_list_change=name_file_product_list_change,
        name_file_product_list_negative=name_file_product_list_negative,
        name_file_product_list_positive=name_file_product_list_positive,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        column_effect_identifier=column_effect_identifier,
        column_effect_estimate=column_effect_estimate,
        column_effect_error=column_effect_error,
        column_effect_p=column_effect_p,
        column_feature_identifier=column_feature_identifier,
        column_feature_name=column_feature_name,
        columns_extra_keep=columns_extra_keep,
        rate_false_discovery=rate_false_discovery,
        report=report,
    )

    pass



#
