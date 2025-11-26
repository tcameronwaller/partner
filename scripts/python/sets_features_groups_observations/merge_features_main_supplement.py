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
# Date, initialization: 21 November 2025
# Date, review or revision: 21 November 2025
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
    path_directory_source=None,
    path_directory_source_features_supplement=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_main=None,
    path_file_source_table_reference_first=None,
    path_file_source_table_reference_second=None,
    path_file_source_table_supplement_first=None,
    path_file_source_table_supplement_second=None,
    column_main_identifier=None,
    column_main_identifier_supplement=None,
    column_main_name=None,
    column_reference_identifier=None,
    column_reference_name=None,
    column_supplement_feature=None,
    column_supplement_observation=None,
    transpose_table_supplement=None,
    report=None,
):
    """
    Parse parameters from text.

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files

        ...

        report (str): whether to print reports

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
    pail["path_directory_source_features_supplement"] = str(
        path_directory_source_features_supplement
    ).strip()

    # Paths to files.
    pail["path_file_source_table_main"] = str(
        path_file_source_table_main
    ).strip()
    pail["path_file_source_table_reference_first"] = str(
        path_file_source_table_reference_first
    ).strip()
    pail["path_file_source_table_reference_second"] = str(
        path_file_source_table_reference_second
    ).strip()
    pail["path_file_source_table_supplement_first"] = str(
        path_file_source_table_supplement_first
    ).strip()
    pail["path_file_source_table_supplement_second"] = str(
        path_file_source_table_supplement_second
    ).strip()

    # Names of columns.
    pail["column_main_identifier"] = str(
        column_main_identifier
    ).strip()
    pail["column_main_identifier_supplement"] = str(
        column_main_identifier_supplement
    ).strip()

    pail["column_main_name"] = str(
        column_main_name
    ).strip()
    pail["column_reference_identifier"] = str(
        column_reference_identifier
    ).strip()
    pail["column_reference_name"] = str(
        column_reference_name
    ).strip()
    pail["column_supplement_feature"] = str(
        column_supplement_feature
    ).strip()
    pail["column_supplement_observation"] = str(
        column_supplement_observation
    ).strip()

    # Boolean, true or false.
    # Iterate on individual of Boolean designations.
    designations = {
        "transpose_table_supplement": transpose_table_supplement,
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
            "merge_features_main_supplement.py"
        )
        print(str("module: " + module))
        function = str(
            "parse_text_parameters()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        print("parameters:")
        print(pail)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source_directory_files_sets_features(
    path_directory_parent=None,
    name_file_child_prefix=None,
    name_file_child_suffix=None,
    name_file_child_not=None,
    report=None,
):
    """
    Read and organize source information.

    Date, review or revision: TCW; 24 November 2025

    arguments:
        path_directory_parent (str): path to parent directory in which to find
            child files
        name_file_child_prefix (str): prefix in name by which to recognize
            relevant child files within parent directory
        name_file_child_suffix (str): suffix in name by which to recognize
            relevant child files within parent directory
        name_file_child_not (str): character string in names of files to exclude
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Extract and filter complete paths to child files within parent directory.
    paths = putly.extract_filter_child_file_names_paths(
        path_directory=path_directory_parent,
        name_file_prefix=name_file_child_prefix,
        name_file_suffix=name_file_child_suffix,
        name_file_not=name_file_child_not,
        report=report,
    )
    # Iterate on paths to files.
    # Read files and collect information.
    pail = dict()
    for path_file in paths:
        # Organize information for write.
        # Extract information from path to directory and file.
        path_directory = os.path.dirname(path_file)
        name_suffix_file = os.path.basename(path_file)
        #name_file = name_suffix_file.replace(str(name_file_child_suffix), "")
        name_file, suffix_file = os.path.splitext(name_suffix_file)
        # Read information from file.
        # Collect information.
        pail[name_file] = putly.read_file_text_list(
            path_file=path_file,
            delimiter="\n",
            unique=True,
        )
        pass

    # Create the union of all sets of features.
    # Explicit method.
    features_union = list()
    for key in pail.keys():
        features_union.extend(copy.deepcopy(pail[key]))
        pass
    pail["union_not_unique"] = features_union
    # Set operation, concise.
    pail["union"] = putly.combine_sets_items_union_unique(
        sets_items=list(pail.values()),
        report=None,
    )

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: merge_features_main_supplement.py")
        print("function: read_source_directory_files_sets_features()")
        putly.print_terminal_partition(level=5)
        print("sets of features")
        print("(set name: set size)")
        putly.print_terminal_partition(level=5)
        for key in pail.keys():
            count = len(pail[key])
            print(str(
                "set " + str(key) + ": " + str(count)
            ))
            pass
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source(
    path_directory_source=None,
    path_directory_source_features_supplement=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_main=None,
    path_file_source_table_reference_first=None,
    path_file_source_table_reference_second=None,
    path_file_source_table_supplement_first=None,
    path_file_source_table_supplement_second=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Date, review or revision: TCW; 24 November 2025

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files

        ...

        report (str): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters

    """

    # Read information from file.
    # Read sets of features in lists.
    pail_lists = read_source_directory_files_sets_features(
        path_directory_parent=path_directory_source_features_supplement,
        name_file_child_prefix="",
        name_file_child_suffix=".txt",
        name_file_child_not="blargharium_291367",
        report=report,
    )

    # Define paths to parent directories.
    path_directory_source_data = os.path.join(
        path_directory_source, "data",
    )
    path_directory_source_parameters = os.path.join(
        path_directory_source, "parameters",
    )

    # Define paths to child files.
    # Tables.
    #path_file_table_sample = os.path.join(
    #    path_directory_source_data, "table_sample.tsv",
    #)
    # Lists.

    # Determine whether paths point to a directory or file that exist.
    #existence_directory = os.path.exists(path_directory)
    existence_file_table_reference_first = os.path.exists(
        path_file_source_table_reference_first
    )
    existence_file_table_reference_second = os.path.exists(
        path_file_source_table_reference_second
    )
    existence_file_table_supplement_first = os.path.exists(
        path_file_source_table_supplement_first
    )
    existence_file_table_supplement_second = os.path.exists(
        path_file_source_table_supplement_second
    )

    # Read and collect information from files.
    pail_tables = dict()

    # Read information from file.
    # Tables.
    pail_tables["table_main"] = pandas.read_csv(
        path_file_source_table_main,
        sep="\t",
        header=0,
        #dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    if (existence_file_table_reference_first):
        pail_tables["table_reference_first"] = pandas.read_csv(
            path_file_source_table_reference_first,
            sep="\t",
            header=0,
            #dtype=types_columns,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        pail["table_reference_first"] = None
        pass
    if (existence_file_table_reference_second):
        pail_tables["table_reference_second"] = pandas.read_csv(
            path_file_source_table_reference_second,
            sep="\t",
            header=0,
            #dtype=types_columns,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        pail["table_reference_second"] = None
        pass

    if (existence_file_table_supplement_first):
        pail_tables["table_supplement_first"] = pandas.read_csv(
            path_file_source_table_supplement_first,
            sep="\t",
            header=0,
            #dtype=types_columns,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        pail["table_supplement_first"] = None
        pass
    if (existence_file_table_supplement_second):
        pail_tables["table_supplement_second"] = pandas.read_csv(
            path_file_source_table_supplement_second,
            sep="\t",
            header=0,
            #dtype=types_columns,
            na_values=[
                "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
            ],
            encoding="utf-8",
        )
    else:
        pail["table_supplement_second"] = None
        pass

    # Bundle information.
    pail = dict()
    pail["lists"] = pail_lists
    pail["tables"] = pail_tables

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: merge_features_main_supplement.py")
        print("function: read_source()")
        putly.print_terminal_partition(level=5)
        print("table_main")
        print(pail["tables"]["table_main"])
        pass
    # Return information.
    return pail


def organize_table_signal_first(
    table=None,
    name_index_features=None,
    name_index_observations=None,
    transpose_table=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 24 November 2025

    arguments:
        table (object): Pandas data-frame table
        ...
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): bundle of information

    """

    # Copy information.
    table = table.copy(deep=True)

    # Determine whether there is a table of supplemental features (signals) to
    # merge into the main table of features and observations.
    if (table is not None):
        # Optional preliminary transposition.
        # Copy information.
        table_format = table.copy(deep=True)
        # Determine whether to apply optional transposition.
        if (transpose_table):
            # Organize indices in table.
            table_format = (
                porg.explicate_table_indices_columns_rows_single_level(
                    table=table_format,
                    index_columns=name_index_observations,
                    index_rows=name_index_features,
                    explicate_indices=True,
                    report=report,
            ))
            # Transpose table.
            table_format = table_format.transpose(
                copy=True
            )
            # Organize indices in table.
            table_format.reset_index(
                level=None,
                inplace=True,
                drop=False, # remove index; do not move to regular columns
            )
            table_format.columns.rename(
                None,
                inplace=True,
            ) # single-dimensional index
        else:
            # Copy information.
            table_format = table.copy(deep=True)
            pass
        # Extract identifiers of features.
        identifiers_features = copy.deepcopy(
            table_format.columns.unique().tolist()
        )
        identifiers_features = list(filter(
            lambda identifier: (identifier != name_index_observations),
            identifiers_features
        ))
    else:
        # Fill null information.
        table_format = None
        identifiers_features = list()
        pass

    # Bundle information.
    pail = dict()
    pail["features_available"] = identifiers_features
    pail["table"] = table_format

    # Report.
    if report:
        # Organize.
        count_columns = table_format.shape[1]
        count_rows = table_format.shape[0]
        # Report.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: merge_features_main_supplement.py")
        print("function: organize_table_signal_first()")
        putly.print_terminal_partition(level=5)
        print("table after any transposition: ")
        print(table_format)
        putly.print_terminal_partition(level=5)
        print(str("count columns: " + str(count_columns)))
        print(str("count rows: " + str(count_rows)))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def filter_features_by_available_supplemental_information(
    sets_features=None,
    features_availability=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 24 November 2025

    arguments:

    TODO: update documentation

    raises:

    returns:
        (dict<object>): bundle of information
    """

    # Copy information.
    sets_features_original = copy.deepcopy(sets_features)
    features_availability = copy.deepcopy(features_availability)

    # Filter lists for sets of features by their availability in main or
    # supplemental tables.
    sets_features_novel = dict()
    for key in sets_features.keys():
        features_filter_unique = list(filter(
            lambda item: item in features_availability,
            copy.deepcopy(sets_features_original[key])
        ))
        features_filter_unique = putly.collect_unique_items(
            items=features_filter_unique,
        )
        sets_features_novel[key] = features_filter_unique
        pass

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        module = str(
            "merge_features_main_supplement.py"
        )
        print(str("module: " + module))
        function = str(
            "filter_features_by_available_supplemental_information()"
        )
        print("function: " + function)
        putly.print_terminal_partition(level=5)
        print("sets of features")
        print("(set name: set size)")
        putly.print_terminal_partition(level=5)
        for name in sets_features_novel.keys():
            count = len(sets_features_novel[name])
            print(str(
                "set " + str(name) + ": " + str(count)
            ))
            pass
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return sets_features_novel







def organize_table_signal_second(
    table=None,
    name_index_features=None,
    name_index_observations=None,
    features_selection=None,
    translations_features=None,
    report=None,
):
    """
    Blank.

    Review: TCW; 26 October 2025

    arguments:
        table_signal (object): Pandas data-frame table
        ...
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): bundle of information

    """

    # Copy information.
    table = table.copy(deep=True)
    features_selection = copy.deepcopy(features_selection)
    translations_features = copy.deepcopy(translations_features)

    # Extract identifiers of observations.
    identifiers_observations = copy.deepcopy(
        table.columns.unique().tolist()
    )
    identifiers_observations = list(filter(
        lambda identifier: (identifier != name_index_features),
        identifiers_observations
    ))
    # Filter columns and rows in table, corresponding to observations and
    # features, respectively.
    table_filter = porg.filter_select_table_columns_rows_by_identifiers(
        table=table,
        index_rows=name_index_features,
        identifiers_columns=identifiers_observations,
        identifiers_rows=features_selection,
        report=report,
    )
    # Copy information.
    table_format = table_filter.copy(deep=True)
    # Organize indices in table.
    table_format = (
        porg.explicate_table_indices_columns_rows_single_level(
            table=table_format,
            index_columns=name_index_observations,
            index_rows=name_index_features,
            explicate_indices=True,
            report=report,
    ))
    # Transpose table.
    table_format = table_format.transpose(
        copy=True
    )
    # Organize indices in table.
    table_format.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_format.columns.rename(
        None,
        inplace=True,
    ) # single-dimensional index
    # Translate names of columns.
    table_format.rename(
        columns=translations_features,
        inplace=True,
    )

    # Report.
    if report:
        # Organize.
        # Report.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: merge_features.py")
        print("function: organize_table_signal()")
        putly.print_terminal_partition(level=5)
        print(table_format)
        pass
    # Return information.
    return table_format




################################################################################
# Procedure


# TODO: TCW; 21 November 2025
# read in features from all files in source directory
# keep the features separate
# filter each feature set against the features available in the supplement table (RNAseq signals)
#   1. adipose RNAseq
#   2. muscle RNAseq
# write out the lists of features (Ensemble, symbols)
# take union of features from all sets
# Merge supplement signals for RNAseq features to the main table


##########
# Call main procedure.


def execute_procedure(
    path_directory_source=None,
    path_directory_source_features_supplement=None,
    path_directory_product=None,
    path_directory_dock=None,
    path_file_source_table_main=None,
    path_file_source_table_reference_first=None,
    path_file_source_table_reference_second=None,
    path_file_source_table_supplement_first=None,
    path_file_source_table_supplement_second=None,
    column_main_identifier=None,
    column_main_identifier_supplement=None,
    column_main_name=None,
    column_reference_identifier=None,
    column_reference_name=None,
    column_supplement_feature=None,
    column_supplement_observation=None,
    transpose_table_supplement=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review: TCW; 24 November 2025

    arguments:
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files

        ...

        report (str): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parse parameters.
    pail_parameters = parse_text_parameters(
        path_directory_source=path_directory_source,
        path_directory_source_features_supplement=(
            path_directory_source_features_supplement
        ),
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_main=path_file_source_table_main,
        path_file_source_table_reference_first=(
            path_file_source_table_reference_first
        ),
        path_file_source_table_reference_second=(
            path_file_source_table_reference_second
        ),
        path_file_source_table_supplement_first=(
            path_file_source_table_supplement_first
        ),
        path_file_source_table_supplement_second=(
            path_file_source_table_supplement_second
        ),
        column_main_identifier=column_main_identifier,
        column_main_identifier_supplement=column_main_identifier_supplement,
        column_main_name=column_main_name,
        column_reference_identifier=column_reference_identifier,
        column_reference_name=column_reference_name,
        column_supplement_feature=column_supplement_feature,
        column_supplement_observation=column_supplement_observation,
        transpose_table_supplement=transpose_table_supplement,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: merge_features_main_supplement.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("path_directory_source: " + str(path_directory_source))
        print("path_directory_product: " + str(path_directory_product))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_source_features_supplement=(
            pail_parameters["path_directory_source_features_supplement"]
        ),
        path_directory_product=pail_parameters["path_directory_product"],
        path_directory_dock=pail_parameters["path_directory_dock"],
        path_file_source_table_main=(
            pail_parameters["path_file_source_table_main"]
        ),
        path_file_source_table_reference_first=(
            pail_parameters["path_file_source_table_reference_first"]
        ),
        path_file_source_table_reference_second=(
            pail_parameters["path_file_source_table_reference_second"]
        ),
        path_file_source_table_supplement_first=(
            pail_parameters["path_file_source_table_supplement_first"]
        ),
        path_file_source_table_supplement_second=(
            pail_parameters["path_file_source_table_supplement_second"]
        ),
        report=pail_parameters["report"],
    )
    pail_source["lists"]
    pail_source["tables"]
    #pail_source["tables"]["table_main"]
    #pail_source["tables"]["table_reference_first"]
    #pail_source["tables"]["table_reference_second"]
    #pail_source["tables"]["table_supplement_first"]
    #pail_source["tables"]["table_supplement_second"]

    ##########
    # Organize table of supplemental features across observations.
    # Extract identifiers of available supplemental features.
    pail_supplement_first = organize_table_signal_first(
        table=pail_source["tables"]["table_supplement_first"],
        name_index_features=(
            pail_parameters["column_supplement_feature"]
        ),
        name_index_observations=(
            pail_parameters["column_supplement_observation"]
        ),
        transpose_table=pail_parameters["transpose_table_supplement"],
        report=pail_parameters["report"],
    )
    #pail_supplement_first["features_available"]
    #pail_supplement_first["table"]
    pail_supplement_second = organize_table_signal_first(
        table=pail_source["tables"]["table_supplement_second"],
        name_index_features=(
            pail_parameters["column_supplement_feature"]
        ),
        name_index_observations=(
            pail_parameters["column_supplement_observation"]
        ),
        transpose_table=pail_parameters["transpose_table_supplement"],
        report=pail_parameters["report"],
    )

    ##########
    # Filter sets of features by availability of supplemental information in
    # the respective tables.
    # Extract identifiers of features.
    features_main = copy.deepcopy(
        pail_source["tables"]["table_main"].columns.unique().tolist()
    )
    # Combine available features from main and supplemental tables.
    features_available_first = list()
    features_available_first.extend(copy.deepcopy(features_main))
    features_available_first.extend(copy.deepcopy(
        pail_supplement_first["features_available"]
    ))
    features_available_second = list()
    features_available_second.extend(copy.deepcopy(features_main))
    features_available_second.extend(copy.deepcopy(
        pail_supplement_second["features_available"]
    ))
    # Filter features in sets by their availability in the respective
    # supplemental table.
    pail_features_first = (
        filter_features_by_available_supplemental_information(
            sets_features=pail_source["lists"],
            features_availability=features_available_first,
            report=pail_parameters["report"],
    ))
    pail_features_second = (
        filter_features_by_available_supplemental_information(
            sets_features=pail_source["lists"],
            features_availability=features_available_second,
            report=pail_parameters["report"],
    ))


    # TODO: TCW; 24 November 2025
    # 1. write out the lists for sets of features
    # 2. prepare translations for names of features (with conditionals)
    # 3. finish organizing the supplemental tables
    # 4. merge supplements to the mains


    sys.exit()


    ##########
    # Organize table of signals for genes from RNAseq.

    # Extract information for translation of names of genes.
    # Adipose.
    table_gene_adipose = pail_source["table_gene_adipose"].copy(deep=True)
    table_gene_adipose["gene_name_prefix"] = table_gene_adipose.apply(
        lambda row: str(
            str("rnaseq_adipose_") + str(row["gene_name"])
        ),
        axis="columns", # apply function to each row
    )
    translations_genes_adipose = (
        porg.extract_translation_keys_values_from_table_columns(
            table=table_gene_adipose,
            column_keys="gene_identifier_base",
            column_values="gene_name_prefix",
            report=True,
    ))
    # Muscle.
    table_gene_muscle = pail_source["table_gene_muscle"].copy(deep=True)
    table_gene_muscle["gene_name_prefix"] = table_gene_muscle.apply(
        lambda row: str(
            str("rnaseq_muscle_") + str(row["gene_name"])
        ),
        axis="columns", # apply function to each row
    )
    translations_genes_muscle = (
        porg.extract_translation_keys_values_from_table_columns(
            table=table_gene_muscle,
            column_keys="gene_identifier_base",
            column_values="gene_name_prefix",
            report=True,
    ))

    # Organize table of signals.
    table_signal_adipose_format = organize_table_signal(
        table=pail_source["table_signal_adipose"],
        name_index_features="identifier_gene",
        name_index_observations="identifier_signal",
        features_selection=pail_source["features"]["rnaseq"],
        translations_features=translations_genes_adipose,
        report=pail_parameters["report"],
    )
    table_signal_muscle_format = organize_table_signal(
        table=pail_source["table_signal_muscle"],
        name_index_features="identifier_gene",
        name_index_observations="identifier_signal",
        features_selection=pail_source["features"]["rnaseq"],
        translations_features=translations_genes_muscle,
        report=pail_parameters["report"],
    )

    # Merge together features from table of signals for genes from RNAseq to
    # the main table of identifications, categories, and features of samples
    # from the clinic and laboratory.
    table_merge = porg.merge_columns_two_tables(
        identifier_first="identifier_signal",
        identifier_second="identifier_signal",
        table_first=table_sample,
        table_second=table_signal_adipose_format,
        preserve_index=False,
        report=report,
    )
    table_merge = porg.merge_columns_two_tables(
        identifier_first="identifier_signal",
        identifier_second="identifier_signal",
        table_first=table_merge,
        table_second=table_signal_muscle_format,
        preserve_index=False,
        report=report,
    )

    ##########
    # Bundle information.
    # Bundles of information for files.
    # Lists.
    pail_write_lists = dict()
    # Tables.
    pail_write_tables = dict()
    pail_write_tables["table_merge"] = table_merge

    ##########
    # Write product information to file.
    # Lists.
    # Tables.
    putly.write_tables_to_file(
        pail_write=pail_write_tables,
        path_directory=pail_parameters["path_directory_product"],
        reset_index_rows=False,
        write_index_rows=False,
        write_index_columns=True,
        type="text",
        delimiter="\t",
        suffix=".tsv",
    )

    pass


# Execute program process in Python.

if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_source = sys.argv[1]
    path_directory_source_features_supplement = sys.argv[2]
    path_directory_product = sys.argv[3]
    path_directory_dock = sys.argv[4]
    path_file_source_table_main = sys.argv[5]
    path_file_source_table_reference_first = sys.argv[6]
    path_file_source_table_reference_second = sys.argv[7]
    path_file_source_table_supplement_first = sys.argv[8]
    path_file_source_table_supplement_second = sys.argv[9]
    column_main_identifier = sys.argv[10]
    column_main_identifier_supplement = sys.argv[11]
    column_main_name = sys.argv[12]
    column_reference_identifier = sys.argv[13]
    column_reference_name = sys.argv[14]
    column_supplement_feature = sys.argv[15]
    column_supplement_observation = sys.argv[16]
    transpose_table_supplement = sys.argv[17]
    report = sys.argv[18]

    # Call function for procedure.
    execute_procedure(
        path_directory_source=path_directory_source,
        path_directory_source_features_supplement=(
            path_directory_source_features_supplement
        ),
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        path_file_source_table_main=path_file_source_table_main,
        path_file_source_table_reference_first=(
            path_file_source_table_reference_first
        ),
        path_file_source_table_reference_second=(
            path_file_source_table_reference_second
        ),
        path_file_source_table_supplement_first=(
            path_file_source_table_supplement_first
        ),
        path_file_source_table_supplement_second=(
            path_file_source_table_supplement_second
        ),
        column_main_identifier=column_main_identifier,
        column_main_identifier_supplement=column_main_identifier_supplement,
        column_main_name=column_main_name,
        column_reference_identifier=column_reference_identifier,
        column_reference_name=column_reference_name,
        column_supplement_feature=column_supplement_feature,
        column_supplement_observation=column_supplement_observation,
        transpose_table_supplement=transpose_table_supplement,
        report=report,
    )

    pass



#
