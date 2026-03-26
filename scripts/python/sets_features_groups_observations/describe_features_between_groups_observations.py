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
    Copyright (C) 2026 Thomas Cameron Waller

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
# Date, initialization: 22 January 2026
# Date, review or revision: 22 January 2026
################################################################################
# Note

# The specialty of this Python script is to describe quantitative features
# and compare them between groups of observations.

# TODO: TCW; 22 January 2026
# Consider including "list_features_category" along with "list_features_quantity",
# and support different types of descriptions for each.

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
    path_file_source_table_features_observations=None,
    path_file_source_table_features=None,
    path_file_source_list_features=None,
    path_file_source_table_groups_observations=None,
    column_identifier_observation=None,
    column_name_observation=None,
    column_identifier_feature=None,
    column_name_feature=None,
    proportion_nonmissing_observations=None,
    allow_replicate_observations=None,
    report=None,
):
    """
    Parse parameters from text.

    arguments:
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_directory_dock_pail (str): path to pail directory for procedure's
            source and product directories and files
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_file_source_table_features_observations (str): path to source file
        path_file_source_table_features (str): path to source file
        path_file_source_list_features (str): path to source file
        path_file_source_table_groups_observations (str): path to source file
        column_identifier_observation (str): name of column in source table
        column_name_observation (str): name of column in source table
        column_identifier_feature (str): name of column in source table
        column_name_feature (str): name of column in source table
        proportion_nonmissing_observations (float): threshold by proportion of
            observations that must have nonmissing values for each feature
        allow_replicate_observations (bool): whether to allow replicate
            observations or to require groups to be mutually exclusive, such
            that any individual observation can only belong to one group
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information
    """

    # Bundle information.
    pail = dict()

    # Parse information.

    # Paths to directories and files.
    pail["path_directory_dock"] = str(path_directory_dock).strip()
    pail["path_directory_dock_pail"] = str(path_directory_dock_pail).strip()
    pail["path_directory_source"] = str(path_directory_source).strip()
    pail["path_directory_product"] = str(path_directory_product).strip()
    pail["path_file_source_table_features_observations"] = str(
        path_file_source_table_features_observations
    ).strip()
    pail["path_file_source_table_features"] = str(
        path_file_source_table_features
    ).strip()
    pail["path_file_source_list_features"] = str(
        path_file_source_list_features
    ).strip()
    pail["path_file_source_table_groups_observations"] = str(
        path_file_source_table_groups_observations
    ).strip()

    # Names and categories.
    # It is problematic to pass any white space in parameters from a script in
    # Bash. Designate the hash symbol "#" as a substitute for white space.
    # It is also problematic to pass an empty string in parameters from a
    # script in Bash. Designate the word "none" as a substitute for missing or
    # empty.
    # Iterate on individual parameters for names and categories.
    names = {
        "column_identifier_observation": column_identifier_observation,
        "column_name_observation": column_name_observation,
        "column_identifier_feature": column_identifier_feature,
        "column_name_feature": column_name_feature,
    }
    for key_name in names.keys():
        # Determine whether parameter has a valid value that is not none.
        if (
            (str(names[key_name]).strip().lower() != "none")
        ):
            # Parse value.
            pail[key_name] = str(
                names[key_name]
            ).strip().replace("#", " ")
        else:
            pail[key_name] = ""
            pass
        pass

    # Numbers.
    # Iterate on individual parameters for numbers.
    numbers = {
        "proportion_nonmissing_observations": (
            proportion_nonmissing_observations
        ),
    }
    nulls = [
        "na", "nan", "NA", "NAN",
    ]
    for key_number in numbers.keys():
        # Determine whether parameter has a valid value.
        if (
            (numbers[key_number] is not None) and
            (len(str(numbers[key_number])) > 0) and
            (str(numbers[key_number]) != "") and
            (str(numbers[key_number]).strip().lower() != "none") and
            (str(numbers[key_number]).strip().lower() not in nulls)
        ):
            # Number is valid.
            pail[key_number] = float(str(numbers[key_number]).strip())
        else:
            # Number is missing or null.
            pail[key_number] = None
            pass
        pass

    # Lists.

    # Booleans, true or false.
    # Iterate on individual of Boolean designations.
    designations = {
        "allow_replicate_observations": allow_replicate_observations,
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
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("describe_features_between_groups_observations.py")
        print("module: " + name_module)
        name_function = str("parse_text_parameters()")
        print("function: " + name_function)
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
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_features_observations=None,
    path_file_source_table_features=None,
    path_file_source_list_features=None,
    path_file_source_table_groups_observations=None,
    report=None,
):
    """
    Read and organize source information.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Date, review or revision: 22 January 2026

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
    existence_file_features = os.path.exists(path_file_source_table_features)

    # Bundle information.
    pail = dict()

    # Read information from file.
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
    if (existence_file_features):
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
    else:
        pail["table_features"] = None
        pass
    pail["list_features"] = putly.read_file_text_list(
        path_file=path_file_source_list_features,
        delimiter="\n",
        unique=True,
    )
    pail["table_groups_observations"] = pandas.read_csv(
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
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("describe_features_between_groups_observations.py")
        print("module: " + name_module)
        name_function = str("read_source()")
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Organize parameters.


def organize_parameters_features_observations(
    table_features_observations=None,
    table_features=None,
    list_features=None,
    table_groups_observations=None,
    column_identifier_feature=None,
    column_name_feature=None,
    column_identifier_observation=None,
    column_name_observation=None,
    report=None,
):
    """
    Organize parameters and information about features.

    Date, review or revision: 22 January 2026

    arguments:
        table_features_observations (object): Pandas data-frame table
        table_features (object): Pandas data-frame table
        list_features (list<str>): names of features
        ...
        column_identifier_feature (str): name of column in source table
        column_name_feature (str): name of column in source table
        ...
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): bundle of information

    """

    # Copy information.
    table_features_observations.copy(deep=True)
    if (table_features is not None):
        table_features = table_features.copy(deep=True)
        pass
    list_features = copy.deepcopy(list_features)
    table_groups_observations.copy(deep=True)

    # Available features.
    # Extract identifiers of features.
    features_observations = copy.deepcopy(
        table_features_observations.columns.unique().tolist()
    )
    # Filter to features with available signals.
    features_selection = list(filter(
        lambda feature: (feature in features_observations),
        copy.deepcopy(list_features)
    ))
    # Collect unique features.
    features_selection = putly.collect_unique_items(
        items=features_selection,
    )
    # Information about features and translations of their names.
    pail_features = sutly.organize_parameters_table_features_translations(
        table_features=table_features,
        column_identifier_feature=column_identifier_feature,
        column_name_feature=column_name_feature,
        features_selection=features_selection,
        prefix_name_feature="",
        report=report,
    )
    #pail_features["table_features_selection"]
    #pail_features["features_selection"]
    #pail_features["features_selection_translation"]
    #pail_features["translations_features"]
    # Organize parameters about groups of observations.
    pail_groups = sutly.organize_parameters_groups_observations(
        table=table_groups_observations,
        column_name="abbreviation",
        report=report,
    )
    #pail_groups["table"]
    #pail_groups["names_groups_observations_sequence"]
    #pail_groups["categories_groups"]
    #pail_groups["records"]
    pail_observations = sutly.organize_parameters_further_groups_observations(
        table_observations=table_features_observations,
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        column_identifier_groups_observations=column_identifier_observation,
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

    # Bundle information.
    pail = dict()
    pail["pail_features"] = pail_features
    pail["pail_groups"] = pail_groups
    pail["pail_observations"] = pail_observations

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("describe_features_between_groups_observations.py")
        print("module: " + name_module)
        name_function = str("organize_parameters_features_observations()")
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Stratify groups of observations.


def collect_entries_tables_features_groups_observations(
    table=None,
    name_index_columns=None,
    name_index_rows=None,
    names_groups_observations_sequence=None,
    groups_observations=None,
    column_group=None,
    replicate_groups=None,
    report=None,
):
    """
    Define and stratify groups of observations in a table.

    From a source table of columns for features and rows for observations, this
    function organizes the product information both as partial, separate,
    stratified tables within a dictionary and as a whole table with a new
    column that designates groups of observations.

    For the collection of stratified tables within the dictionary
    'entries_tables', any individual observation can belong to multiple groups
    and can appear in multiple stratified tables accordingly.

    For the table 'table_group', it is optional to replicate rows for
    individual observations that occur in multiple groups. It is important to
    be aware and cautious when using this option of replication of records for
    observations. This option can be convenient for describing features in
    groups of observations that overlap.

    ----------
    Format of source data table (name: 'table')
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. For versatility, this table does not have
    explicitly defined indices across rows or columns.
    ----------
    identifiers     feature_1 feature_2 feature_3 feature_4 feature_5 ...

    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    ----------
    Format of product data tables, partial (name: 'table')
    ----------
    Format of product, partial data tables organized within entries of a
    dictionary is in wide format with features across columns and values
    corresponding to their observations across rows. A special
    header row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. The table has explicitly named indices across
    columns and rows.
    ----------
    features        feature_1 feature_2 feature_3 feature_4 feature_5 ...
    observations
    observation_1   0.001     0.001     0.001     0.001     0.001     ...
    observation_2   0.001     0.001     0.001     0.001     0.001     ...
    observation_3   0.001     0.001     0.001     0.001     0.001     ...
    observation_4   0.001     0.001     0.001     0.001     0.001     ...
    observation_5   0.001     0.001     0.001     0.001     0.001     ...
    ----------

    ----------
    Format of product data table, whole (name: 'table_group')
    ----------
    Format of product, whole data table is in wide format with features across
    columns and values corresponding to their observations across rows. A
    special header row gives identifiers or names corresponding to each feature
    across columns, and a special column gives identifiers or names
    corresponding to each observation across rows. Another special column
    provides names of categorical groups for these observations. The table has
    explicitly named indices across columns and rows.
    ----------
    features        group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_2   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_3   0.001     0.001     0.001     0.001     0.001
    ----------

    Date, review or revision: 22 January 2026
    Date, review or revision: 25 July 2025
    Date, review or revision: 14 April 2025

    arguments:
        table (object): Pandas data-frame table with columns for features and
            rows for observations
        name_index_columns (str): name for single-level index across columns
        name_index_rows (str): name for single-level index across rows
        names_groups_observations_sequence (list<str>): names of groups for
            observation in specific sequence
        groups_observations (dict<list<str>>): sets of observations in groups
        column_group (str): name for column to designate groups of observations
        replicate_groups (bool): whether to replicate records or rows in table
            for individual observations that belong to multiple groups
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information
            groups_sequence (list<str>): unique names of groups of observations
                corresponding to the separate tables in the original sequence
                of instances of parameters
            entries_tables (dict<object>): entries with names as keys and as
                values, Pandas data-frame tables of data with features and
                observations for analysis
            table_group (object): Pandas data-frame table of data with features
                and observations for analysis

    """

    # Copy information.
    table_source = table.copy(deep=True)
    names_groups_observations_sequence = copy.deepcopy(
        names_groups_observations_sequence
    )
    groups_observations = copy.deepcopy(groups_observations)

    # Available features.
    # Extract identifiers of features.
    columns_all = copy.deepcopy(
        table_source.columns.unique().tolist()
    )

    # Collect information.
    groups_sequence = list()
    entries_tables = dict()
    observations_selection = list()
    sequence_groups_observations = dict()
    index = 0
    for name_group in names_groups_observations_sequence:
        if (
            (name_group is not None) and
            (len(name_group) > 0) and
            (name_group in groups_observations.keys())
        ):
            # Filter columns and rows in table.
            table_part = (
                porg.filter_select_table_columns_rows_by_identifiers(
                    table=table_source,
                    index_rows=name_index_rows,
                    identifiers_columns=columns_all,
                    identifiers_rows=groups_observations[name_group],
                    report=False,
            ))
            # Collect information.
            groups_sequence.append(name_group)
            entries_tables[name_group] = table_part
            observations_selection.extend(groups_observations[name_group])
            sequence_groups_observations[name_group] = index
            index += 1
            pass
        pass
    # Filter rows in table.
    table_selection = table_source.loc[
        table_source[name_index_rows].isin(
            observations_selection
        ), :
    ].copy(deep=True)

    # Determine and fill groups of observations.
    # These next two functions create a new column named "group" and assign
    # categorical names corresponding to specific groups of observations.
    if (replicate_groups):
        # Each observation has potential to belong to multiple groups.
        table_group = porg.determine_fill_table_groups_rows_with_replicates(
            table=table_selection,
            index_rows=name_index_rows,
            column_group=column_group,
            groups_rows=groups_observations,
            report=False,
        )
    else:
        # Each observation can only belong to a single group.
        table_group = porg.determine_fill_table_groups_rows(
            table=table_selection,
            index_rows=name_index_rows,
            column_group=column_group,
            groups_rows=groups_observations,
            report=False,
        )
        pass

    # Sort rows in table by groups.
    table_group = porg.sort_table_rows_by_single_column_reference(
        table=table_group,
        index_rows=name_index_rows,
        column_reference=column_group,
        column_sort_temporary="sort_temporary",
        reference_sort=sequence_groups_observations,
    )

    # Bundle information.
    pail = dict()
    pail["groups_sequence"] = groups_sequence
    pail["entries_tables"] = entries_tables
    pail["table_group"] = table_group

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("describe_features_between_groups_observations.py")
        print("module: " + name_module)
        name_function = str(
            "collect_entries_tables_features_groups_observations()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print(str(
            "whole table after basic preparation and introduction of groups:"
        ))
        print(table_group)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


# Describe quantitative features in groups of observations.


def describe_quantitative_features_by_groups_observations(
    entries_tables=None,
    groups_sequence=None,
    index_columns=None,
    index_rows=None,
    index_features=None,
    columns_features=None,
    translations_feature=None,
    table_group=None,
    column_group=None,
    ttest_one=None,
    ttest_two=None,
    ttest_three=None,
    report=None,
):
    """
    Describe features on quantitative, continuous, interval or ratio scale of
    measurement in terms of their values within groups of observations. For
    each feature, this function prepares a table of summary, descriptive,
    statistical measurements.

    ----------
    Format of source table
    ----------
    Format of source table is in wide format with features across columns and
    values corresponding to their observations across rows. A special header
    row gives identifiers or names corresponding to each feature across
    columns, and a special column gives identifiers or names corresponding to
    each observation across rows. Another special column provides names of
    categorical groups for these observations. The table has explicitly named
    indices across columns and rows.
    ----------
    features        group     feature_1 feature_2 feature_3 feature_4 feature_5
    observation
    observation_1   group_1   0.001     0.001     0.001     0.001     0.001
    observation_2   group_1   0.001     0.001     0.001     0.001     0.001
    observation_3   group_2   0.001     0.001     0.001     0.001     0.001
    observation_4   group_2   0.001     0.001     0.001     0.001     0.001
    observation_5   group_3   0.001     0.001     0.001     0.001     0.001
    ----------

    ----------
    Format of product table
    ----------
    Format of product table is in partial long format with summary statistics
    and measures across columns and features across rows. A special column
    gives identifiers corresponding to each feature across rows. Another
    special column provides names of categorical groups of observations. For
    versatility, this table does not have explicity defined indices across
    columns or rows.
    ----------
    detail    group   mean standard_error standard_deviation median interqua...
    feature
    feature_1 group_1 0.01 0.001          0.001              0.015  0.5
    feature_1 group_2 0.01 0.001          0.001              0.015  0.5
    feature_1 group_3 0.01 0.001          0.001              0.015  0.5
    feature_1 group_4 0.01 0.001          0.001              0.015  0.5
    feature_2 group_1 0.01 0.001          0.001              0.015  0.5
    feature_2 group_2 0.01 0.001          0.001              0.015  0.5
    feature_2 group_3 0.01 0.001          0.001              0.015  0.5
    feature_2 group_4 0.01 0.001          0.001              0.015  0.5
    ----------

    ----------
    An alternative format for a product table.
    ----------
    A disadvantage of this format is that values within columns would have
    different variable types.
    ----------
    group                group_1   group_2   group_3   group_4   ...
    measure
    mean                 0.001     0.001     0.001     0.001     ...
    standard_error       0.001     0.001     0.001     0.001     ...
    standard_deviation   0.001     0.001     0.001     0.001     ...
    95_confidence_low    0.001     0.001     0.001     0.001     ...
    95_confidence_high   0.001     0.001     0.001     0.001     ...
    minimum              0.001     0.001     0.001     0.001     ...
    maximum              0.001     0.001     0.001     0.001     ...
    median               0.001     0.001     0.001     0.001     ...
    count_total          30        30        30        30        ...
    count_valid          25        25        25        25        ...
    ----------

    Date, review or revision: 22 January 2026
    Date, review or revision: 14 April 2025

    arguments:
        entries_tables (dict<object>): entries with names as keys and as
            values, Pandas data-frame tables with features across columns and
            observations across rows
        groups_sequence (list<str>): unique names of groups of observations
            corresponding to the separate tables in the original sequence
            of instances of parameters
        index_columns (str): name for index corresponding to features across
            columns in the original source table
        index_rows (str): name for index corresponding to observations across
            rows in the original source table
        index_features (str): name for index corresponding to features across
            rows in the novel, product table
        columns_features (list<str>): names of columns in original source
            table for a selection of features on a quantitative, continuous,
            interval or ratio scale of measurement
        translations_feature (dict<str>): translations for names of features
        table_group (object): Pandas data-frame table of features across
            columns and observations across rows with values on a quantitative,
            continuous, interval or ratio scale of measurement
        column_group (str): name of column in original source table that
            designates groups of observations across rows
        ttest_one (dict): collection of parameters for T-test
            name (str): name for a column in the summary table that reports the
                p-value of the T-test
            groups (list<str>): names of groups of observations between which
                to perform a T-test, only use names of two groups
            equal_variances (bool): whether to assume that values in both groups
                have equal variances
            independent_groups (bool): whether the groups are independent
                ('True'), or whether the groups are paired, related, or
                otherwise dependent ('False')
            hypothesis_alternative (str): description of the alternative
                hypothesis, either a 'one-sided' or 'one-tailed' analysis in
                which only 'less' or 'greater' is relevant or a 'two-sided' or
                'two-tailed' analysis in which both lesser and greater are
                relevant
        ttest_two (dict): collection of parameters for T-test
        ttest_three (dict): collection of parameters for T-test
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Copy information.
    table_group = table_group.copy(deep=True)
    columns_features = copy.deepcopy(columns_features)
    translations_feature = copy.deepcopy(translations_feature)
    entries_tables = copy.deepcopy(entries_tables)
    groups_sequence = copy.deepcopy(groups_sequence)
    ttest_one = copy.deepcopy(ttest_one)
    ttest_two = copy.deepcopy(ttest_two)

    ##########
    # Describe features in groups of observations.
    table_description_priority = (
        pdesc.describe_features_from_columns_by_separate_tables_rows(
            entries_tables=entries_tables,
            groups_sequence=groups_sequence,
            index_columns=index_columns,
            index_rows=index_rows,
            index_features=index_features,
            columns_features=columns_features,
            translations_feature=translations_feature,
            key_group=column_group,
            threshold_observations=5,
            digits_round=4,
            ttest_one=ttest_one,
            ttest_two=ttest_two,
            ttest_three=ttest_three,
            report=False,
    ))
    table_description_check = (
        pdesc.describe_features_from_table_columns_by_groups_rows(
            table_group=table_group,
            index_columns=index_columns,
            index_rows=index_rows,
            index_features=index_features,
            column_group=column_group,
            groups_sequence=groups_sequence,
            columns_features=columns_features,
            translations_feature=translations_feature,
            key_group=column_group,
            threshold_observations=5,
            digits_round=4,
            ttest_one=ttest_one,
            ttest_two=ttest_two,
            ttest_three=ttest_three,
            report=False,
    ))

    # Collect information.
    pail = dict()
    pail["table_priority"] = table_description_priority
    pail["table_check"] = table_description_check

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("describe_features_between_groups_observations.py")
        print("module: " + name_module)
        name_function = str(
            "describe_quantitative_features_by_groups_observations()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("Table of descriptive statistics and T-tests")
        print(pail["table_priority"])
        print(pail["table_check"])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def prepare_clean_table_description(
    table=None,
    features_sequence=None,
    names_pvalue=None,
    set_threshold_pvalue=None,
    threshold_low_pvalue=None,
    report=None,
):
    """
    Prepare a clean format of the table with descriptive statistics for
    features across groups of observations.

    ----------
    Format of product table
    ----------
    Format of product table is in partial long format with summary statistics
    and measures across columns and features across rows. A special column
    gives identifiers corresponding to each feature across rows. Another
    special column provides names of categorical groups of observations. For
    versatility, this table does not have explicity defined indices across
    columns or rows.
    ----------
    detail    group   mean standard_error standard_deviation median interqua...
    feature
    feature_1 group_1 0.01 0.001          0.001              0.015  0.5
    feature_1 group_2 0.01 0.001          0.001              0.015  0.5
    feature_1 group_3 0.01 0.001          0.001              0.015  0.5
    feature_1 group_4 0.01 0.001          0.001              0.015  0.5
    feature_2 group_1 0.01 0.001          0.001              0.015  0.5
    feature_2 group_2 0.01 0.001          0.001              0.015  0.5
    feature_2 group_3 0.01 0.001          0.001              0.015  0.5
    feature_2 group_4 0.01 0.001          0.001              0.015  0.5
    ----------

    Date, review or revision: 22 January 2026
    Date, review or revision: 23 May 2025

    arguments:
        table (object): Pandas data-frame table of features and groups of
            observations across rows with statistical measures across columns
            to describe the features within these groups of observations
        features_sequence (list<str>): unique names of features in their proper
            sequence across rows in table
        names_pvalue (list<str>): names of columns in table for p-values for
            which to manage custom representation, such as '< 0.0001'
        set_threshold_pvalue (bool): whether to apply threshold and text
            replacement of p-values
        threshold_low_pvalue (float): minimal value of p-value to represent as
            a number, with a textual representation (i.e. '< 0.0001') of any
            value below threshold
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    ##########
    # Define subordinate functions for internal use.
    def determine_pvalue_below_threshold(pvalue=None, threshold_low=None):
        if (pvalue < threshold_low):
            pvalue_new = str("< " + str(threshold_low))
        else:
            pvalue_new = pvalue
        return pvalue_new

    ##########
    # Copy information.
    table = table.copy(deep=True)
    features_sequence = copy.deepcopy(features_sequence)
    # Restore or reset indices to generic default.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table.columns.rename(
        None,
        inplace=True,
    ) # single-dimensional index
    # Copy names of columns in table in their original sequence.
    columns_sequence_source = copy.deepcopy(table.columns.to_list())
    # Create reference of sequence of features.
    sequence_features = dict(zip(
        features_sequence,
        range(len(features_sequence))
    ))
    # Asign sort sequences to features.
    table["sort_feature"] = table.apply(
        lambda row: sequence_features[row["features"]],
        axis="columns", # apply function to each row
    )
    # Sort rows in table.
    table.sort_values(
        by=["group", "sort_feature",],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Organize indices in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Filter and sort columns in table.
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=columns_sequence_source,
        report=False,
    )

    # Manage representation of p-values.
    if (set_threshold_pvalue):
        for name_pvalue in names_pvalue:
            table[name_pvalue] = table.apply(
                lambda row: determine_pvalue_below_threshold(
                    pvalue=row[name_pvalue],
                    threshold_low=threshold_low_pvalue,
                ),
                axis="columns", # apply function to each row
            )
            pass
        pass

    # Report.
    if report:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("describe_features_between_groups_observations.py")
        print("module: " + name_module)
        name_function = str(
            "prepare_clean_table_description()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("Table of descriptive statistics and T-tests")
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


# Plot charts.


def manage_create_write_plot_box_violin_feature_groups(
    column_feature=None,
    entries_tables=None,
    index_columns=None,
    index_rows=None,
    index_features=None,
    groups_sequence=None,
    translations_feature=None,
    translations_group=None,
    path_directory_parent=None,
    report=None,
):
    """
    Manage creation and write to file of box or violin plots to represent
    visually features on quantitative, continuous scales of measurement
    within and between groups of observations.

    Review: TCW; 1 May 2025

    arguments:
        column_feature (str): name of column in original source table for a
            feature on a quantitative, continuous, interval or ratio scale of
            measurement
        entries_tables (dict<object>): entries with names as keys and as
            values, Pandas data-frame tables with features across columns and
            observations across rows
        index_columns (str): name for index corresponding to features across
            columns in the original source table
        index_rows (str): name for index corresponding to observations across
            rows in the original source table
        index_features (str): name for index corresponding to features across
            rows in the novel, product table
        groups_sequence (list<str>): unique names of groups of observations
            corresponding to the separate tables in the original sequence
            of instances of parameters
        translations_feature (dict<str>): translations for names of features
        translations_group (dict<str>): translations for names of groups
        path_directory_parent ((str): path to parent directory within which to
            write files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Extract values for current feature across groups of observations in
    # separate, stratified tables.
    pail_extract = (
        porg.extract_array_values_from_column_by_separate_tables_rows(
            entries_tables=entries_tables,
            column_feature=column_feature,
            groups_sequence=groups_sequence,
            report=False,
    ))
    groups_values_nonmissing = pail_extract["groups_values_nonmissing"]

    # Determine title.
    if (
        (translations_feature is not None) and
        (column_feature in translations_feature.keys())
    ):
        title = translations_feature[column_feature]
    else:
        title = column_feature
        pass

    # Determine names of groups for visual representation on plot chart.
    if (
        (translations_group is not None) and
        (putly.compare_lists_by_inclusion(
            items_dominant=list(translations_group.keys()),
            items_subordinate=pail_extract["names_groups"],
        ))
    ):
        names_groups = list(map(
            lambda name: translations_group[name],
            copy.deepcopy(pail_extract["names_groups"])
        ))
    else:
        names_groups = copy.deepcopy(pail_extract["names_groups"])
        pass

    ##########
    # Create plot chart.
    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Extract parameters for colors.
    if False:
        colors_names_groups = [
            #"purple_violet",
            #"blue_sky",
            #"green_mint",
            #"yellow_sunshine",
            #"purple_violet_faint",
            #"blue_sky_faint",
            #"green_mint_faint",
            #"yellow_sunshine_faint",
            "firebrick",
            "tomato",
            "orange",
            "gold",
            "mediumorchid",
            "mediumslateblue",
            "dodgerblue",
            "mediumturquoise",
        ]
    else:
        colors_names_groups = None
        pass
    if (
        (colors_names_groups is not None) and
        (len(colors_names_groups) > 0)
    ):
        #colors_groups = list(map(
        #    lambda color_name: copy.deepcopy(colors[color_name]),
        #    colors_names_groups
        #))
        colors_groups = list(map(
            lambda color_name: matplotlib.colors.to_rgba(color_name),
            colors_names_groups
        ))
    else:
        colors_groups = [
            (1.0, 1.0, 0.588, 1.0), # r: 255; g: 255; b: 150; alpha: 1.0
            (1.0, 0.784, 0.392, 1.0), # r: 255; g: 200; b: 100; alpha: 1.0
            (0.588, 1.0, 0.588, 1.0), # r: 150; g: 255; b: 150; alpha: 1.0
            (0.392, 0.784, 0.392, 1.0), # r: 100; g: 200; b: 100; alpha: 1.0
            (0.588, 1.0, 1.0, 1.0), # r: 150; g: 255; b: 255; alpha: 1.0
            (0.392, 0.588, 1.0, 1.0), # r: 100; g: 200; b: 255; alpha: 1.0
        ]
        pass

    # Create figure.
    figure_violin = pplot.plot_boxes_groups(
        values_groups=pail_extract["values_nonmissing_groups"],
        names_groups=names_groups,
        title_ordinate=title,
        title_abscissa="",
        title_chart_top_center="",
        colors_groups=colors_groups,
        label_top_center="",
        label_top_left="",
        label_top_right="",
        aspect="landscape",
        orientation_box="vertical",
        axis_linear_minimum=0.0,
        fonts=fonts,
        colors=colors,
        report=report,
    )

    ##########
    # Collect information.
    pail_write_plot_box = dict()
    #pail_write_plot_box[column_feature] = figure_box
    pail_write_plot_violin = dict()
    pail_write_plot_violin[column_feature] = figure_violin

    ##########
    # Write product information to file.

    # Write figures to file.
    #pplot.write_product_plots_parent_directory(
    #    pail_write=pail_write_plot_box,
    #    format="svg", # jpg, png, svg
    #    resolution=300,
    #    path_directory=path_directory_box,
    #)
    pplot.write_product_plots_parent_directory(
        pail_write=pail_write_plot_violin,
        format="jpg", # jpg, png, svg
        resolution=150,
        path_directory=path_directory_parent,
    )

    # Return information.
    pass


def manage_create_write_plots_box_violin_features_groups(
    entries_tables=None,
    index_columns=None,
    index_rows=None,
    index_features=None,
    columns_features=None,
    groups_sequence=None,
    translations_feature=None,
    translations_group=None,
    path_directory_parent=None,
    report=None,
):
    """
    Manage creation and write to file of box or violin plots to represent
    visually features on quantitative, continuous scales of measurement
    within and between groups of observations.

    Review: TCW; 1 May 2025

    arguments:
        entries_tables (dict<object>): entries with names as keys and as
            values, Pandas data-frame tables with features across columns and
            observations across rows
        index_columns (str): name for index corresponding to features across
            columns in the original source table
        index_rows (str): name for index corresponding to observations across
            rows in the original source table
        index_features (str): name for index corresponding to features across
            rows in the novel, product table
        columns_features (list<str>): names of columns in original source
            table for a selection of features on a quantitative, continuous,
            interval or ratio scale of measurement
        groups_sequence (list<str>): unique names of groups of observations
            corresponding to the separate tables in the original sequence
            of instances of parameters
        translations_feature (dict<str>): translations for names of features
        translations_group (dict<str>): translations for names of groups
        path_directory_parent ((str): path to parent directory within which to
            write files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Copy information.
    entries_tables = copy.deepcopy(entries_tables)
    columns_features = copy.deepcopy(columns_features)
    groups_sequence = copy.deepcopy(groups_sequence)
    translations_feature = copy.deepcopy(translations_feature)
    translations_group = copy.deepcopy(translations_group)

    # Iterate on features, apply operations, and collect information.
    for column_feature in columns_features:
        # Manage creation and write to file of plot charts for current feature.
        manage_create_write_plot_box_violin_feature_groups(
            column_feature=column_feature,
            entries_tables=entries_tables,
            index_columns=index_columns,
            index_rows=index_rows,
            index_features=index_features,
            groups_sequence=groups_sequence,
            translations_feature=translations_feature,
            translations_group=translations_group,
            path_directory_parent=path_directory_parent,
            report=report,
        )
        pass
    # Return information.
    pass





################################################################################
# Procedure


##########
# Call main procedure.


def execute_procedure(
    path_directory_dock=None,
    path_directory_dock_pail=None,
    path_directory_source=None,
    path_directory_product=None,
    path_file_source_table_features_observations=None,
    path_file_source_table_features=None,
    path_file_source_list_features=None,
    path_file_source_table_groups_observations=None,
    column_identifier_observation=None,
    column_name_observation=None,
    column_identifier_feature=None,
    column_name_feature=None,
    proportion_nonmissing_observations=None,
    allow_replicate_observations=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Date, review or revision: 22 January 2026

    arguments:
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        path_directory_dock_pail (str): path to pail directory for procedure's
            source and product directories and files
        path_directory_source (str): path to directory for procedure's source
            directories and files
        path_directory_product (str): path to directory for procedure's product
            directories and files
        path_file_source_table_features_observations (str): path to source file
        path_file_source_table_features (str): path to source file
        path_file_source_list_features (str): path to source file
        path_file_source_table_groups_observations (str): path to source file
        column_identifier_observation (str): name of column in source table
        column_name_observation (str): name of column in source table
        column_identifier_feature (str): name of column in source table
        column_name_feature (str): name of column in source table
        proportion_nonmissing_observations (float): threshold by proportion of
            observations that must have nonmissing values for each feature
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
        path_directory_dock=path_directory_dock,
        path_directory_dock_pail=path_directory_dock_pail,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        path_file_source_table_features=path_file_source_table_features,
        path_file_source_list_features=(
            path_file_source_list_features
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        column_identifier_feature=column_identifier_feature,
        column_name_feature=column_name_feature,
        proportion_nonmissing_observations=proportion_nonmissing_observations,
        allow_replicate_observations=allow_replicate_observations,
        report=report,
    )

    ##########
    # Report.
    if pail_parameters["report"]:
        # Organize information.
        # Print information.
        putly.print_terminal_partition(level=3)
        name_system = str("local")
        print("system: " + name_system)
        name_package = str("partner")
        print("package: " + name_package)
        name_module = str("describe_features_between_groups_observations.py")
        print("module: " + name_module)
        name_function = str("execute_procedure()")
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    pail_source = read_source(
        path_directory_dock=pail_parameters["path_directory_dock"],
        path_directory_dock_pail=pail_parameters["path_directory_dock_pail"],
        path_directory_source=pail_parameters["path_directory_source"],
        path_directory_product=pail_parameters["path_directory_product"],
        path_file_source_table_features_observations=(
            pail_parameters["path_file_source_table_features_observations"]
        ),
        path_file_source_table_features=(
            pail_parameters["path_file_source_table_features"]
        ),
        path_file_source_list_features=(
            pail_parameters["path_file_source_list_features"]
        ),
        path_file_source_table_groups_observations=(
            pail_parameters["path_file_source_table_groups_observations"]
        ),
        report=pail_parameters["report"],
    )

    ##########
    # Organize parameters.
    pail_features_observations = organize_parameters_features_observations(
        table_features_observations=pail_source["table_features_observations"],
        table_features=pail_source["table_features"],
        list_features=pail_source["list_features"],
        table_groups_observations=pail_source["table_groups_observations"],
        column_identifier_feature=pail_parameters["column_identifier_feature"],
        column_name_feature=pail_parameters["column_name_feature"],
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_name_observation=pail_parameters["column_name_observation"],
        report=pail_parameters["report"],
    )
    #pail_features_observations["pail_features"]
    #pail_features_observations["pail_groups"]
    #pail_features_observations["pail_observations"]

    ##########
    # Filter columns for features and rows for observations in main table.
    pail_alias = pail_features_observations
    features_selection = copy.deepcopy(
        pail_alias["pail_features"]["features_selection"]
    )
    observations_selection = copy.deepcopy(
        pail_alias["pail_observations"]["observations_selection"]
    )
    table_main = sutly.filter_table_columns_features_rows_observations(
        table=pail_source["table_features_observations"],
        column_identifier_observation=(
            pail_parameters["column_identifier_observation"]
        ),
        column_name_observation=pail_parameters["column_name_observation"],
        column_identifier_groups_observations=(
            pail_parameters["column_identifier_observation"]
        ),
        columns_categories=pail_alias["pail_groups"]["categories_groups"],
        features_selection=features_selection,
        observations_selection=observations_selection,
        filter_table_main=True,
        report=pail_parameters["report"],
    )

    ##########
    # Define groups and stratify observations of relevant features.
    pail_alias = pail_features_observations
    names_groups_sequence = (
        pail_alias["pail_observations"]["names_groups_observations_sequence"]
    )
    pail_parts = collect_entries_tables_features_groups_observations(
        table=table_main,
        name_index_columns="features",
        name_index_rows=pail_parameters["column_identifier_observation"],
        names_groups_observations_sequence=names_groups_sequence, # or None to select all with 'execution' of 1
        groups_observations=(
            pail_alias["pail_observations"]["groups_observations"]
        ),
        column_group="group",
        replicate_groups=True, # caution; allows an observation to belong to multiple groups
        report=pail_parameters["report"],
    )
    #pail_parts["groups_sequence"]
    #pail_parts["entries_tables"]
    #pail_parts["table_group"]
    #print(pail_parts["table_group"])
    #print(list(pail_parts["entries_tables"].keys()))
    #print(pail_parts["entries_tables"]["younger"])

    ##########
    # 4. Describe quantitative features in groups of observations.
    pail_alias = pail_features_observations
    features_selection = copy.deepcopy(
        pail_alias["pail_features"]["features_selection"]
    )
    translations_features = (
        pail_alias["pail_features"]["translations_features"]
    )
    pail_description = describe_quantitative_features_by_groups_observations(
        entries_tables=pail_parts["entries_tables"],
        groups_sequence=pail_parts["groups_sequence"],
        index_columns="features",
        index_rows=pail_parameters["column_identifier_observation"],
        index_features="features",
        columns_features=features_selection,
        translations_feature=translations_features,
        table_group=pail_parts["table_group"],
        column_group="group",
        ttest_one=None,
        ttest_two=None,
        ttest_three=None,
        #ttest_one={
        #    "name": "p_ttest_placebo",
        #    "groups": [
        #        "male-placebo-1", "male-placebo-2",
        #    ],
        #    "equal_variances": True,
        #    "independent_groups": True,
        #    "hypothesis_alternative": "two-sided",
        #}, # or None
        #ttest_two={
        #    "name": "p_ttest_omega3",
        #    "groups": [
        #        "male-omega3-1", "male-omega3-2",
        #    ],
        #    "equal_variances": True,
        #    "independent_groups": True,
        #    "hypothesis_alternative": "two-sided",
        #}, # or None
        #ttest_one={
        #    "name": "p_ttest_age",
        #    "groups": [
        #        "adipose_age_-_younger", "adipose_age_-_older",
        #    ],
        #    "equal_variances": True,
        #    "independent_groups": True,
        #    "hypothesis_alternative": "two-sided",
        #}, # or None
        #ttest_two={
        #    "name": "p_ttest_sex",
        #    "groups": [
        #        "adipose_age_-_older_female", "adipose_age_-_older_male",
        #    ],
        #    "equal_variances": True,
        #    "independent_groups": True,
        #    "hypothesis_alternative": "two-sided",
        #}, # or None
        #ttest_three={
        #    "name": "p_ttest_omega3",
        #    "groups": [
        #        "adipose_placebo_-_older_placebo_after",
        #        "adipose_omega3_-_older_omega3_after",
        #    ],
        #    "equal_variances": True,
        #    "independent_groups": True,
        #    "hypothesis_alternative": "two-sided",
        #}, # or None
        report=pail_parameters["report"],
    )


    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    #print(table)
    #print(table.columns.to_list())
    print(pail_description["table_priority"])
    print(pail_description["table_priority"].columns.to_list())


    ##########
    # 5. Prepare clean version of the description table for report.
    table_description_clean = prepare_clean_table_description(
        table=pail_description["table_priority"],
        features_sequence=features_selection,
        names_pvalue=None,
        #names_pvalue=[
        #    "p_ttest_age",
        #    "p_ttest_sex",
        #    "p_ttest_omega3",
        #],
        set_threshold_pvalue=False,
        threshold_low_pvalue=0.0001,
        report=pail_parameters["report"],
    )

    ##########
    # Write information to file.
    # Bundle information.
    # Bundles of information for files.
    # Texts.
    #pail_write_texts = dict()
    # Objects.
    #pail_write_objects
    # Lists.
    # Tables.
    pail_write_tables = dict()
    pail_write_tables[str("table_description_clean")] = table_description_clean

    ##########
    # Write product information to file.

    # Define paths to directories.
    path_directory_tables = os.path.join(
        path_directory_product, "tables",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_tables,
    )
    # Text.
    # Objects.
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

    ##########
    # Create charts and write to file.

    # Define paths to directories.
    path_directory_charts_violin = os.path.join(
        path_directory_product, "charts", "box",
    )
    # Create directories.
    putly.create_directories(
        path=path_directory_charts_violin,
    )

    # Create plot charts of type box.
    manage_create_write_plots_box_violin_features_groups(
        entries_tables=pail_parts["entries_tables"],
        index_columns="features",
        index_rows=pail_parameters["column_identifier_observation"],
        index_features="features",
        columns_features=features_selection,
        groups_sequence=pail_parts["groups_sequence"],
        translations_feature=translations_features,
        translations_group=None,
        path_directory_parent=path_directory_charts_violin,
        report=False,
    )


    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_directory_dock = sys.argv[1]
    path_directory_dock_pail = sys.argv[2]
    path_directory_source = sys.argv[3]
    path_directory_product = sys.argv[4]
    path_file_source_table_features_observations = sys.argv[5]
    path_file_source_table_features = sys.argv[6]
    path_file_source_list_features = sys.argv[7]
    path_file_source_table_groups_observations = sys.argv[8]
    column_identifier_observation = sys.argv[9]
    column_name_observation = sys.argv[10]
    column_identifier_feature = sys.argv[11]
    column_name_feature = sys.argv[12]
    proportion_nonmissing_observations = sys.argv[13]
    allow_replicate_observations = sys.argv[14]
    report = sys.argv[15]

    # Call function for procedure.
    execute_procedure(
        path_directory_dock=path_directory_dock,
        path_directory_dock_pail=path_directory_dock_pail,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_file_source_table_features_observations=(
            path_file_source_table_features_observations
        ),
        path_file_source_table_features=path_file_source_table_features,
        path_file_source_list_features=(
            path_file_source_list_features
        ),
        path_file_source_table_groups_observations=(
            path_file_source_table_groups_observations
        ),
        column_identifier_observation=column_identifier_observation,
        column_name_observation=column_name_observation,
        column_identifier_feature=column_identifier_feature,
        column_name_feature=column_name_feature,
        proportion_nonmissing_observations=proportion_nonmissing_observations,
        allow_replicate_observations=allow_replicate_observations,
        report=report,
    )

    pass



#
