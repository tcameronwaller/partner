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
# Date, first execution: 28 March 2025
# Date, last execution: 28 March 2025
# Review: TCW; 28 March 2025
################################################################################
# Note

# The specialty of this Python script is to drive multiple regressions from
# parameters within a single table. This script calls versatile functionality
# from the "regression.py" module within the "partner" Python package.

# Useful functionality for preparing the table of data for regression.
# pandas.get_dummies(groups).values

################################################################################
# Installation and importation

# Standard
import sys
import os
import copy

# Relevant
import pandas
import scipy
import numpy

# Custom
import partner.utility as putly
import partner.parallelization as prall
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc
import partner.regression as preg

#dir()
#importlib.reload()

###############################################################################
# Functionality


##########
# Read source information from file.


def define_column_types_table_parameters():
    """
    Defines the types of variables for columns in table of parameters.

    Review: TCW; 31 March 2025

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify types of variables in columns of table.
    types_columns = dict()
    types_columns["execution"] = "string" # "int32"
    types_columns["sequence"] = "string" # "int32"
    types_columns["group"] = "string"
    types_columns["instance"] = "string"
    types_columns["selection_observations"] = "string"
    types_columns["type_regression"] = "string"
    types_columns["formula_text"] = "string"
    types_columns["feature_response"] = "string"
    types_columns["features_predictor_fixed"] = "string"
    types_columns["features_predictor_random"] = "string"
    types_columns["features_continuity_scale"] = "string"
    types_columns["identifier_observations"] = "string"
    types_columns["method_scale"] = "string"
    types_columns["data_path_directory"] = "string"
    types_columns["data_file"] = "string"
    types_columns["review"] = "string"
    types_columns["note"] = "string"
    # Return information.
    return types_columns


def read_source_table_parameters(
    path_file_table_parameters=None,
    path_directory_dock=None,
    report=None,
):
    """
    Read and organize source information about parameters for regressions.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 31 March 2025

    arguments:
        path_file_table_parameters (str): path to source file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters for selection
            of sets from MSigDB

    """

    # Read information from file.

    # Table of parameters for parallel instances.
    types_columns = define_column_types_table_parameters()
    table = pandas.read_csv(
        path_file_table_parameters,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Organize information.
    table["execution"] = pandas.to_numeric(
        table["execution"],
        downcast="integer",
        errors="coerce",
    )
    table["sequence"] = pandas.to_numeric(
        table["sequence"],
        downcast="integer",
        errors="coerce",
    )

    # Collect information.
    records = list()
    for index, row in table.iterrows():
        if (int(row["execution"]) == 1):
            # Collect information and parameters from current row in table.
            record = dict()
            record["execution"] = row["execution"]
            record["sequence"] = row["sequence"]
            record["group"] = str(row["group"]).strip()
            record["instance"] = str(row["instance"]).strip()
            record["name_instance"] = "_".join([
                str(row["group"]).strip(),
                str(row["sequence"]).strip(),
                str(row["instance"]).strip(),
            ])
            record["selection_observations"] = (
                putly.parse_extract_text_keys_values_semicolon_colon_comma(
                    text=row["selection_observations"],
                )
            )["features_values"]
            record["type_regression"] = str(row["type_regression"]).strip()
            record["formula_text"] = str(row["formula_text"]).strip()
            record["feature_response"] = str(row["feature_response"]).strip()
            record["features_predictor_fixed"] = putly.parse_text_list_values(
                text=row["features_predictor_fixed"],
                delimiter=",",
            )
            record["features_predictor_random"] = putly.parse_text_list_values(
                text=row["features_predictor_random"],
                delimiter=",",
            )
            record["features_continuity_scale"] = putly.parse_text_list_values(
                text=row["features_continuity_scale"],
                delimiter=",",
            )
            if (
                (row["identifier_observations"] is not None) and
                (len(str(row["identifier_observations"])) > 0) and
                (str(row["identifier_observations"]).strip().lower() != "none")
            ):
                record["identifier_observations"] = str(
                    row["identifier_observations"]
                ).strip()
            else:
                record["identifier_observations"] = None
                pass
            record["method_scale"] = str(row["method_scale"]).strip()
            record["data_path_directory"] = putly.parse_text_list_values(
                text=row["data_path_directory"],
                delimiter=",",
            )
            record["data_file"] = str(row["data_file"]).strip()
            record["review"] = str(row["review"]).strip()
            record["note"] = str(row["note"]).strip()

            # Collect unique names of columns relevant to instance of
            # parameters from current row in table.
            features_regression = list()
            features_regression.append(record["feature_response"])
            features_regression.extend(record["features_predictor_fixed"])
            features_regression.extend(record["features_predictor_random"])
            if (record["identifier_observations"] is not None):
                features_regression.insert(
                    0,
                    record["identifier_observations"],
                )
                pass
            features_regression = putly.collect_unique_elements(
                elements=features_regression,
            )
            record["features_regression"] = copy.deepcopy(features_regression)
            features_relevant = list()
            dictionaries = [
                "selection_observations",
            ]
            for dictionary in dictionaries:
                if record[dictionary] is not None:
                    features_relevant.extend(list(record[dictionary].keys()))
                    pass
                pass
            features_relevant.extend(record["features_continuity_scale"])
            features_relevant.extend(features_regression)
            features_relevant = putly.collect_unique_elements(
                elements=features_relevant,
            )
            record["features_relevant"] = copy.deepcopy(features_relevant)
            # Collect information and parameters for current row in table.
            records.append(record)
            pass
        pass

    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["records"] = records
    # Report.
    if report:
        # Organize.
        count_records = len(records)
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: read_source_table_parameters()")
        putly.print_terminal_partition(level=5)
        print("parameter table:")
        print(table)
        putly.print_terminal_partition(level=5)
        print("count of records or instances: " + str(count_records))
        putly.print_terminal_partition(level=5)
        print("instance[0]:")
        print(records[0])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


def read_source_table_data(
    name_file_table_data=None,
    directories_path_data=None,
    path_directory_dock=None,
    report=None,
):
    """
    Read and organize source information about data for regression.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 1 April 2025

    arguments:
        name_file_table_data (str): name of file for the table of data with
            features and observations for regression
        directories_path_data (list<str>): names of directories in path at
            which to find the file for the table of data with features and
            observations for regression
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters for selection
            of sets from MSigDB

    """

    # Define paths to directories and files.
    if ("dock" in directories_path_data):
        directories_path_data = list(filter(
            lambda directory: (directory != "dock"),
            directories_path_data
        ))
        pass
    path_directory_data = os.path.join(
        path_directory_dock,
        *directories_path_data, # 'splat' operator unpacks list items
    )
    path_file_table = os.path.join(
        path_directory_data, name_file_table_data,
    )

    # Read information from file.

    # Table of data with values for observations of features.
    table = pandas.read_csv(
        path_file_table,
        sep="\t",
        header=0,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )

    # Report.
    if report:
        # Organize.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: read_source_table_data()")
        putly.print_terminal_partition(level=5)
        print("data table:")
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


##########
# Organize information in table for features and observations.


def organize_table_data(
    table=None,
    selection_observations=None,
    features_relevant=None,
    features_regression=None,
    features_continuity_scale=None,
    index_columns_source=None,
    index_columns_product=None,
    index_rows_source=None,
    index_rows_product=None,
    adjust_scale=None,
    method_scale=None,
    explicate_indices=None,
    report=None,
):
    """
    Organize table of data with features and observations for regression.

    1. filter columns in table for relevant features
    2. filter rows in table for relevant observations
    3. standardize scale of features with measurement values on quantitative,
       continuous scale of measurement, ratio or interval

    ----------
    Format of source data table (name: "table")
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
    Format of product data table (name: "table")
    ----------
    Format of source data table is in wide format with features across columns
    and values corresponding to their observations across rows. A special
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

    Review: TCW; 2 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        selection_observations (dict<list<str>>): names of columns in data
            table for feature variables and their categorical values by
            which to filter rows for observations in data table
        features_relevant (list<str>): names of columns in data table for
            feature variables that are relevant to the current instance of
            parameters
        features_regression (list<str>): names of columns in data table for
            feature variables that are relevant to the actual regression model
        features_continuity_scale (list<str>): names of columns in data table
            for feature variables with values on quantitative, continuous scale
            of measurement, interval or ratio, for which to standardize the
            scale by z score
        index_columns_source (str): name of single-level index across columns
            in table
        index_columns_product (str): name of single-level index across columns
            in table
        index_rows_source (str): name of single-level index across rows in
            table
        index_rows_product (str): name of single-level index across rows in
            table
        adjust_scale (bool): whether to adjust or standardize the scale of
            values for features across observations
        method_scale (str): name of method to use to adjust the scale of values
            for features across observations, either 'z_score' or 'unit_range'
        explicate_indices (bool): whether to explicate, define, or specify
            explicit indices across columns and rows in table
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table = table.copy(deep=True)
    selection_observations = copy.deepcopy(selection_observations)
    features_relevant = copy.deepcopy(features_relevant)
    features_regression = copy.deepcopy(features_regression)
    features_continuity_scale = copy.deepcopy(features_continuity_scale)

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

    # Filter columns in table for features.
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=features_relevant,
        report=report,
    )
    # Filter rows in table for observations.
    table = porg.filter_select_table_rows_by_columns_categories(
        table=table,
        columns_categories=selection_observations,
        report=report,
    )

    # Remove rows from table for observations with any missing values for
    # relevant features that are part of the model for regression.
    table.dropna(
        axis="index",
        how="any",
        subset=features_regression,
        inplace=True,
    )

    # Standardize scale of values for observations of features.
    if (
        (adjust_scale) and
        (str(method_scale).strip().lower() == "z_score")
    ):
        table = pscl.transform_standard_z_score_by_table_columns(
                table=table,
                columns=features_continuity_scale,
                report=report,
        )
    elif (
        (adjust_scale) and
        (str(method_scale).strip().lower() == "unit_range")
    ):
        table = pscl.transform_unit_range_by_table_columns(
                table=table,
                columns=features_continuity_scale,
                report=report,
        )
        pass

    # Organize indices in table.
    # Determine whether parameters specified a column in the table for
    # identifiers of observations that will become index across rows.
    if (
        (index_rows_source is not None) and
        (index_rows_source in (table.columns.tolist()))
    ):
        # Create index across rows from column that already exists in table.
        table = porg.explicate_table_indices_columns_rows_single_level(
            table=table,
            index_columns=index_columns_source,
            index_rows=index_rows_source,
            explicate_indices=explicate_indices,
            report=False,
        )
    else:
        # Name generic index in table.
        table.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table.index.set_names(index_rows_source, inplace=True)
        pass
    # Standardize names of indices.
    table = porg.translate_names_table_indices_columns_rows(
        table=table,
        index_columns_product=index_columns_product,
        index_rows_source=index_rows_source,
        index_rows_product=index_rows_product,
        report=None,
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        name_function = str(
            "organize_table_data()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


# Evaluate information in table for features and observations.


def evaluate_table_data(
    table=None,
    selection_observations=None,
    features_relevant=None,
    features_regression=None,
    features_continuity_scale=None,
    index_columns_source=None,
    index_columns_product=None,
    index_rows_source=None,
    index_rows_product=None,
    adjust_scale=None,
    method_scale=None,
    explicate_indices=None,
    report=None,
):
    """
    Evaluate table of data with features and observations for regression.

    This function is useful to evaluate the variance of values for features
    across observations after adjusting their scale to the common unit range.
    This evaluation enables a manual, supervised approach to a method of
    feature selection that is called 'variance thresholding'. The rationale is
    that features with very little variance are unlikely to contribute to the
    regression model, at least not with a favorable ratio of signal to noise.

    From basic testing, it seems that a random number with a uniform
    distribution between zero (0) and one (1) has a variance of about 10% or
    0.1. After adjusting the scales of features to have unit range, a
    reasonable threshold on the variance might be approximately 1% or 0.01.

    Review: TCW; 2 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        selection_observations (dict<list<str>>): names of columns in data
            table for feature variables and their categorical values by
            which to filter rows for observations in data table
        features_relevant (list<str>): names of columns in data table for
            feature variables that are relevant to the current instance of
            parameters
        features_regression (list<str>): names of columns in data table for
            feature variables that are relevant to the actual regression model
        features_continuity_scale (list<str>): names of columns in data table
            for feature variables with values on quantitative, continuous scale
            of measurement, interval or ratio, for which to standardize the
            scale by z score
        index_columns_source (str): name of single-level index across columns
            in table
        index_columns_product (str): name of single-level index across columns
            in table
        index_rows_source (str): name of single-level index across rows in
            table
        index_rows_product (str): name of single-level index across rows in
            table
        adjust_scale (bool): whether to adjust or standardize the scale of
            values for features across observations
        method_scale (str): name of method to use to adjust the scale of values
            for features across observations, either 'z_score' or 'unit_range'
        explicate_indices (bool): whether to explicate, define, or specify
            explicit indices across columns and rows in table
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Organize information in table.
    # For all features in the regression model (including the response), adjust
    # the scale of values to have common unit range.
    features_scale = copy.deepcopy(features_regression)
    features_category = copy.deepcopy(list(selection_observations.keys()))
    features_category.insert(0, index_rows_source)
    features_scale = list(filter(
        lambda feature: (feature not in features_category),
        features_scale
    ))
    table = organize_table_data(
        table=table,
        selection_observations=selection_observations,
        features_relevant=features_relevant,
        features_regression=features_regression,
        features_continuity_scale=features_scale,
        index_columns_source=index_columns_source,
        index_columns_product=index_columns_product,
        index_rows_source=index_rows_source,
        index_rows_product=index_rows_product,
        adjust_scale=True,
        method_scale="unit_range", # "z_score" or "unit_range"
        explicate_indices=True,
        report=False,
    )
    # Filter columns in table for features.
    table = porg.filter_sort_table_columns(
        table=table,
        columns_sequence=features_scale,
        report=report,
    )

    # Calculate descriptive statistical measures for features.
    table_mean = table.aggregate(
        lambda series: series.mean(),
        axis="index", # apply function to each column
    )
    table_variance = table.aggregate(
        lambda series: series.var(),
        axis="index", # apply function to each column
    )
    table_deviation = table.aggregate(
        lambda series: series.std(),
        axis="index", # apply function to each column
    )
    table_minimum = table.aggregate(
        lambda series: series.min(),
        axis="index", # apply function to each column
    )
    table_maximum = table.aggregate(
        lambda series: series.max(),
        axis="index", # apply function to each column
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        name_function = str(
            "evaluate_table_data()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=4)
        print("table before aggregation:")
        putly.print_terminal_partition(level=5)
        print(table)
        putly.print_terminal_partition(level=4)
        print("series after aggregation:")
        putly.print_terminal_partition(level=5)
        print("mean")
        putly.print_terminal_partition(level=6)
        print(table_mean.iloc[0:10])
        putly.print_terminal_partition(level=5)
        print("variance")
        putly.print_terminal_partition(level=6)
        print(table_variance.iloc[0:10])
        putly.print_terminal_partition(level=5)
        print("standard deviation")
        putly.print_terminal_partition(level=6)
        print(table_deviation.iloc[0:10])
        putly.print_terminal_partition(level=5)
        print("minimum")
        putly.print_terminal_partition(level=6)
        print(table_minimum.iloc[0:10])
        putly.print_terminal_partition(level=5)
        print("maximum")
        putly.print_terminal_partition(level=6)
        print(table_maximum.iloc[0:10])
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


# Check parameters and table of data for performing regression analysis.


def check_parameters_table_data_regression(
    table=None,
    name_instance=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    threshold_features_variance=None,
    threshold_observations_count=None,
    measure_variance=None,
    report=None,
):
    """
    Check parameters for regression, including the table of data with features
    and observations.

    For the check on the variance of features across observations, it is most
    accurate and informative to compare the variances after transforming the
    values of features to have unit range between zero (0) and one (1). Do not
    compare variances after transforming the values of features to standard
    z scores, as this transformation forces values of all features to have the
    same variance, unit standard deviation.

    Review: TCW; 2 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        name_instance (str): compound name for instance of parameters
        type_regression (str): name of type of regression model to use
        formula_text (str): human readable formula for regression model,
            treated as a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        features_predictor_random (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with random effects
        threshold_features_variance (float): threshold minimal variance for
            values of features across observations
        threshold_observations_count (int): threshold minimal count of
            observations that must have nonmissing values of all features
        measure_variance (str): name of measure to use for variance, either
            'variance', 'standard_deviation', or 'coefficient_variation'
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information

    """

    # Copy information.
    table = table.copy(deep=True)
    features_predictor_fixed = copy.deepcopy(features_predictor_fixed)
    features_predictor_random = copy.deepcopy(features_predictor_random)

    # Check that functionality currently supports the type of regression.
    types = [
        "linear",
        "logistic",
        #"linear_r",
        #"logistic_r",
    ]
    if (
        (type_regression is not None) and
        (len(str(type_regression)) > 0) and
        (str(type_regression).strip().lower() != "none") and
        (str(type_regression) in types)
    ):
        check_type = True
    else:
        check_type = False
        pass

    # Check for features in the regression model that features exist in the
    # table of data.
    columns_available = copy.deepcopy(table.columns.to_list())
    check_fixed = putly.compare_lists_by_inclusion(
        items_dominant=columns_available,
        items_subordinate=features_predictor_fixed,
    )
    check_random = putly.compare_lists_by_inclusion(
        items_dominant=columns_available,
        items_subordinate=features_predictor_random,
    )
    if (
        (len(columns_available) > 1) and
        (feature_response in columns_available) and
        (check_fixed) and
        (check_random)
    ):
        check_features_exist = True
    else:
        check_features_exist = False
        pass

    # Check for features in the regression model that adequate observations
    # with nonmissing values exist in the table of data.
    def check_features_nonmissing_observations(
        table=None,
        features=None,
        threshold_observations_count=None,
    ):
        """
        Define subordinate function for internal use within dominant function.
        """
        for feature in features:
            count = table[feature].dropna().count()
            if (count >= threshold_observations_count):
                return True
            else:
                return False
            pass
        pass
    # Copy information.
    table_count = table.copy(deep=True)
    # Evaluate features separately.
    #check_fixed = check_features_nonmissing_observations(
    #    table=table_count,
    #    features=features_predictor_fixed,
    #    threshold_observations_count=threshold_observations_count,
    #)
    # Evaluate features together.
    features_regression = list()
    features_regression.append(feature_response)
    features_regression.extend(features_predictor_fixed)
    features_regression.extend(features_predictor_random)
    table_count.dropna(
        axis="index",
        how="any",
        subset=features_regression,
        inplace=True,
    )
    count_rows = int(table_count.shape[0])
    if (count_rows >= threshold_observations_count):
        check_observations_count = True
    else:
        check_observations_count = False
        pass

    # Check for features in the regression model that values across
    # observations of each feature have adequate variance.
    def check_features_variance(
        table=None,
        features=None,
        measure_variance=None,
        threshold_features_variance=None,
    ):
        """
        Define subordinate function for internal use within dominant function.
        """
        checks = list()
        for feature in features:
            values_raw = table[feature].to_numpy(
                dtype="float64",
                na_value=numpy.nan,
                copy=True,
            )
            pail = pdesc.calculate_variance_measures(
                array=values_raw,
            )
            if (pail[measure_variance] >= threshold_features_variance):
                checks.append(True)
            else:
                checks.append(False)
                pass
            pass
        return all(checks)
    # Evaluate individual features separately.
    check_fixed = check_features_variance(
        table=table,
        features=features_predictor_fixed,
        measure_variance=measure_variance,
        threshold_features_variance=threshold_features_variance,
    )
    check_random = check_features_variance(
        table=table,
        features=features_predictor_random,
        measure_variance=measure_variance,
        threshold_features_variance=threshold_features_variance,
    )
    if (check_fixed and check_random):
        check_variance = True
    else:
        check_variance = False
        pass

    # Check whether there was failure of any checks overall.
    if (
        (check_type) and
        (check_features_exist) and
        (check_observations_count) and
        (check_variance)
    ):
        check_overall = True
    else:
        check_overall = False
        pass

    # Collect information.
    pail = dict()
    pail["check_type"] = check_type
    pail["check_features_exist"] = check_features_exist
    pail["check_observations_count"] = check_observations_count
    pail["check_variance"] = check_variance
    pail["check_overall"] = check_overall

    # Warn.
    if (not check_overall):
        # Print.
        putly.print_terminal_partition(level=3)
        putly.print_terminal_warning()
        putly.print_terminal_partition(level=5)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        name_function = str(
            "check_parameters_table_data_regression()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print(
            "Parameters and or data for a regression analysis failed checks!"
        )
        print("name_instance: " + name_instance)
        putly.print_terminal_partition(level=5)

        pass

    # Report.
    if report:
        # Organize.
        #count_records = len(records)
        # Print.
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        name_function = str(
            "check_parameters_table_data_regression()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("check_type: " + str(check_type))
        print("check_features_exist: " + str(check_features_exist))
        print("check_observations_count: " + str(check_observations_count))
        print("check_variance: " + str(check_variance))
        print("check_overall: " + str(check_overall))
        pass
    # Return information.
    return pail




################################################################################
# Procedure



##########
# Control procedure within branch for parallelization.


def control_procedure_part_branch(
    execution=None,
    sequence=None,
    group=None,
    instance=None,
    name_instance=None,
    selection_observations=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    features_regression=None,
    features_continuity_scale=None,
    features_relevant=None,
    identifier_observations=None,
    method_scale=None,
    data_path_directory=None,
    data_file=None,
    review=None,
    note=None,
    path_file_table_parameters=None,
    path_directory_dock=None,
    report=None,
):
    """
    Execute module's main behavior for a single parallel branch.

    arguments:
        execution (int): logical binary indicator of whether to execute and
            handle the parameters for the current instance
        sequence (int): sequential index for instance's name and sort order
        group (str): categorical group of instances
        instance (str): name or designator of instance
        name_instance (str): compound name for instance of parameters
        selection_observations (dict<list<str>>): names of columns in data
            table for feature variables and their categorical values by
            which to filter rows for observations in data table
        features_continuity_scale (list<str>): names of columns in data table
            for feature variables with values on quantitative, continuous scale
            of measurement, interval or ratio, for which to standardize the
            scale by z score
        type_regression (str): name of type of regression model to use
        formula_text (str): human readable formula for regression model,
            treated as a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        features_predictor_random (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with random effects
        features_regression (list<str>): names of columns in data table for
            feature variables that are relevant to the actual regression model
        features_relevant (list<str>): names of columns in data table for
            feature variables that are relevant to the current instance of
            parameters
        identifier_observations (str): name of column in data table for unique
            identifiers of observations across rows
        method_scale (str): name of method to use to adjust the scale of values
            for features across observations, either 'z_score' or 'unit_range'
        data_path_directory (list<str>): names of directories in path at which
            to find the file for the table of data with features and
            observations for regression
        data_file (str): name of file for the table of data with features and
            observations for regression
        review (str): notes about review of instance in table of parameters
        note (str): notes about instance in table of parameters
        path_file_table_parameters (str): path to source file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Read source information from file.
    table = read_source_table_data(
        name_file_table_data=data_file,
        directories_path_data=data_path_directory,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    ##########
    # Evaluate information in table.
    # This function is useful to evaluate the variance of values for features
    # across observations after adjusting their scale to the common unit range.
    # This evaluation enables a manual, supervised approach to a method of
    # feature selection that is called 'variance thresholding'.
    if False:
        table = evaluate_table_data(
            table=table,
            selection_observations=selection_observations,
            features_relevant=features_relevant,
            features_regression=features_regression,
            features_continuity_scale=features_continuity_scale,
            index_columns_source="features",
            index_columns_product="features",
            index_rows_source=identifier_observations,
            index_rows_product="observations",
            adjust_scale=True,
            method_scale=method_scale, # 'z_score' or 'unit_range'
            explicate_indices=True,
            report=report,
        )
        pass

    ##########
    # Organize information in table.
    if True:
        table = organize_table_data(
            table=table,
            selection_observations=selection_observations,
            features_relevant=features_relevant,
            features_regression=features_regression,
            features_continuity_scale=features_continuity_scale,
            index_columns_source="features",
            index_columns_product="features",
            index_rows_source=identifier_observations,
            index_rows_product="observations",
            adjust_scale=True,
            method_scale=method_scale, # 'z_score' or 'unit_range'
            explicate_indices=True,
            report=False,
        )
        pass

    ##########
    # Check parameters and table of data for performing regression analysis.
    # It only makes sense to compare the variance of features if there has not
    # been scale standardization by z-score.
    if False:
        pail_check = check_parameters_table_data_regression(
            table=table,
            name_instance=name_instance,
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            features_predictor_random=features_predictor_random,
            threshold_features_variance=0.01,
            threshold_observations_count=5,
            measure_variance="variance",
            report=report,
        )
    else:
        pail_check = dict()
        pail_check["check_overall"] = True
        pass

    ##########
    # Perform regression analysis.
    record_extra = dict()
    record_extra["sequence"] = sequence
    record_extra["group"] = group
    record_extra["instance"] = instance
    record_extra["name_instance"] = name_instance
    pail_regression = preg.collect_record_regression_analysis(
        table=table,
        index_columns="features",
        index_rows="observations",
        check_overall=pail_check["check_overall"],
        type_regression=type_regression,
        formula_text=formula_text,
        feature_response=feature_response,
        features_predictor_fixed=features_predictor_fixed,
        features_predictor_random=features_predictor_random,
        record_extra=record_extra,
        report=report,
    )

    ##########
    # Write information to file.
    # Write for each parallel instance of regression.
    # A subsequent procedure will read the information from file and collect it
    # within a summary for all instances of regression.

    ##########
    # Collect information from all regressions in set.
    # Prepare summary.

    # different instances of regression will have different predictors...
    # that's a challenge for preparing the summary table.
    # 1. for all regression instances, determine max count of predictors
    # 2. summary table needs to accommodate that count of predictors
    # 3. within each group of columns for each predictor, the first should be
    # the name of the predictor variable, as below...
    # predictor_1_name predictor_1_coefficient predictor_1_error predictor_1_p ...
    # sex_y            0.0153                  0.0002            0.01
    # use a for loop to assemble the standardized records (for item in range(max count predictors))
    # access real information from predictors where available, then fill with missing


    # "predictors": ";".join(predictors),



    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: control_procedure_part_branch()")
        putly.print_terminal_partition(level=5)
        print("instance: " + instance)
        putly.print_terminal_partition(level=5)
        pass

    pass


##########
# Manage parallelization.


def control_parallel_instance(
    instance=None,
    parameters=None,
):
    """
    Control procedure for a single instance in parallel with others.

    Review: TCW; 31 March 2025

    arguments:
        instance (dict): parameters specific to current instance
            execution (int): logical binary indicator of whether to execute
                and handle the parameters for the current instance
            sequence (int): sequential index for instance's name and sort
                order
            group (str): categorical group of instances
            instance (str): name or designator of instance
            name_instance (str): compound name for instance of parameters
            selection_observations (dict<list<str>>): names of columns in data
                table for feature variables and their categorical values by
                which to filter rows for observations in data table
            type_regression (str): name of type of regression model to use
            formula_text (str): human readable formula for regression model,
                treated as a note for clarification
            feature_response (str): name of column in data table for feature
                variable to include in regression model as response dependent
                variable
            features_predictor_fixed (list<str>): names of columns in data
                table for feature variables to include in regression model as
                predictor independent variables with fixed effects
            features_predictor_random (list<str>): names of columns in data
                table for feature variables to include in regression model as
                predictor independent variables with random effects
            features_regression (list<str>): names of columns in data table for
                feature variables that are relevant to the actual regression
                model
            features_continuity_scale (list<str>): names of columns in data
                table for feature variables with values on quantitative,
                continuous scale of measurement, interval or ratio, for which
                to standardize the scale by z score
            features_relevant (list<str>): names of columns in data table for
                feature variables that are relevant to the current instance of
                parameters
            identifier_observations (str): name of column in data table for
                unique identifiers of observations across rows
            method_scale (str): name of method to use to adjust the scale of
                values for features across observations, either 'z_score' or
                'unit_range'
            data_path_directory (list<str>): names of directories in path at
                which to find the file for the table of data with features and
                observations for regression
            data_file (str): name of file for the data table of information
                about features and observations
            review (str): notes about review of instance in table of parameters
            note (str): notes about instance in table of parameters

        parameters (dict): parameters common to all instances
            path_file_table_parameters (str): path to source file in text format as
                a table with tab delimiters between columns and newline delimiters
                between rows, with parameters for multiple regressions
            path_directory_dock (str): path to dock directory for procedure's
                source and product directories and files
            report (bool): whether to print reports


    raises:

    returns:

    """

    ##########
    # Copy information.
    instance_record = copy.deepcopy(instance)

    ##########
    # Extract parameters.

    # Extract parameters specific to each instance.
    execution = instance_record["execution"]
    sequence = instance_record["sequence"]
    group = instance_record["group"]
    instance = instance_record["instance"]
    name_instance = instance_record["name_instance"]
    selection_observations = instance_record["selection_observations"]
    type_regression = instance_record["type_regression"]
    formula_text = instance_record["formula_text"]
    feature_response = instance_record["feature_response"]
    features_predictor_fixed = instance_record["features_predictor_fixed"]
    features_predictor_random = instance_record["features_predictor_random"]
    features_regression = instance_record["features_regression"]
    features_continuity_scale = instance_record["features_continuity_scale"]
    features_relevant = instance_record["features_relevant"]
    identifier_observations = instance_record["identifier_observations"]
    method_scale = instance_record["method_scale"]
    data_path_directory = instance_record["data_path_directory"]
    data_file = instance_record["data_file"]
    review = instance_record["review"]
    note = instance_record["note"]

    # Extract parameters common across all instances.
    path_file_table_parameters = parameters["path_file_table_parameters"]
    path_directory_dock = parameters["path_directory_dock"]
    report = parameters["report"]

    ##########
    # Control procedure with split branch for parallelization.
    control_procedure_part_branch(
        execution=execution,
        sequence=sequence,
        group=group,
        instance=instance,
        name_instance=name_instance,
        selection_observations=selection_observations,
        type_regression=type_regression,
        formula_text=formula_text,
        feature_response=feature_response,
        features_predictor_fixed=features_predictor_fixed,
        features_predictor_random=features_predictor_random,
        features_regression=features_regression,
        features_continuity_scale=features_continuity_scale,
        features_relevant=features_relevant,
        identifier_observations=identifier_observations,
        method_scale=method_scale,
        data_path_directory=data_path_directory,
        data_file=data_file,
        review=review,
        note=note,
        path_file_table_parameters=path_file_table_parameters,
        path_directory_dock=path_directory_dock,
        report=report,
    )
    pass


def control_parallel_instances(
    instances=None,
    path_file_table_parameters=None,
    path_directory_dock=None,
    report=None,
):
    """
    Control procedure for parallel instances.

    Review: TCW; 31 March 2025

    arguments:
        instances (list<dict>): parameters to control individual instances in
            parallel
        path_file_table_parameters (str): path to source file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports


    raises:

    returns:

    """

    # Collect parameters common across all instances.
    parameters = dict()
    parameters["path_file_table_parameters"] = path_file_table_parameters
    parameters["path_directory_dock"] = path_directory_dock
    parameters["report"] = report

    # Execute procedure iteratively with parallelization across instances.
    if False:
        prall.drive_procedure_parallel(
            function_control=(
                control_parallel_instance
            ),
            instances=instances,
            parameters=parameters,
            cores=5,
            report=True,
        )
    else:
        # Execute procedure directly for testing.
        control_parallel_instance(
            instance=instances[0],
            parameters=parameters,
        )
    pass


##########
# Call main procedure.


def execute_procedure(
    path_file_table_parameters=None,
    path_directory_dock=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    1. read source information from file
       - source information at this level consists of a table of parameters for
          regressions
       - it is necessary to parse information with hierarchical structure from
          the flat text
    2. drive multiple, concurrent procedures in parallel branches
    3. within each procedural branch, read from file and organize the data
       table consisting of features across observations
       - filter features and observations in table
       - standardize scale of values within specific feature variables
    4. within each procedural branch, execute regression analysis using type of
       model, features, and other parameters from the source table of
       parameters
    5. organize report summary information from the regression analysis and
       write this information to file

    arguments:
        path_file_table_parameters (str): path to source file in text format as
            a table with tab delimiters between columns and newline delimiters
            between rows, with parameters for multiple regressions
        path_directory_dock (str): path to dock directory for procedure's
            source and product directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    ##########
    # Parameters.
    if (
        (report is not None) and
        (str(report) != "") and
        (str(report) == "true")
    ):
        report = True
    else:
        report = False
        pass

    ##########
    # Read source information from file.
    pail_source = read_source_table_parameters(
        path_file_table_parameters=path_file_table_parameters,
        path_directory_dock=path_directory_dock,
        report=report,
    )
    for record in pail_source["records"]:
        print(record)
        pass


    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        print("module: script_drive_regressions_from_table_parameters.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("path_file_table_parameters: " + str(path_file_table_parameters))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass


    ##########
    # Control procedure for parallel instances.
    control_parallel_instances(
        instances=pail_source["records"],
        path_file_table_parameters=path_file_table_parameters,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_table_parameters = sys.argv[1]
    path_directory_dock = sys.argv[2]
    report = sys.argv[3]

    # Call function for procedure.
    execute_procedure(
        path_file_table_parameters=path_file_table_parameters,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass



#
