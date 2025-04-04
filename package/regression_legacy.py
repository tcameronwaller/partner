"""
Supply functionality for regression analysis.

This module 'regression' is part of the 'partner' package.

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

###############################################################################
# Notes

# TODO: TCW; 26 June 2024
# Support more regression models
# Ordinary Least Squares
#   - https://www.statsmodels.org/dev/examples/notebooks/generated/ols.html
# Generalized Linear Models
#   - https://www.statsmodels.org/stable/glm.html
# Mixed effects models
# - support pairs of samples from the same subject
# ANOVA
# ANCOVA


###############################################################################
# Installation and importation

# Standard
import os
import csv
import copy
import textwrap
import string
import gzip
import shutil
import textwrap
import itertools
import math

# Relevant

import pandas
import scipy
import numpy
import statsmodels.api
import statsmodels.stats.outliers_influence

# Custom
import partner.utility as utility # this import path for subpackage
import partner.scale as pscale
import partner.description as pdesc

#dir()
#importlib.reload()

###############################################################################
# Functionality


# TODO: TCW; 2 April 2025
# obsolete and not interesting or useful.
def organize_table_cohort_model_variables_for_regression(
    dependence=None,
    independence=None,
    standard_scale=None,
    table=None,
    report=None,
):
    """
    Organizes a table and relevant information for regression analysis.

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    This function accepts pre-determined dependent and independent variables.

    1. select relevant columns in table
    2. drop records with any missing variables
    3. transform dependent and independent variables to z-score standard scale

    arguments:
        dependence (str): name of table's column for dependent variable
        independence (list<str>): names of table's columns for independent
            variables
        standard_scale (bool): whether to transform all independent variables to
            z-score standard scale (mean of zero and standard deviation of one)
        table (object): Pandas data frame of dependent and independent variables
            for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information for regression

    """

    # Copy information.
    table = table.copy(deep=True)
    # Determine whether variables exist in the table.
    #independence_match = all(
    #    lambda variable: (str(variable) in table.columns.to_list()),
    #    independence
    #)
    if (dependence not in table.columns.to_list()):
        dependence_product = ""
    else:
        dependence_product = dependence
    # Select table's columns for relevant variables.
    columns = copy.deepcopy(independence)
    if (dependence in table.columns.to_list()):
        columns.insert(0, dependence)
    table = table.loc[:, table.columns.isin(columns)]
    table = table[[*columns]]

    # Remove columns with inadequate non-missing values or inadequate variance
    # (standard deviation) across rows.
    # This preliminary filter on columns avoids losing table rows for missing
    # values in an unnecessary column.
    # An unnecessary column is any column with minimal non-missing values or
    # with minimal variance across those values.
    table = utility.filter_table_columns_by_nonmissing_variance(
        threshold_valid_proportion_per_column=0.001,
        threshold_column_variance=0.001,
        type_variance="standard_deviation",
        table=table,
        report=report,
    )
    # Drop any rows with missing values in any columns.
    table.dropna(
        axis="index", # drop rows
        how="any",
        subset=None,
        inplace=True,
    )
    # Remove columns with inadequate non-missing values or inadequate variance
    # (standard deviation) across rows.
    # Filter columns again since the removal of missing values changes the
    # variance of columns.
    table = utility.filter_table_columns_by_nonmissing_variance(
        threshold_valid_proportion_per_column=0.05,
        threshold_column_variance=0.001,
        type_variance="standard_deviation",
        table=table,
        report=report,
    )

    # Determine whether to transform all dependent and independent variables to
    # z-score standard scale.
    # Standardization introduces missing values if standard deviation is zero.
    if (standard_scale):
        table = pscale.transform_standard_z_score_by_table_columns(
            table=table,
            columns=independence,
            report=report,
        )
        table.dropna(
            axis="index",
            how="any",
            subset=None,
            inplace=True,
        )
        pass
    # Determine whether the variables have changed.
    # Expect errors if there is removal of the column for the dependent
    # variable.
    independence_product = list(filter(
        lambda column: (str(column) in table.columns.to_list()),
        independence
    ))

    # Determine count of valid samples (cases, observations).
    count_samples = int(table.shape[0])
    # Collect information for regression.
    pail = dict()
    pail["table"] = table
    pail["count_samples"] = count_samples
    pail["dependence"] = dependence_product
    pail["independence"] = independence_product
    # Return information.
    return pail


###############################################
###
##
#
# Interesting and useful...


def create_missing_regression_independent_variable(
    variable=None,
):
    """
    Creates missing values for a single independent variable in a linear
    regression.

    arguments:
        variable (str): name of an independent variable

    raises:

    returns:
        (dict): collection of missing values for a linear regression on the
            independent variable
    """

    # Collect information for independent variable.
    pail = dict()
    # Create missing values for the independent variable in linear regression.
    pail["report_b95ci"] = str("b: NAN; 95% CI: NAN ... NAN")
    pail["report_bep"] = str("b: NAN (NAN); p: NAN")
    pail["parameter"] = float("nan")
    pail["error"] = float("nan")
    pail["interval_99"] = float("nan")
    pail["interval_95"] = float("nan")
    pail["range_95"] = str("NAN ... NAN")
    pail["range_95_below"] = float("nan")
    pail["range_95_above"] = float("nan")
    pail["probability"] = float("nan")
    pail["inflation"] = float("nan")
    # Return information.
    return pail


def create_missing_values_regression_linear(
    dependence=None,
    independence=None,
):
    """
    Creates missing values for a regression.

    arguments:
        dependence (str): name of table's column for dependent variable
        independence (list<str>): names of table's columns for independent
            variables

    raises:

    returns:
        (dict): collection of regression's residuals and statistics
    """

    # Create missing values for statistics on the whole regression model.
    summary = {
        "dependence_actual": "",
        "independence_actual": "",
        "freedom": float("nan"),
        "observations": float("nan"),
        "samples": float("nan"),
        "r_square": float("nan"),
        "r_square_adjust": float("nan"),
        "log_likelihood": float("nan"),
        "akaike": float("nan"),
        "bayes": float("nan"),
        "condition": float("nan"),
    }
    # Include intercept with independent variables.
    independences_names = copy.deepcopy(independence)
    independences_names.append("intercept")

    # Create missing values for regression.
    #probabilities = list(
    #    itertools.repeat(float("nan"), len(values_independence_intercept))
    #)
    pail_tree = dict()
    for variable in independences_names:
        # Create missing values for the independent variable.
        pail_tree[variable] = create_missing_regression_independent_variable(
            variable=variable,
        )
        pass
    # Organize missing residuals.
    residuals = numpy.empty(0)
    # Collect information.
    summary["independence_tree"] = pail_tree
    pail = dict()
    pail["summary"] = summary
    pail["residuals"] = residuals
    # Return information.
    return pail


def create_missing_values_regression_logistic(
    dependence=None,
    independence=None,
):
    """
    Creates missing values for a regression.

    arguments:
        dependence (str): name of table's column for dependent variable
        independence (list<str>): names of table's columns for independent
            variables

    raises:

    returns:
        (dict): collection of regression's residuals and statistics
    """

    # Create missing values for statistics on the whole regression model.
    summary = {
        "dependence_actual": "",
        "independence_actual": "",
        "freedom": float("nan"),
        "observations": float("nan"),
        "samples": float("nan"),
        "r_square_pseudo": float("nan"),
        "log_likelihood": float("nan"),
        "akaike": float("nan"),
        "bayes": float("nan"),
    }
    # Include intercept with independent variables.
    independences_names = copy.deepcopy(independence)
    independences_names.append("intercept")

    # Create missing values for regression.
    pail_tree = dict()
    for variable in independences_names:
        # Create missing values for the independent variable.
        pail_tree[variable] = create_missing_regression_independent_variable(
            variable=variable,
        )
        pass
    # Organize missing residuals.
    residuals = numpy.empty(0)
    # Collect information.
    summary["independence_tree"] = pail_tree
    pail = dict()
    pail["summary"] = summary
    pail["residuals"] = residuals
    # Return information.
    return pail


# inspect
def create_regression_missing_values(
    dependence=None,
    independence=None,
    type=None,
):
    """
    Creates missing values for a regression.

    arguments:
        dependence (str): name of table's column for dependent variable
        independence (list<str>): names of table's columns for independent
            variables
        type (str): type of regression analysis, either 'linear' or 'logistic'

    raises:

    returns:
        (dict): collection of missing values for regression
    """

    # Determine type of regression.
    if (type == "linear"):
        pail = create_missing_values_regression_linear(
            dependence=dependence,
            independence=independence,
        )
    elif (type == "logistic"):
        pail = create_missing_values_regression_logistic(
            dependence=dependence,
            independence=independence,
        )
    else:
        pail = dict()
        pass
    # Return information.
    return pail


# Organize regression summary table for Forest Plot with two groups


def organize_regression_summary_table_for_forest_plots(
    type=None,
    model_contexts=None,
    model_adjustments=None,
    column_variable=None,
    variables=None,
    columns_translations=None,
    labels_categories=None,
    sorts_categories=None,
    column_stratification=None,
    table=None,
    report=None,
):
    """
    This function organizes information from a regression summary table for
    plots.

    The regression summary table ought to represent regressions for a single
    dependent variable across multiple cohorts and regression models.

    arguments:
        type (str): type of regression analysis, either 'linear' or 'logistic'
        model_contexts (list<str>): type of regression model context, 'joint' or
            'marginal'
        model_adjustments (list<str>): type of regression model adjustments,
            'adjust' or 'unadjust'
        column_variable (str): name of table's column for variables for
            selection
        variables (list<str>): names of variables for plots
        columns_translations (dict<str>): translations of names of columns
        labels_categories (dict<str>): labels for categories from independent
            variables
        sorts_categories (dict<str>): sort orders for categories from
            independent variables
        column_stratification (str): name of column for use in stratification of
            records into separate tables
        table (object): Pandas data frame of summary information from multiple
            regressions
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables for plots

    """

    # Copy information.
    table = table.copy(deep=True)
    # Filter records by type of regression, 'linear' or 'logistic'.
    table = table.loc[(table["dependence_type"] == type), :]
    # Filter records by context of regression model, 'joint' or 'marginal'.
    table = table.loc[(table["model_context"].isin(model_contexts)), :]
    # Filter records by context of regression model, 'joint' or 'marginal'.
    table = table.loc[(table["model_adjustment"].isin(model_adjustments)), :]
    # Filter records by values of the independent variable.
    table = table.loc[(table[column_variable].isin(variables)), :]
    # Change names of columns.
    table.rename(
        columns=columns_translations,
        inplace=True,
    )
    table["interval_above"] = table["interval_below"] # only if symmetrical
    table["category_sort"] = table.apply(
        lambda row: sorts_categories[row["category"]],
        axis="columns", # apply function to each row
    )
    table["category_label"] = table.apply(
        lambda row: labels_categories[row["category"]],
        axis="columns", # apply function to each row
    )
    # Stratify records in the table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    table.set_index(
        [column_stratification],
        append=False,
        drop=True,
        inplace=True
    )
    tables_stratification = table.groupby(
        level=column_stratification,
        axis="index",
    )
    # Collect stratification tables for each group.
    pail_tables = dict()
    for name, table_stratification in tables_stratification:
        # Select relevant columns.
        columns = [
            "group", "category", "category_sort", "category_label",
            "value", "interval_below", "interval_above",
        ]
        table_stratification = table_stratification.loc[
            :, table_stratification.columns.isin(columns)
        ]
        # Sort columns.
        table_stratification = table_stratification[[*columns]]
        table_stratification.reset_index(
            level=None,
            inplace=True,
            drop=True,
        )
        # Collect table.
        pail_tables[name] = table_stratification
        pass
        # Report.
        if report:
            utility.print_terminal_partition(level=5)
            print("report: ")
            function_name = str(
                "organize_regression_table_for_forest_plots()"
            )
            print(function_name)
            utility.print_terminal_partition(level=5)
            print(name)
            print(table_stratification)
            pass
    # Return information.
    return pail_tables

#
##
###
########################################################3

# Drivers


# inspect
# order: 2
# TODO: TCW; 2 April 2025
# obsolete and not interesting or useful.
def organize_check_table_information_for_regression(
    dependence=None,
    independence=None,
    standard_scale=None,
    threshold_samples=None,
    table=None,
    type=None,
    report=None,
):
    """
    Drive the organization of a table and regression.

    Table format must have samples (cases, observations) across rows and
    dependent and independent variables (features) across columns.

    arguments:
        dependence (str): name of table's column for dependent variable
        independence (list<str>): names of table's columns for independent
            variables
        standard_scale (bool): whether to transform all independent variables to
            z-score standard scale (mean of zero and standard deviation of one)
        threshold_samples (float): minimal count of samples with non-missing
            values of dependent and independent variables to perform regression
        table (object): Pandas data frame of dependent and independent variables
            for regression
        type (str): type of regression analysis, either 'linear' or 'logistic'
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of regression's residuals and statistics
    """

    # Report.
    if report:
        utility.print_terminal_partition(level=5)
        print("report: ")
        function_name = str(
            "organize_check_table_information_for_regression()"
        )
        print(function_name)
        pass

    # Organize table for regression.
    pail_organization = organize_table_cohort_model_variables_for_regression(
        dependence=dependence,
        independence=independence,
        standard_scale=standard_scale,
        table=table,
        report=False,
    )

    # Determine whether the table includes sufficient information for
    # regression.
    if (
        (len(pail_organization["dependence"]) > 0) and
        (len(pail_organization["independence"]) > 0) and
        (pail_organization["count_samples"] >= threshold_samples) and
        ((type == "linear") or (type == "logistic"))
    ):
        # There is sufficient information for regression.
        pail_organization["validity"] = True
        # Report.
        if report:
            print("Sufficient information for regression.")
            utility.print_terminal_partition(level=5)
    else:
        # There is insufficient information for regression.
        pail_organization["validity"] = False
        # Report.
        if report:
            print("Insufficient information for regression.")
            utility.print_terminal_partition(level=5)
            # Determine and report whether there is insufficient information.
            if (
                (len(pail_organization["dependence"]) == 0) or
                (len(pail_organization["independence"]) == 0)
            ):
                print(
                    "Dependent or independent variables do not exist in table."
                )
                pass
            if (
                (not (pail_organization["count_samples"] >= threshold_samples))
            ):
                print("Insufficient samples in the table.")
                pass
            if (not ((type == "linear") or (type == "logistic"))):
                print("Regression type is neither 'linear' nor 'logistic'.")
                pass
            pass
        pass
    # Return information.
    return pail_organization


# inspect
# order: 1
def organize_check_table_drive_regression(
    table=None,
    type=None,
    dependence=None,
    independence=None,
    report=None,
):
    """
    Organize regressions.

    Table format must have samples (cases, observations) across rows and
    dependent and independent variables (features) across columns.

    arguments:
        table (object): Pandas data frame of dependent and independent variables
            (features) across columns and samples (cases, observations) within
            a specific cohort across rows
        type (str): type of regression analysis, either 'linear' or 'logistic'
        dependence (str): name of table's column for dependent variable
        independence (list<str>): names of table's columns for independent
            variables
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regressions

    """

    # Organize and check information for regression on specific cohort and
    # model.
    pail_check = organize_check_table_information_for_regression(
        dependence=dependence,
        independence=independence,
        standard_scale=True,
        threshold_samples=10,
        table=table,
        type=type,
        report=report,
    )
    # Determine whether to execute the regression.
    if (pail_check["validity"]):
        # Adequate information for regression.
        # Execute regression analysis.
        if (type == "linear"):
            pail_regression = regress_linear_ordinary_least_squares(
                dependence=dependence,
                independence=pail_check["independence"],
                table=pail_check["table"],
                report=report,
            )
        elif (type == "logistic"):
            pail_regression = regress_discrete_logit(
                dependence=dependence,
                independence=pail_check["independence"],
                table=pail_check["table"],
                report=report,
            )
            pass
        # Determine whether any independent variables need missing values.
        if (len(pail_check["independence"]) < len(independence)):
            # Some independent variables were excluded from the regression.
            # Create missing values for these independent variables.
            independence_missing = list(filter(
                lambda variable: (variable not in pail_check["independence"]),
                independence
            ))
            for variable in independence_missing:
                # Create missing values for the independent variable.
                pail_regression["summary"]["independence_tree"][variable] = (
                    create_missing_regression_independent_variable(
                        variable=variable,
                ))
                pass
            pass
        pass
    else:
        # Inadequate information for regression.
        # Create missing values.
        pail_regression = create_regression_missing_values(
            dependence=dependence,
            independence=independence,
            type=type,
        )
    # Return information.
    return pail_regression


# inspect
# order: last, I think
def organize_regressions_summary_table_long(
    records_regressions=None,
    independences_summary=None,
    type=None,
    report=None,
):
    """
    This function organizes the summary table with information for multiple
    regressions.

    A primary function is to filter the table.

    arguments:
        records_regressions (list<dict>): information about regressions
        independences_summary (list<str>): names of independent variables for
            which to include information in the summary table, or "None" to
            include information for all original independent variables
        type (str): type of regression analysis, either 'linear' or 'logistic'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of summary information from multiple
            regressions

    """

    # Collect information.
    records = list()
    # Iterate on regressions.
    for record_regression in records_regressions:
        # Extract names of entries for statistics on the whole model.
        keys_general = copy.deepcopy(list(
            record_regression.keys()
        ))
        entries_general = list(filter(
            lambda entry: ("independence_tree" not in str(entry)),
            keys_general
        ))
        # Extract names of independent variables to include in summary.
        # If 'independences_summary' is not 'None', then only include
        # information in summary table about independent variables that are in
        # the list 'independences_summary'.
        keys_tree = copy.deepcopy(list(
            record_regression["independence_tree"].keys()
        ))
        if (
            (independences_summary is not None) and
            (len(independences_summary) > 0)
        ):
            independences_inclusion = list(filter(
                lambda variable: (str(variable) in independences_summary),
                keys_tree
            ))
        else:
            independences_inclusion = (
                keys_tree
            )
        # Iterate on independent variables.
        for variable in independences_inclusion:
            # Collection information.
            # There is one record for each independent variable in the model.
            record = dict()
            # Collect entries for statistics on the whole model.
            for entry in entries_general:
                record[entry] = copy.deepcopy(record_regression[entry])
                pass
            # Collect name of variable.
            record["variable_key"] = variable
            # Extract names of entries for statistics on the independent
            # variable.
            entries_variable = copy.deepcopy(list(
                record_regression["independence_tree"][variable].keys()
            ))
            # Collect entries for statistics on the independent variable.
            for entry in entries_variable:
                record[entry] = copy.deepcopy(
                    record_regression["independence_tree"][variable][entry]
                )
                pass
            # Collect information.
            records.append(record)
        pass
    pass
    # Organize table.
    table = pandas.DataFrame(data=records)
    table = organize_table_regression_summary(
        type=type,
        table=table,
        report=report,
    )
    # Return information.
    return table


# inspect
def organize_table_regression_summary(
    type=None,
    table=None,
    report=None,
):
    """
    This function organizes the summary table with information for multiple
    regressions.

    A primary function is to filter and sort columns within the table.

    arguments:
        type (str): type of regression analysis, either 'linear' or 'logistic'
        table (object): Pandas data frame of summary information from multiple
            regressions
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of summary information from multiple
            regressions

    """

    # Copy information.
    table = table.copy(deep=True)
    # Define columns according to regression type.
    columns = [
        "cohort", "dependence",
        "model_adjustment", "model_context",
        "variable",
        "report_bep", "report_b95ci",
        "dependence_type",
        "model", "model_note",
        "independence",
        "dependence_actual", "independence_actual",
        "freedom", "observations", "samples",
        "log_likelihood", "akaike", "bayes",
        "variable_key",
        "parameter", "error", "interval_95", "interval_99",
        "range_95", "range_95_low", "range_95_high",
        "range_99", "range_99_low", "range_99_high",
        "probability", "inflation",
    ]
    if (type == "linear"):
        columns.extend(["r_square", "r_square_adjust","condition",])
    elif (type == "logistic"):
        columns.extend(["r_square_pseudo",])
    else:
        columns = list()
        pass
    # Select columns.
    # columns.insert(0, dependence)
    table = table.loc[:, table.columns.isin(columns)]
    # Sort columns.
    table = table[[*columns]]
    # Report.
    if report:
        utility.print_terminal_partition(level=5)
        print("report: ")
        function_name = str(
            "organize_table_regression_summary()"
        )
        print(function_name)
        utility.print_terminal_partition(level=5)
        print(table)
        pass
    # Return information.
    return table



# TODO: TCW; 2 April 2025
# Obsolete, I think. See "script_drive_regressions_from_table_parameters".
def drive_linear_logistic_regression_cohort_model(
    record_cohort_model=None,
    entries_cohorts=None,
    report=None,
):
    """
    Organize and drive linear or logistic regression on a specific cohort and
    model.

    arguments:
        record_cohort_model (dict): information about a regression on a specific
            cohort with specific dependent and independent variables
        entries_cohorts (dict<dict<object>>): information about phenotype
            variables within Pandas data frame tables that have stratification
            for relevant cohorts with names that match "table_cohorts_models"
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regression on a specific cohort and model

    """

    # Organize summary record for regression.
    record_regression = dict()
    record_regression["cohort"] = record_cohort_model["cohort"]
    record_regression["dependence"] = record_cohort_model["dependence"]
    record_regression["dependence_type"] = (
        record_cohort_model["dependence_type"]
    )
    record_regression["model"] = record_cohort_model["model"]
    record_regression["model_context"] = record_cohort_model["model_context"]
    record_regression["model_adjustment"] = (
        record_cohort_model["model_adjustment"]
    )
    record_regression["model_note"] = record_cohort_model["model_note"]
    record_regression["independence"] = record_cohort_model["independence"]
    # Report.
    if report:
        utility.print_terminal_partition(level=3)
        print("report: ")
        report_function_name = str(
            "drive_linear_logistic_regression_cohort_model()"
        )
        print(report_function_name)
        utility.print_terminal_partition(level=5)
        print("cohort: " + str(record_cohort_model["cohort"]))
        print("dependence: " + str(record_cohort_model["dependence"]))
        print("independent variables: ")
        print(record_cohort_model["independence"])
        utility.print_terminal_partition(level=5)
        pass
    # Determine whether there is sufficient information for regression.
    # 1. Variable table exists for specific cohort.
    if (record_cohort_model["cohort"] in entries_cohorts.keys()):
        # Determine variable table for stratification cohort.
        table_cohort = (
            entries_cohorts[record_cohort_model["cohort"]]["table"]
        )
        # Extract names of independent variables.
        independence_string = str(record_cohort_model["independence"])
        independence_split = independence_string.split(";")
        independence_strip = list(map(
            lambda value: value.strip(), independence_split
        ))
        # Drive regression.
        pail_regression = organize_check_table_drive_regression(
            table=table_cohort,
            type=record_cohort_model["dependence_type"],
            dependence=record_cohort_model["dependence"],
            independence=independence_strip,
            report=report,
        )
    else:
        # Variable table does not exist for current cohort.
        # Report.
        if report:
            utility.print_terminal_partition(level=5)
            print("Missing regression information for cohort...")
            utility.print_terminal_partition(level=5)
        # Create missing values for regression.
        pail_regression = create_regression_missing_values(
            dependence=record_cohort_model["dependence"],
            independence=record_cohort_model["independence"],
            type=record_cohort_model["dependence_type"],
        )
    # Collect information from regression.
    record_regression.update(pail_regression["summary"])

    # Return information.
    return record_regression


# TODO: TCW; 2 April 2025
# Obsolete, I think. See "script_drive_regressions_from_table_parameters".
def drive_linear_logistic_regressions_cohorts_models(
    entries_cohorts=None,
    table_cohorts_models=None,
    independences_summary=None,
    filter_execution=None,
    type=None,
    report=None,
):
    """
    Organize and drive linear or logistic regressions.

    Format of entry "table" within "entries_cohorts"...
    Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows, with an explicit index

    Format of "table_cohorts_models"...
    columns: "execution", "cohort", "cohort_sort", "dependence",
    "dependence_sort", "dependence_type", "independence",
    "model", "model_sort", "model_context", "model_adjustment", "model_note"
    Rows in Column "dependence" have actual names of variable columns in the
    main source table.
    Rows in Column "independence" have list of independent variables with
    semicolon (";") delimiters.

    arguments:
        entries_cohorts (dict<dict<object>>): information about phenotype
            variables within Pandas data frame tables that have stratification
            for relevant cohorts with names that match "table_cohorts_models"
        table_cohorts_models (object): Pandas data frame that specifies cohorts,
            dependent variables, and independent variables (models) for
            regression
        independences_summary (list<str>): names of independent variables for
            which to include information in the summary table, or "None" to
            include information for all original independent variables
        filter_execution (bool): whether to filter records in cohort-model table
            by logical binary "execution" variable
        type (str): type of regression analysis, either 'linear' or 'logistic'
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regressions

    """

    # Filter relevant cohorts-dependences-models by "execution" flag.
    if filter_execution:
        table_cohorts_models = table_cohorts_models.loc[
            (
                (table_cohorts_models["execution"] == 1)
            ), :
        ]
        pass
    # Filter relevant cohorts-dependences-models by "dependence_type".
    table_cohorts_models = table_cohorts_models.loc[
        (
            (table_cohorts_models["dependence_type"] == type)
        ), :
    ]
    # Iterate on records that specify cohors, dependent variables, and
    # independent variables for regressions.
    records_cohorts_models = table_cohorts_models.to_dict(
        orient="records",
    )

    # Collect summary records for each regression.
    records_regressions = list()
    # Iterate across cohorts.
    for record_cohort_model in records_cohorts_models:
        # Drive regression on single cohort and model.
        record_regression = (
            drive_linear_logistic_regression_cohort_model(
                record_cohort_model=record_cohort_model,
                entries_cohorts=entries_cohorts,
                report=report,
        ))
        # Collect records.
        records_regressions.append(record_regression)
        pass

    # Organize table.
    #table_regressions_raw = pandas.DataFrame(data=records_regressions)
    table_regressions_long = organize_regressions_summary_table_long(
        records_regressions=records_regressions,
        independences_summary=independences_summary,
        type=type,
        report=report,
    )
    # Compile information.
    pail = dict()
    pail["table"] = table_regressions_long
    # Return information.
    return pail



##########
# Write


def write_product_table(
    name=None,
    table=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        name (str): base name for file
        table (object): Pandas data-frame table to write to file
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    # Reset index.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    # Specify directories and files.
    path_table = os.path.join(
        path_directory, str(name + ".tsv")
    )
    # Write information to file.
    table.to_csv(
        path_or_buf=path_table,
        sep="\t",
        header=True,
        index=False,
    )
    pass


def write_product_tables(
    pail_write=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict<object>): collection of information to write to file
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    for name in pail_write.keys():
        write_product_table(
            name=name,
            table=pail_write[name],
            path_directory=path_directory,
        )
    pass


def write_product(
    pail_write=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict<dict<object>>): collection of information to write to
            file
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    # Export information.
    write_product_tables(
        pail_write=pail_write["tables"],
        path_directory=path_directory,
    )
    pass




#
