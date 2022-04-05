"""
Supply functionality for Regression Analysis.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of Promiscuity
    (https://github.com/tcameronwaller/promiscuity/).

    Promiscuity supports data analysis in multiple other projects.
    Copyright (C) 2021 Thomas Cameron Waller

    Promiscuity is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    Promiscuity is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with Promiscuity. If not, see <http://www.gnu.org/licenses/>.
"""

###############################################################################
# Notes

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
import sklearn
import sklearn.decomposition
import scipy
import numpy
import statsmodels.api
import statsmodels.stats.outliers_influence

# Custom

import promiscuity.utility as utility # this import path for subpackage


#dir()
#importlib.reload()

###############################################################################
# Functionality


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
        threshold_valid_proportion_per_column=0.05,
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
        table = utility.standardize_scale_values_specific_table_columns(
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


def determine_confidence_interval_range_text(
    estimate=None,
    interval_low=None,
    interval_high=None,
):
    """
    Prepares a textual representation of the confidence interval range about an
    estimate.

    arguments:
        estimate (float): value of an estimate
        interval_low (float): value of lower confidence interval
        interval_high (float): value of higher confidence interval

    raises:

    returns:
        (str): textual representation of confidence interval range
    """

    range = str(
        str(round((float(estimate) - float(interval_low)), 5)) +
        " ... " +
        str(round((float(estimate) + float(interval_high)), 5))
    )
    return range


def organize_linear_logistic_regression_independence_tree(
    independence=None,
    model_parameters=None,
    model_parameter_errors=None,
    model_probabilities=None,
    table_independence_intercept=None,
):
    """
    Organizes a dictionary tree of information about each independent variable
    from a logistic or linear regression.

    arguments:
        independence (list<str>): names of table's columns for independent
            variables
        model_parameters (dict): parameter coefficients
        model_parameter_errors (dict): standard errors of estimates of parameter
            coefficients
        model_probabilities (dict): probabilities of estimates of parameter
            coefficients
        table_independence_intercept (object): Pandas data frame of independent
            variables with a constant intercept, the same source from the
            regression

    raises:

    returns:
        (dict<dict>): dictionary tree of information about each independent
            variable from a logistic or linear regression
    """

    # Copy information.
    table_independence_intercept = table_independence_intercept.copy(deep=True)

    # Collect information about independent variables.
    pail_tree = dict()

    # Intercept.
    if (
        ("const" in model_parameters.index) and
        ("const" in model_parameter_errors.index) and
        ("const" in model_probabilities.index)
    ):
        # Collect information for intercept.
        pail_tree["intercept"] = dict()
        pail_tree["intercept"]["variable"] = "intercept"
        #pail_tree["intercept"]["parameter"] = report.params[0]
        pail_tree["intercept"]["parameter"] = model_parameters["const"]
        pail_tree["intercept"]["error"] = model_parameter_errors["const"]
        pail_tree["intercept"]["interval_95"] = float(
            1.96 * pail_tree["intercept"]["error"]
        )
        pail_tree["intercept"]["range_95"] = (
            determine_confidence_interval_range_text(
                estimate=pail_tree["intercept"]["parameter"],
                interval_low=pail_tree["intercept"]["interval_95"],
                interval_high=pail_tree["intercept"]["interval_95"],
        ))
        pail_tree["intercept"]["inflation"] = float("nan")
        pail_tree["intercept"]["probability"] = model_probabilities["const"]
    else:
        # Report.
        if report:
            print("Warning: regression data does not have constant intercept.")
        # Create missing values for intercept.
        pail_tree["intercept"] = create_missing_regression_independent_variable(
            variable="intercept",
        )
    # Independent variables.
    # Initiate counter at 1 to assume that intercept is at index 0.
    counter = 1
    # Accommodate index for intercept.
    for variable in independence:
        # Collect information for intercept.
        pail_tree[variable] = dict()
        pail_tree[variable]["variable"] = str(variable)
        # Coefficient or parameter.
        #pail_tree[variable]["parameter"] = report.params[counter]
        pail_tree[variable]["parameter"] = model_parameters[variable]
        # Parameter standard error
        pail_tree[variable]["error"] = model_parameter_errors[variable]
        pail_tree[variable]["interval_95"] = float(
            1.96 * pail_tree[variable]["error"]
        )
        pail_tree[variable]["range_95"] = (
            determine_confidence_interval_range_text(
                estimate=pail_tree[variable]["parameter"],
                interval_low=pail_tree[variable]["interval_95"],
                interval_high=pail_tree[variable]["interval_95"],
        ))
        # Probability.
        pail_tree[variable]["probability"] = model_probabilities[variable]
        # Variance Inflation Factor (VIF).
        inflation_value = (
            statsmodels.stats.outliers_influence.variance_inflation_factor(
                table_independence_intercept.to_numpy(),
                counter
            )
        )
        pail_tree[variable]["inflation"] = round(inflation_value, 3)
        # Increment index.
        counter += 1
        pass
    # Return information.
    return pail_tree


def regress_discrete_logit(
    dependence=None,
    independence=None,
    table=None,
    report=None,
):
    """
    Regresses a binary dependent variable against multiple independent variables
    and returns relevant parameters and statistics.

    Table format must have samples (cases, observations) across rows and
    dependent and independent variables (features) across columns.

    Description of formats for StatsModels...

    Format of dependent variable is a vector of scalar values.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4]

    Format of independent variable(s) is a matrix: a first-dimension vector of
    samples (observations) and for each sample a second-dimension vector of
    variables' (features') scalar values.
    StatsModels also requires a constant for the intercept.
    [
        [1.3, 5.2, 1.0],
        [1.5, 5.1, 1.0],
        [1.2, 5.5, 1.0],
        ...
    ]

    arguments:
        dependence (str): name of table's column for dependent variable
        independence (list<str>): names of table's columns for independent
            variables
        table (object): Pandas data frame of dependent and independent variables
            for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of regression's residuals and statistics
    """

    # Determine count of valid samples (cases, observations).
    count_samples = int(table.shape[0])
    # Extract values of dependent and independent variables.
    # Can pass values of dependent and independent variables as NumPy arrays or
    # as Pandas Series and Dataframe respectively.
    # Passing variables as Pandas Series and Dataframe preserves variable names.
    values_dependence = table[dependence].to_numpy()

    # Keep independent variables in Pandas dataframe to preserve variables'
    # names.
    #values_independence = data.loc[ :, independence].to_numpy()
    table_independence = table.loc[ :, independence]
    # Introduce constant value for intercept.
    # If any column in the independent variables already has constant
    # values, then the function skips it by default.
    # It is necessary to change parameter "has_constant" to avoid this
    # conditional behavior.
    table_independence_intercept = statsmodels.api.add_constant(
        table_independence,
        prepend=True, # insert intercept constant first
        has_constant="add", # introduce new intercept constant regardless
    )
    columns_independence = copy.deepcopy(
        table_independence_intercept.columns.to_list()
    )
    #matrix_independence = table.to_numpy()
    # Define model.
    # statsmodels.discrete.discrete_model.Logit
    model = statsmodels.api.Logit(
        table[dependence],
        table_independence_intercept,
        missing="drop",
    )
    pail_raw = model.fit()
    # Report.
    if report:
        print("--------------------------------------------------")
        print(
            "Report source: " +
            "regress_dependent_independent_variables_linear_ordinary()"
        )
        print("--------------------------------------------------")
        print("Version check: TCW 28 September 2021")
        print("Information from regression:")
        print(pail_raw.summary())
        #utility.print_terminal_partition(level=3)
        #print(dir(pail_raw))
        #print(pail_raw.params)
        #print(pail_raw.pvalues)
        pass

    # Organize residuals.
    residuals = pail_raw.resid_generalized

    ##########
    # Collect parameters, errors, probabilities, and statistics.
    model_parameters = pandas.Series(data=pail_raw.params)
    model_parameter_errors = pandas.Series(data=pail_raw.bse)
    model_probabilities = pandas.Series(data=pail_raw.pvalues)
    parameters = dict()
    parameter_errors = dict()
    parameter_intervals = dict()
    parameter_ranges = dict()
    probabilities = dict()
    inflations = dict()
    if ("const" in model_parameters.index):
        #parameters["intercept_parameter"] = report.params[0]
        parameters["intercept_parameter"] = model_parameters["const"]
    else:
        parameters["intercept_parameter"] = float("nan")
        # Report.
        if report:
            utility.print_terminal_partition(level=4)
            print("Warning: regression data does not have constant intercept.")
            print(independence)
    if ("const" in model_parameter_errors.index):
        parameter_errors["intercept_error"] = model_parameter_errors["const"]
        parameter_intervals["intercept_interval_95"] = float(
            1.96 * parameter_errors["intercept_error"]
        )
        parameter_ranges["intercept_range_95"] = (
            determine_confidence_interval_range_text(
                estimate=parameters["intercept_parameter"],
                interval_low=parameter_intervals["intercept_interval_95"],
                interval_high=parameter_intervals["intercept_interval_95"],
        ))
    else:
        parameter_errors["intercept_error"] = float("nan")
        parameter_intervals["intercept_interval_95"] = float("nan")
        parameter_ranges["intercept_range_95"] = str("nan ... nan")
        # Report.
        if report:
            utility.print_terminal_partition(level=4)
            print("Warning: regression data does not have constant intercept.")
            print(independence)
    if ("const" in model_probabilities.index):
        #probabilities["intercept_probability"] = report.pvalues[0]
        probabilities["intercept_probability"] = (
            model_probabilities["const"]
        )
    else:
        probabilities["intercept_probability"] = float("nan")
        # Report.
        if report:
            utility.print_terminal_partition(level=4)
            print("Warning: regression data does not have constant intercept.")
            print(independence)
    inflations["intercept_inflation"] = float("nan")
    # Iterate on each independent variable.
    # Initiate counter at 1 to assume that intercept is at index 0.
    counter = 1
    # Accommodate index for intercept.
    for variable in independence:
        # Coefficient or parameter.
        parameter = str(variable + ("_parameter"))
        #parameters[parameter] = report.params[counter]
        parameters[parameter] = model_parameters[variable]
        # Parameter standard error
        parameter_error = str(variable + ("_error"))
        parameter_errors[parameter_error] = model_parameter_errors[variable]
        parameter_interval = str(variable + ("_interval_95"))
        parameter_intervals[parameter_interval] = float(
            1.96 * parameter_errors[parameter_error]
        )
        parameter_range = str(variable + ("_range_95"))
        parameter_ranges[parameter_range] = (
            determine_confidence_interval_range_text(
                estimate=parameters[parameter],
                interval_low=parameter_intervals[parameter_interval],
                interval_high=parameter_intervals[parameter_interval],
        ))
        # Probability.
        probability = str(variable + ("_probability"))
        #probabilities[probability] = report.pvalues[counter]
        probabilities[probability] = model_probabilities[variable]
        # Variance Inflation Factor (VIF).
        inflation = str(variable + ("_inflation"))
        inflation_value = (
            statsmodels.stats.outliers_influence.variance_inflation_factor(
                table_independence_intercept.to_numpy(),
                counter
            )
        )
        inflations[inflation] = round(inflation_value, 3)
        # Increment index.
        counter += 1
        pass
    summary = {
        "dependence_actual": str(dependence),
        "independence_actual": ";".join(independence),
        "freedom": pail_raw.df_model,
        "observations": pail_raw.nobs,
        "samples": count_samples,
        "r_square_pseudo": pail_raw.prsquared,
        "log_likelihood": pail_raw.llf,
        "akaike": pail_raw.aic,
        "bayes": pail_raw.bic,
    }
    summary.update(parameters)
    summary.update(parameter_errors)
    summary.update(parameter_intervals)
    summary.update(parameter_ranges)
    summary.update(probabilities)
    summary.update(inflations)

    # Compile information.
    pail = dict()
    pail["summary"] = summary
    pail["residuals"] = residuals
    # Return information.
    return pail


def regress_linear_ordinary_least_squares(
    dependence=None,
    independence=None,
    table=None,
    report=None,
):
    """
    Regresses a quantitative continuous dependent variable against multiple
    independent variables and returns relevant parameters and statistics.

    Table format must have samples (cases, observations) across rows and
    dependent and independent variables (features) across columns.

    Description of formats for StatsModels...

    Format of dependent variable is a vector of scalar values.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4]

    Format of independent variable(s) is a matrix: a first-dimension vector of
    samples (observations) and for each sample a second-dimension vector of
    variables' (features') scalar values.
    StatsModels also requires a constant for the intercept.
    [
        [1.3, 5.2, 1.0],
        [1.5, 5.1, 1.0],
        [1.2, 5.5, 1.0],
        ...
    ]

    arguments:
        dependence (str): name of table's column for dependent variable
        independence (list<str>): names of table's columns for independent
            variables
        table (object): Pandas data frame of dependent and independent variables
            for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of regression's residuals and statistics
    """

    # Determine count of valid samples (cases, observations).
    count_samples = int(table.shape[0])
    # Extract values of dependent and independent variables.
    # Can pass values of dependent and independent variables as NumPy arrays or
    # as Pandas Series and Dataframe respectively.
    # Passing variables as Pandas Series and Dataframe preserves variable names.
    values_dependence = table[dependence].to_numpy()

    # Keep independent variables in Pandas dataframe to preserve variables'
    # names.
    #values_independence = data.loc[ :, independence].to_numpy()
    table_independence = table.loc[ :, independence]
    # Introduce constant value for intercept.
    # If any column in the independent variables already has constant
    # values, then the function skips it by default.
    # It is necessary to change parameter "has_constant" to avoid this
    # conditional behavior.
    table_independence_intercept = statsmodels.api.add_constant(
        table_independence,
        prepend=True, # insert intercept constant first
        has_constant="add", # introduce new intercept constant regardless
    )
    columns_independence = copy.deepcopy(
        table_independence_intercept.columns.to_list()
    )
    #matrix_independence = table.to_numpy()
    # Define model.
    model = statsmodels.api.OLS(
        table[dependence],
        table_independence_intercept,
        missing="drop",
    )
    pail_raw = model.fit()
    # Report.
    if report:
        print("--------------------------------------------------")
        print(
            "Report source: " +
            "regress_dependent_independent_variables_linear_ordinary()"
        )
        print("--------------------------------------------------")
        print("Version check: TCW 28 September 2021")
        print("Information from regression:")
        print(pail_raw.summary())
        #utility.print_terminal_partition(level=3)
        #print(dir(pail_raw))
        #print(pail_raw.params)
        #print(pail_raw.pvalues)
        pass

    # Organize residuals.
    residuals = pail_raw.resid

    ##########
    # Collect parameters, errors, probabilities, and statistics.
    model_parameters = pandas.Series(data=pail_raw.params)
    model_parameter_errors = pandas.Series(data=pail_raw.bse)
    model_probabilities = pandas.Series(data=pail_raw.pvalues)
    parameters = dict()
    parameter_errors = dict()
    parameter_intervals = dict()
    parameter_ranges = dict()
    probabilities = dict()
    inflations = dict()
    if ("const" in model_parameters.index):
        #parameters["intercept_parameter"] = report.params[0]
        parameters["intercept_parameter"] = model_parameters["const"]
    else:
        parameters["intercept_parameter"] = float("nan")
        # Report.
        if report:
            utility.print_terminal_partition(level=4)
            print("Warning: regression data does not have constant intercept.")
            print(independence)
    if ("const" in model_parameter_errors.index):
        parameter_errors["intercept_error"] = model_parameter_errors["const"]
        parameter_intervals["intercept_interval_95"] = float(
            1.96 * parameter_errors["intercept_error"]
        )
        parameter_ranges["intercept_range_95"] = (
            determine_confidence_interval_range_text(
                estimate=parameters["intercept_parameter"],
                interval_low=parameter_intervals["intercept_interval_95"],
                interval_high=parameter_intervals["intercept_interval_95"],
        ))
    else:
        parameter_errors["intercept_error"] = float("nan")
        parameter_intervals["intercept_interval_95"] = float("nan")
        parameter_ranges["intercept_range_95"] = str("nan ... nan")
        # Report.
        if report:
            utility.print_terminal_partition(level=4)
            print("Warning: regression data does not have constant intercept.")
            print(independence)
    if ("const" in model_probabilities.index):
        #probabilities["intercept_probability"] = report.pvalues[0]
        probabilities["intercept_probability"] = (
            model_probabilities["const"]
        )
    else:
        probabilities["intercept_probability"] = float("nan")
        # Report.
        if report:
            utility.print_terminal_partition(level=4)
            print("Warning: regression data does not have constant intercept.")
            print(independence)
    inflations["intercept_inflation"] = float("nan")
    # Iterate on each independent variable.
    # Initiate counter at 1 to assume that intercept is at index 0.
    counter = 1
    # Accommodate index for intercept.
    for variable in independence:
        # Coefficient or parameter.
        parameter = str(variable + ("_parameter"))
        #parameters[parameter] = report.params[counter]
        parameters[parameter] = model_parameters[variable]
        # Parameter standard error
        parameter_error = str(variable + ("_error"))
        parameter_errors[parameter_error] = model_parameter_errors[variable]
        parameter_interval = str(variable + ("_interval_95"))
        parameter_intervals[parameter_interval] = float(
            1.96 * parameter_errors[parameter_error]
        )
        parameter_range = str(variable + ("_range_95"))
        parameter_ranges[parameter_range] = (
            determine_confidence_interval_range_text(
                estimate=parameters[parameter],
                interval_low=parameter_intervals[parameter_interval],
                interval_high=parameter_intervals[parameter_interval],
        ))
        # Probability.
        probability = str(variable + ("_probability"))
        #probabilities[probability] = report.pvalues[counter]
        probabilities[probability] = model_probabilities[variable]
        # Variance Inflation Factor (VIF).
        inflation = str(variable + ("_inflation"))
        inflation_value = (
            statsmodels.stats.outliers_influence.variance_inflation_factor(
                table_independence_intercept.to_numpy(),
                counter
            )
        )
        inflations[inflation] = round(inflation_value, 3)
        # Increment index.
        counter += 1
        pass
    summary = {
        "dependence_actual": str(dependence),
        "independence_actual": ";".join(independence),
        "freedom": pail_raw.df_model,
        "observations": pail_raw.nobs,
        "samples": count_samples,
        "r_square": pail_raw.rsquared,
        "r_square_adjust": pail_raw.rsquared_adj,
        "log_likelihood": pail_raw.llf,
        "akaike": pail_raw.aic,
        "bayes": pail_raw.bic,
        "condition": pail_raw.condition_number,
    }
    summary.update(parameters)
    summary.update(parameter_errors)
    summary.update(parameter_intervals)
    summary.update(parameter_ranges)
    summary.update(probabilities)
    summary.update(inflations)

    # Compile information.
    pail = dict()
    pail["summary"] = summary
    pail["residuals"] = residuals
    # Return information.
    return pail


def regress_tree_discrete_logit(
    dependence=None,
    independence=None,
    table=None,
    report=None,
):
    """
    Regresses a binary dependent variable against multiple independent variables
    and returns relevant parameters and statistics.

    Table format must have samples (cases, observations) across rows and
    dependent and independent variables (features) across columns.

    Description of formats for StatsModels...

    Format of dependent variable is a vector of scalar values.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4]

    Format of independent variable(s) is a matrix: a first-dimension vector of
    samples (observations) and for each sample a second-dimension vector of
    variables' (features') scalar values.
    StatsModels also requires a constant for the intercept.
    [
        [1.3, 5.2, 1.0],
        [1.5, 5.1, 1.0],
        [1.2, 5.5, 1.0],
        ...
    ]

    arguments:
        dependence (str): name of table's column for dependent variable
        independence (list<str>): names of table's columns for independent
            variables
        table (object): Pandas data frame of dependent and independent variables
            for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of regression's residuals and statistics
    """

    # Copy information.
    table = table.copy(deep=True)
    # Determine count of valid samples (cases, observations).
    count_samples = int(table.shape[0])
    # Extract values of dependent and independent variables.
    # Can pass values of dependent and independent variables as NumPy arrays or
    # as Pandas Series and Dataframe respectively.
    # Passing variables as Pandas Series and Dataframe preserves variable names.
    values_dependence = table[dependence].to_numpy()

    # Keep independent variables in Pandas dataframe to preserve variables'
    # names.
    #values_independence = data.loc[ :, independence].to_numpy()
    table_independence = table.loc[ :, independence]
    # Introduce constant value for intercept.
    # If any column in the independent variables already has constant
    # values, then the function skips it by default.
    # It is necessary to change parameter "has_constant" to avoid this
    # conditional behavior.
    table_independence_intercept = statsmodels.api.add_constant(
        table_independence,
        prepend=True, # insert intercept constant first
        has_constant="add", # introduce new intercept constant regardless
    )
    columns_independence = copy.deepcopy(
        table_independence_intercept.columns.to_list()
    )
    #matrix_independence = table.to_numpy()
    # Define model.
    # statsmodels.discrete.discrete_model.Logit
    model = statsmodels.api.Logit(
        table[dependence],
        table_independence_intercept,
        missing="drop",
    )
    pail_raw = model.fit()
    # Report.
    if report:
        print("--------------------------------------------------")
        print(
            "Report source: " +
            "regress_dependent_independent_variables_linear_ordinary()"
        )
        print("--------------------------------------------------")
        print("Version check: TCW 28 September 2021")
        print("Information from regression:")
        print(pail_raw.summary())
        #utility.print_terminal_partition(level=3)
        #print(dir(pail_raw))
        #print(pail_raw.params)
        #print(pail_raw.pvalues)
        pass

    ##########
    # Organize residuals.
    residuals = pail_raw.resid_generalized

    ##########
    # Collect parameters, errors, probabilities, and statistics for whole model.
    summary = {
        "dependence_actual": str(dependence),
        "independence_actual": ";".join(independence),
        "freedom": pail_raw.df_model,
        "observations": pail_raw.nobs,
        "samples": count_samples,
        "r_square_pseudo": pail_raw.prsquared,
        "log_likelihood": pail_raw.llf,
        "akaike": pail_raw.aic,
        "bayes": pail_raw.bic,
    }

    ##########
    # Collect parameters, errors, probabilities, and statistics for independent
    # variables.
    model_parameters = pandas.Series(data=pail_raw.params)
    model_parameter_errors = pandas.Series(data=pail_raw.bse)
    model_probabilities = pandas.Series(data=pail_raw.pvalues)
    pail_tree = organize_linear_logistic_regression_independence_tree(
        independence=independence,
        model_parameters=model_parameters,
        model_parameter_errors=model_parameter_errors,
        model_probabilities=model_probabilities,
        table_independence_intercept=table_independence_intercept,
    )
    # Collect information.
    summary["independence_tree"] = pail_tree
    pail = dict()
    pail["summary"] = summary
    pail["residuals"] = residuals
    # Return information.
    return pail


def regress_tree_linear_ordinary_least_squares(
    dependence=None,
    independence=None,
    table=None,
    report=None,
):
    """
    Regresses a quantitative continuous dependent variable against multiple
    independent variables and returns relevant parameters and statistics.

    Table format must have samples (cases, observations) across rows and
    dependent and independent variables (features) across columns.

    Description of formats for StatsModels...

    Format of dependent variable is a vector of scalar values.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4]

    Format of independent variable(s) is a matrix: a first-dimension vector of
    samples (observations) and for each sample a second-dimension vector of
    variables' (features') scalar values.
    StatsModels also requires a constant for the intercept.
    [
        [1.3, 5.2, 1.0],
        [1.5, 5.1, 1.0],
        [1.2, 5.5, 1.0],
        ...
    ]

    arguments:
        dependence (str): name of table's column for dependent variable
        independence (list<str>): names of table's columns for independent
            variables
        table (object): Pandas data frame of dependent and independent variables
            for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of regression's residuals and statistics
    """

    # Copy information.
    table = table.copy(deep=True)
    # Determine count of valid samples (cases, observations).
    count_samples = int(table.shape[0])
    # Extract values of dependent and independent variables.
    # Can pass values of dependent and independent variables as NumPy arrays or
    # as Pandas Series and Dataframe respectively.
    # Passing variables as Pandas Series and Dataframe preserves variable names.
    values_dependence = table[dependence].to_numpy()

    # Keep independent variables in Pandas dataframe to preserve variables'
    # names.
    #values_independence = data.loc[ :, independence].to_numpy()
    table_independence = table.loc[ :, independence]
    # Introduce constant value for intercept.
    # If any column in the independent variables already has constant
    # values, then the function skips it by default.
    # It is necessary to change parameter "has_constant" to avoid this
    # conditional behavior.
    table_independence_intercept = statsmodels.api.add_constant(
        table_independence,
        prepend=True, # insert intercept constant first
        has_constant="add", # introduce new intercept constant regardless
    )
    columns_independence = copy.deepcopy(
        table_independence_intercept.columns.to_list()
    )
    #matrix_independence = table.to_numpy()
    # Define model.
    model = statsmodels.api.OLS(
        table[dependence],
        table_independence_intercept,
        missing="drop",
    )
    pail_raw = model.fit()
    # Report.
    if report:
        print("--------------------------------------------------")
        print(
            "Report source: " +
            "regress_dependent_independent_variables_linear_ordinary()"
        )
        print("--------------------------------------------------")
        print("Version check: TCW 28 September 2021")
        print("Information from regression:")
        print(pail_raw.summary())
        #utility.print_terminal_partition(level=3)
        #print(dir(pail_raw))
        #print(pail_raw.params)
        #print(pail_raw.pvalues)
        pass

    ##########
    # Organize residuals.
    residuals = pail_raw.resid

    ##########
    # Collect parameters, errors, probabilities, and statistics for whole model.
    summary = {
        "dependence_actual": str(dependence),
        "independence_actual": ";".join(independence),
        "freedom": pail_raw.df_model,
        "observations": pail_raw.nobs,
        "samples": count_samples,
        "r_square": pail_raw.rsquared,
        "r_square_adjust": pail_raw.rsquared_adj,
        "log_likelihood": pail_raw.llf,
        "akaike": pail_raw.aic,
        "bayes": pail_raw.bic,
        "condition": pail_raw.condition_number,
    }

    ##########
    # Collect parameters, errors, probabilities, and statistics for independent
    # variables.
    model_parameters = pandas.Series(data=pail_raw.params)
    model_parameter_errors = pandas.Series(data=pail_raw.bse)
    model_probabilities = pandas.Series(data=pail_raw.pvalues)
    pail_tree = organize_linear_logistic_regression_independence_tree(
        independence=independence,
        model_parameters=model_parameters,
        model_parameter_errors=model_parameter_errors,
        model_probabilities=model_probabilities,
        table_independence_intercept=table_independence_intercept,
    )
    # Collect information.
    summary["independence_tree"] = pail_tree
    pail = dict()
    pail["summary"] = summary
    pail["residuals"] = residuals

    # TEMPORARY!!!

    print(pail["summary"])
    # Return information.
    return pail


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
    pail["parameter"] = float("nan")
    pail["error"] = float("nan")
    pail["interval_95"] = float("nan")
    pail["range_95"] = str("nan ... nan")
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


# Drivers


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
    else:
        # There is insufficient information for regression.
        pail_organization["validity"] = False
        # Report.
        if report:
            print("Insufficient information for regression.")
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
        threshold_samples=50,
        table=table,
        type=type,
        report=report,
    )
    # Determine whether to execute the regression.
    if (pail_check["validity"]):
        # Adequate information for regression.
        # Execute regression analysis.
        if (type == "linear"):
            pail_regression = regress_tree_linear_ordinary_least_squares(
                dependence=dependence,
                independence=pail_check["independence"],
                table=pail_check["table"],
                report=report,
            )
        elif (type == "logistic"):
            pail_regression = regress_tree_discrete_logit(
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


# TODO: TCW, 01 April 2022
# give this function more information... for each individual regression, need to know...
# name of cohort, dependent variable, model identifier, model note, list of all independent variables
# iterate on the independent variables and make long-format table

# TODO: TCW, 01 April 2022
# TODO: I need to reorganize the summary report tables
# TODO: 1. LONG format 2. WIDE format

# Create both a "long" table and a "wide" table for versatility... same info... different formats


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
        # Collection information.
        record = dict()
        # Extract names of entries for statistics on the whole model.
        entries_general = list(filter(
            lambda entry: ("independence_tree" not in str(entry)),
            record_regression.keys()
        ))
        # Extract names of independent variables to include in summary.
        # If 'independences_summary' is not 'None', then only include
        # information in summary table about independent variables that are in
        # the list 'independences_summary'.
        if (
            (independences_summary is not None) and
            (len(independences_summary) > 0)
        ):
            independences_inclusion = list(filter(
                lambda variable: (str(variable) in independences_summary),
                record_regression["independence_tree"].keys()
            ))
        else:
            independences_inclusion = (
                record_regression["independence_tree"].keys()
            )
        print("check: independence tree keys...")
        print(record_regression["independence_tree"].keys())
        # Iterate on independent variables.
        for variable in independences_inclusion:
            # Collect entries for statistics on the whole model.
            for entry in entries_general:
                record[entry] = record_regression[entry]
                pass
            # Collect name of variable.
            record["variable_key"] = variable
            # Extract names of entries for statistics on the independent
            # variable.
            entries_variable = (
                record_regression["independence_tree"][variable].keys()
            )
            # Collect entries for statistics on the independent variable.
            for entry in entries_variable:
                record[entry] = (
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
    if (type == "linear"):
        columns = [
            "cohort",
            #"dependence", "dependence_type",
            #"model",
            "model_note",
            "variable", #"variable_key",
            "parameter", #"error", "interval_95", "range_95",
            #"probability", "inflation",
            #"freedom", "observations", "samples",
            #"r_square", "r_square_adjust", "log_likelihood", "akaike", "bayes",
            #"condition",
            #"independence",
            #"dependence_actual", "independence_actual",
        ]
    elif (type == "logistic"):
        columns = [
            "cohort",
            "dependence", "dependence_type",
            "model", "model_note",
            "variable", "variable_key",
            "parameter", "error", "interval_95", "range_95",
            "probability", "inflation",
            "freedom", "observations", "samples",
            "r_square_pseudo", "log_likelihood", "akaike", "bayes", "condition",
            "independence",
            "dependence_actual", "independence_actual",
        ]
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
    "dependence_sort", "dependence_type", "model", "model_sort", "independence",
    "model_note"
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





#
