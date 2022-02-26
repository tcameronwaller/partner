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
        standard_scale (bool): whether to transform all dependent and
            independent variables to z-score standard scale
        table (object): Pandas data frame of dependent and independent variables
            for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information for regression

    """

    # Copy information.
    table = table.copy(deep=True)
    # Determine whether the variables exist.
    independence_source = list(filter(
        lambda column: (str(column) in table.columns.to_list()),
        independence
    ))
    # Select table's columns for relevant variables.
    columns = copy.deepcopy(independence_source)
    if (dependence in table.columns.to_list()):
        columns.insert(0, dependence)
    table = table.loc[:, table.columns.isin(columns)]
    table = table[[*columns]]

    # Remove columns with inadequate non-missing values or inadequate variance
    # (standard deviation) across rows.
    # This preliminary filter on columns avoids losing table rows for missing
    # values in an unnecessary column.
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
        "independence": ";".join(independence),
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
        "independence": ";".join(independence),
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

    # Create missing values for regression.
    #probabilities = list(
    #    itertools.repeat(float("nan"), len(values_independence_intercept))
    #)
    parameters = dict()
    parameter_errors = dict()
    parameter_intervals = dict()
    parameter_ranges = dict()
    probabilities = dict()
    inflations = dict()
    parameters["intercept_parameter"] = float("nan")
    parameter_errors["intercept_error"] = float("nan")
    parameter_intervals["intercept_interval_95"] = float("nan")
    parameter_ranges["intercept_range_95"] = str("nan ... nan")
    probabilities["intercept_probability"] = float("nan")
    inflations["intercept_inflation"] = float("nan")
    for variable in independence:
        parameter = str(variable + ("_parameter"))
        parameters[parameter] = float("nan")
        parameter_error = str(variable + ("_error"))
        parameter_errors[parameter_error] = float("nan")
        parameter_interval = str(variable + ("_interval_95"))
        parameter_intervals[parameter_interval] = float("nan")
        parameter_range = str(variable + ("_range_95"))
        parameter_ranges[parameter_range] = str("nan ... nan")
        probability = str(variable + ("_probability"))
        probabilities[probability] = float("nan")
        inflation = str(variable + ("_inflation"))
        inflations[inflation] = float("nan")
        pass
    summary = {
        "independence": "",
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
    summary.update(parameters)
    summary.update(parameter_errors)
    summary.update(parameter_intervals)
    summary.update(parameter_ranges)
    summary.update(probabilities)
    summary.update(inflations)
    residuals = numpy.empty(0)
    # Compile information.
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

    # Create missing values for regression.
    #probabilities = list(
    #    itertools.repeat(float("nan"), len(values_independence_intercept))
    #)
    parameters = dict()
    parameter_errors = dict()
    parameter_intervals = dict()
    parameter_ranges = dict()
    probabilities = dict()
    inflations = dict()
    parameters["intercept_parameter"] = float("nan")
    parameter_errors["intercept_error"] = float("nan")
    parameter_intervals["intercept_interval_95"] = float("nan")
    parameter_ranges["intercept_range_95"] = str("nan ... nan")
    probabilities["intercept_probability"] = float("nan")
    inflations["intercept_inflation"] = float("nan")
    for variable in independence:
        parameter = str(variable + ("_parameter"))
        parameters[parameter] = float("nan")
        parameter_error = str(variable + ("_error"))
        parameter_errors[parameter_error] = float("nan")
        parameter_interval = str(variable + ("_interval_95"))
        parameter_intervals[parameter_interval] = float("nan")
        parameter_range = str(variable + ("_range_95"))
        parameter_ranges[parameter_range] = str("nan ... nan")
        probability = str(variable + ("_probability"))
        probabilities[probability] = float("nan")
        inflation = str(variable + ("_inflation"))
        inflations[inflation] = float("nan")
        pass
    summary = {
        "independence": "",
        "freedom": float("nan"),
        "observations": float("nan"),
        "samples": float("nan"),
        "r_square_pseudo": float("nan"),
        "log_likelihood": float("nan"),
        "akaike": float("nan"),
        "bayes": float("nan"),
    }
    summary.update(parameters)
    summary.update(parameter_errors)
    summary.update(parameter_intervals)
    summary.update(parameter_ranges)
    summary.update(probabilities)
    summary.update(inflations)
    residuals = numpy.empty(0)
    # Compile information.
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


def drive_organize_table_regress_linear_logistic(
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
        standard_scale (bool): whether to transform all dependent and
            independent variables to z-score standard scale
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
        utility.print_terminal_partition(level=2)
        print("report: ")
        function_name = str(
            "drive_organize_table_regress_linear_logistic()"
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
    # Determine and report whether the table does not include information for
    # all requested independent variables.
    if (
        report and
        (len(pail_organization["independence"]) < len(independence))
    ):
            utility.print_terminal_partition(level=5)
            print(
                "Warning: Some requested independent variables are " +
                "uavailable in the table for regression."
            )
            print("requested independent variables:")
            print(independence)
            print("available independent variables:")
            print(pail_organization["independence"])
            pass
    # Determine whether the table has adequate information for regression.
    if (
        (dependence in pail_organization["table"].columns.to_list()) and
        (len(pail_organization["independence"]) > 0) and
        (pail_organization["count_samples"] >= threshold_samples)
    ):
        # Adequate information for regression.
        # Execute regression analysis.
        if (type == "linear"):
            pail_regression = regress_linear_ordinary_least_squares(
                dependence=dependence,
                independence=pail_organization["independence"],
                table=pail_organization["table"],
                report=report,
            )
        elif (type == "logistic"):
            pail_regression = regress_discrete_logit(
                dependence=dependence,
                independence=pail_organization["independence"],
                table=pail_organization["table"],
                report=report,
            )
        else:
            pail_regression = create_regression_missing_values(
                dependence=dependence,
                independence=independence,
                type=type,
            )
            pass
    else:
        # Inadequate information for regression.
        # Create missing values.
        pail_regression = create_regression_missing_values(
            dependence=dependence,
            independence=independence,
            type=type,
        )
        # Determine and report cause of error.
        if (
            report and
            (dependence not in pail_organization["table"].columns.to_list())
        ):
            utility.print_terminal_partition(level=5)
            print("Error: Inadequate information for regression.")
            print(
                "Dependent variable is unavailable in the table."
            )
            pass
        if (
            report and
            (len(pail_organization["independence"]) == 0)
        ):
            utility.print_terminal_partition(level=5)
            print("Error: Inadequate information for regression.")
            print(
                "There are not any independent variables available " +
                "in the table."
            )
            pass
        if (
            report and
            (not (pail_organization["count_samples"] >= threshold_samples))
        ):
            utility.print_terminal_partition(level=5)
            print("Error: Inadequate information for regression.")
            print(
                "There are inadequate samples in the table."
            )
            print(pail_organization["count_samples"])
            pass
        pass
    # Return information.
    return pail_regression


# TODO: TCW, 24 February 2022
# TODO: This function is unnecessary.
# TODO: instead pass the independent variables directly from the table row record on top-level iteration...

def determine_cohort_model_variables_from_reference_table(
    cohort=None,
    model=None,
    dependence=None,
    table=None,
    report=None,
):
    """
    Determine the independent variables for a regression model that match a
    cohort and dependent variable in a reference table.

    Table format must have columns "cohort", "dependence", and "independence".
    Rows in Column "dependence" have actual names of variable columns in the
    main source table (not visible to this function).
    Rows in Column "independence" have list of independent variables with
    semicolon (";") delimiters.

    arguments:
        cohort (str): name of a stratification cohort for regression analysis
        model (str): name of a model for regression analysis, normally
            "complex", "simple", or "unadjust"
        dependence (str): name of table's column for dependent variable
        table (object): Pandas data frame of cohorts, models, dependent
            variables, and independent variables for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of names of cohort, dependent variable, and
            independent variables for regression
    """

    # Copy information.
    table = table.copy(deep=True)
    # Select table records for cohort.
    table_cohort = table.loc[
        table["cohort"] == str(cohort), :
    ]
    # Select table records for dependent variable.
    table_dependence = table_cohort.loc[
        table_cohort["dependence"] == str(dependence), :
    ]
    # Select table records for model.
    match_model = (str(model) in table_dependence["model"].to_list())
    if (match_model):
        table_model = table_dependence.loc[
            table_dependence["model"] == str(model), :
        ]
        # Determine whether cohort-model-dependence defines a single set of
        # independent variables.
        if (int(table_model.shape[0]) == 1):
            independence_string = copy.deepcopy(
                table_model.iloc[0].at["independence"]
            )
        elif (int(table_dependence.shape[0]) > 1):
            print(
                "WARNING: multiple sets of independent variables for given " +
                "cohort, model, and dependent variable."
            )
            print(
                "Returning first set of independent variables."
            )
            independence_string = copy.deepcopy(
                table_model["independence"].to_list()[0]
            )
    else:
        print(
            "WARNING: no information for given " +
            "cohort, model, and dependent variable."
        )
        print(
            "Returning empty set of independent variables."
        )
        independence_string = ""

    # Extract names of independent variables.
    independence_split = independence_string.split(";")
    independence_strip = list(map(
        lambda value: value.strip(), independence_split
    ))
    # Collect information.
    pail = dict()
    pail["cohort"] = cohort
    pail["model"] = model
    pail["match"] = match_model
    pail["dependence"] = dependence
    pail["independence"] = independence_strip
    # Return information.
    return pail


def drive_cohort_model_linear_logistic_regression(
    table=None,
    table_cohorts_models=None,
    cohort=None,
    model=None,
    dependence=None,
    type=None,
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
        table_cohorts_models (object): Pandas data frame of cohorts, models,
            dependent variables, and independent variables for regression
        cohort (str): name of a stratification cohort for regression analysis
        model (str): name of a model for regression analysis, normally
            "complex", "simple", or "unadjust"
        dependence (str): name of table's column for dependent variable
        type (str): type of regression analysis, either 'linear' or 'logistic'
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regressions

    """

    # Extract details for the regression analysis.
    pail_model = determine_cohort_model_variables_from_reference_table(
        cohort=cohort,
        model=model,
        dependence=dependence,
        table=table_cohorts_models,
        report=report,
    )
    if (pail_model["match"]):
        # Report.
        if report:
            utility.print_terminal_partition(level=2)
            print("report: ")
            print("drive_cohort_model_linear_regression()")
            utility.print_terminal_partition(level=5)
            print("cohort: " + str(cohort))
            print("model: " + str(model))
            print("dependent variable: " + str(dependence))
            print("independent variables: ")
            print(pail_model["independence"])
            utility.print_terminal_partition(level=5)
            pass
        # Determine type of regression analysis.
        pail_regression = drive_organize_table_regress_linear_logistic(
            dependence=dependence,
            independence=pail_model["independence"],
            standard_scale=True,
            threshold_samples=50,
            table=table,
            type=type,
            report=report,
        )
    else:
        # Report.
        if report:
            utility.print_terminal_partition(level=2)
            print("report: ")
            print("drive_cohort_model_linear_regression()")
            utility.print_terminal_partition(level=5)
            print("cohort: " + str(cohort))
            print("model: " + str(model))
            print("dependent variable: " + str(dependence))
            print("independent variables: ")
            print(pail_model["independence"])
            utility.print_terminal_partition(level=5)
            print("Missing information for model...")
            utility.print_terminal_partition(level=5)
        pail_regression = create_regression_missing_values(
            dependence=dependence,
            independence=pail_model["independence"],
            type=type,
        )
    return pail_regression


def organize_table_regression_summaries(
    independence=None,
    table=None,
    type=type,
    report=None,
):
    """
    This function organizes the summary table with information for multiple
    regressions.

    A primary function is to filter the table.

    arguments:
        independence (list<str>): name of independent variables of interest
        table (object): Pandas data frame of summary information from multiple
            regressions
        type (str): type of regression analysis, either 'linear' or 'logistic'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of summary information from multiple
            regressions

    """

    # Copy information.
    table = table.copy(deep=True)
    # Determine type of regression.
    if (type == "linear"):
        # Define columns of interest.
        columns = list()
        columns.append("name")
        columns.append("cohort")
        columns.append("dependence")
        columns.append("model")
        columns.append("independence")
        columns.append("freedom")
        columns.append("observations")
        columns.append("samples")
        columns.append("r_square")
        columns.append("r_square_adjust")
        columns.append("log_likelihood")
        columns.append("akaike")
        columns.append("bayes")
        columns.append("condition")
        for variable in independence:
            parameter = str(variable + ("_parameter"))
            parameter_error = str(variable + ("_error"))
            parameter_interval = str(variable + ("_interval_95"))
            parameter_range = str(variable + ("_range_95"))
            probability = str(variable + ("_probability"))
            inflation = str(variable + ("_inflation"))
            columns.append(parameter)
            columns.append(parameter_error)
            columns.append(parameter_interval)
            columns.append(parameter_range)
            columns.append(probability)
            columns.append(inflation)
            pass
    elif (type == "logistic"):
        # Define columns of interest.
        columns = list()
        columns.append("name")
        columns.append("cohort")
        columns.append("dependence")
        columns.append("model")
        columns.append("independence")
        columns.append("freedom")
        columns.append("observations")
        columns.append("samples")
        columns.append("r_square_pseudo")
        columns.append("log_likelihood")
        columns.append("akaike")
        columns.append("bayes")
        for variable in independence:
            parameter = str(variable + ("_parameter"))
            parameter_error = str(variable + ("_error"))
            parameter_interval = str(variable + ("_interval_95"))
            parameter_range = str(variable + ("_range_95"))
            probability = str(variable + ("_probability"))
            inflation = str(variable + ("_inflation"))
            columns.append(parameter)
            columns.append(parameter_error)
            columns.append(parameter_interval)
            columns.append(parameter_range)
            columns.append(probability)
            columns.append(inflation)
            pass
    else:
        columns = list()
        pass
    # Select columns.
    # columns.insert(0, dependence)
    table = table.loc[:, table.columns.isin(columns)]
    # Sort columns.
    table = table[[*columns]]
    # Return information.
    return table




#
