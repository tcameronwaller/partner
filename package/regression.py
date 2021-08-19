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

# Custom

import promiscuity.utility as utility # this import path for subpackage


#dir()
#importlib.reload()

###############################################################################
# Functionality




# TODO: TCW 19 August 2021
# TODO: I need to filter and standardize the table columns... BEFORE

def regress_linear_ordinary_least_squares(
    dependence=None,
    independence=None,
    threshold_samples=None,
    table=None,
    report=None,
):
    """
    Regresses a quantitative continuous dependent variable against multiple
    independent variables and returns relevant parameters and statistics.

    Table format must have samples (observations) across rows and dependent and
    independent variables (features) across columns.

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
        threshold_samples (float): minimal count of samples with non-missing
            values of dependent and independent variables to perform regression
        table (object): Pandas data frame of dependent and independent variables
            for regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of regression's residuals and statistics
    """

    # Copy information.
    table = table.copy(deep=True)
    # Remove observations with any missing values.
    columns = copy.deepcopy(independence)
    columns.insert(0, dependence)
    table = table.loc[:, table.columns.isin(columns)]
    table.dropna(
        axis="index",
        how="any",
        subset=None,
        inplace=True,
    )
    table = table[[*columns]]
    # Determine count of valid samples (observations).
    count_samples = table.shape[0]

    # Note
    # It is very important to keep track of the order of independent variables
    # in order to match parameters and probabilities to the correct variables.

    # Determine whether data have sufficient observations for regression.
    if count_samples >= threshold_samples:
        # Extract values of dependent and independent variables.
        values_dependence = table[dependence].to_numpy()
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
            has_constant="add", # Introduce new intercept constant regardless
        )
        columns_independence = copy.deepcopy(
            table_independence_intercept.columns.to_list()
        )
        #matrix_independence = table.to_numpy()
        # Define model.
        model = statsmodels.api.OLS(
            values_dependence,
            table_independence_intercept,
            missing="drop",
        )
        pail_raw = model.fit()
        # Report.
        if report:
            utility.print_terminal_partition(level=2)
            print(
                "Report source: " +
                "regress_dependent_independent_variables_linear_ordinary()"
            )
            utility.print_terminal_partition(level=3)
            print("Information from regression:")
            print(pail_raw.summary())
            #utility.print_terminal_partition(level=3)
            #print(dir(pail_raw))
            #print(pail_raw.params)
            #print(pail_raw.pvalues)
            pass

        # Organize residuals.
        residuals = pail_raw.resid
        # Collect parameters, probabilities, and statistics.
        model_parameters = pandas.Series(data=pail_raw.params)
        model_probabilities = pandas.Series(data=pail_raw.pvalues)
        parameters = dict()
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
            inflations[inflation] = round(inflation_value, 2)
            # Increment index.
            counter += 1
            pass
        summary = {
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
        summary.update(probabilities)
        summary.update(inflations)
    else:
        # Compile information.
        #probabilities = list(
        #    itertools.repeat(float("nan"), len(values_independence_intercept))
        #)
        parameters = dict()
        probabilities = dict()
        inflations = dict()
        parameters["intercept_parameter"] = float("nan")
        probabilities["intercept_probability"] = float("nan")
        inflations["intercept_inflation"] = float("nan")
        for variable in independence:
            parameter = str(variable + ("_parameter"))
            parameters[parameter] = float("nan")
            probability = str(variable + ("_probability"))
            probabilities[probability] = float("nan")
            inflation = str(variable + ("_inflation"))
            inflations[inflation] = float("nan")
            pass
        summary = {
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
        summary.update(probabilities)
        summary.update(inflations)
        residuals = numpy.empty(0)
    # Compile information.
    pail = dict()
    pail["summary"] = summary
    pail["residuals"] = residuals
    # Return information.
    return pail




#
