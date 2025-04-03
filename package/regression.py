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
import partner.utility as putly
import partner.organization as porg
import partner.scale as pscl
import partner.description as pdesc

#dir()
#importlib.reload()

###############################################################################
# Functionality

# Note: TCW; 4 October 2022
# There are many different methods to "standardize", "scale", "transform", or
# "standardize" the distributions of multiple variables to simplify their
# comparisons to each other.




# TODO: TCW; 2 April 2025
# I don't actually need to keep the "table_independence_intercept"
# keep this function and the regression tree records as simple as possible
# later on it'll be more convenient or cleaner to derive extra info for summary table

# flatten
# predictor_1_name: "intercept"
# predictor_1_parameter:

# or... keep the dictionary tree, but then extract using a translation
# 1: intercept
# 2: sex_y


##########
# Organize information from regresion.


def define_sequence_features(
    features=None,
):
    """
    Define specific sequence of features.

    Review: TCW; 2 April 2025

    arguments:
        features (list<str>): names of features for which to define specific
            sequence

    raises:

    returns:
        (dict<str>): collection of information

    """

    # Copy information.
    features = copy.deepcopy(features)

    # Option 1.
    if True:
        sequence_features = dict()
        index = 0
        for name in features:
            sequence_features[index] = name
            index += 1
            pass
    # Option 2.
    if False:
        sequence_features = dict(zip(
            range(len(features)),
            features
        ))
    # Option 3.
    if False:
        sequence_features = {
            key: value for key, value in enumerate(features)
        }

    # Return information.
    return sequence_features


def define_features_sequence(
    features=None,
):
    """
    Define specific sequence of features.

    Review: TCW; 2 April 2025

    arguments:
        features (list<str>): names of features for which to define specific
            sequence

    raises:

    returns:
        (dict<int>): collection of information

    """

    # Copy information.
    features = copy.deepcopy(features)

    # Option 1.
    if True:
        features_sequence = dict()
        index = 0
        for name in features:
            features_sequence[name] = index
            index += 1
            pass
    # Option 2.
    if False:
        features_sequence = dict(zip(
            features,
            range(len(features))
        ))
    # Option 3.
    if False:
        features_sequence = {
            key: value for value, key in enumerate(features)
        }

    # Return information.
    return features_sequence


def parse_regression_model_predictor(
    table_predictor_intercept=None,
    indicator=None,
    name_source=None,
    name_product=None,
    parameters=None,
    standard_errors=None,
    pvalues=None,
    variance_inflation=None,
):
    """
    Parse raw information about predictor independent variables from a model
    of a linear or logistic regression.

    Review: TCW; 2 April 2025

    arguments:
        table_predictor_intercept (object): Pandas data-frame table of data
            with features and observations for regression
        indicator (int): numeric sequence and predictor feature
        name_source (str): name of predictor feature
        name_product (str): name of predictor feature
        parameters (dict): parameters for predictors from regression model
        standard_errors (dict): standard errors for predictors from regression
            model
        pvalues (dict): pvalues for predictors from regression model
        variance_inflation (bool): whether to determine variance inflation
            factor for the predictor feature

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Determine code.
    code = str("predictor_" + str(indicator))
    # Collect information.
    pail = dict()
    pail["predictor"] = dict()
    pail["sequence"] = [
        str(code + "_name"),
        str(code + "_parameter"),
        str(code + "_error"),
        str(code + "_pvalue"),
        str(code + "_variance_inflation"),
    ]
    pail["predictor"][str(code + "_name")] = name_product
    pail["predictor"][str(code + "_parameter")] = float(
        parameters[name_source]
    )
    pail["predictor"][str(code + "_error")] = float(
        standard_errors[name_source]
    )
    pail["predictor"][str(code + "_pvalue")] = float(pvalues[name_source])
    if (variance_inflation):
        inflation_value = float(
            statsmodels.stats.outliers_influence.variance_inflation_factor(
                table_predictor_intercept.to_numpy(),
                indicator
            )
        )
        pail["predictor"][str(code + "_variance_inflation")] = float(
            inflation_value
        )
    else:
        pail["predictor"][str(code + "_variance_inflation")] = float("nan")
        pass

    # Return information.
    return pail


def parse_regression_model_predictors(
    table_predictor_intercept=None,
    features_predictor=None,
    features_sequence=None,
    name_intercept_source=None,
    name_intercept_product=None,
    parameters=None,
    standard_errors=None,
    pvalues=None,
):
    """
    Parse raw information about predictor independent variables from a model
    of a linear or logistic regression.

    Review: TCW; 2 April 2025

    arguments:
        table_predictor_intercept (object): Pandas data-frame table of data
            with features and observations for regression
        features_predictor (list<str>): names of features to include in
            regression model as predictor independent variables
        features_sequence (dict<int>): numeric sequence and names of features
        name_intercept_source (str): name of parameters for intercept
        name_intercept_product (str): name of parameters for intercept
        parameters (dict): parameters for predictors from regression model
        standard_errors (dict): standard errors for predictors from regression
            model
        pvalues (dict): pvalues for predictors from regression model

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Copy information.
    table_predictor_intercept = table_predictor_intercept.copy(deep=True)
    features_predictor = copy.deepcopy(features_predictor)
    features_sequence = copy.deepcopy(features_sequence)

    # Collect information.
    pail = dict()
    pail["predictors"] = dict()
    pail["sequence"] = list()

    # Intercept.
    # Invert dictionary for convenience.
    #sequence_features = {v : k for k, v in features_sequence.items()}
    indicator = features_sequence[name_intercept_product]
    if (
        (name_intercept_source in parameters.index) and
        (name_intercept_source in standard_errors.index) and
        (name_intercept_source in pvalues.index)
    ):
        # Collect information.
        pail_predictor = parse_regression_model_predictor(
            table_predictor_intercept=table_predictor_intercept,
            indicator=indicator,
            name_source=name_intercept_source,
            name_product=name_intercept_product,
            parameters=parameters,
            standard_errors=standard_errors,
            pvalues=pvalues,
            variance_inflation=False,
        )
        pail["predictors"].update(pail_predictor["predictor"])
        pail["sequence"].extend(pail_predictor["sequence"])
        pass

    # Main predictor features.
    for feature in features_predictor:
        indicator = features_sequence[feature]
        # Collect information.
        pail_predictor = parse_regression_model_predictor(
            table_predictor_intercept=table_predictor_intercept,
            indicator=indicator,
            name_source=feature,
            name_product=feature,
            parameters=parameters,
            standard_errors=standard_errors,
            pvalues=pvalues,
            variance_inflation=True,
        )
        pail["predictors"].update(pail_predictor["predictor"])
        pail["sequence"].extend(pail_predictor["sequence"])
        pass

    # Return information.
    return pail


# TODO: TCW; 2 April
# obsolete, I think, but check again and test the new method
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
        model_parameters (dict): estimate values of the model's parameter
            coefficients
        model_parameter_errors (dict): standard errors for estimates of the
            model's parameter coefficients
        model_probabilities (dict): probabilities for estimates of the model's
            parameter coefficients
        table_independence_intercept (object): Pandas data-frame table of
            independent variables with a constant intercept, the same source
            from the regression

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
        # Coefficient or parameter.
        #pail_tree["intercept"]["parameter"] = report.params[0]
        pail_tree["intercept"]["parameter"] = float(model_parameters["const"])
        # Standard error of parameter.
        pail_tree["intercept"]["error"] = float(model_parameter_errors["const"])
        # Confidence intervals and ranges.
        pail_confidence = pdesc.determine_95_99_confidence_intervals_ranges(
            estimate=pail_tree["intercept"]["parameter"],
            standard_error=pail_tree["intercept"]["error"],
        )
        pail_tree["intercept"].update(pail_confidence)
        # Probability.
        pail_tree["intercept"]["probability"] = float(
            model_probabilities["const"]
        )
        # Variance Inflation Factor (VIF).
        # Missing or undefined for intercept?
        pail_tree["intercept"]["inflation"] = float("nan")
        # Report summaries.
        pail_tree["intercept"]["report_b95ci"] = str(
            "b: " + str(round(pail_tree["intercept"]["parameter"], 5)) +
            "; 95% CI: " + str(pail_tree["intercept"]["range_95"])
        )
        pail_tree["intercept"]["report_bep"] = str(
            "b: " + str(round(pail_tree["intercept"]["parameter"], 5)) +
            " (" + str(round(pail_tree["intercept"]["error"], 5)) +
            "); p: " + str(round(pail_tree["intercept"]["probability"], 5))
        )
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
        pail_tree[variable]["parameter"] = float(model_parameters[variable])
        # Standard error parameter.
        pail_tree[variable]["error"] = float(model_parameter_errors[variable])
        # Confidence intervals and ranges.
        pail_confidence = pdesc.determine_95_99_confidence_intervals_ranges(
            estimate=pail_tree[variable]["parameter"],
            standard_error=pail_tree[variable]["error"],
        )
        pail_tree[variable].update(pail_confidence)
        # Probability.
        pail_tree[variable]["probability"] = float(
            model_probabilities[variable]
        )
        # Variance Inflation Factor (VIF).
        inflation_value = float(
            statsmodels.stats.outliers_influence.variance_inflation_factor(
                table_independence_intercept.to_numpy(),
                counter
            )
        )
        pail_tree[variable]["inflation"] = round(inflation_value, 5)
        # Report summaries.
        pail_tree[variable]["report_b95ci"] = str(
            "b: " + str(round(pail_tree[variable]["parameter"], 5)) +
            "; 95% CI: " + str(pail_tree[variable]["range_95"])
        )
        pail_tree[variable]["report_bep"] = str(
            "b: " + str(round(pail_tree[variable]["parameter"], 5)) +
            " (" + str(round(pail_tree[variable]["error"], 5)) +
            "); p: " + str(round(pail_tree[variable]["probability"], 5))
        )
        # Increment index.
        counter += 1
        pass
    # Return information.
    return pail_tree


def parse_organize_linear_logistic_predictors(
    table_predictor_intercept=None,
    features_predictor=None,
    parameters=None,
    standard_errors=None,
    pvalues=None,
    report=None,
):
    """
    Parse and organize information about predictor independent variables from a
    model of a linear or logistic regression.

    Review: TCW; 2 April 2025

    arguments:
        table_predictor_intercept (object): Pandas data-frame table of data
            with features and observations for regression
        features_predictor (list<str>): names of features to include in
            regression model as predictor independent variables
        parameters (dict): parameters for predictors from regression model
        standard_errors (dict): standard errors for predictors from regression
            model
        pvalues (dict): pvalues for predictors from regression model
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Copy information.
    table_predictor_intercept = table_predictor_intercept.copy(deep=True)
    features_predictor = copy.deepcopy(features_predictor)
    features_predictor_intercept = copy.deepcopy(features_predictor)

    # Determine count and sequence of predictors.
    features_predictor_intercept.insert(0, "intercept")
    count_predictors = len(features_predictor_intercept)
    features_sequence = define_features_sequence(
        features=features_predictor_intercept,
    )

    # Parse information for predictor features from regression model.
    pail_predictors = parse_regression_model_predictors(
        table_predictor_intercept=table_predictor_intercept,
        features_predictor=features_predictor,
        features_sequence=features_sequence,
        name_intercept_source="const",
        name_intercept_product="intercept",
        parameters=parameters,
        standard_errors=standard_errors,
        pvalues=pvalues,
    )
    pail_predictors["predictors"]["count_predictors"] = count_predictors
    pail_predictors["sequence"].insert(0, "count_predictors")

    print("!!!!!!!!!!!!!!!!!!!!!!")
    print(pail_predictors["predictors"])
    print(pail_predictors["sequence"]) # for filter and sort on table columns

    ##########
    # Collect parameters, errors, probabilities, and statistics for independent
    # variables.
    if False:
        pail_tree = organize_linear_logistic_regression_independence_tree(
            table_independence_intercept=table_predictor_intercept,
            independence=features_predictor_fixed,
            model_parameters=model_parameters,
            model_parameter_errors=model_parameter_errors,
            pvalues=pvalues,
        )


    # Return information.
    #return pail
    pass




##########
# Regression models







##########
# Manage regression analysis.



def regress_linear_ordinary_least_squares(
    table=None,
    feature_response=None,
    features_predictor_fixed=None,
    report=None,
):
    """
    Regress by linear ordinary least squares, applying implementation from
    the StatsModels package in Python.

    ----------
    Format of source data table (name: "table")
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

    Format of information for response dependent variable is a one-dimensional
    vector of scalar values on a quantitative, continuous scale of measurement,
    interval or ratio.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4]

    Format of information for predictor independent variable(s) is a
    two-dimensional matrix: a first-dimension vector corresponding to
    observations or samples and for each observation a second-dimension vector
    of scalar values corresponding to each feature or variable.
    StatsModels also requires the intentional addition of a constant for the
    intercept.
    [
        [1.3, 5.2, 1.0],
        [1.5, 5.1, 1.0],
        [1.2, 5.5, 1.0],
        ...
    ]

    implementation: 'statsmodels.regression.linear_model.OLS()'

    Review: TCW; 2 April 2025

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Extract values of response and predictor variables.
    # It is possible to pass values of response and predictor features either
    # as NumPy arrays or as Pandas Series and data-frame table respectively.
    # Passing variables as Pandas Series and Dataframe preserves variable names.
    #values_response = table[feature_response].to_numpy()

    # Keep information for predictor features in Pandas data-frame table to
    # preserve names of variables.
    #values_predictor = table.loc[ :, features_predictor_fixed].to_numpy()
    table_predictor = table.loc[ :, features_predictor_fixed].copy(deep=True)

    # Introduce constant intercept to predictor features.
    # If any column in the table of information for predictor features already
    # has constant values across observations, then the function skips it by
    # default. It is necessary to change parameter "has_constant" to avoid this
    # conditional behavior.
    table_predictor_intercept = statsmodels.api.add_constant(
        table_predictor,
        prepend=True, # insert intercept constant before other features
        has_constant="add", # force introduction of new intercept constant
    )
    columns_predictor = copy.deepcopy(
        table_predictor_intercept.columns.to_list()
    )
    #matrix_predictor = table.to_numpy()
    # Define model.
    model = statsmodels.api.OLS(
        table[feature_response],
        table_predictor_intercept,
        missing="drop",
    )
    # Fit the model.
    handle_model = model.fit()

    # Collect information.
    pail = dict()
    pail["model"] = handle_model
    pail["table_predictor_intercept"] = table_predictor_intercept

    # Return information.
    return pail


def manage_regression_linear_ordinary_least_squares(
    table=None,
    index_columns=None,
    index_rows=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    report=None,
):
    """
    Manage regression analysis according to parameters.

    ----------
    Format of source data table (name: "table")
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
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        formula_text (str): human readable formula for regression model,
            treated as a note for clarification
        feature_response (str): name of column in data table for feature
            variable to include in regression model as response dependent
            variable
        features_predictor_fixed (list<str>): names of columns in data table
            for feature variables to include in regression model as predictor
            independent variables with fixed effects
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Copy information.
    table = table.copy(deep=True)
    features_predictor_fixed = copy.deepcopy(features_predictor_fixed)

    # Regress.
    pail_regression = regress_linear_ordinary_least_squares(
        table=table,
        feature_response=feature_response,
        features_predictor_fixed=features_predictor_fixed,
        report=report,
    )

    # Organize information from regression for predictor features.
    parameters = pandas.Series(data=pail_regression["model"].params)
    standard_errors = pandas.Series(data=pail_regression["model"].bse)
    pvalues = pandas.Series(data=pail_regression["model"].pvalues)
    pail_predictors = parse_organize_linear_logistic_predictors(
        table_predictor_intercept=pail_regression["table_predictor_intercept"],
        features_predictor=features_predictor_fixed,
        parameters=parameters,
        standard_errors=standard_errors,
        pvalues=pvalues,
        report=report,
    )


    if False:

        # Collect information.
        pail = dict()
        pail["implementation"] = "statsmodels.regression.linear_model.OLS()"
        pail["count_observations_table"] = int(table.shape[0])
        pail["count_observations_model"] = pail_regression["model"].nobs
        pail["degrees_freedom_model"] = pail_regression["model"].df_model
        pail["degrees_freedom_residual"] = pail_regression["model"].df_resid
        pail["residuals"] = pail_regression["model"].resid
        pail["r_square"] = pail_regression["model"].rsquared
        pail["r_square_adjust"] = pail_regression["model"].rsquared_adj
        pail["log_likelihood"] = pail_regression["model"].llf
        pail["akaike"] = pail_regression["model"].aic
        pail["bayes"] = pail_regression["model"].bic
        pail["condition"] = pail_regression["model"].condition_number

        #pail.update(pail_predictors) <--- needs to be flat

        # Collect information.
        #pail["pail_tree"] = pail_tree
        pail["model"] = pail_regression["model"]

        # Report.
        if report:
            putly.print_terminal_partition(level=4)
            print("package: partner")
            print("module: regression.py")
            name_function = str(
                "perform_regression_analysis()"
            )
            print("function: " + name_function)
            putly.print_terminal_partition(level=5)
            print("text formula for regression model:")
            print(formula_text)
            putly.print_terminal_partition(level=5)
            print("summary from regression model:")
            print(pail_regression["model"].summary())
            #print(dir(handle_regression))
            #print(handle_regression.params)
            #print(handle_regression.pvalues)
            pass

        print("*********************")
        print(model_parameters)

        print("!!!!!!!!!!!!!!!!!!!!!!")
        print(pail_tree)
        # Return information.
        return pail
    pass


def determine_type_regression_analysis(
    table=None,
    index_columns=None,
    index_rows=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    report=None,
):
    """
    Determine the type of regression analysis according to parameters.

    ----------
    Format of source data table (name: "table")
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
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
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
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Determine type of regression.
    if (type_regression == "linear"):
        # Perform linear ordinary least squares regression.
        pail = manage_regression_linear_ordinary_least_squares(
            table=table,
            index_columns="features",
            index_rows="observations",
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            report=report,
        )
        pass
    elif (type_regression == "logistic"):
        # Perform logistic regression.
        pail = None
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "perform_regression_analysis()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return pail


##########
# Collect record of information for regression analysis.


def collect_record_regression_analysis(
    table=None,
    index_columns=None,
    index_rows=None,
    check_overall=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    record_extra=None,
    report=None,
):
    """
    Collect information within a record to describe a regression analysis.

    For an example on the application of this function in the context of
    extensive ancillary functionality to organize parameters and data and to
    drive multiple instances of regression concurrently in parallel, refer to
    the Python script 'script_drive_regressions_from_table_parameters.py'
    within the 'partner' repository. Within its subdirectory 'demonstration',
    the 'partner' repository also includes an example of the table of
    parameters that controls this Python script.

    ----------
    Format of source data table (name: "table")
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
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        check_overall (bool): whether to execute the regression or create
            missing values
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
        record_extra (dict): collection of secondary information to join to the
            record of primary information from the regression
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # Collect information.
    record = dict()
    record["check_overall"] = check_overall
    record["type_regression"] = type_regression
    record["formula_text"] = formula_text
    record["feature_response"] = feature_response
    record["features_predictor_fixed"] = features_predictor_fixed
    record["features_predictor_random"] = features_predictor_random
    record.update(record_extra)
    # Determine whether to perform regression.
    if (check_overall):
        # Perform regression.
        record["regression"] = determine_type_regression_analysis(
            table=table,
            index_columns="features",
            index_rows="observations",
            type_regression=type_regression,
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            features_predictor_random=features_predictor_random,
            report=report,
        )
        pass
    else:
        # Create a null record.
        record["regression"] = None
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "collect_record_regression_analysis()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return record




###############################################################################
# End
