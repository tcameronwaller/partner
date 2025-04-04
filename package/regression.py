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


##########
# Prepare table of parameters for many different response features, all with
# identical predictor features and other parameters.

# Refer to Bash script "template_create_parameter_table.sh".



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

    Review: TCW; 3 April 2025

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
    pail["entries"] = [
        str(code + "_name"),
        str(code + "_parameter"),
        str(code + "_error"),
        str(code + "_ci99_low"),
        str(code + "_ci99_high"),
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
    pail_confidence = pdesc.determine_95_99_confidence_intervals_ranges(
        estimate=parameters[name_source],
        standard_error=standard_errors[name_source],
    )
    pail["predictor"][str(code + "_ci99_low")] = float(
        pail_confidence["range_99_low"]
    )
    pail["predictor"][str(code + "_ci99_high")] = float(
        pail_confidence["range_99_high"]
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


def parse_regression_intercept_predictors(
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

    Review: TCW; 3 April 2025

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
    parameters = copy.deepcopy(parameters)
    standard_errors = copy.deepcopy(standard_errors)
    pvalues = copy.deepcopy(pvalues)

    # Collect information.
    pail = dict()
    pail["intercept_predictors"] = dict()
    pail["entries"] = list()

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
        pail["intercept_predictors"].update(pail_predictor["predictor"])
        pail["entries"].extend(pail_predictor["entries"])
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
        pail["intercept_predictors"].update(pail_predictor["predictor"])
        pail["entries"].extend(pail_predictor["entries"])
        pass

    # Return information.
    return pail


def parse_organize_regression_intercept_predictors(
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
    count_intercept_predictors = len(features_predictor_intercept)
    features_sequence = define_features_sequence(
        features=features_predictor_intercept,
    )

    # Parse information for predictor features from regression model.
    pail = parse_regression_intercept_predictors(
        table_predictor_intercept=table_predictor_intercept,
        features_predictor=features_predictor,
        features_sequence=features_sequence,
        name_intercept_source="const",
        name_intercept_product="intercept",
        parameters=parameters,
        standard_errors=standard_errors,
        pvalues=pvalues,
    )

    # Collect information.
    pail["intercept_predictors"]["count_intercept_predictors"] = (
        count_intercept_predictors
    )
    pail["entries"].insert(0, "count_intercept_predictors")

    # Return information.
    return pail


##########
# Regression models


def regress_continuous_linear_ordinary_least_squares(
    table=None,
    feature_response=None,
    features_predictor_fixed=None,
    report=None,
):
    """
    Regress a response feature on a quantitative, continuous, ratio or interval
    scale of measurement against either a single or multiple predictor features
    by a linear ordinary least squares (OLS) model, applying an implementation
    from the StatsModels package in Python.

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
    vector of scalar values on a quantitative, continuous, ratio or interval
    scale of measurement.
    [1.3, 1.5, 1.2, 1.0, 1.7, 1.5, 1.9, 1.1, 1.3, 1.4,]

    Format of information for predictor independent variable(s) is a
    two-dimensional matrix: a first-dimension vector corresponding to
    observations or samples and for each observation a second-dimension vector
    of scalar values corresponding to each feature or variable.
    StatsModels also requires the intentional addition of a constant for the
    intercept.
    [
        [1.3, 5.2, 1.0,],
        [1.5, 5.1, 1.0,],
        [1.2, 5.5, 1.0,],
        ...
    ]

    Implementation: 'statsmodels.regression.linear_model.OLS()'

    References:

    1. Documentation for implementation in StatsModels
       - title: 'statsmodels.regression.linear_model.OLS'
       - site: https://www.statsmodels.org/stable/generated/statsmodels.
                regression.linear_model.OLS.html

    Review: TCW; 3 April 2025

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
    # Passing variables as Pandas Series and Dataframe preserves variable
    # names.
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


def regress_discrete_logistic_logit(
    table=None,
    feature_response=None,
    features_predictor_fixed=None,
    report=None,
):
    """
    Regress a response feature on a discrete, binary scale of measurement
    against either a single or multiple predictor features by a generalized
    linear model (GLM) with the logit link function, applying an implementation
    from the StatsModels package in Python.

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
    vector of scalar values on a discrete, binary scale of measurement.
    [0, 1, 1, 0, 1, 1, 0, 0, 0, 1,]

    Format of information for predictor independent variable(s) is a
    two-dimensional matrix: a first-dimension vector corresponding to
    observations or samples and for each observation a second-dimension vector
    of scalar values corresponding to each feature or variable.
    StatsModels also requires the intentional addition of a constant for the
    intercept.
    [
        [1.3, 5.2, 1.0,],
        [1.5, 5.1, 1.0,],
        [1.2, 5.5, 1.0,],
        ...
    ]

    Complete or Quasi-Complete Separation of discrete (or binary) values of any
    independent variable between discrete values of the dependent variable can
    prevent convergence of the logistic regression either by maximum likelihood
    or regularized maximum likelihood methods. Evaluate the relationships
    between dependent and independent variables to avoid any Complete or
    Quasi-Complete Separation.

    Implementation: 'statsmodels.discrete.discrete_model.Logit()'

    References:

    1. Documentation for implementation in StatsModels
       - title: 'statsmodels.discrete.discrete_model.Logit'
       - site: https://www.statsmodels.org/stable/generated/statsmodels.
                discrete.discrete_model.Logit.html

    Review: TCW; 3 April 2025

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
    # Passing variables as Pandas Series and Dataframe preserves variable
    # names.
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
    model = statsmodels.api.Logit(
        table[feature_response],
        table_predictor_intercept,
        missing="drop",
    )
    # Fit the model.
    handle_model = model.fit(
        maxiter=100, # maximal count of iterations to convergence
    )

    # Regularization of the fit method can help to avoid overfitting and promote
    # convergence.
    # Fit the model using regularized maximum likelihood.
    #handle_model = model.fit_regularized(
    #    method="l1", # solver method
    #    maxiter=100,
    #)

    # Collect information.
    pail = dict()
    pail["model"] = handle_model
    pail["table_predictor_intercept"] = table_predictor_intercept

    # Return information.
    return pail


##########
# Manage regression analysis.


def manage_regression_continuous_linear_ordinary_least_squares(
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

    Review: TCW; 3 April 2025

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
    pail_regression = regress_continuous_linear_ordinary_least_squares(
        table=table,
        feature_response=feature_response,
        features_predictor_fixed=features_predictor_fixed,
        report=report,
    )
    # Extract residuals.
    residuals = pail_regression["model"].resid

    # Organize information from regression for predictor features.
    parameters = pandas.Series(data=pail_regression["model"].params)
    standard_errors = pandas.Series(data=pail_regression["model"].bse)
    pvalues = pandas.Series(data=pail_regression["model"].pvalues)
    pail_predictors = parse_organize_regression_intercept_predictors(
        table_predictor_intercept=pail_regression["table_predictor_intercept"],
        features_predictor=features_predictor_fixed,
        parameters=parameters,
        standard_errors=standard_errors,
        pvalues=pvalues,
        report=report,
    )

    # Collect information.
    pail = dict()
    # Regression model as a whole.
    pail["entries_model"] = [
        "implementation",
        "count_observations_table",
        "count_observations_model",
        "degrees_freedom_model",
        "degrees_freedom_residual",
        "r_square",
        "r_square_adjust",
        "r_square_pseudo",
        "log_likelihood",
        "akaike",
        "bayes",
        "condition",
    ]
    pail["implementation"] = str(
        "statsmodels.regression.linear_model.OLS()"
    )
    pail["count_observations_table"] = int(table.shape[0])
    pail["count_observations_model"] = int(pail_regression["model"].nobs)
    pail["degrees_freedom_model"] = int(
        pail_regression["model"].df_model
    )
    pail["degrees_freedom_residual"] = int(
        pail_regression["model"].df_resid
    )
    pail["r_square"] = pail_regression["model"].rsquared
    pail["r_square_adjust"] = pail_regression["model"].rsquared_adj
    pail["r_square_pseudo"] = float("nan") # hold place for logistic regression
    pail["log_likelihood"] = pail_regression["model"].llf
    pail["akaike"] = pail_regression["model"].aic
    pail["bayes"] = pail_regression["model"].bic
    pail["condition"] = pail_regression["model"].condition_number
    # Intercept and predictors as parts of regression model.
    pail["entries_intercept_predictors"] = pail_predictors["entries"]
    pail.update(
        pail_predictors["intercept_predictors"]
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "manage_regression_continuous_linear_ordinary_least_squares()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("text formula for regression model:")
        print(formula_text)
        putly.print_terminal_partition(level=5)
        print("summary from regression model:")
        print(pail_regression["model"].summary())
        #print(dir(pail_regression))
        #print(pail_regression.params)
        #print(pail_regression.pvalues)
        pass

    # Return information.
    return pail


def manage_regression_discrete_logistic_logit(
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

    Review: TCW; 3 April 2025

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
    pail_regression = regress_discrete_logistic_logit(
        table=table,
        feature_response=feature_response,
        features_predictor_fixed=features_predictor_fixed,
        report=report,
    )
    # Extract residuals.
    residuals = pail_regression["model"].resid_generalized

    # Organize information from regression for predictor features.
    parameters = pandas.Series(data=pail_regression["model"].params)
    standard_errors = pandas.Series(data=pail_regression["model"].bse)
    pvalues = pandas.Series(data=pail_regression["model"].pvalues)
    pail_predictors = parse_organize_regression_intercept_predictors(
        table_predictor_intercept=pail_regression["table_predictor_intercept"],
        features_predictor=features_predictor_fixed,
        parameters=parameters,
        standard_errors=standard_errors,
        pvalues=pvalues,
        report=report,
    )

    # Collect information.
    pail = dict()
    # Regression model as a whole.
    pail["entries_model"] = [
        "implementation",
        "count_observations_table",
        "count_observations_model",
        "degrees_freedom_model",
        "degrees_freedom_residual",
        "r_square",
        "r_square_adjust",
        "r_square_pseudo",
        "log_likelihood",
        "akaike",
        "bayes",
        "condition",
    ]
    pail["implementation"] = str(
        "statsmodels.discrete.discrete_model.Logit()"
    )
    pail["count_observations_table"] = int(table.shape[0])
    pail["count_observations_model"] = int(pail_regression["model"].nobs)
    pail["degrees_freedom_model"] = int(
        pail_regression["model"].df_model
    )
    pail["degrees_freedom_residual"] = int(
        pail_regression["model"].df_resid
    )
    pail["r_square"] = float("nan") # hold place for linear regression
    pail["r_square_adjust"] = float("nan") # hold place for linear regression
    pail["r_square_pseudo"] = pail_regression["model"].prsquared
    pail["log_likelihood"] = pail_regression["model"].llf
    pail["akaike"] = pail_regression["model"].aic
    pail["bayes"] = pail_regression["model"].bic
    pail["condition"] = float("nan") # hold place for linear regression
    # Intercept and predictors as parts of regression model.
    pail["entries_intercept_predictors"] = pail_predictors["entries"]
    pail.update(
        pail_predictors["intercept_predictors"]
    )

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "manage_regression_discrete_logistic_logit()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("text formula for regression model:")
        print(formula_text)
        putly.print_terminal_partition(level=5)
        print("summary from regression model:")
        print(pail_regression["model"].summary())
        #print(dir(pail_regression))
        #print(pail_regression.params)
        #print(pail_regression.pvalues)
        pass

    # Return information.
    return pail


# See document string for informative notes about regressions.
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

    ----------
    Notes

    Scale transformation of predictor independent variables
    - In multiple linear regression, it is common to standardize the scale of
       of predictor independent variables that are on a quantitative,
       continuous scale of measurement, such as ratio or interval.
    - Standardizing the scale of these predictor independent variables makes it
       more reasonable and convenient to compare the magnitudes of estimates
       for parameter coefficients between these different predictor independent
       variables.
    - When the regression model also includes additional predictor independent
       variables that are on a binary scale of measurement, such as dummy
       variables for categories, it is most reasonable to leave these on the
       binary scale.
    - The rough interpretation of the slope parameter coefficient for a binary
       predictor independent variable is the difference in means of the
       response dependent variable between the two groups of observations
       corresponding to the two values of the binary predictor.
    - In regression with multiple predictor independent variables, the rough
       interpretation of the slope for a binary predictor is the same, while
       holding the other predictor covariates constant, perhaps at their
       respective means.

    References:

    1. Forum discussion about standardizing the scales of dummy variables
       title: 'Standardizing dummy variable in multiple linear regression?'
       site: https://stats.stackexchange.com/questions/414053/standardizing-
              dummy-variable-in-multiple-linear-regression

    2. Forum discussion about interpretation of regression for binary variable
       title: 'Binary variable in multiple linear regression?'
       site: https://www.reddit.com/r/AskStatistics/comments/kaiwk2/
              binary_variable_in_multiple_linear_regression/

    ----------
    Notes

    Interpretation
    - The parameter coefficient for the intercept is the model's prediction of
       the response when all predictors have values of zero (0). If the model
       includes dummy predictors to represent features with categorical values,
       then the value of the intercept corresponds to the reference category
       that corresponds to values of zero (0) for dummy predictors that
       represent all the other categories. To represent "n" categories", "n-1"
       dummy predictors are necessary.
    - The parameter coefficient for a continuous predictor indicates that, when
       holding all other covariate predictors constant at any of their
       respective values, a one unit increase in the predictor of interest
       corresponds on average to the coefficient's unit change in the response.

    ----------
    Notes

    Types of Regression Models
    - Continuous response feature
     - Linear Ordinary Least Squares (OLS)
      - requires homoscedasticity of variance in response across ranges of all
         predictors
     - Linear Generalized Least Squares (GLS)
      - allows heteroscedasticity of variance in response across ranges of all
         predictors
    - Discrete response feature
     - Generalized Linear Model with Logit link function
      - response is on discrete binary scale of measurement
      - link function is the logarithm of the odds ratio
     - Generalized Linear Model with Probit (probability unit) link function
      - response is on discrete binary scale of measurement
      - link function bases on the cumulative density function of the standard
         normal distribution
      - less accurate and less interpretable for odds ratios
     - Generalized Linear Model with Poisson link function
      - response is on discrete ordinal scale of measurement
     - Generalized Linear Model with Negative Binomial link function
      - response is on discrete ordinal scale of measurement

    References:

    1. Types of linear regression models in StatsModels
       title: 'Linear Regression'
       site: https://www.statsmodels.org/stable/regression.html

    2. Types of discrete response regression models in StatsModels
       title: 'Regression with Discrete Dependent Variable'
       site: https://www.statsmodels.org/stable/discretemod.html

    3. Types of regression with mixed effects in StatsModels
       title: 'Linear Mixed Effects Models'
       site: https://www.statsmodels.org/stable/mixed_linear.html

    4. Forum discussion of differences between logit, probit, and other link
       functions for regression on discrete responses
       title: 'Difference between logit and probit models'
       site: https://stats.stackexchange.com/questions/20523/difference-
              between-logit-and-probit-models/30909#30909

    ----------
    Notes

    Linear regression against predictors with mixed fixed and random effects
    - ...

    References:

    1. Article about linear regression analysis with mixed fixed and random effects
       title: 'Linear Mixed Effect Models in python using mtcars'
       site: https://medium.com/@josef.waples/linear-mixed-effect-models-in-
              python-using-mtcars-064c6c3e5b32
       note:
       - This article also provides a clear and insightful interpretation of
          the results of regression analysis.

    ----------
    Review: TCW; 3 April 2025
    - On 3 April 2025, TCW confirmed that extracted values from linear OLS
       and discrete generalized Logit regression models, intercept, and
       predictors in the summary table match those reported directly in the
       summary from the implementations of the respective regression models on
       demonstration data.

    arguments:
        table (object): Pandas data-frame table of data with features
            and observations for regression
        index_columns (str): name of single-level index across columns in table
        index_rows (str): name of single-level index across rows in table
        type_regression (str): name of type of regression model to use, either
            'continuous_ols', 'continuous_gls', 'discrete_logit',
            'discrete_probit', 'discrete_poisson', 'discrete_negbinom'
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
    if (type_regression == "continuous_ols"):
        # Perform continuous linear ordinary least squares regression.
        pail = manage_regression_continuous_linear_ordinary_least_squares(
            table=table,
            index_columns="features",
            index_rows="observations",
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            report=report,
        )
        pass
    elif (type_regression == "discrete_logit"):
        # Perform logistic regression.
        pail = manage_regression_discrete_logistic_logit(
            table=table,
            index_columns="features",
            index_rows="observations",
            formula_text=formula_text,
            feature_response=feature_response,
            features_predictor_fixed=features_predictor_fixed,
            report=report,
        )
        pass
    else:
        pail = dict()
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


def collect_organize_record_regression_analysis(
    table=None,
    index_columns=None,
    index_rows=None,
    sequence=None,
    group=None,
    instance=None,
    name_instance=None,
    check_overall=None,
    type_regression=None,
    formula_text=None,
    feature_response=None,
    features_predictor_fixed=None,
    features_predictor_random=None,
    record_extra=None,
    delimiter_list_items=None,
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
        sequence (int): sequential index for instance's name and sort order
        group (str): categorical group of instances
        instance (str): name or designator of instance
        name_instance (str): compound name for instance of parameters
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
        delimiter_list_items (str): delimiter to place between items when
            collapsing or flattening lists
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of information from regression

    """

    # pail (dict): collection of information about a regression analysis
    #     entries_parameter_instance (list<str>): names and sequence of entries
    #         with information about the instance of parameters for the
    #         regression
    #     entries_model (list<str>): names and sequence of entries with
    #         information about the regression model as a whole
    #     entries_intercept_predictors (list<str>): names and sequence of
    #         entries with with information about model's intercept and
    #         predictors
    #     record (dict<str>): record for collection and assembly as rows within
    #         a Pandas data-frame table, using the lists of entries to sort the
    #         columns in the table

    # Collect information.
    pail = dict()
    pail["entries_parameter_instance"] = [
        "sequence",
        "group",
        "instance",
        "name_instance",
        "check_overall",
        "type_regression",
        "formula_text",
        "feature_response",
        "features_predictor_fixed",
        "features_predictor_random",
    ]
    pail["entries_parameter_instance"].extend(list(record_extra.keys()))
    pail["record"] = dict()
    pail["record"]["sequence"] = sequence
    pail["record"]["group"] = group
    pail["record"]["instance"] = instance
    pail["record"]["name_instance"] = name_instance
    pail["record"]["check_overall"] = check_overall
    pail["record"]["type_regression"] = type_regression
    pail["record"]["formula_text"] = formula_text
    pail["record"]["feature_response"] = feature_response
    pail["record"]["features_predictor_fixed"] = delimiter_list_items.join(
        features_predictor_fixed
    )
    pail["record"]["features_predictor_random"] = delimiter_list_items.join(
        features_predictor_random
    )
    pail["record"].update(record_extra)
    # Determine whether to perform regression.
    if (check_overall):
        # Perform regression.
        pail_regression = determine_type_regression_analysis(
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
        # Collect information.
        pail["entries_model"] = copy.deepcopy(pail_regression["entries_model"])
        pail["entries_intercept_predictors"] = copy.deepcopy(
            pail_regression["entries_intercept_predictors"]
        )
        pail["record"].update(pail_regression)
        del pail["record"]["entries_model"]
        del pail["record"]["entries_intercept_predictors"]
    else:
        # Collect information.
        pail["entries_model"] = None
        pail["entries_intercept_predictors"] = None
        pass

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("package: partner")
        print("module: regression.py")
        name_function = str(
            "collect_organize_record_regression_analysis()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=4)
        pass
    # Return information.
    return pail




###############################################################################
# End
