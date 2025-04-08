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
# Date, first execution: 7 April 2025
# Date, last execution: 7 April 2025
# Review: TCW; 7 April 2025
################################################################################
# Note

# title: 'Comparing ANOVA to Mixed Effects Linear Regression for analysis of
#         placebo-controlled repeat measures'

# The purpose of this Python script is to demonstrate a comparison between
# the methods, analysis of variance (ANOVA) with repeated measures and linear
# regression with mixed effects.
# This demonstration applies repeated measures ANOVA and mixed effects linear
# regression to the same data from the same study.
# The study design is


# References:

# 1. Scientific literature article with peer review comparing analysis of
#     variance (ANOVA) with repeat measures to linear regression with mixed
#     effects to analyze the effect of treatment between groups that receive
#     intervention or placebo control
#  PubMed:
#  title: 'Mixed-model regression analysis and dealing with interindividual
#          differences'
#  journal: 'Methods in Enzymology'
#  year: 2004
#  authors: Hans P.A. Van Dongen; ... Greg Maislin
#  note:
#  - the study design included measurement at baseline before treatment and at
#    3 time points after treatment by either placebo or intervention
#  - TCW transferred the raw data from Table 1 to file "table_data.tsv" within
#    the directory 'partner/repository/demonstration/15081686_vandongen_2004'.
#  - TCW encoded the placebo group as value '0' for feature 'condition'.
#  - TCW encoded the intervention group as value '1' for feature 'condition'.
#  - TCW endoced the 'day 0', 'day 1', 'day 2', 'day 3' time points of
#    measurement as values '0', '1', '2', and '3' respectively for feature
#    'time_point'.


################################################################################
# Installation and importation

# Standard
import sys
import os
import copy
import textwrap

# Relevant
import pandas
import scipy
import numpy
#import statsmodels.api
#import pingouin

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


def define_type_table_columns():
    """
    Defines the types of variables for columns in table.

    Review: TCW; 7 April 2025

    arguments:

    raises:

    returns:
        (dict<str>): variable types of columns within table

    """

    # Specify types of variables in columns of table.
    types_columns = dict()
    types_columns["identifier_subject"] = "string"
    types_columns["condition"] = "int32"
    types_columns["time_point"] = "int32"
    types_columns["measurement"] = "float32"
    # Return information.
    return types_columns


def read_source_table_data(
    name_file_table_data=None,
    path_directory_source=None,
    report=None,
):
    """
    Read and organize source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    Review: TCW; 7 april 2025

    arguments:
        name_file_table_data (str): name of source file for table of data in
            text format as a table with tab delimiters between columns and
            newline delimiters between rows
        path_directory_source (str): path to directory from which to read
            source files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of data with features and
            observations for regression

    """

    # Define paths to child files.
    path_file_table_data = os.path.join(
        path_directory_source, name_file_table_data,
    )

    # Read information from file.

    # Table of parameters for parallel instances.
    types_columns = define_type_table_columns()
    table = pandas.read_csv(
        path_file_table_data,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Organize information.
    table["condition"] = pandas.to_numeric(
        table["condition"],
        downcast="integer",
        errors="coerce",
    )
    table["time_point"] = pandas.to_numeric(
        table["time_point"],
        downcast="integer",
        errors="coerce",
    )
    table["measurement"] = pandas.to_numeric(
        table["measurement"],
        downcast="float",
        errors="coerce",
    )

    # Report.
    if report:
        # Organize.
        # Print.
        putly.print_terminal_partition(level=3)
        print("package: partner")
        name_module = str(
            "script_demonstration_compare_anova_regression_placebo.py"
        )
        print("module: " + name_module)
        name_function = str(
            "read_source_table_data()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("table:")
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return table


################################################################################
# Procedure


##########
# Call main procedure.


def execute_procedure(
    name_file_table_data=None,
    path_directory_source=None,
    path_directory_product=None,
    path_directory_dock=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        name_file_table_data (str): name of source file for table of data in
            text format as a table with tab delimiters between columns and
            newline delimiters between rows
        path_directory_source (str): path to directory from which to read
            source files
        path_directory_product (str): path to directory for procedure's product
            directories and files
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
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("package: partner")
        name_module = str(
            "script_demonstration_compare_anova_regression_placebo.py"
        )
        print("module: " + name_module)
        name_function = str(
            "execute_procedure()"
        )
        print("function: " + name_function)
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("name_file_table_data: " + str(name_file_table_data))
        print("path_directory_source: " + str(path_directory_source))
        print("path_directory_product: " + str(path_directory_product))
        print("path_directory_dock: " + str(path_directory_dock))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read source information from file.
    table_source = read_source_table_data(
        name_file_table_data=name_file_table_data,
        path_directory_source=path_directory_source,
        report=report,
    )

    ##########
    # Organize information in table.
    table_source["observations"] = table_source["identifier_subject"]
    table_source["identifier_subject_raw"] = table_source["identifier_subject"]
    table_source["identifier_subject"] = table_source.apply(
        lambda row: str(
            "subject_" + row["identifier_subject_raw"]
        ),
        axis="columns", # apply function to each row
    )
    selection_observations = dict()
    selection_observations["condition"] = [0, 1,]
    selection_observations["time_point"] = [0, 1,]
    features_relevant = [
        "observations",
        "identifier_subject",
        "condition",
        "time_point",
        "time_point_0",
        "time_point_1",
        "condition_by_time_point_1",
        "measurement",
    ]
    table = preg.organize_table_data(
        table=table_source,
        selection_observations=selection_observations,
        features_relevant=features_relevant,
        features_regression=features_relevant,
        features_continuity_scale=list(),
        index_columns_source="features",
        index_columns_product="features",
        index_rows_source="observations",
        index_rows_product="observations",
        adjust_scale=False,
        method_scale=None, # 'z_score' or 'unit_range'
        explicate_indices=True,
        report=report,
    )

    # Two-Way ANOVA with Repeated Measures.
    # This version of ANOVA does not fit the study design and data, because the
    # predictor feature for 'condition' distinguishes between subjects and not
    # within subjects.
    #pail_anova = preg.determine_type_analysis_variance_anova(
    #    table=table,
    #    index_columns="features",
    #    index_rows="observations",
    #    type_anova="2way_repeat",
    #    formula_text="(measurement) ~ (condition) + (time_point)",
    #    feature_response="measurement",
    #    features_predictor_between=["condition",],
    #    features_predictor_within=["time_point_1"],
    #    subject="identifier_subject",
    #    report=report,
    #)

    # Two-Way ANOVA with Repeated Measures and Mixed Effects.
    pail_anova_mix = preg.determine_type_analysis_variance_anova(
        table=table,
        index_columns="features",
        index_rows="observations",
        type_anova="2way_repeat_mixed",
        formula_text="(measurement) ~ (condition) + (time_point_1)",
        feature_response="measurement",
        features_predictor_between=["condition",],
        features_predictor_within=["time_point_1"],
        subject="identifier_subject",
        report=report,
    )
    summary_text_anova = preg.prepare_summary_text_anova(
        title="Two-Way ANOVA with Repeated Measures and Mixed Effects",
        formula_text="(measurement) ~ (condition) + (time_point_1)",
        type_anova="2way_repeat_mixed",
        text_anova=str(pail_anova_mix["anova"].round(3)),
        text_t_1=str(pail_anova_mix["t_test_within"]),
        text_t_2=str(pail_anova_mix["t_test_between"]),
        report=report,
    )

    print(table)

    # Linear Regression with Mixed Effects.
    formula_text = str(
        "(measurement) ~ (condition) + (time_point_1) + " +
        "(condition_by_time_point_1)"
    )
    pail_regression = preg.determine_type_regression_analysis(
        table=table,
        index_columns="features",
        index_rows="observations",
        type_regression="continuous_mixed",
        formula_text=formula_text,
        feature_response="measurement",
        features_predictor_fixed=[
            "condition", "time_point_1", "condition_by_time_point_1",
        ],
        features_predictor_random=["time_point_1",],
        groups_random="identifier_subject",
        report=report,
    )
    print(pail_regression["table_summary"])


    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    name_file_table_data = sys.argv[1]
    path_directory_source = sys.argv[2]
    path_directory_product = sys.argv[3]
    path_directory_dock = sys.argv[4]
    report = sys.argv[5]

    # Call function for procedure.
    execute_procedure(
        name_file_table_data=name_file_table_data,
        path_directory_source=path_directory_source,
        path_directory_product=path_directory_product,
        path_directory_dock=path_directory_dock,
        report=report,
    )

    pass



#
