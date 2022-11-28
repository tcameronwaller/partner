"""
Supply functionality for extraction and organization of information from
reports of other tools.

This module is not directly executable.

This subpackage 'promiscuity' provides executable functionality under the
management of a higher level package. Importation paths represent this
hierarchy.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of Promiscuity
    (https://github.com/tcameronwaller/promiscuity/).

    Promiscuity supports data analysis in multiple other projects.
    Copyright (C) 2022 Thomas Cameron Waller

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
import sklearn.preprocessing
import scipy
import numpy
import statsmodels.api

# Custom
import promiscuity.utility as utility # this import path for subpackage

#dir()
#importlib.reload()

###############################################################################
# Functionality


# parameters:
# 1. path to parent directory
# 2. file name suffix by which to recognize relevant files in parent directory
# 3. type of information to extract: "heritability" or "correlation"
# functionality:
# read names of all files within the parent directory
# iterate on all files, calling the extraction function for each
# collect information from the extraction record for each file
# append the file name to the extraction record
# organize the extracted information in a table
# return table


def read_extract_from_all_ldsc_files_in_directory(
    path_directory=None,
    file_name_pattern=None,
    file_name_pattern_not=None,
    analysis=None,
    report=None,
):
    """
    Reads and extracts information about LDSC SNP heritability or genetic
        correlation analyses from all files within a parent directory.

        Linkage Disequilibrium Score Regression (LDSC)

    arguments:
        path_directory (str): full path to parent directory that contains
            relevant files
        file_name_pattern (str): character string in names of relevant files
        file_name_pattern_not (str): character string not in names of relevant
            files
        analysis (str): type of analysis in LDSC, either 'heritability' or
            'correlation'
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table with variables across columns and
            samples or instances across rows

    """

    names_files = utility.extract_directory_file_names_filter_by_name(
        path=path_directory,
        name=file_name_pattern,
        name_not=file_name_pattern_not,
    )

    print("here is the list of file names")
    print(names_files)

    # Collect information from analysis on each study.
    records = list()
    # Iterate on files for each study.
    for name_file in names_files:
        # Define full path to file.
        path_file = os.path.join(
            path_directory, name_file,
        )

        # print file path
        print(name_file)
        print(path_file)

        #
        #record = read_extract_heritability_design_study_detail(
        #    file_name=file_name,
        #    file_prefix=file_prefix,
        #    file_suffix=file_suffix,
        #    design=design,
        #    study=child_directory,
        #    path_parent_directory=path_child_directory,
        #)
        #records.append(record)
        pass

    # Organize table.
    if False:
        table = utility.convert_records_to_dataframe(
            records=records
        )
        table.reset_index(
            level=None,
            inplace=True,
            drop=True,
        )
        columns = [
            "design",
            "study",
            "summary_error", "summary_interval",
            "variants",
            "heritability", "standard_error", "confidence_95_range",
            "confidence_95_low", "confidence_95_high",
            "ratio",
            "ratio_standard_error",
        ]
        table = table.loc[
            :, table.columns.isin(columns)
        ]
        table = table[[*columns]]
        table.sort_values(
            by=["design", "study",],
            axis="index",
            ascending=True,
            inplace=True,
        )
    # Return information.
    #return table
    return dict()




###############################################################################
# Procedure
# Currently, this module is not executable.
