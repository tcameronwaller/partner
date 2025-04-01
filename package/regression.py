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

# Note: TCW; 4 October 2022
# There are many different methods to "standardize", "scale", "transform", or
# "standardize" the distributions of multiple variables to simplify their
# comparisons to each other.




#
