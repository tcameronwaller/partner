"""
Supply functionality for parallelization in Python.

This module is not directly executable.

This module within subpackage 'partner' provides executable functionality under
the management of a higher level package. Importation paths must represent this
hierarchy.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Monroe, North Carolina 28110
    United States of America

License:

    This file is part of project 'partner'
    (https://github.com/tcameronwaller/partner/).

    Project 'partner' supports data analysis in multiple other projects.
    Copyright (C) 2024 Thomas Cameron Waller

    Project 'partner' is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    Project 'partner' is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'partner'.
    If not, see <http://www.gnu.org/licenses/>.
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
import math
import itertools
import functools
import multiprocessing
import datetime
import time


# Relevant

import pandas
import sklearn
import sklearn.preprocessing
import scipy
import numpy
import statsmodels.api

# Custom
import partner.utility as putly
import partner.extraction as pextr
import partner.organization as porg
import partner.description as pdesc
import partner.regression as preg
import partner.plot as plot
#import partner.parallelization as prallel


#dir()
#importlib.reload()

###############################################################################
# Functionality


def drive_procedure_parallel(
    function_control=None,
    instances=None,
    parameters=None,
    cores=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    Review: TCW; 28 March 2024

    arguments:
        function_control (object): handle to function to control iterative
            execution of a process
        instances (list): list of instances to pass to each iterative execution
            of the control function, with each instance potentially being a
            simple varialbe or a dictionary collection of information
        parameters (dict): additional parameters that are common for all
            iterative instances, such as file paths for read and write
        cores (int): count of corse to use for parallel iteration
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Determine begin date and time.
    time_begin = datetime.datetime.now()

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print(
            "Count of instances for iterative execution: " +
            str(len(instances))
        )
        print("Time at begin: " + str(time_begin))
        pass

    # Set up partial function for iterative execution.
    # Each iteration uses the same values of the collection of additional
    # parameters.
    execute_instance_parallel = functools.partial(
        function_control,
        parameters=parameters,
    )

    # Initialize multiprocessing pool.
    #pool = multiprocessing.Pool(processes=os.cpu_count())
    pool = multiprocessing.Pool(processes=cores)

    # Initialize test instance.
    instances_check = []

    #log = pool.map(execute_instance_parallel, instances_check)
    log = pool.map(execute_instance_parallel, instances)

    # Pause procedure.
    time.sleep(5.0)

    # Determine end date and time.
    time_end = datetime.datetime.now()
    # Determine duration.
    time_duration = (time_end - time_begin)

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Time at end: " + str(time_end))
        print("Process duration: " + str(time_duration))
        pass
    pass






###############################################################################
# Procedure
# Currently, this module is not executable.

##########
