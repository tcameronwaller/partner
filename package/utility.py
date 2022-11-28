"""
Supply basic utility functionality.

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
import scipy
import numpy
import statsmodels.api

# Custom

#dir()
#importlib.reload()

###############################################################################
# Functionality


# General.


def convert_string_low_alpha_num(characters):
    """
    Converts string of characters to lower case with only alphabetical or
    numerical characters.

    arguments:
        characters (str): characters in a string

    raises:

    returns:
        (str): characters in a string

    """

    # Convert all characters to lower case.
    characters_lower = characters.lower()
    # Remove all characters other than alphabetical or numerical characters.
    characters_novel = characters_lower
    for character in characters_lower:
        if (
            (character not in string.ascii_letters) and
            (character not in string.digits)
        ):
            characters_novel = characters_novel.replace(character, "")
    return characters_novel


def extract_subdirectory_names(path=None):
    """
    Extracts names of subdirectories.

    Reads all contents of a directory and returns only those of subdirectories.

    arguments:
        path (str): path to directory

    raises:

    returns:
        (list<str>): names of subdirectories

    """

    contents = os.listdir(path=path)
    names = list(filter(
        lambda content: os.path.isdir(os.path.join(path, content)),
        contents
    ))
    return names


def extract_directory_file_names(path=None):
    """
    Extracts names of files within a directory.

    Reads all contents of a directory and returns only those of files (not
    directories).

    arguments:
        path (str): path to directory

    raises:

    returns:
        (list<str>): names of files within the parent directory

    """

    contents = os.listdir(path=path)
    files = list(filter(
        lambda content: (
            (os.path.isfile(os.path.join(path, content))) and
            (not os.path.isdir(os.path.join(path, content)))
        ), contents
    ))
    return files


def extract_directory_file_names_filter_by_name(
    path=None,
    name=None,
    name_not=None,
):
    """
    Extracts the names of all files within a parent directory and then filters
    these file names by whether they include a specific string.

    arguments:
        path (str): path to directory
        name (str): characters within names of files to keep
        name_not (str): characters not within names of files to keep

    raises:

    returns:
        (list<str>): names of files within the parent directory

    """

    names_files = extract_directory_file_names(path=path)
    names_files_keep = list(filter(
        lambda name_file: (
            (str(name) in str(name_file)) and
            (str(name_not) not in str(name_file))
        ), names_files
    ))
    return names_files_keep


def create_directories(path=None):
    """
    Creates a directory if it does not already exist.

    arguments:
        path (str): path to directory

    raises:

    returns:

    """

    if not os.path.exists(path):
        os.makedirs(path)


def remove_file(path=None):
    """
    Removes a file if it exists.

    arguments:
        path (str): path to file

    raises:

    returns:

    """

    if os.path.exists(path):
        os.remove(path)


def remove_directory(path=None):
    """
    Removes a directory if it exists.

    arguments:
        path (str): path to directory

    raises:

    returns:

    """

    if os.path.exists(path):
        # Only for empty directory.
        #os.rmdir(path)
        # For non empty directory.
        shutil.rmtree(path)


def remove_empty_directory(path=None):
    """
    Removes a directory if it is empty.

    arguments:
        path (str): path to directory

    raises:

    returns:

    """

    if (os.path.exists(path)) and (len(os.listdir(path)) < 1):
        os.rmdir(path)


def create_directory(path=None):
    """
    Creates a directory if it is does not already exist.

    arguments:
        path (str): path to directory

    raises:

    returns:

    """

    if not (os.path.exists(path)):
        os.mkdir(path)


def compress_file_gzip(path=None):
    """
    Copies and compresses a file by gzip.

    arguments:
        path (str): path to file

    raises:

    returns:

    """

    with open(path, "rb") as file_source:
        with gzip.open((path + ".gz"), "wb") as file_product:
            shutil.copyfileobj(file_source, file_product)
    pass


def decompress_file_gzip(path=None):
    """
    Copies and decompresses a file by gzip.

    Removes any ".gz" suffix from the file name.

    arguments:
        path (str): path to file

    raises:

    returns:

    """

    split_strings = path.split(".")
    if split_strings[-1] == "gz":
        path_novel = ".".join(split_strings[0:-1])
    else:
        path_novel = path
    with gzip.open(path, "rb") as file_source:
        with open(path_novel, "wb") as file_product:
            shutil.copyfileobj(file_source, file_product)
    pass


def print_terminal_partition(level=None):
    """
    Prints string to terminal to partition reports for clarity.

    arguments:
        level (int): level or type of partition to create and print.

    raises:

    returns:

    """

    if level == 1:
        partition = textwrap.dedent("""\



            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------



            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------



            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------



        """)
    elif level == 2:
        partition = textwrap.dedent("""\



            --------------------------------------------------
            --------------------------------------------------
            --------------------------------------------------



        """)
    elif level == 3:
        partition = textwrap.dedent("""\



            ----------
            ----------
            ----------



        """)
    elif level == 4:
        partition = textwrap.dedent("""\

            ----------

        """)
    elif level == 5:
        partition = textwrap.dedent("""----------""")
    else:
        partition = ""
    print(partition)
    pass


def determine_logical_or_combination_binary_missing(
    first=None,
    second=None,
    single_false_sufficient=None,
):
    """
    Determines the logical "or" combination of two binary logical variables,
    each of which can be "true" (1), "false" (0), or "missing" ("nan").

    Current implementation evaluates to "true" (1) even if only one of the
    variables is true and the other is false or missing. This interpretation is
    consistent with "or" logic for the "true" scenario, but not for the "false"
    scenario. Hence in this implementation, a single "true" is sufficient.

    Option "single_false_sufficient" is most appropriate in cases in which both
    variables have similar meaning and their combination has potential to reduce
    missing information and loss of samples.

    arguments:
        first (float): first binary logical variable that can be true (1),
            false (0), or missing ("nan")
        second (float): second binary logical variable that can be true (1),
            false (0), or missing ("nan")
        single_false_sufficient (bool): whether to assign a false value if one
            variable is valid and false while the other variable is missing

    raises:

    returns:
        (float): interpretation value

    """

    # pandas.isna() or math.isnan()

    # Comparison.
    if (
        ((not math.isnan(first)) and (first == 1)) or
        ((not math.isnan(second)) and (second == 1))
    ):
        # Either one of the two variables is valid and true.
        # Interpret as true even if one variable is true and the other is
        # false.
        # Also interpret as true even if only one variable is valid and true and
        # the other is missing.
        value = 1
    elif (
        ((not math.isnan(first)) and (first == 0)) and
        ((not math.isnan(second)) and (second == 0))
    ):
        # Both variables are valid and false.
        value = 0
    elif (
        ((not math.isnan(first)) and (first == 0)) and
        ((math.isnan(second)))
    ):
        # Only first variable is valid and false.
        # Second variable is missing and could either be true or false.
        if (single_false_sufficient):
            value = 0
        else:
            value = float("nan")
    elif (
        ((math.isnan(first))) and
        ((not math.isnan(second)) and (second == 0))
    ):
        # Only second variable is valid and false.
        # First variable is missing and could either be true or false.
        if (single_false_sufficient):
            value = 0
        else:
            value = float("nan")
    elif (
        (math.isnan(first)) and
        (math.isnan(second))
    ):
        # Both variables are missing.
        value = float("nan")
    else:
        value = float("nan")
    # Return information.
    return value


def determine_binary_categorical_product_of_two_binary_variables(
    product=None,
    first=None,
    second=None,
):
    """
    Determines any one of four binary categorical variables that represent the
    product of two binary variables.

    Here are the combinations in which each product variable has a binary true
    value (1).

    [prefix]_1: ("first" variable = 1) and ("second" variable = 1)
    [prefix]_2: ("first" variable = 1) and ("second" variable = 0)
    [prefix]_3: ("first" variable = 0) and ("second" variable = 1)
    [prefix]_4: ("first" variable = 0) and ("second" variable = 0)

    If either "first" or "second" variables have null, missing values then all
    product variables also have null, missing values.

    arguments:
        product (int): count of product definition, 1 through 4
        first (str): name of first column with binary logical representation of
            a variable
        second (str): name of second column with binary logical representation
            of a variable

    raises:

    returns:
        (bool): interpretation value

    """

    # Determine product value.
    if (
        (math.isnan(first)) or
        (math.isnan(second))
    ):
        value = float("nan")
    else:
        # The relevant variables have valid values.
        if (product == 1):
            if (
                (0.5 <= first and first < 1.5) and
                (0.5 <= second and second < 1.5)
            ):
                # first = 1
                # second = 1
                value = 1
            else:
                value = 0
        elif (product == 2):
            if (
                (0.5 <= first and first < 1.5) and
                (-0.5 <= second and second < 0.5)
            ):
                # first = 1
                # second = 0
                value = 1
            else:
                value = 0
        elif (product == 3):
            if (
                (-0.5 <= first and first < 0.5) and
                (0.5 <= second and second < 1.5)
            ):
                # first = 0
                # second = 1
                value = 1
            else:
                value = 0
        elif (product == 4):
            if (
                (-0.5 <= first and first < 0.5) and
                (-0.5 <= second and second < 0.5)
            ):
                # first = 0
                # second = 0
                value = 1
            else:
                value = 0
            pass
        pass
    # Return.
    return value


def determine_binary_categorical_products_of_two_binary_variables(
    table=None,
    first=None,
    second=None,
    prefix=None,
    report=None,
):
    """
    Determines four binary categorical variables that represent the product of
    two binary variables.
    Uses a prefix in the names of the four product variables.

    arguments:
        table (object): Pandas data frame of feature variables across columns
            and observation records across rows
        first (str): name of first column with binary logical representation of
            a variable
        second (str): name of second column with binary logical representation
            of a variable
        prefix (str): prefix for name of the four product columns
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Copy data.
    table = table.copy(deep=True)
    # Define product variables.
    table[str(prefix + "_1")] = table.apply(
        lambda row:
            determine_binary_categorical_product_of_two_binary_variables(
                product=1,
                first=row[first],
                second=row[second],
            ),
        axis="columns", # apply across rows
    )
    table[str(prefix + "_2")] = table.apply(
        lambda row:
            determine_binary_categorical_product_of_two_binary_variables(
                product=2,
                first=row[first],
                second=row[second],
            ),
        axis="columns", # apply across rows
    )
    table[str(prefix + "_3")] = table.apply(
        lambda row:
            determine_binary_categorical_product_of_two_binary_variables(
                product=3,
                first=row[first],
                second=row[second],
            ),
        axis="columns", # apply across rows
    )
    table[str(prefix + "_4")] = table.apply(
        lambda row:
            determine_binary_categorical_product_of_two_binary_variables(
                product=4,
                first=row[first],
                second=row[second],
            ),
        axis="columns", # apply across rows
    )

    # Organize information for report.
    table_report = table.copy(deep=True)
    columns_report = [
        first,
        second,
        str(prefix + "_1"),
        str(prefix + "_2"),
        str(prefix + "_3"),
        str(prefix + "_4"),
    ]
    table_report = table_report.loc[
        :, table_report.columns.isin(columns_report)
    ]
    table_report = table_report[[*columns_report]]
    table_report.dropna(
        axis="index",
        how="any",
        subset=columns_report,
        inplace=True,
    )
    table_1 = table_report.loc[
        (table_report[str(prefix + "_1")] > 0.5), :
    ]
    table_2 = table_report.loc[
        (table_report[str(prefix + "_2")] > 0.5), :
    ]
    table_3 = table_report.loc[
        (table_report[str(prefix + "_3")] > 0.5), :
    ]
    table_4 = table_report.loc[
        (table_report[str(prefix + "_4")] > 0.5), :
    ]
    # Report.
    if report:
        print_terminal_partition(level=2)
        print(
            "report: " +
            "determine_binary_categorical_products_of_two_binary_variables()"
        )
        #print_terminal_partition(level=3)
        #print(table_report)
        print_terminal_partition(level=3)
        print("Counts of records (rows) in each product category...")
        print("Category 1: " + str(table_1.shape[0]))
        print("Category 2: " + str(table_2.shape[0]))
        print("Category 3: " + str(table_3.shape[0]))
        print("Category 4: " + str(table_4.shape[0]))
        pass
    # Return information.
    return table


def interpret_raw_string_value_missingness_convert_to_float(
    value_raw=None,
):
    """
    Inteprets the missingness of a raw string value and converts this value to
    a float.

    Note:

    arguments:
        value_raw (str): raw string value for a floating point number

    raises:

    returns:
        (float): interpretation value

    """

    # Alternative.
    #table["variable"] = table["variable_raw"].astype("string").copy(
    #    deep=True,
    #)
    #table["variable"].replace(
    #    "",
    #    numpy.nan,
    #    inplace=True,
    #)
    #table["variable"] = table["variable"].astype("float")

    # Define string values indicative of missing values.
    missingness = [
        "nan", "naN", "nAN", "NAN", "NaN", "NAn", "Nan", "nAn",
        "NA", "N/A", "N/a", "n/A", "n/a",
        "<NAN>", "<nan>", "<NA>", "<na>",
        "null",
    ]

    # Interpret value.
    if (
        (not pandas.isna(value_raw)) and
        (len(str(value_raw)) > 0) and
        (str(value_raw) not in missingness)
    ):
        # The variable has a valid, non-missing value.
        # Convert the value to a float.
        interpretation = float(value_raw)
    else:
        # The variable has a missing or uninterpretable value
        interpretation = float("nan")
    # Return.
    return interpretation


def prioritize_combination_values_string(
    value_priority=None,
    value_spare=None,
):
    """
    Determines the combination or consensus value between a clear priority and
    spare.

    arguments:
        value_priority (str): priority value
        value_spare (str): spare value that is only relevant if the priority
            value is missing

    raises:

    returns:
        (str): choice value

    """

    # Determine identifier of priority genotype.
    if (
        (not pandas.isna(value_priority)) and
        (len(str(value_priority)) > 0)
    ):
        # Priority value is not missing.
        choice = str(copy.deepcopy(value_priority))
    elif (
        (not pandas.isna(value_spare)) and
        (len(str(value_spare)) > 0)
    ):
        # Priority value is missing.
        # Spare value is not missing.
        choice = str(copy.deepcopy(value_spare))
        pass
    else:
        # Both priority and spare values are missing.
        # There is not a value available.
        choice = ""
    # Return information.
    return choice


def prioritize_combination_values_float(
    value_priority=None,
    value_spare=None,
):
    """
    Determines the combination or consensus value between a clear priority and
    spare.

    arguments:
        value_priority (float): priority value
        value_spare (float): spare value that is only relevant if the priority
            value is missing

    raises:

    returns:
        (float): choice value

    """

    # Determine identifier of priority genotype.
    if (
        (not pandas.isna(value_priority))
    ):
        # Priority value is not missing.
        choice = float(copy.deepcopy(value_priority))
    elif (
        (not pandas.isna(value_spare))
    ):
        # Priority value is missing.
        # Spare value is not missing.
        choice = float(copy.deepcopy(value_spare))
        pass
    else:
        # Both priority and spare values are missing.
        # There is not a value available.
        choice = float("nan")
    # Return information.
    return choice


##########
# Variable interpretations relevant and specific to human physiology


def determine_human_physiology_age(
    value_raw=None,
):
    """
    Determines value of age with consideration of relevant range for human
    physiology.

    arguments:
        value_raw (str): raw string value for an age

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret string value's missingness and convert to float.
    value = interpret_raw_string_value_missingness_convert_to_float(
        value_raw=value_raw,
    )
    # Comparison.
    if (
        (not pandas.isna(value)) and
        (0.0 <= value and value < 150.0)
    ):
        # Value has a valid, non-missing value.
        # Value is within a human physiologically relevant range.
        age = value
    else:
        # Value has a missing value or is not within a human physiologically
        # relevant range.
        age = float("nan")
    # Return information.
    return age


def determine_human_physiology_body_mass_index(
    value_raw=None,
):
    """
    Determines value of body mass index (BMI) with consideration of relevant
    range for human physiology.

    Maximum realistic body mass index (BMI) is less than 185, which would be the
    BMI of a person with height 6 feet 1 inches and body weight 1,400 pounds.

    https://www.nhlbi.nih.gov/health/educational/lose_wt/BMI/bmicalc.htm

    arguments:
        value_raw (str): raw string value for a body mass index

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret string value's missingness and convert to float.
    value = interpret_raw_string_value_missingness_convert_to_float(
        value_raw=value_raw,
    )
    # Comparison.
    if (
        (not pandas.isna(value)) and
        (5.0 <= value and value < 190.0)
    ):
        # Value has a valid, non-missing value.
        # Value is within a human physiologically relevant range.
        body = value
    else:
        # Value has a missing value or is not within a human physiologically
        # relevant range.
        body = float("nan")
    # Return information.
    return body


##########
# Read text from file


def read_file_text_table(path_file=None, names=None, delimiter=None):
    """
    Reads and organizes source information from file

    This function reads and organizes relevant information from file.

    arguments:
        path_file (str): path to directory and file
        names (list<str>): names for values in each row of table
        delimiter (str): delimiter between values in the table

    returns:
        (list<dict>): tabular information from file

    raises:

    """

    # Read information from file
    #with open(path_file_source, "r") as file_source:
    #    content = file_source.read()
    with open(path_file, "r") as file_source:
        reader = csv.DictReader(
            file_source, fieldnames=names, delimiter=delimiter
        )
        information = list(map(lambda row: dict(row), list(reader)))
    # Return information
    return information


def read_file_text(
    path_file=None
):
    """
    Reads information from file as a text string.

    arguments:
        path_file (str): path to directory and file

    returns:
        (str): information from file

    raises:

    """

    # Read information from file
    #with open(path_file_source, "r") as file_source:
    #    content = file_source.read()
    with open(path_file, "r") as file_source:
        content = file_source.read()
    # Return information
    return content


def read_file_text_list(
    delimiter=None,
    path_file=None
):
    """
    Reads and organizes source information from file.

    Delimiters include "\n", "\t", ",", " ".

    arguments:
        delimiter (str): delimiter between elements in list
        path_file (str): path to directory and file

    returns:
        (list<str>): information from file

    raises:

    """

    # Read information from file
    content = read_file_text(path_file=path_file)
    # Split content by line delimiters.
    values_split = content.split(delimiter)
    values_strip = list(map(lambda value: value.strip(), values_split))
    # Return information
    return values_strip


def read_file_text_lines(
    path_file=None,
    start=None,
    stop=None,
):
    """
    Reads a range of lines from a text file.

    arguments:
        path_file (str): path to directory and file
        start (int): index of line to start reading
        stop (int): index of line to stop reading

    returns:
        (list<str>): information from file

    raises:

    """

    # Read information from file
    #with open(path_file_source, "r") as file_source:
    #    content = file_source.read()
    lines = list()
    with open(path_file, "r") as file_source:
        line = file_source.readline()
        count = 0
        while line:
            if (start <= count and count < stop):
                lines.append(line.strip())
            line = file_source.readline()
            count += 1
    # Return information.
    return lines


def read_file_text_lines_elements(
    path_file=None,
    delimiter=None,
    index=None,
    start=None,
    stop=None,
    report=None,
):
    """
    Reads all lines from a text file and collects elements from a single
    index within each line.

    arguments:
        path_file (str): path to directory and file
        delimiter (str): delimiter between elements in list
        index (int): index of element to collect from each line
        start (int): index of line to start reading
        stop (int): index of line to stop reading
        report (bool): whether to print element from each line

    returns:
        (list<str>): information from file

    raises:

    """

    # Read information from file.
    elements = list()
    with open(path_file, "r") as file_source:
        line = file_source.readline()
        count = 0
        while line:
            if (start <= count and count < stop):
                row = line.strip().split(delimiter)
                element = row[index]
                elements.append(element)
                if report:
                    print("line " + str(count) + ": " + str(element))
            line = file_source.readline()
            count += 1
    # Return information.
    return elements


def read_all_pandas_tables_files_within_parent_directory(
    path_directory_parent=None,
    types_pandas_table_read=None,
    report=None,
):
    """
    Reads all files within parent directory as Pandas data-frame tables and
    organizes these within a dictionary collection with entry names from the
    original file names.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        path_directory_parent (str): path to parent directory for files
        types_pandas_table_read (dict<str>): variable types for read of table in
            Pandas
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables with entry names
            (keys) derived from original names of files

    """

    # Read names of files from parent directory.
    contents = os.listdir(path=path_directory_parent)
    names_files = list(filter(
        lambda content: (os.path.isfile(os.path.join(
            path_directory_parent, content
        ))),
        contents
    ))
    # Report.
    if report:
        print("Names of files within parent directory:")
        print(names_files)
    # Iterate on names of files to read and organize tables.
    # Collect tables.
    pail = dict()
    for name_file in names_files:
        # Determine name of table.
        name_table = name_file.replace(".tsv", "")
        # Specify directories and files.
        path_table = os.path.join(
            path_directory_parent,
            name_file,
        )
        # Read information from file.
        table = pandas.read_csv(
            path_table,
            sep="\t",
            header=0,
            dtype=types_pandas_table_read,
            na_values=["nan", "na", "NAN", "NA",],
        )
        table.reset_index(
            level=None,
            inplace=True,
            drop=True,
        )
        # Collect table.
        pail[name_table] = table.copy(deep=True)
        pass
    # Return information.
    return pail


def write_file_text_table(
    information=None,
    path_file=None,
    names=None,
    delimiter=None,
    header=None,
):
    """
    Writes to file in text format information from a list of dictionaries.

    arguments:
        information (list<str>): information
        path_file (str): path to directory and file
        names (list<str>): names for values in each row of table
        delimiter (str): delimiter between values in the table
        header (bool): whether to write column headers in file

    returns:

    raises:

    """

    # Write information to file
    #with open(out_file_path_model, "w") as out_file:
    #    out_file.write(content_identifier)
    with open(path_file, "w") as file_product:
        writer = csv.DictWriter(
            file_product, fieldnames=names, delimiter=delimiter
        )
        if header:
            writer.writeheader()
        writer.writerows(information)


def write_file_text_list(
    elements=None,
    delimiter=None,
    path_file=None
):
    """
    Writes to file in text format information from an array of strings.

    Delimiters include "\n", "\t", ",", " ".

    arguments:
        elements (list<str>): character elements
        delimiter (str): delimiter between elements in list
        path_file (str): path to directory and file

    returns:

    raises:

    """

    # Write information to file
    with open(path_file, "w") as file_product:
        string = delimiter.join(elements)
        file_product.write(string)
    pass


def print_file_lines(path_file=None, start=None, stop=None):
    """
    Reads and prints a specific range of lines from a file.

    arguments:
        path_file (str): path to directory and file
        start (int): index of line in file at which to start reading
        stop (int): index of line in file at which to stop reading

    returns:

    raises:

    """

    # Read information from file
    #with open(path_file_source, "r") as file_source:
    #    content = file_source.read()
    with open(path_file, "r") as file_source:
        line = file_source.readline()
        count = 0
        while line and (start <= count and count < stop):
            print("***line {}***: {}".format(count, line.strip()))
            line = file_source.readline()
            count += 1


def parse_text_list_values(
    collection=None,
    delimiter=None,
):
    """
    Parse a textual representation of a list or array of values.

    arguments:
        collection (str): textual string of values in a list or array
        delimiter (str): delimiter between values

    raises:

    returns:
        (list<str>): values

    """

    if (
        (len(str(collection)) > 0) and (delimiter in str(collection))
    ):
        values_split = str(collection).split(delimiter)
        #values = list()
        #for value_raw in values_split:
        #    value = str(value_raw).strip()
        #    if (len(value) > 0):
        #        values.append(str(value))
        #    pass
        values_strip = list(map(lambda value: str(value).strip(), values_split))
        values = list(filter(lambda value: (len(str(value)) > 0), values_strip))
    else:
        values = list()
    return values


def determine_any_actual_values_match_comparisons(
    values_actual=None,
    values_comparison=None,
):
    """
    Determines whether any actual values match a list of comparison values.

    arguments:
        values_actual (list<str>): actual values to test for any matches to
            comparison values
        values_comparison (list<str>): set of values for comparison

    raises:

    returns:
        (bool): whether any actual values match the comparison values

    """

    match = False
    if (len(values_actual) > 0):
        for value in values_actual:
            if (
                (len(str(value)) > 0) and
                (value in values_comparison)
            ):
                match = True
                pass
            pass
        pass
    return match


def find(match=None, sequence=None):
    """
    Finds the first element in a sequence to match a condition, otherwise none

    arguments:
        match (function): condition for elements to match
        sequence (list): sequence of elements

    returns:
        (object | NoneType): first element from sequence to match condition or
            none

    raises:

    """

    for element in sequence:
        if match(element):
            return element
    return None


def find_index(match=None, sequence=None):
    """
    Finds index of first element in sequence to match a condition, otherwise -1

    arguments:
        match (function): condition for elements to match
        sequence (list): sequence of elements

    returns:
        (int): index of element if it exists or -1 if it does not exist

    raises:

    """

    for index, element in enumerate(sequence):
        if match(element):
            # Element matches condition
            # Return element's index
            return index
    # Not any elements match condition
    # Return -1
    return -1


def find_all(match=None, sequence=None):
    """
    Finds all elements in a sequence to match a condition, otherwise none.

    arguments:
        match (function): condition for elements to match
        sequence (list): sequence of elements

    returns:
        (list<dict> | NoneType): elements from sequence to match condition or
            none

    raises:

    """

    matches = []
    for element in sequence:
        if match(element):
            matches.append(element)
    if len(matches) > 0:
        return matches
    else:
        return None


def select_elements_by_sets(
    names=None,
    sets=None,
    count=None,
):
    """
    Selects unique elements that belong to at a minimal count of specific sets.

    This algorithm assumes that each element can belong to multiple sets but
    that each set's elements are unique or nonredundant.

    Data format
    - elements (dict)
    -- element (list)
    --- set (str)

    arguments:
        names (list<str>): names of sets to consider
        sets (dict<list<str>>): sets of elements
        count (int): minimal count of sets to which an element must belong

    returns:
        (list): unique elements that belong to minimal count of specific sets

    raises:

    """

    # Collect sets to which each element belongs.
    elements = dict()
    for name in names:
        for element in collect_unique_elements(elements_original=sets[name]):
            if element not in elements:
                elements[element] = list()
                elements[element].append(name)
            else:
                elements[element].append(name)
        pass
    pass

    # Filter elements by count of sets.
    passes = list()
    for element in elements:
        if len(elements[element]) >= count:
            passes.append(element)

    # Return elements that pass filter.
    return passes


def collect_unique_elements(elements=None):
    """
    Collects unique elements

    arguments:
        elements (list): sequence of elements

    returns:
        (list): unique elements

    raises:

    """

    elements_unique = []
    for element in elements:
        if element not in elements_unique:
            elements_unique.append(element)
    return elements_unique


def collect_value_from_records(key=None, records=None):
    """
    Collects a single value from multiple records

    arguments:
        key (str): key of value in each record
        records (list<dict>): sequence of records

    raises:

    returns:
        (list): values from records

    """

    def access(record):
        return record[key]
    return list(map(access, records))


def collect_values_from_records(key=None, records=None):
    """
    Collects values from multiple records.

    arguments:
        key (str): key of value in each record
        records (list<dict>): sequence of records

    raises:

    returns:
        (list): values from records

    """

    collection = []
    for record in records:
        collection.extend(record[key])
    return collection


def compare_lists_by_inclusion(list_one=None, list_two=None):
    """
    Compares lists by inclusion.

    Returns True if all elements in list_two are in list_one.

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (bool): whether first list includes all elements from second

    raises:

    """

    def match(element_two=None):
        return element_two in list_one
    matches = list(map(match, list_two))
    return all(matches)


def compare_lists_by_mutual_inclusion(list_one=None, list_two=None):
    """
    Compares lists by mutual inclusion

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (bool): whether each list includes all elements from the other

    raises:

    """

    forward = compare_lists_by_inclusion(
        list_one=list_one,
        list_two=list_two
    )
    reverse = compare_lists_by_inclusion(
        list_one=list_two,
        list_two=list_one
    )
    return forward and reverse


def filter_common_elements(list_minor=None, list_major=None):
    """
    Filters elements in minor list by whether the major list also includes them.

    arguments:
        list_minor (list): list of elements
        list_major (list): list of elements

    returns:
        (list): elements from minor list that major list also includes

    raises:

    """

    def match(element=None):
        return element in list_major
    return list(filter(match, list_minor))


def filter_unique_common_elements(list_one=None, list_two=None):
    """
    Filters elements by whether both of two lists include them

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (list): elements that both of two lists include

    raises:

    """

    elements_common = filter_common_elements(
        list_one=list_one,
        list_two=list_two,
    )
    elements_unique = collect_unique_elements(
        elements_original=elements_common
    )
    return elements_unique


def filter_unique_union_elements(list_one=None, list_two=None):
    """
    Filters unique elements from union of two lists.

    arguments:
        list_one (list): list of elements
        list_two (list): list of elements

    returns:
        (list): elements that both of two lists include

    raises:

    """

    union = list_one + list_two
    unique = collect_unique_elements(elements_original=union)
    return unique


def filter_unique_exclusion_elements(
    elements_exclusion=None,
    elements_total=None
):
    """
    Filters unique elements by exclusion.

    arguments:
        elements_exclusion (list): list of elements
        elements_total (list): list of elements

    returns:
        (list): elements from total list not in exclusion list

    raises:

    """

    elements_novel = []
    for element in elements_total:
        if element not in elements_exclusion:
            elements_novel.append(element)
    elements_unique = collect_unique_elements(elements_original=elements_novel)
    return elements_unique


def collect_records_targets_by_categories(
    target=None,
    category=None,
    records=None
):
    """
    Collects values of a target attribute that occur together in records with
    each value of another category attribute.
    Each record has a single value of the target attribute.
    Each record can have either a single value or multiple values of the
    category attribute.
    These collections do not necessarily include only unique values of the
    target attribute.

    arguments:
        target (str): name of attribute in records to collect for each category
        category (str): name of attribute in records to define categories
        records (list<dict>): records with target and category attributes

    raises:

    returns:
        (dict<list<str>>): values of the target attribute that occur together
            in records with each value of the category attribute

    """

    def collect_record_target_by_category(
        target_value=None,
        category_value=None,
        collection_original=None,
    ):
        collection_novel = copy.deepcopy(collection_original)
        # Determine whether collection includes the category's value.
        if category_value in collection_novel.keys():
            # Collection includes the category's value.
            target_values = collection_novel[category_value]
            target_values.append(target_value)
            # Include target's value in collection.
            collection_novel[category_value] = target_values
        else:
            # Collection does not include the category's value.
            # Include category's value and target's value in collection.
            collection_novel[category_value] = [target_value]
        return collection_novel

    collection = {}
    for record in records:
        target_value = record[target]
        category_values = record[category]
        if isinstance(category_values, list):
            for category_value in category_values:
                collection = collect_record_target_by_category(
                    target_value=target_value,
                    category_value=category_value,
                    collection_original=collection
                )
        else:
            category_value = category_values
            collection = collect_record_target_by_category(
                target_value=target_value,
                category_value=category_value,
                collection_original=collection
            )
    return collection


def collect_values_from_records_in_reference(
    key=None, identifiers=None, reference=None
):
    """
    Collects a single value from a specific record in a reference.

    arguments:
        key (str): key of value in record
        identifiers (list<str>): identifiers of records in reference
        reference (dict<dict<str>>): reference of records

    raises:

    returns:
        (list<str>): values from records

    """

    values = []
    for identifier in identifiers:
        record = reference[identifier]
        value = record[key]
        values.append(value)
    return values


def filter_nonempty_elements(elements_original=None):
    """
    Filters nonempty elements.

    arguments:
        elements_original (list<str>): sequence of elements

    returns:
        (list<str>): non-empty elements

    raises:

    """

    elements_novel = []
    for element in elements_original:
        if len(str(element)) > 0:
            elements_novel.append(element)
    return elements_novel


def filter_entries_identifiers(
    identifiers=None,
    entries_original=None
):
    """
    Filters nodes and links by identifiers.

    arguments:
        identifiers (list<str>): identifiers of elements to keep
        entries_original (dict<dict>): entries

    raises:

    returns:
        (dict<dict>): entries

    """

    entries_novel = {}
    for entry in entries_original.values():
        if entry["identifier"] in identifiers:
            entries_novel[entry["identifier"]] = entry
    return entries_novel


def calculate_standard_score(
    value=None,
    mean=None,
    deviation=None,
):
    """
    Calculates the standard score, z-score, of a value.

    arguments:
        value (float): value to transform to standard score space
        mean (float): mean of distribution from which value comes
        deviation (float): standard deviation of distribution from which value
            comes

    raises:

    returns:
        (float): value in standard score space

    """

    if deviation != 0:
        return ((value - mean) / deviation)
    else:
        return float("nan")


def calculate_standard_scores(
    values=None,
    mean=None,
    deviation=None,
):
    """
    Calculates the standard scores, z-scores, of values.

    arguments:
        value (list<float>): values to transform to standard score space
        mean (float): mean of distribution from which value comes
        deviation (float): standard deviation of distribution from which value
            comes

    raises:

    returns:
        (list<float>): value in standard score space

    """

    values_standard = list()
    for value in values:
        value_standard = calculate_standard_score(
            value=value,
            mean=mean,
            deviation=deviation,
        )
        values_standard.append(value_standard)
    return values_standard


def combine_unique_elements_pairwise_orderless(
    elements=None,
):
    """
    Combines elements in orderless pairs.

    ABCD: AB AC AD BC BD CD

    arguments:
        elements (list): elements of interest

    returns:
        (list<tuple>): unique pairs of elements

    raises:

    """

    # Select unique elements.
    elements_unique = collect_unique_elements(elements_original=elements)
    # Combine elements in pairs.
    pairs = list(itertools.combinations(elements_unique, 2))
    # Return information.
    return pairs


def combine_unique_elements_pairwise_order(
    elements=None,
):
    """
    Combines elements in ordered pairs.

    ABCD: AB AC AD BA BC BD CA CB CD DA DB DC

    arguments:
        elements (list): elements of interest

    returns:
        (list<tuple>): unique pairs of elements

    raises:

    """

    # Select unique elements.
    elements_unique = collect_unique_elements(elements_original=elements)
    # Combine elements in pairs.
    pairs = list(itertools.permutations(elements_unique, 2))
    # Return information.
    return pairs


# Pandas


def convert_records_to_dataframe(records=None):
    """
    Converts information from list of dictionaries to Pandas data frame.

    arguments:
        records (list<dict>): records.

    raises:

    returns:
        (object): Pandas data frame.

    """

    return pandas.DataFrame(data=records)


def convert_dataframe_to_records(data=None):
    """
    Converts information from Pandas data frame to list of dictionaries.

    arguments:
        data (object): Pandas data frame.

    raises:

    returns:
        (list<dict>): records.

    """

    data_copy = data.copy(deep=True)
    data_copy.reset_index(
        level=None,
        inplace=True,
    )
    return data_copy.to_dict(orient="records")


def convert_table_columns_variables_types_float(
    columns=None,
    table=None,
):
    """
    Converts data variable types.

    arguments:
        columns (list<str>): names of columns
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # table[column] = pandas.DataFrame[column].astype(
    #    "float32",
    #    copy=True,
    #    errors="raise",
    #)

    # Copy information in table.
    table = table.copy(deep=True)
    # Convert data variable types.
    for column in columns:
        table[column] = pandas.to_numeric(
            table[column],
            errors="coerce", # force any parse error values to missing "NaN"
            downcast="float", # cast type to smallest float type
        )
    # Return information.
    return table


def segregate_data_two_thresholds(
    data=None,
    abscissa=None,
    ordinate=None,
    threshold_abscissa=None,
    selection_abscissa=None,
    threshold_ordinate=None,
    selection_ordinate=None,
):
    """
    Segregates data by values against thresholds on two dimensions.

    arguments:
        data (object): Pandas data frame of groups, series, and values
        abscissa (str): name of data column with independent variable
        ordinate (str): name of data column with dependent variable
        threshold_abscissa (float): threshold for abscissa
        selection_abscissa (str): selection criterion for abscissa's values
            against threshold
        threshold_ordinate (float): threshold for ordinate
        selection_ordinate (str): selection criterion for ordinate's values
            against threshold

    raises:

    returns:
        (dict): collection of data that pass and fail against thresholds

    """

    # Copy data.
    data_pass = data.copy(deep=True)
    data_fail = data.copy(deep=True)
    # Abscissa.
    if selection_abscissa == ">=":
        data_pass = data_pass.loc[
            data_pass[abscissa] >= threshold_abscissa,
            :
        ]
    elif selection_abscissa == "<=":
        data_pass = data_pass.loc[
            data_pass[abscissa] <= threshold_abscissa,
            :
        ]
    # Ordinate.
    if selection_ordinate == ">=":
        data_pass = data_pass.loc[
            data_pass[ordinate] >= threshold_ordinate,
            :
        ]
    elif selection_ordinate == "<=":
        data_pass = data_pass.loc[
            data_pass[ordinate] <= threshold_ordinate,
            :
        ]

    # Fail.
    data_fail = data_fail.loc[
        ~data_fail.index.isin(data_pass.index),
        :
    ]

    # Return information.
    collection = dict()
    collection["pass"] = data_pass
    collection["fail"] = data_fail
    return collection


def filter_rows_columns_by_threshold_proportion(
    data=None,
    dimension=None,
    threshold=None,
    proportion=None,
):
    """
    Filters either rows or columns by whether a certain proportion of values
    are greater than a minimal thresholds.

    Persistence of a row or column requires at least a specific proportion of
    values beyond a specific threshold.

    Filter rows by consideration of values across columns in each row.
    Filter columns by consideration of values across rows in each column.

    arguments:
        data (object): Pandas data frame of values
        dimension (str): dimension to filter, either "row" or "column"
        threshold (float): minimal signal
        proportion (float): minimal proportion of rows or columns that must
            pass threshold

    raises:

    returns:
        (object): Pandas data frame of values

    """

    def count_true(slice=None, count=None):
        values = slice.tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Determine count from proportion.
    if dimension == "row":
        # Filter rows by consideration of columns for each row.
        columns = data.shape[1]
        count = round(proportion * columns)
    elif dimension == "column":
        # Filter columns by consideration of rows for each columns.
        rows = data.shape[0]
        count = round(proportion * rows)

    # Determine whether values exceed threshold.
    data_threshold = (data >= threshold)
    # Determine whether count of values exceed threshold.
    if dimension == "row":
        axis = "columns"
    elif dimension == "column":
        axis = "index"
    # This aggregation operation produces a series.
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis=axis,
    )

    # Select rows and columns with appropriate values.
    if dimension == "row":
        data_pass = data.loc[data_count, : ]
    elif dimension == "column":
        data_pass = data.loc[:, data_count]
    return data_pass


def filter_rows_columns_by_threshold_outer_proportion(
    data=None,
    dimension=None,
    threshold_high=None,
    threshold_low=None,
    proportion=None,
):
    """
    Filters either rows or columns by whether a certain proportion of values
    are outside of low and high thresholds.

    Persistence of a row or column requires at least a specific proportion of
    values beyond a specific threshold.

    Filter rows by consideration of values across columns in each row.
    Filter columns by consideration of values across rows in each column.

    arguments:
        data (object): Pandas data frame of values
        dimension (str): dimension to filter, either "row" or "column"
        threshold_high (float): value must be greater than this threshold
        threshold_low (float): value must be less than this threshold
        proportion (float): minimal proportion of rows or columns that must
            pass threshold

    raises:

    returns:
        (object): Pandas data frame of values

    """

    def count_true(slice=None, count=None):
        values = slice.values.tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Copy data.
    data = data.copy(deep=True)

    # Determine count from proportion.
    if dimension == "row":
        # Filter rows by consideration of columns for each row.
        columns = data.shape[1]
        count = round(proportion * columns)
    elif dimension == "column":
        # Filter columns by consideration of rows for each columns.
        rows = data.shape[0]
        count = round(proportion * rows)

    # Determine whether values exceed threshold.
    data_threshold = data.applymap(
        lambda value: ((value <= threshold_low) or (threshold_high <= value))
    )
    # Determine whether count of values exceed threshold.
    if dimension == "row":
        axis = "columns"
    elif dimension == "column":
        axis = "index"
    # This aggregation operation produces a series.
    data_count = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis=axis,
    )

    # Select rows and columns with appropriate values.
    if dimension == "row":
        data_pass = data.loc[data_count, : ]
    elif dimension == "column":
        data_pass = data.loc[:, data_count]
    return data_pass


def calculate_table_false_discovery_rates(
    threshold=None,
    probability=None,
    discovery=None,
    significance=None,
    table=None,
):
    """
    Calculates false discovery rates (FDRs) from probabilities.

    arguments:
        threshold (float): value of alpha, or family-wise error rate of false
            discoveries
        probability (str): name of table's column of probabilities
        discovery (str): name for table's column of false discovery rates
        significance (str): name for table's column of FDR significance
        table (object): Pandas data frame with column of probabilities across
            observations in rows

    raises:

    returns:
        (object): Pandas data frame of probabilities and false discovery rates

    """

    # Copy information.
    table = table.copy(deep=True)

    # False discovery rate method cannot accommodate missing values.
    # Remove null values.
    table_null_boolean = pandas.isna(table[probability])
    table_null = table.loc[table_null_boolean]
    table_valid = table.dropna(
        axis="index",
        how="any",
        subset=[probability],
        inplace=False,
    )
    # Calculate false discovery rates from probabilities.
    probabilities = table_valid[probability].to_numpy()
    if len(probabilities) > 3:
        report = statsmodels.stats.multitest.multipletests(
            probabilities,
            alpha=threshold,
            method="fdr_bh", # use Benjamini-Hochberg False Discovery Rate (FDR)
            is_sorted=False,
        )
        significances = report[0] # valid to reject null hypothesis
        #significances = numpy.invert(rejects)
        discoveries = report[1]
        table_valid[significance] = significances
        table_valid[discovery] = discoveries
    else:
        table_valid[significance] = float("nan")
        table_valid[discovery] = float("nan")
        pass
    table_null[significance] = False
    table_null[discovery] = float("nan")

    # Combine null and valid portions of data.
    table_discoveries = table_valid.append(
        table_null,
        ignore_index=False,
    )
    # Return information.
    return table_discoveries


def impute_continuous_variables_missing_values_half_minimum(
    columns=None,
    table=None,
    report=None,
):
    """
    Imputes missing values to half of their minimum value.

    arguments:
        columns (list<str>): names of columns for continuous, ratio-scale
            variables for which to impute missing values
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame

    """

    # Copy table.
    table = table.copy(deep=True)
    # Convert data variable types.
    for column in columns:
        # Determine name for new column.
        column_imputation = str(column + "_imputation")
        # Determine half minimum value.
        minimum = numpy.nanmin(table[column].to_numpy())
        half_minimum = float(minimum / 2)
        # Report.
        if report:
            print_terminal_partition(level=2)
            print(
                "report: " +
                "impute_continuous_variables_missing_values_half_minimum()"
            )
            print_terminal_partition(level=3)
            print("column: " + str(column))
            print("minimum: " + str(minimum))
            print("half minimum: " + str(half_minimum))
        # Impute missing values to half minimum.
        table[column_imputation] = table[column].fillna(
            value=half_minimum,
            inplace=False,
        )
    # Return information.
    return table


def calculate_relative_variance(
    array=None,
):
    """
    Calculates the relative variance.

    Variance is mathematically never negative.
    Use the absolute value of the mean.

    arguments:
        array (object): NumPy array of ratio-scale values

    raises:

    returns:
        (float): value of relative variance

    """

    mean = numpy.nanmean(array)
    mean_absolute = math.fabs(mean)
    variance = numpy.nanvar(array)
    if (mean_absolute > 0.0000001):
        relative_variance = (variance / mean_absolute)
    else:
        relative_variance = float("nan")
    # Return information.
    return relative_variance


def count_data_factors_groups_elements(
    factors=None,
    element=None,
    count=None,
    data=None,
):
    """
    Counts elements in groups by factor.

    arguments:
        factors (list<str>): names of factor columns in data
        element (str): name of column in data for elements to count
        count (str): name for column of counts in count data
        data (object): Pandas data frame of elements in groups by factors

    raises:

    returns:
        (object): Pandas data frame of counts of elements in each factor group
    """

    # Copy data.
    data_copy = data.copy(deep=True)

    # Organize data.
    data_copy.reset_index(
        level=None,
        inplace=True,
    )
    columns = factors
    columns = copy.deepcopy(factors)
    columns.append(element)
    data_columns = data_copy.loc[
        :, data_copy.columns.isin(columns)
    ]
    data_columns.drop_duplicates(
        subset=None,
        keep="first",
        inplace=True,
    )
    data_columns.set_index(
        factors,
        append=False,
        drop=True,
        inplace=True
    )
    # Count rows (elements) in each factor group of the data.
    # This process is similar to iteration on groups and collection of the
    # counts of elements in each group.
    # Notice the use of ".size()" instead of ".shape[0]" to count rows.
    # As the groups are series, ".size()" counts the elements properly.
    data_counts = data_columns.groupby(
        level=factors,
        sort=True,
        as_index=False
    ).size().to_frame(
        name=count
    )
    data_counts.reset_index(
        level=None,
        inplace=True,
    )
    return data_counts


def report_contingency_table_stratification_by_missingness(
    column_stratification=None,
    stratifications=None,
    column_missingness=None,
    table=None,
    report=None,
):
    """
    Organizes and reports information about a contingency table between two
    values of two variables.

    The primary variable is the priority for stratification.
    The secondary variable has stratification by whether values are missing or
    non-missing.

    arguments:
        column_stratification (str): name of column for primary stratification
        stratifications (list): two values of the column by which to
            stratify the table's rows
        column_missingness (str): name of column for which to stratify valid
            and null or missing values
        table (object): Pandas data frame of features (columns) across
            observations (rows)
        report (bool): whether to print reports

    raises:

    returns:


    """

    # Copy information.
    table = table.copy(deep=True)
    # Reset index.
    table.reset_index(
        level=None,
        inplace=True
    )
    # Reduce the table to the relevant columns.
    table = table.loc[
        :, table.columns.isin([
            column_stratification, column_missingness,
        ])
    ]
    # Define binary representation of the two relevant values of the primary
    # stratification column.
    column_primary = str(column_stratification + "_relevant")
    table[column_primary] = table[column_stratification].apply(
        lambda value:
            0 if (value == stratifications[0]) else
            (1 if (value == stratifications[1]) else
            (float("nan")))
    )
    # Define binary representation of missingness in the secondary column.
    column_secondary = str(column_missingness + "_missing")
    table[column_secondary] = table[column_missingness].apply(
        lambda value:
            0 if (not pandas.isna(value)) else 1
    )
    # Remove any rows with missing values in the primary or secondary columns.
    table.dropna(
        axis="index",
        how="any",
        subset=[column_primary, column_secondary],
        inplace=True,
    )
    # Contingency table.
    table_contingency = pandas.crosstab(
        table[column_primary],
        table[column_secondary],
        rownames=[column_primary],
        colnames=[column_secondary],
    )
    # Chi-square test.
    (chi2, probability, freedom, expectation) = scipy.stats.chi2_contingency(
        table_contingency.to_numpy(),
        correction=True,
    )
    # Percentages.
    count_total = table.shape[0]
    count_0_0 = table_contingency.to_numpy()[0][0]
    count_0_1 = table_contingency.to_numpy()[0][1]
    count_1_0 = table_contingency.to_numpy()[1][0]
    count_1_1 = table_contingency.to_numpy()[1][1]
    percentage_0_0 = round(float((count_0_0 / count_total) * 100), 2)
    percentage_0_1 = round(float((count_0_1 / count_total) * 100), 2)
    percentage_1_0 = round(float((count_1_0 / count_total) * 100), 2)
    percentage_1_1 = round(float((count_1_1 / count_total) * 100), 2)
    entry_0_0 = str(str(count_0_0) + " (" + str(percentage_0_0) + "%)")
    entry_0_1 = str(str(count_0_1) + " (" + str(percentage_0_1) + "%)")
    entry_1_0 = str(str(count_1_0) + " (" + str(percentage_1_0) + "%)")
    entry_1_1 = str(str(count_1_1) + " (" + str(percentage_1_1) + "%)")
    name_missing_false = str(column_missingness + "-valid")
    name_missing_true = str(column_missingness + "-missing")
    entries = dict()
    entries[column_stratification] = stratifications
    entries[name_missing_false] = [entry_0_0, entry_1_0]
    entries[name_missing_true] = [entry_0_1, entry_1_1]
    table_report = pandas.DataFrame(data=entries)
    table_report.set_index(
        column_stratification,
        append=False,
        drop=True,
        inplace=True
    )
    # Report.
    if report:
        print_terminal_partition(level=2)
        print(
            "Contingency table and Chi2 test for independence."
        )
        print(str(
            column_stratification + " values: " +
            str(stratifications[0]) + ", " + str(stratifications[1])
        ))
        print("versus")
        print(str(column_missingness + " missingness"))
        print_terminal_partition(level=4)
        print(table_report)
        #print(table_contingency)
        #print(table_contingency.to_numpy())
        print_terminal_partition(level=5)
        print("chi2: " + str(chi2))
        print("probability: " + str(probability))
    pass


def filter_table_columns_by_nonmissing_variance(
    threshold_valid_proportion_per_column=None,
    threshold_column_variance=None,
    type_variance=None,
    table=None,
    report=None,
):
    """
    Filters columns in a table by their proportions of non-missing values and by
    their relative varances across rows.

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    arguments:
        threshold_valid_proportion_per_column (float): minimal proportion of
            a column's rows that must have a valid value
        threshold_column_variance (float): minimal value of variance measure in
            each column's values
        type_variance (str): type of measurement of variance, either
            "standard_deviation" or "relative_variance"
        table (object): Pandas data frame with variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information for the singular value decomposition

    """

    # Nested function.
    def match_column_variance(
        name=None,
        threshold_column_variance=None,
        variances=None,
    ):
        if (str(name) in variances.keys()):
            variance = variances[name]
            if (
                (not math.isnan(variance)) and
                (threshold_column_variance <= variance)
            ):
                match = True
            else:
                match = False
        else:
            match = False
        return match

    # Copy information.
    table_source = table.copy(deep=True)
    table_product = table.copy(deep=True)
    # Drop any columns with inadequate valid (non-missing) values across rows.
    if (0.001 <= threshold_valid_proportion_per_column):
        rows = table_product.shape[0]
        threshold = round(rows*threshold_valid_proportion_per_column)
        table_product.dropna(
            axis="columns", # drop columns
            thresh=threshold,
            subset=None,
            inplace=True,
        )
    # Drop any columns with minimal relative variance.
    if (0.00001 <= threshold_column_variance):
        if (type_variance == "standard_deviation"):
            series_variance = table_product.aggregate(
                lambda column: numpy.nanstd(column.to_numpy()),
                axis="index", # apply function to each column
            )
            #series_variance = table_product.aggregate(
            #    lambda column: column.std(),
            #    axis="index", # apply function to each column
            #)
        elif (type_variance == "relative_variance"):
            series_variance = table_product.aggregate(
                lambda column: calculate_relative_variance(
                    array=column.to_numpy()
                ),
                axis="index", # apply function to each column
            )
            pass
        variances = series_variance.to_dict()
        columns = copy.deepcopy(table.columns.to_list())
        columns_variance = list(filter(
            lambda column_trial: match_column_variance(
                name=column_trial,
                threshold_column_variance=(
                    threshold_column_variance
                ),
                variances=variances,
            ),
            columns
        ))
        table_product = table_product.loc[
            :, table_product.columns.isin(columns_variance)
        ]
    # Determine any columns removed.
    columns_exclusion = list(filter(
        lambda column: (str(column) not in table_product.columns.tolist()),
        table_source.columns.tolist()
    ))

    # Report.
    if report:
        print_terminal_partition(level=2)
        print(
            "Report from: " +
            "filter_table_columns_by_nonmissing_relative_variance()"
        )
        print_terminal_partition(level=3)
        print("count rows in source table: " + str(table_source.shape[0]))
        print("count columns in source table: " + str(table_source.shape[1]))
        print("series variance: " + str(type_variance))
        print(series_variance.iloc[0:25])
        print("count rows in product table: " + str(table_product.shape[0]))
        print("count columns in product table: " + str(table_product.shape[1]))
        print("any columns removed from table: ")
        print(columns_exclusion)
    # Return.
    return table_product


def extract_first_search_string_from_table_column_main_string(
    string_source=None,
    search_strings_1=None,
):
    """
    Searches a longer main character string for any matches to shorter search
    strings.

    Example main string: 'female_premenopause_unadjust_oestradiol'
    Example tier 1 search string: 'female_premenopause'
    Example tier 2 search string: 'female'

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    arguments:
        string_source (str): a longer character main string from which to extract
            any of the shorter character search strings
        search_strings_1 (list<str>): first tier of shorter, search character
            strings for which to search before extraction of second tier strings

    raises:

    returns:
        (str): match string from first search

    """

    # Initialize match string.
    match_string_1 = ""
    # Iterate on search strings.
    for search_string in search_strings_1:
        if (search_string in str(string_source)):
            match_string_1 = copy.deepcopy(search_string)
            pass
        pass
    # Return information.
    return match_string_1


def extract_second_search_string_from_table_column_main_string(
    string_source=None,
    match_string_1=None,
    search_strings_2=None,
):
    """
    Searches a longer main character string for any matches to shorter search
    strings.

    Example main string: 'female_premenopause_unadjust_oestradiol'
    Example tier 1 search string: 'female_premenopause'
    Example tier 2 search string: 'female'

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    arguments:
        string_source (str): a longer character main string from which to extract
            any of the shorter character search strings
        match_string_1 (str): any match from previous search on first tier
            search strings
        search_strings_2 (list<str>): second tier of shorter, search character
            strings for which to search after extraction of first tier strings

    raises:

    returns:
        (str): match string from second search

    """

    # Copy information.
    string_source = copy.deepcopy(string_source)
    match_string_1 = copy.deepcopy(match_string_1)
    # Find the remainder of the longer character main string after extraction of
    # any matches from previous search on first tier search strings.
    # Initialize remainder string.
    remainder_string = copy.deepcopy(string_source)
    if (len(match_string_1) > 0):
        remainder_string = remainder_string.replace(str(match_string_1), "")
    # Initialize match string.
    match_string_2 = ""
    # Iterate on search strings.
    for search_string in search_strings_2:
        if (search_string in str(remainder_string)):
            match_string_2 = copy.deepcopy(search_string)
            pass
        pass
    # Return information.
    return match_string_2


def extract_third_search_string_from_table_column_main_string(
    string_source=None,
    match_string_1=None,
    match_string_2=None,
    search_strings_3=None,
):
    """
    Searches a longer main character string for any matches to shorter search
    strings.

    Example main string: 'female_premenopause_unadjust_oestradiol'
    Example tier 1 search string: 'female_premenopause'
    Example tier 2 search string: 'female'

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    arguments:
        string_source (str): a longer character main string from which to extract
            any of the shorter character search strings
        match_string_1 (str): any match from previous search on first tier
            search strings
        match_string_2 (str): any match from previous search on first tier
            search strings
        search_strings_3 (list<str>): second tier of shorter, search character
            strings for which to search after extraction of first tier strings

    raises:

    returns:
        (str): match string from second search

    """

    # Copy information.
    string_source = copy.deepcopy(string_source)
    match_string_1 = copy.deepcopy(match_string_1)
    match_string_2 = copy.deepcopy(match_string_2)
    # Find the remainder of the longer character main string after extraction of
    # any matches from previous search on first tier search strings.
    # Initialize remainder string.
    remainder_string = copy.deepcopy(string_source)
    if (len(match_string_1) > 0):
        remainder_string = remainder_string.replace(str(match_string_1), "")
    if (len(match_string_2) > 0):
        remainder_string = remainder_string.replace(str(match_string_2), "")
    # Initialize match string.
    match_string_3 = ""
    # Iterate on search strings.
    for search_string in search_strings_3:
        if (search_string in str(remainder_string)):
            match_string_3 = copy.deepcopy(search_string)
            pass
        pass
    # Return information.
    return match_string_3


def combine_first_second_third_search_string_matches(
    match_string_1=None,
    match_string_2=None,
    match_string_3=None,
):
    """
    Combines matches from first, second, and third searches on a main character
    string.

    arguments:
        match_string_1 (str): any match from previous search on first tier
            search strings
        match_string_2 (str): any match from previous search on second tier
            search strings
        match_string_3 (str): any match from previous search on third tier
            search strings

    raises:

    returns:
        (str): combination of match strings from first, second, and third
            searches

    """

    # Copy information.
    match_string_1 = copy.deepcopy(match_string_1)
    match_string_2 = copy.deepcopy(match_string_2)
    match_string_3 = copy.deepcopy(match_string_3)
    # Determine how to combine matches from first and second searches.
    if (
        (len(match_string_1) > 0) and
        (len(match_string_2) == 0) and
        (len(match_string_3) == 0)
    ):
        combination_string = match_string_1
    elif (
        (len(match_string_1) == 0) and
        (len(match_string_2) > 0) and
        (len(match_string_3) == 0)
    ):
        combination_string = match_string_2
    elif (
        (len(match_string_1) == 0) and
        (len(match_string_2) == 0) and
        (len(match_string_3) > 0)
    ):
        combination_string = match_string_3
    elif (
        (len(match_string_1) > 0) and
        (len(match_string_2) > 0) and
        (len(match_string_3) == 0)
    ):
        combination_string = str(match_string_1 + "_" + match_string_2)
    elif (
        (len(match_string_1) > 0) and
        (len(match_string_2) == 0) and
        (len(match_string_3) > 0)
    ):
        combination_string = str(match_string_1 + "_" + match_string_3)
    elif (
        (len(match_string_1) == 0) and
        (len(match_string_2) > 0) and
        (len(match_string_3) > 0)
    ):
        combination_string = str(match_string_2 + "_" + match_string_3)
    elif (
        (len(match_string_1) > 0) and
        (len(match_string_2) > 0) and
        (len(match_string_3) > 0)
    ):
        combination_string = str(
            match_string_1 + "_" + match_string_2 + "_" + match_string_3
        )
    else:
        combination_string = ""
    # Return information.
    return combination_string


def extract_overlap_search_strings_from_table_column_main_string(
    column_source=None,
    column_target=None,
    temporary_column_prefix=None,
    search_strings_1=None,
    search_strings_2=None,
    search_strings_3=None,
    table=None,
    report=None,
):
    """
    Searches longer main character strings for shorter search strings in two
    separate steps to avoid problems with search strings that contain each
    other.

    Example main string: 'female_premenopause_unadjust_oestradiol_bioavailable'
    Example 'cohort' tier 1 search string: 'female_premenopause'
    Example 'cohort' tier 2 search string: 'female'

    For the search and extraction to be effective, the longer, main character
    string must contain only a single instance of a single shorter, search
    character string in either the first or second search tiers.

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    arguments:
        column_source (str): name of the original table column that contains the
            longer, main character strings from which to extract the shorter
            character search strings
        column_target (str): name for novel table column in which to collect
            shorter strings after extraction from the longer strings in the
            original table column
        temporary_column_prefix (str): character string prefix for names of
            temporary table columns
        type_variance (str): type of measurement of variance, either
            "standard_deviation" or "relative_variance"
        search_strings_1 (list<str>): first tier of shorter, search character
            strings for which to search before extraction of second tier strings
        search_strings_2 (list<str>): second tier of shorter, search character
            strings for which to search after extraction of first tier strings
        search_strings_3 (list<str>): third tier of shorter, search character
            strings for which to search after extraction of first tier strings
        table (object): Pandas data-frame table with a column of longer
            character strings from which to extract multiple shorter character
            strings
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Copy information.
    table = table.copy(deep=True)
    # Define names of temporary columns for use in first and second searches.
    column_match_1 = str(temporary_column_prefix + "_match_1")
    column_match_2 = str(temporary_column_prefix + "_match_2")
    column_match_3 = str(temporary_column_prefix + "_match_3")
    # Search for and extract first tier strings.
    table[column_match_1] = table.apply(
        lambda row:
            extract_first_search_string_from_table_column_main_string(
                string_source=row[column_source],
                search_strings_1=search_strings_1,
            ),
        axis="columns", # apply function to each row
    )
    # Search for and extract second tier strings.
    table[column_match_2] = table.apply(
        lambda row:
            extract_second_search_string_from_table_column_main_string(
                string_source=row[column_source],
                match_string_1=row[column_match_1],
                search_strings_2=search_strings_2,
            ),
        axis="columns", # apply function to each row
    )
    # Search for and extract third tier strings.
    table[column_match_3] = table.apply(
        lambda row:
            extract_third_search_string_from_table_column_main_string(
                string_source=row[column_source],
                match_string_1=row[column_match_1],
                match_string_2=row[column_match_2],
                search_strings_3=search_strings_3,
            ),
        axis="columns", # apply function to each row
    )
    # Collect information from first and second searches.
    # For the search and extraction to be effective, there should only be a
    # single match from either the first or second search.
    # If there are matches from both the first and second searches, then it is
    # necessary to combine them.
    table[column_target] = table.apply(
        lambda row:
            combine_first_second_third_search_string_matches(
                match_string_1=row[column_match_1],
                match_string_2=row[column_match_2],
                match_string_3=row[column_match_3],
            ),
        axis="columns", # apply function to each row
    )
    # Remove temporary columns from first and second searches.
    table.drop(
        labels=[column_match_1, column_match_2, column_match_3,],
        axis="columns",
        inplace=True
    )
    # Report.
    if report:
        print_terminal_partition(level=2)
        print(
            "Report from: " +
            "extract_overlap_search_strings_from_table_column_main_string()"
        )
        print_terminal_partition(level=5)
    # Return.
    return table


def drive_extract_search_strings_from_table_columns_main_strings(
    pail_extractions=None,
    temporary_column_prefix=None,
    column_source=None,
    table=None,
    report=None,
):
    """
    Drives two-tier searches and extractions on longer main character strings
    within a single column of table.

    Example main string: 'female_premenopause_unadjust_oestradiol_bioavailable'
    Example 'cohort' tier 1 search string: 'female_premenopause'
    Example 'cohort' tier 2 search string: 'female'
    Example 'cohort' tier 3 search string: 'male'
    Example 'phenotype' tier 1 search string: 'oestradiol_bioavailable'
    Example 'phenotype' tier 2 search string: 'oestradiol'

    Table format: Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows with an explicit index

    arguments:
        pail_extractions (dict<dict<list<str>>>): collection of tier 1 and tier
            2 shorter, search character strings across multiple variables of
            interest
        temporary_column_prefix (str): character string prefix for names of
            temporary table columns
        column_source (str): name of the original table column that contains the
            longer, main character strings from which to extract the shorter
            character search strings
        table (object): Pandas data-frame table with a column of longer
            character strings from which to extract multiple shorter character
            strings
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Iterate on variables for extraction.
    for variable in list(pail_extractions.keys()):
        table = extract_overlap_search_strings_from_table_column_main_string(
            column_source=column_source,
            column_target=str(variable),
            temporary_column_prefix=temporary_column_prefix,
            search_strings_1=pail_extractions[variable]["search_1"],
            search_strings_2=pail_extractions[variable]["search_2"],
            search_strings_3=pail_extractions[variable]["search_3"],
            table=table,
            report=report,
        )
        pass
    # Report.
    if report:
        print_terminal_partition(level=2)
        print(
            "Report from: " +
            "drive_extract_search_strings_from_table_columns_" +
            "main_strings()"
        )
        print_terminal_partition(level=5)
    # Return.
    return table


def calculate_table_column_pair_correlations(
    column_one=None,
    column_two=None,
    table=None,
    report=None,
):
    """
    Calculates correlations between a pair of columns in a table.

    arguments:
        column_one (str): name of first column
        column_two (str): name of second column
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about correlations



    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Convert columns to float.
    table = convert_table_columns_variables_types_float(
        columns=[column_one, column_two],
        table=table,
    )
    # Remove table rows with missing values in relevant columns.
    table.dropna(
        axis="index",
        how="any",
        subset=[column_one, column_two],
        inplace=True,
    )
    # Calculate correlations.
    correlation_pearson, probability_pearson = scipy.stats.pearsonr(
        table[column_one].to_numpy(),
        table[column_two].to_numpy(),
    )
    correlation_spearman, probability_spearman = scipy.stats.spearmanr(
        table[column_one].to_numpy(),
        table[column_two].to_numpy(),
    )
    correlation_kendall, probability_kendall = scipy.stats.kendalltau(
        table[column_one].to_numpy(),
        table[column_two].to_numpy(),
    )
    # Report.
    if report:
        print_terminal_partition(level=4)
        print("Correlations between pair of columns")
        print("column one: " + str(column_one))
        print("column_two: " + str(column_two))
        print("valid value pairs for correlation: " + str(table.shape[0]))
        print("Pearson correlation: " + str(correlation_pearson))
        print("Pearson probability: " + str(probability_pearson))
        print("Spearman correlation: " + str(correlation_spearman))
        print("Spearman probability: " + str(probability_spearman))
        print("Kendall correlation: " + str(correlation_kendall))
        print("Kendall probability: " + str(probability_kendall))
        pass
    # Collect information.
    pail = dict()
    pail["pairs"] = copy.deepcopy(int(table.shape[0]))
    pail["correlation_pearson"] = correlation_pearson
    pail["probability_pearson"] = probability_pearson
    pail["correlation_spearman"] = correlation_spearman
    pail["probability_spearman"] = probability_spearman
    pail["correlation_kendall"] = correlation_kendall
    pail["probability_kendall"] = probability_kendall
    # Return information.
    return pail


def drive_calculate_table_column_pair_correlations(
    entries_cohorts=None,
    name_one=None,
    name_two=None,
    records_comparisons=None,
    report=None,
):
    """
    Drives the calculation of Pearson, Spearman, and Kendall correlations
    between pairs of variables (columns) within tables representing
    stratification cohorts. Organizes information from these correlations within
    a summary table.

    arguments:
        entries_cohorts (dict<dict<object>>): information about variables within
            Pandas data frame tables that represent stratification cohorts
        name_one (str): common name for the first variable in comparisons
        name_two (str): common name for the second variable in comparisons
        records_comparisons (list<dic<str>>): information for comparisons
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about correlations

    """

    # Collect records of information about correlations between variables within
    # cohorts.
    records_correlations = list()
    # Calculate and report correlations between variables within cohorts.
    for record_comparison in records_comparisons:
        record = dict()
        record["cohort"] = record_comparison["cohort"]
        record[name_one] = record_comparison["one"]
        record[name_two] = record_comparison["two"]
        pail = calculate_table_column_pair_correlations(
            column_one=record[name_one],
            column_two=record[name_two],
            table=entries_cohorts[record["cohort"]]["table"],
            report=report,
        )
        record.update(pail)
        records_correlations.append(record)
        pass

    # Organize table.
    table = pandas.DataFrame(data=records_correlations)
    # Select columns.
    # columns.insert(0, dependence)
    columns = [
        "cohort",
        name_one, name_two,
        "pairs",
        "correlation_pearson", "probability_pearson",
        "correlation_spearman", "probability_spearman",
        "correlation_kendall", "probability_kendall",
    ]
    table = table.loc[:, table.columns.isin(columns)]
    # Sort columns.
    table = table[[*columns]]
    # Report.
    if report:
        print_terminal_partition(level=4)
        print("report: ")
        print("drive_calculate_table_column_pair_correlations()")
        print_terminal_partition(level=5)
        print(table)
        pass
    # Return information.
    return table



##########################
# Basic operations on Pandas dataframe tables
###############################


def merge_columns_tables_supplements_to_main(
    identifier_main=None,
    identifier_supplement=None,
    table_main=None,
    tables_supplements=None,
    report=None,
):
    """
    Merge columns from multiple source tables to columns of a single main table.

    This function does not discard any existing indices of the main and
    supplement tables. This behavior might produce redundant or unnecessary
    columns.

    arguments:
        identifier_main (str): name of column in main table on which
            to merge
        identifier_supplement (str): name of column in supplement tables on
            which to merge
        table_main (object): Pandas data frame of information
        tables_supplements (list<object>): collection of Pandas data frames of
            information
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information

    """

    # Copy information in table.
    table_main = table_main.copy(deep=True)

    # Reset index in original main table.
    table_main.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )

    # Collect columns in original main table.
    # This collection must include the index.
    columns_collection = copy.deepcopy(table_main.columns.to_list())

    # Organize main table's index.
    table_main[identifier_main] = table_main[identifier_main].astype(
        "string",
        copy=True,
        errors="raise",
    )
    table_main.set_index(
        identifier_main,
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True,
    )

    # Copy information in table.
    table_merge = table_main.copy(deep=True)

    # Iterate on tables for polygenic scores.
    for table_supplement in tables_supplements:
        # Reset index in original supplement table.
        table_supplement.reset_index(
            level=None,
            inplace=True,
            drop=False, # remove index; do not move to regular columns
        )
        # Collect columns in original supplement table.
        # This collection must include the index.
        columns_supplement = copy.deepcopy(table_supplement.columns.to_list())
        columns_collection.extend(columns_supplement)
        # Organize supplement table's index.
        table_supplement[identifier_supplement] = (
            table_supplement[identifier_supplement].astype(
                "string",
                copy=True,
                errors="raise",
        ))
        table_supplement.set_index(
            identifier_supplement,
            append=False,
            drop=True, # move regular column to index; remove original column
            inplace=True
        )
        # Merge data tables using database-style join.
        # Alternative is to use DataFrame.join().
        table_merge = pandas.merge(
            table_merge, # left table
            table_supplement, # right table
            left_on=None, # "identifier_main"
            right_on=None, # "identifier_supplement"
            left_index=True, # "identifier_main"
            right_index=True, # "identifier_supplement"
            how="outer", # keep union of keys from both tables
            #suffixes=("_main", "_identifiers"), # deprecated?
        )
        pass
    # Organize table's index.
    table_merge.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    #table_merge.set_index(
    #    "identifier_genotype",
    #    append=False,
    #    drop=True, # move regular column to index; remove original column
    #    inplace=True
    #)
    # Determine unique columns from collection.
    columns_collection = list(set(columns_collection))
    # Remove unnecessary columns from transformations on tables.
    # The goal is to avoid any columns from the merge operations themselves.
    #table_merge.drop(
    #    labels=["index_x", "index_y", "index",],
    #    axis="columns",
    #    inplace=True
    #)
    table_merge = table_merge.loc[
        :, table_merge.columns.isin(columns_collection)
    ]
    # Report.
    if report:
        print_terminal_partition(level=2)
        print("report: ")
        print("merge_columns_tables_supplements_to_main()")
        print_terminal_partition(level=3)
        print("table columns: " + str(int(table_merge.shape[1])))
        print("table rows: " + str(int(table_merge.shape[0])))
        print("columns")
        print(table_merge.columns.to_list())
        print(table_merge)
        pass
    # Return information.
    return table_merge


def merge_columns_two_tables(
    identifier_first=None,
    identifier_second=None,
    table_first=None,
    table_second=None,
    report=None,
):
    """
    Merge columns from two tables.

    This function does not discard any existing indices of the main and
    supplement tables. This behavior might produce redundant or unnecessary
    columns.

    arguments:
        identifier_first (str): name of column in first table on which to merge
        identifier_second (str): name of column in second table on which to
            merge
        table_first (object): Pandas data frame of information
        tables_second (object): Pandas data frame of information
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information

    """

    # Copy information in table.
    table_first = table_first.copy(deep=True)
    table_second = table_second.copy(deep=True)

    # Reset indices in original tables.
    table_first.reset_index(
        level=None,
        inplace=True,
        drop=False, # do not remove index; move to regular columns
    )
    table_second.reset_index(
        level=None,
        inplace=True,
        drop=False, # do not remove index; move to regular columns
    )

    # Collect columns in original tables.
    # This collection must include the indices.
    columns_collection = copy.deepcopy(table_first.columns.to_list())
    columns_second = copy.deepcopy(table_second.columns.to_list())
    columns_collection.extend(columns_second)

    # Organize first table's index.
    table_first[identifier_first] = table_first[identifier_first].astype(
        "string",
        copy=True,
        errors="raise",
    )
    table_first.set_index(
        identifier_first,
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True,
    )
    # Organize second table's index.
    table_second[identifier_second] = table_second[identifier_second].astype(
        "string",
        copy=True,
        errors="raise",
    )
    table_second.set_index(
        identifier_second,
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True,
    )

    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table = pandas.merge(
        table_first, # left table
        table_second, # right table
        left_on=None, # "identifier_first"
        right_on=None, # "identifier_second"
        left_index=True, # "identifier_first"
        right_index=True, # "identifier_second"
        how="outer", # keep union of keys from both tables
        #suffixes=("_main", "_identifiers"), # deprecated?
    )

    # Organize table's index.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    #table.set_index(
    #    "identifier_genotype",
    #    append=False,
    #    drop=True, # move regular column to index; remove original column
    #    inplace=True
    #)
    # Determine unique columns from collection.
    columns_collection = list(set(columns_collection))
    # Remove unnecessary columns from transformations on tables.
    # The goal is to avoid any columns from the merge operations themselves.
    #table_merge.drop(
    #    labels=["index_x", "index_y", "index",],
    #    axis="columns",
    #    inplace=True
    #)
    table = table.loc[
        :, table.columns.isin(columns_collection)
    ]

    # Report.
    if report:
        print_terminal_partition(level=2)
        print("report: ")
        print("merge_columns_two_tables()")
        print_terminal_partition(level=3)
        print("table columns: " + str(int(table.shape[1])))
        print("table rows: " + str(int(table.shape[0])))
        print("columns")
        print(table.columns.to_list())
        print(table)
        pass
    # Return information.
    return table


def report_table_column_categories_rows(
    column=None,
    table=None,
):
    """
    Reports the total count of unique, categorical values in a single column of
    a table and also reports counts of rows (samples) with each value.

    These functions do not ignore missing values.

    arguments:
        column (str): name of column in table with categorical values
        table (object): Pandas data frame table of variables (features) across
            columns and samples (records) across rows

    raises:

    returns:

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Count unique, nonmissing categorical values in the table's column.
    # table[column].nunique(dropna=False)
    count_unique = table[column].unique().size
    values_unique = table[column].unique().tolist()
    # Count rows (samples) with each categorical value of the column (variable).
    # This represents the frequency of each categorical value.
    frequencies = table[column].value_counts(
        normalize=False,
        sort=True,
        ascending=False,
        bins=None,
        dropna=False, # whether to ignore missing values
    )
    # Report.
    print_terminal_partition(level=3)
    print("report: ")
    name_function = (
        "report_table_column_categories_rows()"
    )
    print(name_function)
    print_terminal_partition(level=4)
    print("column: " + str(column))
    print("count of unique values in column: " + str(count_unique))
    print("unique values:")
    print(values_unique)
    print("frequencies of values in column:")
    print(frequencies)
    pass


##########
# Processes on simple Records (lists of dictionaries)
# Useful for management of multiple stratification cohorts.


def report_records_name_table_size(
    records=None,
):
    """
    Reports basic information from records of stratification cohorts, including
    the cohort 'name' and the counts of columns and rows in the cohort 'table'.

    At a minimum, records include a name (key: "name") and a table
    (key: "table").

    arguments:
        records (list<dict>): records with information about cohorts

    raises:

    returns:

    """

    # Report.
    print_terminal_partition(level=3)
    print("report: ")
    print("report_records_name_table_size()")
    print_terminal_partition(level=4)
    # Copy information.
    records = copy.deepcopy(records)
    # Iterate on records for stratification cohorts.
    for record in records:
        print_terminal_partition(level=5)
        print("table name: " + str(record["name"]))
        print("table columns: " + str(int(record["table"].shape[1])))
        print("table rows: " + str(int(record["table"].shape[0])))
        pass
    pass


def structure_from_records_to_entries(
    records=None,
):
    """
    Structures information from an original list of dictionary records to a
    novel dictionary with an entry for each original record. The name
    (key: "name") from each original record becomes the key for the record's
    novel entry.

    At a minimum, original records include a name (key: "name").

    arguments:
        records (list<dict>): records with information about cohorts

    raises:

    returns:
        (dict<dict>): entries with information about cohorts

    """

    # Organize dictionary entries for cohorts.
    entries = dict()
    for record in records:
        # Copy information.
        entries[str(record["name"])] = copy.deepcopy(record)
        pass
    # Return information
    return entries


def filter_records_by_name(
    names=None,
    records=None,
    report=None,
):
    """
    Filters records by name (key: "name").

    At a minimum, records include a name (key: "name").

    arguments:
        names (list<str>): names of records to keep
        records (list<dict>): records
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): records

    """

    # Copy information.
    records = copy.deepcopy(records)
    # Filter records by name.
    records_filter = list(filter(
        lambda record: (str(record["name"]) in names),
        records
    ))
    # Report.
    if report:
        print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "filter_records_by_name()"
        )
        print(name_function)
        print_terminal_partition(level=3)
        count_records = int(len(records))
        count_records_filter = int(len(records_filter))
        print("count of original records: " + str(count_records))
        print("count of novel records: " + str(count_records_filter))
        pass
    # Return information
    return records_filter





# Stratifications by continuous variables


################### WORK ZONE ###########################


# Principal components

# TODO: switch to statsmodels PCA
# https://www.statsmodels.org/devel/generated/statsmodels.multivariate.pca.pca.html


def calculate_principal_components_statsmodels(
    table=None,
    report=None,
):
    """
    Calculates the principal components using the statsmodels package.

    Format of data should have features across columns and observations across
    rows.

    This function then transforms the data frame to a matrix.

    This matrix has observations across dimension zero and features across
    dimension one.

    arguments:
        table (object): Pandas data frame of variables (features) across
            columns and samples (cases, observations) across rows with an
            explicit index
        report (bool): whether to print reports to terminal

    raises:

    returns:
        (dict): information about data's principal components

    """

    # Copy information.
    table = table.copy(deep=True)
    index = copy.deepcopy(table.index)
    # Organize matrix.
    # Matrix format has samples (cases, observations) across dimension 0 (rows)
    # and variables (features) across dimension 1 (columns).
    matrix = table.to_numpy()
    # Calculate Principal Components.
    # If there is specification of "ncomp", then function only returns
    # Eigenvalues, loadings, Eigenvectors, and principal components to this
    # count.
    # Function sorts Eigenvalues in decreasing order.
    # Sort order of Eigenvectors mush match the sort order of Eigenvalues.
    # Statsmodels erroneously returns "loadings" that have identical values and
    # dimensions as the Eigenvectors.
    pail_components = statsmodels.multivariate.pca.PCA(
        matrix,
        ncomp=3, # None # temporarily reduce count to 3
        standardize=True,
        gls=False,
        weights=None,
        method="eig", # "svd", "eig", "nipals"
        missing=None, # None or "drop-row"
    )



    pass



def calculate_principal_components_sklearn(
    data=None,
    components=None,
    report=None,
):
    """
    Calculates the principal components.

    Format of data should have features across columns and observations across
    rows.

    This function then transforms the data frame to a matrix.

    This matrix has observations across dimension zero and features across
    dimension one.

    TODO: (TCW 16 August 2021) I think the data format description is incorrect
    ... I think features are across columns and observations across rows

    arguments:
        data (object): Pandas data frame of signals with features across rows
            and observations across columns
        components (int): count of principle components
        report (bool): whether to print reports to terminal

    raises:

    returns:
        (dict): information about data's principal components

    """

    # Organize data.
    data_copy = data.copy(deep=True)

    # Organize data for principle component analysis.
    # Organize data as an array of arrays.
    matrix = data_copy.to_numpy()
    # Execute principle component analysis.
    pca = sklearn.decomposition.PCA(n_components=components)
    #report = pca.fit_transform(matrix)
    pca.fit(matrix)
    # Report extent of variance that each principal component explains.
    variance_ratios = pca.explained_variance_ratio_
    component_numbers = range(len(variance_ratios) + 1)
    variance_series = {
        "component": component_numbers[1:],
        "variance": variance_ratios
    }
    data_component_variance = pandas.DataFrame(data=variance_series)
    # Transform data by principal components.
    matrix_component = pca.transform(matrix)
    # Match components to samples.
    observations = data_copy.index.to_list()
    records = []
    for index, row in enumerate(matrix_component):
        record = {}
        record["observation"] = observations[index]
        for count in range(components):
            name = "component_" + str(count + 1)
            record[name] = row[count]
        records.append(record)
    data_component = (
        convert_records_to_dataframe(records=records)
    )
    data_component.set_index(
        ["observation"],
        append=False,
        drop=True,
        inplace=True
    )
    # Report
    if report:
        print_terminal_partition(level=3)
        print("dimensions of original data: " + str(data.shape))
        print("dimensions of novel data: " + str(data_component.shape))
        print("proportional variance:")
        print(data_component_variance)
        print_terminal_partition(level=3)
    # Compile information.
    information = {
        "data_observations_components": data_component,
        "data_components_variances": data_component_variance,
    }
    # Return information.
    return information



##########################################################


# Pairwise correlations in symmetric adjacency matrix


def collect_pairwise_correlations_pobabilities(
    method=None,
    features=None,
    data=None,
):
    """
    Calculates and collects correlation coefficients and probabilities across
    observations between pairs of features. These values are symmetrical and
    not specific to order in pairs of features.

    arguments:
        method (str): method for correlation, pearson, spearman, or kendall
        features (list<str>): names of features
        data (object): Pandas data frame of features (columns) across
            observations (rows)

    raises:

    returns:
        (dict): correlation coefficients and probabilities for pairs of
            features

    """

    # Organize data.
    data_copy = data.copy(deep=True)

    # Determine ordered, pairwise combinations of genes.
    # ABCD: AB AC AD BC BD CD
    pairs = combine_unique_elements_pairwise_orderless(
        elements=features,
    )

    # Collect counts and correlations across pairs of features.
    collection = dict()
    # Iterate on features.
    # Generate appropriate values for self-pairs.
    for feature in features:
        collection[feature] = dict()
        collection[feature][feature] = dict()
        collection[feature][feature]["count"] = float("nan")
        collection[feature][feature]["correlation"] = 1.0
        collection[feature][feature]["probability"] = 0.0
        collection[feature][feature]["discovery"] = float("nan")
        collection[feature][feature]["significance"] = True

    # Iterate on pairs of features.
    for pair in pairs:
        # Select data for pair of features.
        data_pair = data_copy.loc[:, list(pair)]
        # Remove observations with missing values for either feature.
        data_pair.dropna(
            axis="index",
            how="any",
            inplace=True,
        )
        # Determine observations for which pair of features has matching
        # values.
        count = data_pair.shape[0]
        # Calculate correlation.
        if count > 1:
            if method == "pearson":
                correlation, probability = scipy.stats.pearsonr(
                    data_pair[pair[0]].to_numpy(),
                    data_pair[pair[1]].to_numpy(),
                )
            elif method == "spearman":
                correlation, probability = scipy.stats.spearmanr(
                    data_pair[pair[0]].to_numpy(),
                    data_pair[pair[1]].to_numpy(),
                )
            elif method == "kendall":
                correlation, probability = scipy.stats.kendalltau(
                    data_pair[pair[0]].to_numpy(),
                    data_pair[pair[1]].to_numpy(),
                )
            pass
        else:
            correlation = float("nan")
            probability = float("nan")
            pass

        # Collect information.
        if not pair[0] in collection:
            collection[pair[0]] = dict()
            pass
        collection[pair[0]][pair[1]] = dict()
        collection[pair[0]][pair[1]]["count"] = count
        collection[pair[0]][pair[1]]["correlation"] = correlation
        collection[pair[0]][pair[1]]["probability"] = probability

        pass

    # Return information.
    return collection


def calculate_pairwise_correlations_discoveries(
    collection=None,
    threshold=None,
):
    """
    Calculates and collects Benjamini-Hochberg false discovery rates from
    probabilities of correlations across observations between pairs of
    features.

    arguments:
        collection (dict): correlation coefficients and probabilities for pairs
            of features
        threshold (float): value of alpha, or family-wise error rate of false
            discoveries

    raises:

    returns:
        (dict): correlation coefficients, probabilities, and discoveries for
            pairs of features

    """

    # Extract probabilities.
    data_probability = extract_pairs_values_long(
        collection=collection,
        key="probability",
        name="probability",
    )

    # Calculate false discovery rates.
    data_discovery = calculate_false_discovery_rate(
        threshold=threshold,
        probability="probability",
        discovery="discovery",
        significance="significance",
        data_probabilities=data_probability,
    )

    # Introduce false discovery rates to collection.
    collection_discovery = copy.deepcopy(collection)
    records = convert_dataframe_to_records(data=data_discovery)
    for record in records:
        one = record["feature_one"]
        two = record["feature_two"]
        collection_discovery[one][two]["discovery"] = record["discovery"]
        collection_discovery[one][two]["significance"] = record["significance"]
        pass
    # Return information.
    return collection_discovery


def extract_pairs_values_long(
    collection=None,
    key=None,
    name=None,
):
    """
    Extracts  a matrix of correlation coefficients.

    arguments:
        collection (dict): collection of correlation coefficients and
            probabilities across observations between pairs of features
        key (str): key to extract values, either correlation or probability
        name (str): name of column for values

    raises:

    returns:
        (object): Pandas data frame of correlation coefficients

    """

    # Copy.
    collection = copy.deepcopy(collection)

    # Collect counts and correlations across pairs of features.
    records = list()
    # Iterate on dimension one of entities.
    for feature_one in collection:
        # Iterate on dimension two of entities.
        for feature_two in collection[feature_one]:
            # Collect information.
            record = dict()
            record["feature_one"] = feature_one
            record["feature_two"] = feature_two
            record[name] = collection[feature_one][feature_two][key]
            # Collect information.
            records.append(record)
            pass
        pass

    # Organize data.
    data = convert_records_to_dataframe(records=records)
    #data.set_index(
    #    ["feature_one", "feature_two"],
    #    append=False,
    #    drop=True,
    #    inplace=True
    #)
    # Return information.
    return data


def organize_symmetric_adjacency_matrix(
    features=None,
    collection=None,
    key=None,
    threshold=None,
    fill=None,
):
    """
    Organizes a symmetric adjacency matrix of correlation coefficients.

    arguments:
        features (list<str>): names of features
        collection (dict): collection of correlation coefficients and
            probabilities across observations between pairs of features
        key (str): key to extract values, either correlation or probability
        threshold (bool): whether to include only values that are significant
            by the false discovery rate
        fill (float): value with which to fill missing values or those that are
            insignificant

    raises:

    returns:
        (object): Pandas data frame of correlation coefficients

    """

    # Collect counts and correlations across pairs of features.
    records = list()
    # Iterate on dimension one of entities.
    for feature_one in features:
        # Collect information.
        record = dict()
        record["feature"] = feature_one

        # Iterate on dimension two of entities.
        for feature_two in features:
            # Access value.
            if feature_two in collection[feature_one]:
                value = collection[feature_one][feature_two][key]
                if threshold:
                    significance = (
                        collection[feature_one][feature_two]["significance"]
                    )
            elif feature_one in collection[feature_two]:
                value = (collection[feature_two][feature_one][key])
                if threshold:
                    significance = (
                        collection[feature_two][feature_one]["significance"]
                    )
            else:
                if threshold:
                    significance = False
            pass
            # Collect information.
            if threshold:
                if significance:
                    record[feature_two] = value
                else:
                    record[feature_two] = fill
            else:
                record[feature_two] = value
        # Collect information.
        records.append(record)
        pass

    # Organize data.
    data = convert_records_to_dataframe(records=records)
    data.sort_values(
        by=["feature"],
        axis="index",
        ascending=True,
        inplace=True,
    )
    data.set_index(
        ["feature"],
        append=False,
        drop=True,
        inplace=True
    )
    data.rename_axis(
        columns="features",
        axis="columns",
        copy=False,
        inplace=True
    )
    # Return information.
    return data


def cluster_adjacency_matrix(
    data=None,
):
    """
    Organizes a matrix of correlation coefficients.

    This function is currently specific to a symmetric adjacency matrix.

    arguments:
        data (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of values

    """

    data = data.copy(deep=True)
    # Cluster.
    columns = data.columns.to_numpy()#.tolist()
    rows = data.index.to_numpy()#.tolist()
    matrix = numpy.transpose(data.to_numpy())
    linkage = scipy.cluster.hierarchy.linkage(
        matrix,
        method="average", # "single", "complete", "average"
        metric="euclidean",
        optimal_ordering=True,
    )
    dendrogram = scipy.cluster.hierarchy.dendrogram(
        linkage,
    )
    dimension = len(matrix)
    # Access seriation from dendrogram leaves.
    leaves = dendrogram["leaves"]
    # Generate new matrix with values of 1.0 that will persist on the diagonal.
    #matrix_order = numpy.zeros((dimension, dimension))
    matrix_cluster = numpy.full(
        (dimension, dimension),
        fill_value=1.0,
    )
    # Sort matrix values.
    a,b = numpy.triu_indices(dimension, k=1)
    matrix_cluster[a,b] = (
        matrix[[leaves[i] for i in a], [leaves[j] for j in b]]
    )
    matrix_cluster[b,a] = matrix_cluster[a,b]
    # Sort matrix row and column labels.
    indices = range(0, dimension)
    rows_sort = list(map(
        lambda index: rows[leaves[index]],
        indices
    ))
    columns_sort = list(map(
        lambda index: columns[leaves[index]],
        indices
    ))
    # Organize data.
    data_cluster = pandas.DataFrame(
        data=matrix_cluster,
        index=rows_sort,
        columns=columns_sort,
    )
    # Return information.
    return data_cluster


def cluster_data_columns(
    data=None,
):
    """
    Clusters features on columns by their similarities across instances on
    rows.

    arguments:
        data (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of values

    """

    data = data.copy(deep=True)
    # Cluster.
    columns = data.columns.to_numpy()#.tolist()
    rows = data.index.to_numpy()#.tolist()
    index_name = data.index.name
    # Plan to cluster across columns.
    # Organize columns across dimension zero.
    matrix = numpy.transpose(data.to_numpy())
    linkage = scipy.cluster.hierarchy.linkage(
        matrix,
        method="average", # "single", "complete", "average"
        metric="euclidean",
        optimal_ordering=True,
    )
    dendrogram = scipy.cluster.hierarchy.dendrogram(
        linkage,
    )
    # Access seriation from dendrogram leaves.
    leaves = dendrogram["leaves"]
    # Sort matrix row and column labels.
    indices = range(0, len(matrix))
    matrix_cluster = list(map(
        lambda index: matrix[leaves[index]],
        indices
    ))
    columns_sort = list(map(
        lambda index: columns[leaves[index]],
        indices
    ))
    # Organize data.
    data_cluster = pandas.DataFrame(
        data=numpy.transpose(matrix_cluster),
        index=rows,
        columns=columns_sort,
    )
    data_cluster.rename_axis(
        index_name,
        axis="index",
        inplace=True,
    )
    # Return information.
    return data_cluster


def cluster_data_rows(
    data=None,
):
    """
    Clusters instances on rows by their similarities across features on
    columns.

    arguments:
        data (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of values

    """

    data = data.copy(deep=True)
    # Cluster.
    columns = data.columns.to_numpy()#.tolist()
    rows = data.index.to_numpy()#.tolist()
    index_name = data.index.name
    # Plan to cluster across columns.
    # Organize rows across dimension zero.
    matrix = data.to_numpy()
    linkage = scipy.cluster.hierarchy.linkage(
        matrix,
        method="average", # "single", "complete", "average"
        metric="euclidean",
        optimal_ordering=True,
    )
    dendrogram = scipy.cluster.hierarchy.dendrogram(
        linkage,
    )
    # Access seriation from dendrogram leaves.
    leaves = dendrogram["leaves"]
    # Sort matrix row and column labels.
    indices = range(0, len(matrix))
    matrix_cluster = list(map(
        lambda index: matrix[leaves[index]],
        indices
    ))
    rows_sort = list(map(
        lambda index: rows[leaves[index]],
        indices
    ))
    # Organize data.
    data_cluster = pandas.DataFrame(
        data=matrix_cluster,
        index=rows_sort,
        columns=columns,
    )
    data_cluster.rename_axis(
        index_name,
        axis="index",
        inplace=True,
    )
    # Return information.
    return data_cluster


def cluster_data_rows_by_group(
    group=None,
    index=None,
    data=None,
):
    """
    Clusters instances on rows by their similarities across features on
    columns.

    arguments:
        group (str): name of column to use for groups
        index (str): name of column to use for index during cluster
        data (object): Pandas data frame of values

    raises:

    returns:
        (object): Pandas data frame of values

    """

    data = data.copy(deep=True)
    groups = data.groupby(
        level=[group],
    )
    data_collection = pandas.DataFrame()
    for name, data_group in groups:
        data_group = data_group.copy(deep=True)
        data_group.reset_index(
            level=None,
            inplace=True
        )
        data_group.set_index(
            [index],
            append=False,
            drop=True,
            inplace=True
        )
        data_cluster = cluster_data_rows(
            data=data_group,
        )
        data_collection = data_collection.append(
            data_cluster,
            ignore_index=False,
        )
    # Return information.
    return data_collection


def filter_features_by_threshold_outer_count(
    data=None,
    threshold_high=None,
    threshold_low=None,
    count=None,
):
    """
    Filters features in a symmetric adjacency matrix by whether a certain count
    of values are outside of low and high thresholds.

    arguments:
        data (object): Pandas data frame of values
        threshold_high (float): value must be greater than this threshold
        threshold_low (float): value must be less than this threshold
        count (int): minimal count of rows or columns that must
            pass threshold

    raises:

    returns:
        (object): Pandas data frame of values

    """

    def count_true(slice=None, count=None):
        values = slice.to_numpy().tolist()
        values_true = list(itertools.compress(values, values))
        return (len(values_true) >= count)

    # Copy data.
    data = data.copy(deep=True)
    # Determine whether values exceed threshold.
    data_threshold = data.applymap(
        lambda value: ((value <= threshold_low) or (threshold_high <= value))
    )
    # This aggregation operation produces a series.
    data_count_row = data_threshold.aggregate(
        lambda slice: count_true(slice=slice, count=count),
        axis="columns",
    )
    # Select rows with appropriate values.
    data_row = data.loc[data_count_row, : ]
    # Extract identifiers of features.
    features = collect_unique_elements(
        elements_original=data_row.index.tolist()
    )
    # Select identical features across rows and columns.
    data_feature_row = data.loc[
        data.index.isin(features),
        :
    ]
    data_feature_column = data_feature_row.loc[
        :, data_feature_row.columns.isin(features)
    ]
    # Return information.
    return data_feature_column


def organize_feature_signal_correlations(
    method=None,
    threshold_high=None,
    threshold_low=None,
    count=None,
    discovery=None,
    features=None,
    data_signal=None,
):
    """
    Calculates and organizes correlation coefficients between pairs of features
    across observations.

    arguments:
        method (str): method for correlation, pearson, spearman, or kendall
        threshold_high (float): value must be greater than this threshold
        threshold_low (float): value must be less than this threshold
        count (int): minimal count of rows or columns that must
            pass threshold
        discovery (float): value of false discovery rate for which to include
            correlation coefficients
        features (list<str>): identifiers of features
        data_signal (object): Pandas data frame of features' signals across
            observations

    raises:

    returns:
        (object): Pandas data frame of genes' pairwise correlations

    """

    # Organize data.
    data = data_signal.copy(deep=True)

    # Calculate correlations between gene pairs of their pantissue signals
    # across persons.
    # Only collect genes with correlations beyond thresholds.
    collection = collect_pairwise_correlations_pobabilities(
        method=method,
        features=features,
        data=data,
    )
    # Calculate false discovery rates for each primary feature.
    collection_discovery = calculate_pairwise_correlations_discoveries(
        collection=collection,
        threshold=discovery,
    )

    # Organize data matrix for correlation coefficients.
    # Filter correlation coefficients by value of false discovery rate.
    # Clustering algorithm requires that matrix not include missing values.
    data_correlation = organize_symmetric_adjacency_matrix(
        features=features,
        collection=collection_discovery,
        key="correlation",
        threshold=True,
        fill=0.0,
    )
    # Filter genes by their pairwise correlations.
    data_pass = filter_features_by_threshold_outer_count(
        data=data_correlation,
        threshold_high=threshold_high,
        threshold_low=threshold_low,
        count=count,
    )
    # Cluster data.
    data_cluster = cluster_adjacency_matrix(
        data=data_pass,
    )
    #data_cluster = data_pass
    # Return information.
    return data_cluster


# Metabolic information.


# TODO: this function is actually a nice model for string formatting
def prepare_curation_report(
    compartments=None,
    processes=None,
    reactions=None,
    metabolites=None
):
    """
    Prepares a summary report on curation of metabolic sets and entities.

    arguments:
        compartments (dict<dict>): information about compartments
        processes (dict<dict>): information about processes
        reactions (dict<dict>): information about reactions
        metabolites (dict<dict>): information about metabolites

    returns:
        (str): report of summary information

    raises:

    """

    # Count compartments.
    count_compartments = len(compartments)
    # Count processes.
    count_processes = len(processes)
    # Count reactions.
    count_reactions = len(reactions)
    # Count metabolites.
    count_metabolites = len(metabolites)
    # Count reactions with references to MetaNetX.
    count_one = count_entities_with_references(
        references=["metanetx"],
        entities=reactions
    )
    proportion_one = count_one / count_reactions
    percentage_one = round((proportion_one * 100), 2)
    # Count reactions with references either to genes or enzyme commission.
    count_two = count_entities_with_references(
        references=["gene", "enzyme"],
        entities=reactions
    )
    proportion_two = count_two / count_reactions
    percentage_two = round((proportion_two * 100), 2)
    # Count metabolites with references to MetaNetX.
    count_three = count_entities_with_references(
        references=["metanetx"],
        entities=metabolites
    )
    proportion_three = count_three / count_metabolites
    percentage_three = round((proportion_three * 100), 2)
    # Count metabolites with references to Human Metabolome Database (HMDB) and
    # PubChem.
    count_four = count_entities_with_references(
        references=["hmdb", "pubchem"],
        entities=metabolites
    )
    proportion_four = count_four / count_metabolites
    percentage_four = round((proportion_four * 100), 2)
    # Compile information.
    report = textwrap.dedent("""\

        --------------------------------------------------
        curation report

        compartments: {count_compartments}
        processes: {count_processes}
        reactions: {count_reactions}
        metabolites: {count_metabolites}

        reactions in MetaNetX: {count_one} ({percentage_one} %)
        reactions with gene or enzyme: {count_two} ({percentage_two} %)
        metabolites in MetaNetX: {count_three} ({percentage_three} %)
        metabolites with HMDB or PubChem: {count_four} ({percentage_four} %)

        --------------------------------------------------
    """).format(
        count_compartments=count_compartments,
        count_processes=count_processes,
        count_reactions=count_reactions,
        count_metabolites=count_metabolites,
        count_one=count_one,
        percentage_one=percentage_one,
        count_two=count_two,
        percentage_two=percentage_two,
        count_three=count_three,
        percentage_three=percentage_three,
        count_four=count_four,
        percentage_four=percentage_four
    )
    # Return information.
    return report


###############################################################################
# Procedure
# Currently, this module is not executable.
