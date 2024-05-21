"""
Organize procedural code for extraction of information from raw text report logs
that the LDSC tool creates for estimates of SNP heritability (h2) and
genetic correlation (rg).

This module makes use of the 'extraction' module within the 'partner' package.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This file is part of project 'partner'
    (https://github.com/tcameronwaller/partner/).

    Project 'partner' supports data analysis for multiple projects in
    biomedical research.
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

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 21 May 2024
# Date, last execution: 21 May 2024
# Review: TCW; 21 May 2024
################################################################################
# Note

# An optional enhancement is to use the functions "putly.read_file_text_list()
# and "putly.sort_table_rows_by_list_indices()" to sort rows in the extraction
# tables according to parameter text files of indices.
# As of 21 May 2024, my preferred alternative is to handle filters and sorts all
# together in post-processing.


################################################################################
# Installation and importation

# Standard.
import sys
import os
import copy

# Relevant.
import pandas
import scipy
import numpy

# Custom
import partner.utility as putly # this import path for subpackage
import partner.extraction as pextr
#import partner.scale as pale
#import partner.description as pdesc

#dir()
#importlib.reload()

###############################################################################
# Functionality



################################################################################
# Procedure


def execute_procedure(
    type_analysis=None,
    path_directory_source=None,
    traversal=None,
    name_file_source_prefix=None,
    name_file_source_suffix=None,
    name_file_source_not=None,
    name_file_product=None,
    path_directory_product=None,
    path_directory_temporary=None,
    report=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        type_analysis (str): type of analysis in LDSC, either 'heritability' or
            'correlation', corresponding to information for extraction
        path_directory_source (str): path to parent directory in which to find
            child source files
        traversal (bool): whether to extract from all files in child source
            directories, preserving names of child directories
        name_file_source_prefix (str): prefix in name by which to recognize
            relevant source child files within parent directory
        name_file_source_suffix (str): suffix in name by which to recognize
            relevant source child files within parent directory
        name_file_source_not (str): character string in names of files to
            exclude
        name_file_product (str): either none for traversal of child source
            directories or the name for product file
        path_directory_product (str): path to parent directory in which to write
            child product files
        path_directory_temporary (str): path to parent directory for temporary
            directories and files
        report (bool): whether to print reports

    raises:

    returns:

    """

    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("Extraction from LDSC report logs.")
        print(type_analysis)
        print(path_directory_source)
        print(traversal)
        print(name_file_source_prefix)
        print(name_file_source_suffix)
        print(name_file_source_not)
        print(name_file_product)
        print(path_directory_product)
        print(path_directory_temporary)
        putly.print_terminal_partition(level=4)

    # Collect information.
    pail_write = dict()
    # Determine whether to traverse subdirectories.
    if (traversal):
        # Extract names of child directories within parent directory.
        names_directories = putly.extract_child_directory_names(
            path_directory=path_directory_source,
        )
        names_directories_ldsc = list(filter(
            lambda name: (name != "batch"),
            names_directories
        ))
        # Iterate on child subdirectories within parent directory.
        for name_directory in names_directories_ldsc:
            # Determine name for product file.
            name_file_product = str(
                "table_" + name_directory
            )
            # Extract information from reports of analyses in LDSC
            path_directory_source_child = os.path.join(
                path_directory_source, name_directory,
            )
            pail_write[name_file_product] = (
                pextr.read_extract_from_all_ldsc_files_in_directory(
                    path_directory=path_directory_source_child,
                    name_file_prefix=name_file_source_prefix,
                    name_file_suffix=name_file_source_suffix,
                    name_file_not=name_file_source_not,
                    type_analysis=type_analysis,
                    report=report,
            ))
            pass
        pass
    else:
        if (
            (len(str(name_file_product).strip()) == 0) or
            (str(name_file_product).strip() == "none")
        ):
            # Determine name for product file.
            name_file_product = str(
                "table_" + type_analysis
            )
            pass
        # Extract information from reports of analyses in LDSC.
        pail_write[name_file_product] = (
            pextr.read_extract_from_all_ldsc_files_in_directory(
                path_directory=path_directory_source,
                name_file_prefix=name_file_source_prefix,
                name_file_suffix=name_file_source_suffix,
                name_file_not=name_file_source_not,
                analysis=type_analysis,
                report=report,
        ))
        pass

    # Write product information to file.
    putly.write_product_tables(
        pail_write=pail_write,
        path_directory=path_directory_product,
    )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    name_script = sys.argv[0]
    type_analysis = sys.argv[1]
    path_directory_source = sys.argv[2]
    traversal = sys.argv[3]
    name_file_source_prefix = sys.argv[4]
    name_file_source_suffix = sys.argv[5]
    name_file_source_not = sys.argv[6]
    name_file_product = sys.argv[7]
    path_directory_product = sys.argv[8]
    path_directory_temporary = sys.argv[9]
    report = sys.argv[10]
    # Parse and interpret values.
    traversal = putly.parse_text_boolean(
        string=traversal,
    )
    report = putly.parse_text_boolean(
        string=report,
    )
    # Call function for procedure.
    execute_procedure(
        type_analysis=type_analysis,
        path_directory_source=path_directory_source,
        traversal=traversal,
        name_file_source_prefix=name_file_source_prefix,
        name_file_source_suffix=name_file_source_suffix,
        name_file_source_not=name_file_source_not,
        name_file_product=name_file_product,
        path_directory_product=path_directory_product,
        path_directory_temporary=path_directory_temporary,
        report=report,
    )

    pass



#
