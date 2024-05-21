"""
...
"""

################################################################################
# Notes

# TODO: TCW; 17 July 2023
# TODO: Need to change the name of the column(s) for the genotype identifier.


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
import partner.scale as pale
import partner.description as pdesc


#dir()
#importlib.reload()

###############################################################################
# Functionality



################################################################################
# Procedure


# TODO: TCW; 17 May 2024
# TODO: implement module's main behavior
# TODO: derive functionality from "psychiatry_biomarkers.extraction_ldsc.py"



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
        type_analysis (str): type of analysis, either heritability or
            correlation, corresponding to information for extraction
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

    # Parameters.
    report = True

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
