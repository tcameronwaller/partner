"""
...
"""

################################################################################
# Notes

################################################################################
# Installation and importation

# Standard.
import sys
import os

# Relevant.
import pandas

# Custom
import promiscuity.utility as utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_organize_table_polygenic_scores(
    path_table=None,
    name_column_identifier=None,
    name_column_allele_total=None,
    name_column_allele_dosage=None,
    name_column_score=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_table (str): path to file for table
        name_column_identifier (str): name of column in source table for
            identifier of genotypes
        name_column_allele_total (str): name of column in source table for
            count of total non-missing alleles in genotypes
        name_column_allele_dosage (str): name of column in source table for
            count of alleles considered in calculation of polygenic score
        name_column_score (str): name of column in source table for
            sum (not average or mean) polygenic score across genotypes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of polygenic scores across genotypes

    """

    # Read information from file.
    table = pandas.read_csv(
        path_table,
        sep="\t",
        header=0,
        dtype={
            [name_column_identifier]: "string",
            [name_column_allele_total]: "int",
            [name_column_allele_dosage]: "int",
            [name_column_score]: "float",
        },
        na_values=["nan", "na", "NAN", "NA",],
    )
    # Organize information in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    # Translate names of columns.
    translations = dict()
    translations[name_column_identifier] = "identifier"
    translations[name_column_allele_total] = "count_allele_total"
    translations[name_column_allele_dosage] = "count_allele_dosage"
    translations[name_column_score] = "score"
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Select relevant columns.
    columns = [
        "identifier", "count_allele_total", "count_allele_dosage", "score",
    ]
    table = table.loc[:, table.columns.isin(columns)]
    table = table[[*columns]]
    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("Table of polygenic scores across genotypes:")
        print(table)
        utility.print_terminal_partition(level=4)
    # Return information.
    return table


def read_source_directory_files_polygenic_scores(
    path_directory_parent=None,
    name_file_child_prefix=None,
    name_file_child_suffix=None,
    name_column_identifier=None,
    name_column_allele_total=None,
    name_column_allele_dosage=None,
    name_column_score=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_directory_parent (str): path to parent directory in which to find
            child files
        name_file_child_prefix (str): prefix in name by which to recognize
            relevant child files within parent directory
        name_file_child_suffix (str): suffix in name by which to recognize
            relevant child files within parent directory
        name_column_identifier (str): name of column in source table for
            identifier of genotypes
        name_column_allele_total (str): name of column in source table for
            count of total non-missing alleles in genotypes
        name_column_allele_dosage (str): name of column in source table for
            count of alleles considered in calculation of polygenic score
        name_column_score (str): name of column in source table for
            sum (not average or mean) polygenic score across genotypes
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables with entry names
            (keys) derived from original names of files

    """

    # Read all matching files within parent directory and organize paths to
    # these files.
    paths = utility.read_paths_match_child_files_within_parent_directory(
        path_directory_parent=path_directory_parent,
        name_file_child_prefix=name_file_source_prefix,
        name_file_child_suffix=name_file_source_suffix,
        report=report,
    )
    # Read files as Pandas dataframe tables.
    # Iterate on names of files to read and organize tables.
    # Collect tables.
    pail = dict()
    for path in paths:
        # Extract name of file.
        name_file = os.path.basename(path)
        name_simple = name_file.replace(str(name_file_child_suffix), "")
        # Read file and organize information in table.
        table = read_organize_table_polygenic_scores(
            path_table=path,
            name_column_identifier=name_column_identifier,
            name_column_allele_total=name_column_allele_total,
            name_column_allele_dosage=name_column_allele_dosage,
            name_column_score=name_column_score,
            report=True,
        )
        # Collect table.
        pail[name_simple] = table.copy(deep=True)
        pass
    # Return information.
    return pail












################################################################################
# Procedure


def execute_procedure():
    """
    Function to execute module's main behavior.

    arguments:

    returns:

    raises:

    """

    # Read from source files the tables for polygenic scores.
    pail_score_tables = read_source_directory_files_polygenic_scores(
        path_directory_parent=path_directory_source,
        name_file_child_prefix=name_file_source_prefix,
        name_file_child_suffix=name_file_source_suffix,
        name_column_identifier=name_column_identifier,
        name_column_allele_total=name_column_allele_total,
        name_column_allele_dosage=name_column_allele_dosage,
        name_column_score=name_column_score,
        report=True,
    )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    name_script = sys.argv[0]
    path_directory_source = sys.argv[1]
    name_file_source_prefix = sys.argv[2]
    name_file_source_suffix = sys.argv[3]
    name_column_identifier = sys.argv[4]
    name_column_allele_total = sys.argv[5]
    name_column_allele_dosage = sys.argv[6]
    name_column_score = sys.argv[7]
    path_file_product = sys.argv[8]

    # Read from source files the tables for polygenic scores.
    pail_score_tables = read_source_directory_files_polygenic_scores(
        path_directory_parent=path_directory_source,
        name_file_child_prefix=name_file_source_prefix,
        name_file_child_suffix=name_file_source_suffix,
        name_column_identifier=name_column_identifier,
        name_column_allele_total=name_column_allele_total,
        name_column_allele_dosage=name_column_allele_dosage,
        name_column_score=name_column_score,
        report=True,
    )

    #execute_procedure()



#
