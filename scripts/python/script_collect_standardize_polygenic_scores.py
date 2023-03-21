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
import copy

# Relevant.
import pandas
import scipy
import numpy

# Custom
import promiscuity.utility as utility
import promiscuity.scale as pscale

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_organize_table_polygenic_scores(
    path_table=None,
    name_table=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_table (str): path to file for table
        name_table (str): name of table that distinguishes it from all others
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of polygenic scores across genotypes

    """

    # Read information from file.
    types_columns = dict()
    types_columns["identifier"] = "string"
    types_columns["score"] = "float"
    table = pandas.read_csv(
        path_table,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=["nan", "na", "NAN", "NA",],
    )
    # Organize information in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    # Assign new names to columns.
    identifier = str("identifier")
    score = str("score_" + name_table)
    # Translate names of columns.
    translations = dict()
    translations["identifier"] = identifier
    translations["score"] = score
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Select relevant columns.
    columns = [
        identifier, score,
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
    name_file_child_not=None,
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
        name_file_child_not (str): character string in names of files to exclude
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
        name_file_child_not=name_file_child_not,
        report=report,
    )
    # Read files as Pandas dataframe tables.
    # Iterate on names of files to read and organize tables.
    # Collect tables.
    pail = dict()
    for path in paths:
        # Extract name of file and table that distinguishes it from all others.
        name_file = os.path.basename(path)
        name_table_1 = name_file.replace(str(name_file_child_prefix), "")
        name_table = name_table_1.replace(str(name_file_child_suffix), "")
        # Read file and organize information in table.
        table = read_organize_table_polygenic_scores(
            path_table=path,
            name_table=name_table,
            report=report,
        )
        # Collect table.
        pail[name_table] = table.copy(deep=True)
        pass
    # Return information.
    return pail


def merge_tables_polygenic_scores(
    pail_tables=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        pail_tables (dict<object>): collection of Pandas data-frame tables with
            entry names (keys) derived from original names of files
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table

    """

    # Extract list of tables.
    #names = list(pail_tables.keys())
    tables = list(pail_tables.values())
    # Merge tables.
    table_merge = utility.merge_columns_tables_supplements_to_main(
        identifier_main="identifier",
        identifier_supplement="identifier",
        table_main=tables[0],
        tables_supplements=tables[1:],
        report=True,
    )
    # Organize table after merge.
    table_merge.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table_merge.set_index(
        "identifier",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True,
    )
    # Sort table columns.
    #columns = copy.deepcopy(table_merge.columns.to_list())
    #columns_sort = sorted(columns, reverse=False)
    #table_merge = table_merge[[*columns_sort]]
    table_merge.sort_index(
        axis="columns",
        ascending=True,
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=3)
        print("report: ")
        print("merge_tables_polygenic_scores()")
        utility.print_terminal_partition(level=4)
        print("table columns: " + str(int(table_merge.shape[1])))
        print("table rows: " + str(int(table_merge.shape[0])))
        print("columns")
        print(table_merge.columns.to_list())
        print(table_merge)
        pass
    # Return information.
    return table_merge


def standardize_polygenic_scores(
    table=None,
    report=None,
):
    """
    Standardizes polygenic scores to z-scores with mean of zero and standard
    deviation of one.

    arguments:
        table (object): Pandas data-frame table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table
    """

    # Copy information in table.
    table = table.copy(deep=True)

    # Extract names of columns in table.
    columns = copy.deepcopy(table.columns.to_list())
    columns_score = list(filter(
        lambda column: (str("score_") in str(column)),
        copy.deepcopy(columns)
    ))
    # Standardize scores to z-score scale with mean of zeroa and standard
    # deviation of one.
    # Standardization introduces missing values if standard deviation is zero.
    table = pscale.drive_transform_variables_distribution_scale_z_score(
        table=table,
        columns=columns_score,
        report=report,
    )

    # Select relevant columns.
    # Extract names of columns in table.
    columns = copy.deepcopy(table.columns.to_list())
    columns_keep = list(filter(
        lambda column: (str("score_") in str(column)),
        copy.deepcopy(columns)
    ))
    table = table.loc[:, table.columns.isin(columns_keep)]
    table = table[[*columns_keep]]
    # Report.
    if report:
        utility.print_terminal_partition(level=3)
        print("report: ")
        print("standardize_polygenic_scores()")
        utility.print_terminal_partition(level=4)
        print("table columns: " + str(int(table.shape[1])))
        print("table rows: " + str(int(table.shape[0])))
        print("columns")
        print(table.columns.to_list())
        print(table)
        pass
    # Return information.
    return table


def write_product_table_text_tab(
    table=None,
    path_file=None,
):
    """
    Writes product information to file.

    arguments:
        table (object): table of information to write to file
        path_file (str): path to which to write file

    raises:

    returns:

    """

    # Write information to file.
    table.to_csv(
        path_or_buf=path_file,
        sep="\t",
        header=True,
        index=True,
    )
    pass


################################################################################
# Procedure


def execute_procedure(
    path_directory_source=None,
    name_file_source_prefix=None,
    name_file_source_suffix=None,
    name_file_source_not=None,
    path_file_product=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_directory_source (str): path to parent directory in which to find
            child files
        name_file_source_prefix (str): prefix in name by which to recognize
            relevant child files within parent directory
        name_file_source_suffix (str): suffix in name by which to recognize
            relevant child files within parent directory
        name_file_source_not (str): character string in names of files to
            exclude
        path_file_product (str): path to product file after combination and
            standardization of polygenic scores across chromosomes

    raises:

    returns:

    """

    # Read from source files the tables for polygenic scores.
    pail_tables = read_source_directory_files_polygenic_scores(
        path_directory_parent=path_directory_source,
        name_file_child_prefix=name_file_source_prefix,
        name_file_child_suffix=name_file_source_suffix,
        name_file_child_not=name_file_source_not,
        report=False,
    )

    # Merge together tables for polygenic scores.
    table_merge = merge_tables_polygenic_scores(
        pail_tables=pail_tables,
        report=True,
    )

    # Calculate standardizations of polygenic scores.
    table_standardization = standardize_polygenic_scores(
        table=table_merge,
        report=True,
    )

    # Compare the polygenic score sum against the mean.
    utility.calculate_table_column_pair_correlations(
        column_one="score_bmi_sbayesr_1",
        column_two="score_bmi_sbayesr_2",
        table=table_standardization,
        report=True,
    )
    utility.calculate_table_column_pair_correlations(
        column_one="score_bmi_sbayesr_1",
        column_two="score_bmi_sbayesr_3",
        table=table_standardization,
        report=True,
    )
    utility.calculate_table_column_pair_correlations(
        column_one="score_bmi_sbayesr_2",
        column_two="score_bmi_sbayesr_3",
        table=table_standardization,
        report=True,
    )
    utility.calculate_table_column_pair_correlations(
        column_one="score_bmi_sbayesr_1",
        column_two="score_bmi_ldpred2",
        table=table_standardization,
        report=True,
    )
    utility.calculate_table_column_pair_correlations(
        column_one="score_bmi_sbayesr_1",
        column_two="score_bmi_ldpred2_2",
        table=table_standardization,
        report=True,
    )
    utility.calculate_table_column_pair_correlations(
        column_one="score_bmi_sbayesr_1",
        column_two="score_bmi_prsice2",
        table=table_standardization,
        report=True,
    )
    utility.calculate_table_column_pair_correlations(
        column_one="score_bmi_ldpred2",
        column_two="score_bmi_ldpred2_2",
        table=table_standardization,
        report=True,
    )
    utility.calculate_table_column_pair_correlations(
        column_one="score_bmi_ldpred2",
        column_two="score_bmi_prsice2",
        table=table_standardization,
        report=True,
    )

    # Write table to file.
    write_product_table_text_tab(
        table=table_standardization,
        path_file=path_file_product,
    )


    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    name_script = sys.argv[0]
    path_directory_source = sys.argv[1]
    name_file_source_prefix = sys.argv[2]
    name_file_source_suffix = sys.argv[3]
    name_file_source_not = sys.argv[4]
    path_file_product = sys.argv[5]

    # Call function for procedure.
    execute_procedure(
        path_directory_source=path_directory_source,
        name_file_source_prefix=name_file_source_prefix,
        name_file_source_suffix=name_file_source_suffix,
        name_file_source_not=name_file_source_not,
        path_file_product=path_file_product,
    )

    pass



#
