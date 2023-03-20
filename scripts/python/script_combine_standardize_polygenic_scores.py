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

# Custom
import promiscuity.utility as utility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_organize_table_polygenic_scores(
    path_table=None,
    name_table=None,
    name_column_identifier=None,
    name_column_allele_total=None,
    name_column_allele_dosage=None,
    name_column_score_sum=None,
    name_column_score_mean=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_table (str): path to file for table
        name_table (str): name of table that distinguishes it from all others
        name_column_identifier (str): name of column in source table for
            identifier of genotypes
        name_column_allele_total (str): name of column in source table for
            count of total non-missing alleles in genotypes
        name_column_allele_dosage (str): name of column in source table for
            count of alleles considered in calculation of polygenic score
        name_column_score_sum (str): name of column in source table for
            polygenic score as sum of allelic effects across genotypes
        name_column_score_mean (str): name of column in source table for
            polygenic score as mean of allelic effects (sum divided by count of
            total nonmissing alleles) across genotypes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of polygenic scores across genotypes

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("Parameters:")
        print(name_column_identifier)
        print(name_column_allele_total)
        print(name_column_allele_dosage)
        print(name_column_score_sum)
        print(name_column_score_mean)
        utility.print_terminal_partition(level=4)

    # Read information from file.
    types_columns = dict()
    types_columns[name_column_identifier] = "string"
    types_columns[name_column_allele_total] = "int"
    types_columns[name_column_allele_dosage] = "int"
    types_columns[name_column_score_sum] = "float"
    types_columns[name_column_score_mean] = "float"
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
    allele_total = str("count_allele_total_" + name_table)
    allele_dosage = str("count_allele_dosage_" + name_table)
    score_sum = str("score_sum_" + name_table)
    score_mean = str("score_mean_" + name_table)
    # Simplify identifiers of genotypes.
    # PLINK2 combines "IID" and "FID" identifiers.
    table[identifier] = table.apply(
        lambda row: str(str(row[name_column_identifier]).split("_")[0]),
        axis="columns", # apply function to each row
    )
    # Translate names of columns.
    translations = dict()
    #translations[name_column_identifier] = identifier
    translations[name_column_allele_total] = allele_total
    translations[name_column_allele_dosage] = allele_dosage
    translations[name_column_score_sum] = score_sum
    translations[name_column_score_mean] = score_mean
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Select relevant columns.
    columns = [
        identifier, allele_total, allele_dosage,
        score_sum, score_mean,
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
    name_column_identifier=None,
    name_column_allele_total=None,
    name_column_allele_dosage=None,
    name_column_score_sum=None,
    name_column_score_mean=None,
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
        name_column_score_sum (str): name of column in source table for
            polygenic score as sum of allelic effects across genotypes
        name_column_score_mean (str): name of column in source table for
            polygenic score as mean of allelic effects (sum divided by count of
            total nonmissing alleles) across genotypes
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
            name_column_identifier=name_column_identifier,
            name_column_allele_total=name_column_allele_total,
            name_column_allele_dosage=name_column_allele_dosage,
            name_column_score_sum=name_column_score_sum,
            name_column_score_mean=name_column_score_mean,
            report=True,
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


def combine_standardize_polygenic_scores(
    table=None,
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

    # Copy information in table.
    table = table.copy(deep=True)

    # Extract names of columns in table.
    columns = copy.deepcopy(table.columns.to_list())
    columns_count_allele_total = list(filter(
        lambda column: (str("count_allele_total") in str(column)),
        copy.deepcopy(columns)
    ))
    columns_count_allele_dosage = list(filter(
        lambda column: (str("count_allele_dosage") in str(column)),
        copy.deepcopy(columns)
    ))
    columns_score_sum = list(filter(
        lambda column: (str("score_sum") in str(column)),
        copy.deepcopy(columns)
    ))
    columns_score_mean = list(filter(
        lambda column: (str("score_mean") in str(column)),
        copy.deepcopy(columns)
    ))
    # Combine columns by sum.
    table["count_allele_total"] = table.apply(
        lambda row: utility.calculate_sum_row_column_values(
            columns=columns_count_allele_total,
            row=row.copy(deep=True),
        ),
        axis="columns", # apply function to each row
    )
    table["count_allele_dosage"] = table.apply(
        lambda row: utility.calculate_sum_row_column_values(
            columns=columns_count_allele_dosage,
            row=row.copy(deep=True),
        ),
        axis="columns", # apply function to each row
    )
    table["score_sum"] = table.apply(
        lambda row: utility.calculate_sum_row_column_values(
            columns=columns_score_sum,
            row=row.copy(deep=True),
        ),
        axis="columns", # apply function to each row
    )
    table["score_mean"] = table.apply(
        lambda row: utility.calculate_sum_row_column_values(
            columns=columns_score_mean,
            row=row.copy(deep=True),
        ),
        axis="columns", # apply function to each row
    )

    # Select relevant columns.
    columns_keep = [
        "count_allele_total", "count_allele_dosage",
        "score_sum", "score_mean",
    ]
    table = table.loc[:, table.columns.isin(columns_keep)]
    table = table[[*columns_keep]]
    # Report.
    if report:
        utility.print_terminal_partition(level=3)
        print("report: ")
        print("combine_standardize_polygenic_scores()")
        utility.print_terminal_partition(level=4)
        print("table columns: " + str(int(table.shape[1])))
        print("table rows: " + str(int(table.shape[0])))
        print("columns")
        print(table.columns.to_list())
        print(table)
        pass
    # Return information.
    return table




################################################################################
# Procedure


def execute_procedure(
    path_directory_source=None,
    name_file_source_prefix=None,
    name_file_source_suffix=None,
    name_file_source_not=None,
    name_column_identifier=None,
    name_column_allele_total=None,
    name_column_allele_dosage=None,
    name_column_score_sum=None,
    name_column_score_mean=None,
    path_file_product=None,
):
    """
    Function to execute module's main behavior.

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
        name_column_score_sum (str): name of column in source table for
            polygenic score as sum of allelic effects across genotypes
        name_column_score_mean (str): name of column in source table for
            polygenic score as mean of allelic effects (sum divided by count of
            total nonmissing alleles) across genotypes
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
        name_column_identifier=name_column_identifier,
        name_column_allele_total=name_column_allele_total,
        name_column_allele_dosage=name_column_allele_dosage,
        name_column_score_sum=name_column_score_sum,
        name_column_score_mean=name_column_score_mean,
        report=True,
    )

    # Merge together tables for polygenic scores.
    table_merge = merge_tables_polygenic_scores(
        pail_tables=pail_tables,
        report=True,
    )

    # Calculate combinations and standardizations of polygenic scores.
    table_combination = combine_standardize_polygenic_scores(
        table=table_merge,
        report=True,
    )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    name_script = sys.argv[0]
    path_directory_source = sys.argv[1]
    name_file_source_prefix = sys.argv[2]
    name_file_source_suffix = sys.argv[3]
    name_file_source_not = sys.argv[4]
    name_column_identifier = sys.argv[5]
    name_column_allele_total = sys.argv[6]
    name_column_allele_dosage = sys.argv[7]
    name_column_score_sum = sys.argv[8]
    name_column_score_mean = sys.argv[9]
    path_file_product = sys.argv[10]

    # Call function for procedure.
    execute_procedure(
        path_directory_source=path_directory_source,
        name_file_source_prefix=name_file_source_prefix,
        name_file_source_suffix=name_file_source_suffix,
        name_file_source_not=name_file_source_not,
        name_column_identifier=name_column_identifier,
        name_column_allele_total=name_column_allele_total,
        name_column_allele_dosage=name_column_allele_dosage,
        name_column_score_sum=name_column_score_sum,
        name_column_score_mean=name_column_score_mean,
        path_file_product=path_file_product,
    )

    pass



#
