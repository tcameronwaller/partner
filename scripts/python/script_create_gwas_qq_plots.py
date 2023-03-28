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
import promiscuity.plot as pplot

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_organize_table_gwas(
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
    types_columns["BETA"] = "float"
    types_columns["SE"] = "float"
    types_columns["P"] = "float"
    table = pandas.read_csv(
        path_table,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=["nan", "na", "NAN", "NA",],
        compression="infer",
    )
    # Organize information in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    # Select relevant columns.
    columns = [
        "BETA", "SE", "P"
    ]
    table = table.loc[:, table.columns.isin(columns)]
    table = table[[*columns]]
    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("Table of GWAS summary statistics:")
        print(table)
        utility.print_terminal_partition(level=4)
    # Return information.
    return table


def read_source_directory_files_gwas(
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
        name_table = name_file.replace(str(name_file_child_suffix), "")
        # Read file and organize information in table.
        table = read_organize_table_gwas(
            path_table=path,
            name_table=name_table,
            report=report,
        )
        # Collect table.
        pail[name_table] = table.copy(deep=True)
        pass
    # Return information.
    return pail


# TODO: TCW; 27 March 2023
# TODO: generate the "expected" values by some distribution
# https://www.broadinstitute.org/diabetes-genetics-initiative/plotting-genome-wide-association-results
# https://www.broadinstitute.org/files/shared/diabetes/scandinavs/qqplot.R
# https://physiology.med.cornell.edu/people/banfelder/qbio/resources_2013/2013_1_Mezey.pdf
# expectation = range(1, int(len(probabilities_sort), 1))

def create_qq_plot(
    name=None,
    table=None,
):
    """
    Creates a QQ Plot for a table of GWAS summary statistics.

    arguments:
        name (str): name of table and plot
        table (object): Pandas data-frame table

    raises:

    returns:
        (object): plot figure object

    """

    # Extract and prepare probability values.
    probabilities = table["P"].dropna().to_numpy()
    probabilities_log = numpy.log10(probabilities)
    probabilities_sort = numpy.sort(
        probabilities_log,
        axis=0,
        kind="stable",
    )
    # Prepare expectation values.
    expectations = range(
        1,
        int(len(probabilities_sort),
        1,
    )
    expectations_log = numpy.log10(
        numpy.divide(expectations, (int(len(probabilities_sort)) + 1))
    )

    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = plot.plot_qq_gwas(
        array=probabilities_sort,
        title=name,
        fonts=fonts,
        colors=colors,
        label_title=name, # ""
    )
    # Return information.
    return figure


def drive_create_qq_plots(
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
        (dict<object>): collection of plot objects

    """

    # Collect information for plots.
    pail_plots = dict()
    # Iterate on tables.
    for name in pail_tables.keys():
        table = pail_tables[name].copy(deep=True)
        # Create plot.
        pail_plots[name] = create_qq_plot(
            name=name,
            table=table,
        )
        pass
    # Report.
    if report:
        utility.print_terminal_partition(level=3)
        print("report: ")
        print("drive_create_qq_plots()")
        utility.print_terminal_partition(level=4)
        print("count of tables: " + str(len(pail_tables.keys())))
        print("count of plots: " + str(len(pail_plots.keys())))
        pass
    # Return information.
    return pail_plots


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
    path_directory_product=None,
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
        path_directory_product (str): path to product directory

    raises:

    returns:

    """

    # Read from source files the tables for polygenic scores.
    pail_tables = read_source_directory_files_gwas(
        path_directory_parent=path_directory_source,
        name_file_child_prefix=name_file_source_prefix,
        name_file_child_suffix=name_file_source_suffix,
        name_file_child_not=name_file_source_not,
        report=False,
    )

    # Create QQ Plots.
    pail_plots = drive_create_qq_plots(
        pail_tables=pail_tables,
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
    path_directory_product = sys.argv[5]

    # Call function for procedure.
    execute_procedure(
        path_directory_source=path_directory_source,
        name_file_source_prefix=name_file_source_prefix,
        name_file_source_suffix=name_file_source_suffix,
        name_file_source_not=name_file_source_not,
        path_directory_product=path_directory_product,
    )

    pass



#
