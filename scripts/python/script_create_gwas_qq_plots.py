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


def read_organize_probabilities_from_gwas_table(
    path_table=None,
    name_table=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Use a variable of type "float64" to store probabilities in order to preserve
    very small values.
    Do not read in other columns in order to reduce use of memory.

    arguments:
        path_table (str): path to file for table
        name_table (str): name of table that distinguishes it from all others
        report (bool): whether to print reports

    raises:

    returns:
        (object): NumPy array of probability values (p-values) from summary
            statistics for a GWAS

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("Reading GWAS summary statistics:")
        print(path_table)
        utility.print_terminal_partition(level=4)
    # Count lines in text file.
    path_table_temporary = str(path_table + ".temporary")
    utility.decompress_file_gzip(
        path_file_source=path_table,
        path_file_product=path_table_temporary,
    )
    count_lines = utility.count_file_text_lines(
        path_file=path_table_temporary,
    )
    utility.remove_file(path=path_table_temporary)
    # Determine extent of precision.
    types_columns = dict()
    if (count_lines > 1.3E7):
        # Precision of type float64: 2.2E-308 - 1.7E+308
        types_columns["P"] = "float32"
    else:
        types_columns["P"] = "float64"

    # Read information from file.
    probabilities = pandas.read_csv(
        path_table,
        sep=" ", # white space delimiter
        header=0,
        usecols=["P"],
        dtype=types_columns,
        na_values=["nan", "na", "NAN", "NA",],
        compression="infer",
    )["P"].dropna().to_numpy()
    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("Count of lines in table:")
        print(count_lines)
        print("Count of probabilities from GWAS summary statistics:")
        print(probabilities.size)
        utility.print_terminal_partition(level=4)
    # Return information.
    return probabilities


def create_qq_plot(
    name=None,
    probabilities=None,
):
    """
    Creates a QQ Plot for a table of GWAS summary statistics.

    arguments:
        name (str): name of table and plot
        probabilities (object): NumPy array of probability values (p-values)
            from summary statistics for a GWAS

    raises:

    returns:
        (object): plot figure object

    """

    # Define fonts.
    fonts = pplot.define_font_properties()
    # Define colors.
    colors = pplot.define_color_properties()
    # Create figure.
    figure = pplot.plot_scatter_qq_gwas(
        probabilities=probabilities,
        title=name,
        label="",
        fonts=fonts,
        colors=colors,
    )
    # Return information.
    return figure


def drive_read_gwas_create_write_qq_plots(
    path_directory_source=None,
    name_file_child_prefix=None,
    name_file_child_suffix=None,
    name_file_child_not=None,
    path_directory_product=None,
    report=None,
):
    """
    Reads summary statistics from file, creates QQ Plots, and writes figures to
    file.

    arguments:
        path_directory_source (str): path to parent directory in which to find
            child files
        name_file_child_prefix (str): prefix in name by which to recognize
            relevant child files within parent directory
        name_file_child_suffix (str): suffix in name by which to recognize
            relevant child files within parent directory
        name_file_child_not (str): character string in names of files to exclude
        path_directory_product (str): path to product directory
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): list of base names of product files

    """

    # Read all matching files within parent directory and organize paths to
    # these files.
    paths = utility.read_paths_match_child_files_within_parent_directory(
        path_directory_parent=path_directory_source,
        name_file_child_prefix=name_file_source_prefix,
        name_file_child_suffix=name_file_source_suffix,
        name_file_child_not=name_file_child_not,
        report=report,
    )
    # Read files as Pandas dataframe tables.
    # Iterate on names of files to read and organize tables.
    # Collect names of product files.
    names_product = list()
    for path in paths:
        # Extract name of file and table that distinguishes it from all others.
        name_file = os.path.basename(path)
        name = name_file.replace(str(name_file_child_suffix), "")
        names_product.append(name)
        # Read file and organize information in table.
        probabilities = read_organize_probabilities_from_gwas_table(
            path_table=path,
            name_table=name,
            report=report,
        )
        # Create plot.
        figure = create_qq_plot(
            name=name,
            probabilities=probabilities,
        )
        # Write product information to file.
        pplot.write_product_plot_figure(
            name=name,
            figure=figure,
            path_parent=path_directory_product,
        )
        pass
    # Return information.
    return names_product


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

    # Drive procedure sequentially.
    names_product = drive_read_gwas_create_write_qq_plots(
        path_directory_source=path_directory_source,
        name_file_child_prefix=name_file_source_prefix,
        name_file_child_suffix=name_file_source_suffix,
        name_file_child_not=name_file_source_not,
        path_directory_product=path_directory_product,
        report=True,
    )
    # Report.
    utility.print_terminal_partition(level=4)
    print("Base names of product files:")
    print(names_product)
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
