"""
...
"""

################################################################################
# Notes

# TODO: TCW; 21 May 2024
# TODO: need to update execution structure and context to resemble the cleaner
# convention from "extract_ldsc_heritability_correlation.py"



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
import partner.utility as putility

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source_file_table(
    path_file_source=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_file_source (str): path to file for table
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of polygenic scores across genotypes

    """

    # Read information from file.
    types_columns = dict()
    types_columns["label"] = "string"
    types_columns["heritability_observed"] = "float"
    types_columns["standard_error_observed"] = "float"
    types_columns["prevalence_sample"] = "float"
    types_columns["prevalence_population"] = "float"
    table = pandas.read_csv(
        path_file_source,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=["nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",],
    )
    # Organize information in table.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    # Report.
    if report:
        putility.print_terminal_partition(level=4)
        print("Read source table.")
        print("Source table shape: " + str(table.shape))
        print("Columns:")
        print(table.columns.to_list())
        count_rows = table.shape[0]
        print("Table rows: " + str(count_rows))
        putility.print_terminal_partition(level=4)
    # Return information.
    return table


def calculate_heritability_observed_liability(
    value_observed=None,
    prevalence_sample=None,
    prevalence_population=None,
):
    """
    Calculate the conversion of the SNP heritability and its standard error from
    the observed scale to the liability scale.

    Consider implementing this conversion using a different distribution other
    than the normal distribution. Consider the Gamma or log-normal
    distributions.

    continuous normal random distribution: "scipy.stats.norm()" in Python SciPy
    probability density function (PDF): "dnorm()" in R; "pdf()" in Python SciPy
    cumulative density function (CDF): "pnorm()" in R; "cdf()" in Python SciPy
    quantile inverse CDF: "qnorm()" in R; "ppf()" in Python SciPy

    Converstion of SNP heritabiltiy from observed scale to liability scale:
    1. Equation 23 in Lee et al, AJHG, 2011 (PubMed:21376301)
    2. Neale Lab, "Heritability 201: Types of heritability and how we
        estimate it", 20 September 2017, <http://www.nealelab.is/blog/2017/9/13/
        heritability-201-types-of-heritability-and-how-we-estimate-it>.
    3. Raymond Walters, "Heritability estimate", Google Groups, 27 July 2016,
        <https://groups.google.com/g/ldsc_users/c/tcSMjp7d0Zk?pli=1>.

    arguments:
        value_observed (float): value of heritability or standard error on
            observed scale
        prevalence_sample (float): prevalence of trait in sample or study cohort
        prevalence_population (float): prevalence of trait in population

    raises:

    returns:
        (float): value of heritability or standard error on liability scale

    """

    # Calculate value on liability scale.
    k = prevalence_population
    p = prevalence_sample
    distribution = scipy.stats.norm() # continuous normal random distribution
    z = distribution.pdf(distribution.ppf(k))
    factor = ((k*(1-k))/(z**2)) * ((k*(1-k))/(p*(1-p)))
    value_liability = value_observed * factor
    # Return information.
    return value_liability


def write_product_table(
    table=None,
    path_file=None,
):
    """
    Writes product information in a Pandas data frame table to file.

    arguments:
        table (object): Pandas data frame table
        path_file (str): path to product file

    raises:

    returns:

    """

    # Write information to file.
    table.to_csv(
        path_or_buf=path_file,
        sep="\t",
        header=True,
        index=False,
        na_rep="NA",
    )
    pass


################################################################################
# Procedure


def execute_procedure(
    path_file_source=None,
    path_file_product=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_file_source (str): path to source file
        path_file_product (str): path to product file

    raises:

    returns:

    """

    # Read source file.
    table = read_source_file_table(
        path_file_source=path_file_source,
        report=True,
    )

    # Convert SNP heritabilities from observed scale to liability scale.
    table["heritability_liability"] = table.apply(
        lambda row:
            calculate_heritability_observed_liability(
                value_observed=row["heritability_observed"],
                prevalence_sample=row["prevalence_sample"],
                prevalence_population=row["prevalence_population"],
            ),
        axis="columns", # apply function to each row
    )
    table["standard_error_liability"] = table.apply(
        lambda row:
            calculate_heritability_observed_liability(
                value_observed=row["standard_error_observed"],
                prevalence_sample=row["prevalence_sample"],
                prevalence_population=row["prevalence_population"],
            ),
        axis="columns", # apply function to each row
    )

    # Report.
    putility.print_terminal_partition(level=2)
    print(table)

    # Write.
    write_product_table(
        table=table,
        path_file=path_file_product,
    )
    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    name_script = sys.argv[0]
    path_file_source = sys.argv[1]
    path_file_product = sys.argv[2]

    # Call function for procedure.
    execute_procedure(
        path_file_source=path_file_source,
        path_file_product=path_file_product,
    )

    pass



#
