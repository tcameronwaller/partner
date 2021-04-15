
"""
...

This module collects and organizes information about heritability estimates for
metabolites.

Source of GWAS summary statistics is Shin et al, Nature Genetics, 2014
(PubMed:24816252). Metabolite identifiers are from Metabolon Inc.

Method for estimation of heritability from GWAS summary statistics was linkage
disequilibrium score (LDSC) regression in tool LD SCore
(https://github.com/bulik/ldsc).

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Import modules from specific path without having to install a general package
# I would have to figure out how to pass a path variable...
# https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path


# Standard

import sys
#print(sys.path)
import os
import math
import statistics
import pickle
import copy
import random
import itertools
import time

# Relevant

import numpy
import pandas
import scipy.stats

# Custom
import utility


###############################################################################
# Functionality


##########
# Initialization


def initialize_directories(
    restore=None,
    path_dock=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        restore (bool): whether to remove previous versions of data
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = path_dock
    paths["organization"] = os.path.join(path_dock, "organization")
    paths["coombes_polygene"] = os.path.join(
        path_dock, "coombes_prs_gems_gain_mayo_all_2020-10-13"
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["organization"])
    # Initialize directories.
    utility.create_directories(
        path=paths["organization"]
    )
    # Return information.
    return paths


def read_source(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_table_reference_shin_2014 = os.path.join(
        path_dock, "metabolite_reference",
        "24816252_shin_2014", "table_metabolite_reference.tsv"
    )
    path_table_metabolite_heritabilities = os.path.join(
        path_dock, "heritability_correlation_2021-04-12",
        "table_shin_2014_heritabilities.tsv"
    )

    # Read information from file.
    table_reference_shin_2014 = pandas.read_csv(
        path_table_reference_shin_2014,
        sep="\t",
        header=0,
        #dtype="string",
    )
    table_metabolite_heritabilities = pandas.read_csv(
        path_table_metabolite_heritabilities,
        sep="\t",
        header=0,
        #dtype="string",
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report from read_source()")
        print(table_reference_shin_2014)
        utility.print_terminal_partition(level=2)
    # Compile and return information.
    return {
        "table_reference_shin_2014": table_reference_shin_2014,
        "table_metabolite_heritabilities": table_metabolite_heritabilities,
    }


##########
# Organize tables' formats


def organize_metabolite_reference_table(
    table=None,
    identifier=None,
    name=None,
    identity=None,
):
    """
    Organizes information about general attributes.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        identifier (str): name of column for metabolite identifier
        name (str): name of column for metabolite biochemical name
        identity (str): name of column as binary logical indicator of whether
            the metabolite has a known identity

    raises:

    returns:
        (dict): collection of information about phenotype variables

    """

    # Copy data.
    table = table.copy(deep=True)
    # Translate column names.
    translations = dict()
    translations[identifier] = "identifier"
    translations[name] = "name"
    translations[identity] = "identity"
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Select relevant columns.
    table = table.loc[
        :, table.columns.isin(["identifier", "name", "identity"])
    ]
    # Organize table.
    table.reset_index(
        level=None,
        inplace=True
    )
    #table["identity"].astype("float")
    table["identity"] = pandas.to_numeric(
        table["identity"],
        errors="coerce", # force any invalid values to missing or null
        downcast="float",
    )
    table["identifier"].astype("string")
    #table.set_index(
    #    "identifier",
    #    drop=True,
    #    inplace=True,
    #)
    # Return information.
    return table


def simplify_name_text_as_key(
    name=None,
):
    """
    Simplifies the text in a name for use as a key in table merge.

    arguments:
        name (str): name for simplification

    raises:

    returns:
        (str): simple name for use as key in merge

    """

    # Remove white space.
    name_strip = str(name).strip()
    # Revert custom string changes by Brandon Coombes.
    name_strip = name_strip.replace(
        "docosapentaenoate__n3",
        "docosapentaenoate (n3 DPA; 22:5n3)"
    )
    name_strip = name_strip.replace(
        "glycerol_3-phosphate",
        "glycerol 3-phosphate (G3P)"
    )
    name_strip = name_strip.replace(
        "eicosenoate__20:1n9",
        "eicosenoate (20:1n9 or 11)"
    )
    name_strip = name_strip.replace(
        "15-methylpalmitate__isobar",
        "15-methylpalmitate (isobar with 2-methylpalmitate)"
    )
    name_strip = name_strip.replace(
        "dimethylarginine__SDMA",
        "dimethylarginine (SDMA + ADMA)"
    )
    name_strip = name_strip.replace(
        "dehydroisoandrosterone_sulfate",
        "dehydroisoandrosterone sulfate (DHEA-S)"
    )
    name_strip = name_strip.replace(
        "linolenate_[alpha",
        "linolenate [alpha or gamma; (18:3n3 or 6)]"
    )
    name_strip = name_strip.replace(
        "dihomo-linolenate__20:3n3",
        "dihomo-linolenate (20:3n3 or n6)"
    )
    name_strip = name_strip.replace(
        "ascorbate__Vitamin", "ascorbate (Vitamin C)"
    )
    name_strip = name_strip.replace("bilirubin__E,Z", "bilirubin (E,Z or Z,E)")
    # Define translation table.
    translations = str.maketrans("", "", " ,:()_-*")
    # Translate string characters.
    name_simple = name_strip.translate(translations)
    # Return.
    return name_simple


def select_table_metabolites_valid_identities_heritabilities(
    table=None,
    table_reference=None,
    threshold_metabolite_heritability=None,
    report=None,
):
    """
    Selects identifiers of metabolites from Metabolon with valid identities.

    arguments:
        table (object): Pandas data frame of metabolites' heritability estimates
            and genetic correlation estimates against a phenotype of interest
        table_reference (object): Pandas data frame of metabolites' identifiers
            and names from study
        threshold_metabolite_heritability (float): threshold for metabolite
            heritability
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of information about metabolites, their identifiers,
            and their names

    """

    # Select metabolites with valid identities.
    table_identity = table_reference.loc[
        (table_reference["identity"] > 0.5), :
    ]
    table_identity.reset_index(
        level=None,
        inplace=True
    )
    table_identity.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    metabolites_identity = table_identity.index.to_list()
    # Select table rows for metabolites with valid identities.
    table = table.loc[
        table.index.isin(metabolites_identity), :
    ]
    # Select table rows for metabolites with valid heritability estimates.
    table["heritability"] = pandas.to_numeric(
        table["heritability"],
        errors="coerce", # force any invalid values to missing or null
        downcast="float",
    )
    table = table.loc[
        (table["heritability"] >= threshold_metabolite_heritability), :
    ]
    metabolites_heritability = table.index.to_list()
    # Report.
    if report:
        # Column name translations.
        utility.print_terminal_partition(level=2)
        print("select_table_metabolites_valid_identities_heritabilities()")
        utility.print_terminal_partition(level=3)
        print(
            "Count of identifiable metabolites: " +
            str(len(metabolites_identity))
        )
        print(
            "Count of identifiable metabolites with valid heritability: " +
            str(len(metabolites_heritability))
        )
        utility.print_terminal_partition(level=3)
        print(table)
    # Return information.
    return table


def organize_metabolites_heritabilities_polygenic_associations_table(
    table=None,
    table_reference=None,
    threshold_metabolite_heritability=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        table (object): Pandas data frame of metabolites' heritability estimates
            and genetic correlation estimates against a phenotype of interest
        table_reference (object): Pandas data frame of metabolites' identifiers
            and names from study
        threshold_metabolite_heritability (float): threshold for metabolite
            heritability

    raises:

    returns:
        (object): Pandas data frame of metabolites' heritability estimates and
            genetic correlation estimates against a phenotype of interest

    """

    # Copy information.
    table_reference = table_reference.copy(deep=True)
    table = table.copy(deep=True)
    # Organize table.
    #table = table.loc[
    #    (len(str(table["identifier"])) > 0), :
    #]
    # Rename columns.
    translations = dict()
    translations["name_reference"] = "name"
    translations["identity_reference"] = "identity"
    translations["PRS"] = "name_prs"
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Drop rows with missing values.
    table.dropna(
        axis="index",
        how="any",
        subset=["identifier", "name", "identity", "heritability", "pval"],
        inplace=True,
    )
    # Select and sort relevant table columns.
    columns = [
        "identifier",
        "name", "name_prs", "identity",
        "heritability", "heritability_standard_error",
        "BETA", "SE", "lcl", "ucl", "pval", #"r2.mayo", "r2.gain",
    ]
    table = table.loc[
        :, table.columns.isin(columns)
    ]
    table = table[[*columns]]
    # Sort table rows.
    table.sort_values(
        by=["BETA",],
        axis="index",
        ascending=False,
        na_position="last",
        inplace=True,
    )
    table.sort_values(
        by=["pval",],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    table.reset_index(
        level=None,
        inplace=True
    )
    table.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Filter metabolites.
    table = select_table_metabolites_valid_identities_heritabilities(
        table=table,
        table_reference=table_reference,
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=True,
    )
    # Calculate False Discovery Rates (FDRs).
    table = utility.calculate_table_false_discovery_rates(
        threshold=0.05,
        probability="pval",
        discovery="discovery",
        significance="significance",
        table=table,
    )
    # Return information.
    return table


def read_organize_polygenic_metabolite_phenotype_regression_table(
    file=None,
    path_source_directory=None,
    threshold_metabolite_heritability=None,
    table_reference=None,
    table_heritability=None,
    report=None,
):
    """
    Reads and organize tables from regressions between polygenic estimate
    metabolites and phenotypes.

    arguments:
        file (str): name of a file
        path_source_directory (str): path to source parent directory for files
            of tables summarizing metabolite associations to phenotypes
        threshold_metabolite_heritability (float): threshold for metabolite
            heritability
        table_reference (object): Pandas data frame of metabolites' identifiers
            and names from study
        table_heritability (object): Pandas data frame of metabolites'
            heritability estimates
        report (bool): whether to print reports

    raises:

    returns:
        (object): table summarizing regression associations between polygenic
            estimates of metabolites and phenotype

    """

    # Copy information.
    table_reference_original = table_reference.copy(deep=True)
    table_reference = table_reference.copy(deep=True)
    table_heritability = table_heritability.copy(deep=True)
    # Define path to file.
    path_file = os.path.join(
        path_source_directory, file
    )
    # Read information from file.
    table_association = pandas.read_csv(
        path_file,
        sep=",",
        header=0,
        #dtype="string",
    )
    # Simplify metabolite names for use as keys.
    table_association["name_key"] = table_association.apply(
        lambda row:
            simplify_name_text_as_key(
                name=row["PRS"],
            ),
        axis="columns", # apply across rows
    )
    table_reference["name_key"] = table_reference.apply(
        lambda row:
            simplify_name_text_as_key(
                name=row["name"],
            ),
        axis="columns", # apply across rows
    )
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_identifier = table_reference.merge(
        table_association,
        how="outer", # "outer"
        left_on="name_key",
        right_on="name_key",
        suffixes=("_association", "_reference"),
    )
    table_merge = table_identifier.merge(
        table_heritability,
        how="outer", # "outer"
        left_on="identifier",
        right_on="identifier",
        suffixes=("_reference", "_heritability"),
    )
    # Organize the summary collection table.
    table = organize_metabolites_heritabilities_polygenic_associations_table(
        table=table_merge,
        table_reference=table_reference_original,
        threshold_metabolite_heritability=threshold_metabolite_heritability,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report from...")
        print("read_organize_polygenic_metabolite_phenotype_regression_table()")
        print(table)
        utility.print_terminal_partition(level=2)
    # Return.
    return table


def read_organize_polygenic_metabolite_phenotype_regression_tables(
    threshold_metabolite_heritability=None,
    path_source_directory=None,
    table_reference_shin_2014=None,
    table_metabolite_heritabilities=None,
    report=None,
):
    """
    Reads and organize tables from regressions between polygenic estimate
    metabolites and phenotypes.

    arguments:
        threshold_metabolite_heritability (float): threshold for metabolite
            heritability
        path_source_directory (str): path to source parent directory for files
            of tables summarizing metabolite associations to phenotypes
        table_reference_shin_2014 (object): Pandas data frame of metabolites'
            identifiers and names from study
        table_metabolite_heritabilities (object): Pandas data frame of
            metabolites' heritability estimates
        report (bool): whether to print reports

    raises:

    returns:
        (dict): tables summarizing regression associations between polygenic
            estimates of metabolites and phenotypes

    """

    # Organize metabolite reference table.
    table_reference = organize_metabolite_reference_table(
        table=table_reference_shin_2014,
        identifier="identifier_study",
        name="name",
        identity="identity",
    )
    # Iterate on original tables.
    # Organize tables and collect them in new format.
    pail = dict()
    # Collect names of files for metabolites' heritabilities.
    files = utility.extract_directory_file_names(path=path_source_directory)
    files_relevant = list(filter(
        lambda content: ("GEMS" in content), files
    ))
    for file in files_relevant:
        name = str(file.replace(".csv", ""))
        print(name)
        pail[name] = (
            read_organize_polygenic_metabolite_phenotype_regression_table(
                file=file,
                path_source_directory=path_source_directory,
                threshold_metabolite_heritability=threshold_metabolite_heritability,
                table_reference=table_reference,
                table_heritability=table_metabolite_heritabilities,
                report=report,
        ))
        pass
    # Return information.
    return pail


def write_product_table(
    name=None,
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        name (str): base name for file
        information (object): information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_table = os.path.join(
        path_parent, str(name + ".tsv")
    )
    # Write information to file.
    information.to_csv(
        path_or_buf=path_table,
        sep="\t",
        header=True,
        index=True,
    )
    pass


def write_product_tables(
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        path_parent (str): path to parent directory

    raises:

    returns:

    """

    for name in information.keys():
        write_product_table(
            name=name,
            information=information[name],
            path_parent=path_parent,
        )
    pass


def write_product(
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Summary tables in format.
    write_product_tables(
        information=information["tables"],
        path_parent=paths["organization"],
    )
    pass



###############################################################################
# Procedure


def execute_procedure(
    path_dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # 3. read in Coombes' metabolite regression tables (all)
    # 4. iterate on Coombes' metabolite regression tables
    # 5. match metabolite names to Shin 2014 metabolite identifiers
    # 6. merge Coombes' metabolite regression information with metabolite heritabilities
    # 7. filter metabolites by whether identifiable and SNP-heritability > 0.05
    # 8. calculate Benjamini-Hochberg False-Discovery Rates


    # Report version.
    utility.print_terminal_partition(level=1)
    print(path_dock)
    print("version check: 1")
    # Pause procedure.
    time.sleep(5.0)

    # Initialize directories.
    paths = initialize_directories(
        restore=False,
        path_dock=path_dock,
    )
    # Read source information from file.
    source = read_source(
        path_dock=path_dock,
        report=True,
    )

    #print(source["table_reference_shin_2014"])
    #print(source["table_metabolite_heritabilities"])

    # Read and organize tables for regressions between polygenic estimate
    # metabolites and phenotypes.
    pail = read_organize_polygenic_metabolite_phenotype_regression_tables(
        threshold_metabolite_heritability=0.05, # metabolite heritability
        path_source_directory=paths["coombes_polygene"],
        table_reference_shin_2014=source["table_reference_shin_2014"],
        table_metabolite_heritabilities=(
            source["table_metabolite_heritabilities"]
        ),
        report=True,
    )

    # Collect information.
    information = dict()
    information["tables"] = pail
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )




    pass



if (__name__ == "__main__"):
    execute_procedure()
