
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

# Relevant

import numpy
import pandas
import scipy.stats

# Custom
import promiscuity.utility as utility


###############################################################################
# Functionality


##########
# Initialization


def initialize_heritability_directories(
    heritability_studies=None,
    paths=None,
    path_dock=None,
    restore=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        heritability_studies (list<str>): identifiers of studies with
            heritability reports
        paths (dict<str>): collection of paths to directories for procedure's
            files
        path_dock (str): path to dock directory for source and product
            directories and files
        restore (bool): whether to remove previous versions of data

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    paths = copy.deepcopy(paths)
    paths["heritability"] = dict()
    for study in heritability_studies:
        paths["heritability"][study] = os.path.join(
            path_dock, "heritability", study
        )
    return paths


def initialize_correlation_directories(
    primary_studies=None,
    secondary_studies=None,
    paths=None,
    path_dock=None,
    restore=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        primary_studies (list<str>): identifiers of studies for primary
            phenotypes, first in correlation hierarchy
        secondary_studies (list<str>): identifiers of studies for secondary
            phenotypes, second in correlation hierarchy
        paths (dict<str>): collection of paths to directories for procedure's
            files
        path_dock (str): path to dock directory for source and product
            directories and files
        restore (bool): whether to remove previous versions of data

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    paths = copy.deepcopy(paths)
    paths["genetic_correlation"] = os.path.join(
        path_dock, "genetic_correlation",
    )
    paths["correlation"] = dict()
    for study_first in primary_studies:
        for study_second in secondary_studies:
            paths["correlation"][study_first] = dict()
            paths["correlation"][study_first][study_second] = (
                os.path.join(
                    path_dock, "genetic_correlation",
                    study_first, study_second
                )
            )
    return paths


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
    heritability_studies = [
        "24816252_shin_2014", "33437055_panyard_2021",
        "30124842_yengo_2018", "30239722_pulit_2018",
        "30482948_walters_2018", "30718901_howard_2019", "31043756_stahl_2019",
        "00000000_pgc3_2021_all", "00000000_pgc3_2021_bd1",
        "00000000_pgc3_2021_bd2",
        "29906448_ruderfer_2018_scz_bpd_vs_ctl",
        "29906448_ruderfer_2018_scz_vs_bpd",
        "00000000_ripke_2021",
    ]
    primary_studies = [
        "30124842_yengo_2018", "30239722_pulit_2018",
        "30482948_walters_2018", "30718901_howard_2019", "31043756_stahl_2019",
        "00000000_pgc3_2021_all", "00000000_pgc3_2021_bd1",
        "00000000_pgc3_2021_bd2",
        "29906448_ruderfer_2018_scz_bpd_vs_ctl",
        "29906448_ruderfer_2018_scz_vs_bpd",
        "00000000_ripke_2021",
    ]
    secondary_studies = [
        "24816252_shin_2014",
    ]
    paths = initialize_heritability_directories(
        heritability_studies=heritability_studies,
        paths=paths,
        path_dock=path_dock,
        restore=restore,
    )
    paths = initialize_correlation_directories(
        primary_studies=primary_studies,
        secondary_studies=secondary_studies,
        paths=paths,
        path_dock=path_dock,
        restore=restore,
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
        path_dock, "parameters", "bipolar_metabolism", "metabolite_reference",
        "24816252_shin_2014", "table_metabolite_reference.tsv"
    )
    path_table_reference_panyard_2021 = os.path.join(
        path_dock, "parameters", "bipolar_metabolism", "metabolite_reference",
        "33437055_panyard_2021", "table_metabolite_reference.tsv"
    )
    # Read information from file.
    table_reference_shin_2014 = pandas.read_csv(
        path_table_reference_shin_2014,
        sep="\t",
        header=0,
        #dtype="string",
    )
    table_reference_panyard_2021 = pandas.read_csv(
        path_table_reference_panyard_2021,
        sep="\t",
        header=0,
        #dtype="string",
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_reference_panyard_2021)
        utility.print_terminal_partition(level=2)
    # Compile and return information.
    return {
        "table_reference_shin_2014": table_reference_shin_2014,
        "table_reference_panyard_2021": table_reference_panyard_2021,
    }


# Heritability


def read_collect_organize_metabolites_heritabilities_studies(
    table_reference_panyard_2021=None,
    paths=None,
    report=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        table_reference_panyard_2021 (object): Pandas data frame of metabolites'
            identifiers and names from study
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): heritability estimations for metabolites from multiple GWAS

    """

    # Collect metabolites' heritabilities from each GWAS.
    pail = dict()
    #pail["table_shin_2014"] =
    pail["table_panyard_2021"] = read_collect_metabolites_heritabilities(
        table_reference=table_reference_panyard_2021,
        path_parent=paths["heritability_panyard_2021"],
        report=report,
    )
    # Return information.
    return pail


# Genetic correlation.


def read_collect_organize_metabolites_correlations_studies(
    table_reference_shin_2014=None,
    table_reference_panyard_2021=None,
    paths=None,
    report=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        table_reference_shin_2014 (object): Pandas data frame of metabolites'
            identifiers and names from study
        table_reference_panyard_2021 (object): Pandas data frame of metabolites'
            identifiers and names from study
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): heritability estimations for metabolites from multiple GWAS

    """

    # Collect metabolites' heritabilities from each GWAS.
    pail = dict()
    #pail["table_shin_2014"] =
    pail["table_yengo_2018_shin_2014"] = read_collect_metabolites_genetic_correlations(
        table_reference=table_reference_shin_2014,
        path_source_directory=(
            paths["correlation"]["30124842_yengo_2018"]["24816252_shin_2014"]
        ),
        report=report,
    )
    pail["table_pulit_2018_shin_2014"] = read_collect_metabolites_genetic_correlations(
        table_reference=table_reference_shin_2014,
        path_source_directory=(
            paths["correlation"]["30239722_pulit_2018"]["24816252_shin_2014"]
        ),
        report=report,
    )
    # Return information.
    return pail


# Combination


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
    table["identity"].astype("float")
    table["identifier"].astype("string")
    table.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Return information.
    return table


def read_extract_metabolite_heritability(
    file=None,
    path_source_directory=None,
):
    """
    Reads and extracts information from log of LDSC for heritability estimation
    from GWAS summary statistics.

    arguments:
        file (str): name of a file
        path_source_directory (str): path to source parent directory for files
            with genetic correlation estimations for metabolites

    raises:

    returns:
        (dict): information about estimation of a metabolite's heritability

    """

    # Extract metabolite's identifier.
    identifier = str(
        file.replace("heritability_", "").replace(".log", "")
    )
    # Define path to file.
    path_file = os.path.join(
        path_source_directory, file
    )
    # Initialize variables.
    variants = float("nan")
    heritability = float("nan")
    heritability_error = float("nan")
    ratio = float("nan")
    ratio_error = float("nan")
    # Read relevant lines from file.
    lines = utility.read_file_text_lines(
        path_file=path_file,
        start=22,
        stop=30,
    )
    # Extract information from lines.
    prefix_variants = "After merging with regression SNP LD, "
    suffix_variants = " SNPs remain."
    prefix_heritability = "Total Observed scale h2: "
    prefix_ratio = "Ratio: "
    for line in lines:
        if prefix_variants in line:
            variants = float(
                line.replace(prefix_variants, "").replace(suffix_variants, "")
            )
        elif prefix_heritability in line:
            content = line.replace(prefix_heritability, "")
            contents = content.split(" (")
            heritability_test = contents[0]
            if (not "NA" in heritability_test):
                heritability = float(contents[0])
                heritability_error = float(contents[1].replace(")", ""))
            pass
        elif (
            (not math.isnan(heritability)) and
            (prefix_ratio in line)
        ):
            content = line.replace(prefix_ratio, "")
            contents = content.split(" (")
            ratio_test = contents[0]
            if (not "NA" in ratio_test):
                ratio = float(contents[0])
                ratio_error = float(
                    contents[1].replace(")", "")
                )
            pass
        pass
    # Collect information.
    record = dict()
    record["identifier"] = identifier
    record["heritability_variants"] = variants
    record["heritability"] = heritability
    record["heritability_standard_error"] = heritability_error
    record["heritability_ratio"] = ratio
    record["heritability_ratio_standard_error"] = ratio_error
    # Return information.
    return record


def read_collect_metabolites_heritabilities(
    path_source_directory=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        path_source_directory (str): path to source parent directory for files
            with genetic correlation estimations for metabolites

    raises:

    returns:
        (object): Pandas data frame of metabolites' heritability estimates

    """

    # Collect names of files for metabolites' heritabilities.
    files = utility.extract_directory_file_names(path=path_source_directory)
    files_relevant = list(filter(
        lambda content: ("heritability" in content), files
    ))
    records = list()
    for file in files_relevant:
        record = read_extract_metabolite_heritability(
            file=file,
            path_source_directory=path_source_directory,
        )
        records.append(record)
        pass
    # Organize heritability table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    table.sort_values(
        by=["heritability"],
        axis="index",
        ascending=False,
        inplace=True,
    )
    table.reset_index(
        level=None,
        inplace=True
    )
    table["identifier"].astype("string")
    table.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Return information.
    return table


def read_extract_phenotype_metabolite_genetic_correlation(
    file=None,
    path_source_directory=None,
):
    """
    Reads and extracts information from log of LDSC for estimation of
    genetic correlation between a phenotype of interest and a metabolite.

    phenotype 1: phenotype of interest compared accross all metabolites
    phenotype 2: single metabolite of interest

    arguments:
        file (str): name of a file
        path_source_directory (str): path to source parent directory for files
            with genetic correlation estimations for metabolites

    raises:

    returns:
        (dict): information about estimation of a metabolite's genetic
            correlation to phenotype

    """

    # Extract metabolite's identifier.
    identifier = str(
        file.replace("correlation_", "").replace(".log", "")
    )
    # Define path to file.
    path_file = os.path.join(
        path_source_directory, file
    )
    # Initialize variables.
    variants = float("nan")
    correlation = float("nan")
    correlation_error = float("nan")
    correlation_absolute = float("nan")
    probability = float("nan")
    # Read relevant lines from file.
    lines = utility.read_file_text_lines(
        path_file=path_file,
        start=25,
        stop=57,
    )
    # Extract information from lines.
    prefix_variants = ""
    suffix_variants = " SNPs with valid alleles."
    prefix_correlation = "Genetic Correlation: "
    prefix_probability = "P: "
    for line in lines:
        if suffix_variants in line:
            variants = float(
                line.replace(prefix_variants, "").replace(suffix_variants, "")
            )
        elif prefix_correlation in line:
            content = line.replace(prefix_correlation, "")
            contents = content.split(" (")
            correlation_test = contents[0]
            if (not "nan" in correlation_test):
                correlation = float(contents[0])
                correlation_absolute = math.fabs(correlation)
                correlation_error = float(contents[1].replace(")", ""))
            pass
        elif (
            (not math.isnan(correlation)) and
            (prefix_probability in line)
        ):
            probability = float(line.replace(prefix_probability, ""))
            pass
        pass
    # Collect information.
    record = dict()
    record["identifier"] = identifier
    record["correlation_variants"] = variants
    record["correlation"] = correlation
    record["correlation_standard_error"] = correlation_error
    record["correlation_absolute"] = correlation_absolute
    record["correlation_probability"] = probability
    # Return information.
    return record


def read_collect_metabolites_genetic_correlations(
    path_source_directory=None,
):
    """
    Reads and collects estimations of genetic correlation between phenotype and
        metabolites.

    arguments:
        path_source_directory (str): path to source parent directory for files
            with genetic correlation estimations for metabolites

    raises:

    returns:
        (object): Pandas data frame of metabolites' heritability estimates

    """

    # Collect names of files for metabolites' heritabilities.
    files = utility.extract_directory_file_names(path=path_source_directory)
    files_relevant = list(filter(
        lambda content: ("correlation" in content), files
    ))
    records = list()
    for file in files_relevant:
        record = read_extract_phenotype_metabolite_genetic_correlation(
            file=file,
            path_source_directory=path_source_directory,
        )
        records.append(record)
        pass
    # Organize heritability table.
    table = utility.convert_records_to_dataframe(
        records=records
    )
    table.sort_values(
        by=["correlation_absolute"],
        axis="index",
        ascending=False,
        na_position="last",
        inplace=True,
    )
    table.reset_index(
        level=None,
        inplace=True
    )
    table["identifier"].astype("string")
    table.set_index(
        "identifier",
        drop=True,
        inplace=True,
    )
    # Return information.
    return table


def determine_metabolite_valid_identity(
    name=None,
):
    """
    Determine whether a single metabolite has a valid identity from Metabolon.

    arguments:
        name (str): name of metabolite from Metabolon reference

    raises:

    returns:
        (float): ordinal representation of person's frequency of alcohol
            consumption

    """

    # Determine whether the variable has a valid (non-missing) value.
    if (len(str(name)) > 2):
        # The variable has a valid value.
        if (str(name).strip().lower().startswith("x-")):
            # Metabolite has an indefinite identity.
            identity = 0
        else:
            # Metabolite has a definite identity.
            identity = 1
    else:
        # Name is empty.
        #identity = float("nan")
        identity = 0
    # Return information.
    return identity


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
    metabolites_identity = table_identity.index.to_list()
    # Select table rows for metabolites with valid identities.
    table = table.loc[
        table.index.isin(metabolites_identity), :
    ]
    # Select table rows for metabolites with valid heritability estimates.
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


def organize_metabolites_heritabilities_correlations_table(
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
        probability="correlation_probability",
        discovery="correlation_discovery",
        significance="correlation_significance",
        table=table,
    )
    # Sort table rows.
    table.sort_values(
        by=["correlation_absolute"],
        axis="index",
        ascending=False,
        na_position="last",
        inplace=True,
    )
    table.sort_values(
        by=["correlation_probability",],
        axis="index",
        ascending=True,
        na_position="last",
        inplace=True,
    )
    # Sort table columns.
    columns_sequence = [
        #"identifier",
        "name",
        "phenotype_heritability",
        "phenotype_heritability_error",
        "correlation", "correlation_standard_error",
        "correlation_absolute",
        "correlation_probability",
        "correlation_discovery",
        "correlation_significance",
        "correlation_variants",
        "heritability",
        "heritability_standard_error",
        "heritability_ratio",
        "heritability_ratio_standard_error",
        "heritability_variants",
    ]
    table = table[[*columns_sequence]]
    # Return information.
    return table


def read_collect_combine_study(
    table_reference=None,
    file_phenotype_heritability=None,
    path_phenotype_heritability=None,
    path_metabolite_heritabilities=None,
    path_correlations=None,
    threshold_metabolite_heritability=None,
    report=None,
):
    """
    Reads, collects, and organizes metabolite heritability estimates.

    arguments:
        table_reference (object): Pandas data frame of metabolites' identifiers
            and names from study
        file_phenotype_heritability (str): name of file for heritability
            estimation for phenotype
        path_phenotype_heritability (str): path to source parent directory for
            heritability estimation for phenotype
        path_metabolite_heritabilities (str): path to source parent directory
            for files with heritability estimations for metabolites
        path_correlations (str): path to source parent directory for files with
            genetic correlation estimations for phenotype and metabolites
        threshold_metabolite_heritability (float): threshold for metabolite
            heritability
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of metabolites' heritability estimates and
            genetic correlation estimates against a phenotype of interest

    """

    # Organize metabolite reference table.
    table_reference = organize_metabolite_reference_table(
        table=table_reference,
        identifier="identifier_study",
        name="name",
        identity="identity",
    )

    pail_phenotype = read_extract_metabolite_heritability(
        file=file_phenotype_heritability,
        path_source_directory=path_phenotype_heritability,
    )

    table_correlations = read_collect_metabolites_genetic_correlations(
        path_source_directory=path_correlations,
    )

    table_metabolite_heritabilities = read_collect_metabolites_heritabilities(
        path_source_directory=path_metabolite_heritabilities,
    )

    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge = table_reference.merge(
        table_correlations,
        how="outer",
        left_on="identifier",
        right_on="identifier",
        suffixes=("_reference", "_correlation"),
    )
    table_merge = table_merge.merge(
        table_metabolite_heritabilities,
        how="outer",
        left_on="identifier",
        right_on="identifier",
        suffixes=("_reference", "_heritability"),
    )
    # Introduce columns for phenotype heritability.
    table_merge["phenotype_heritability"] = pail_phenotype["heritability"]
    table_merge["phenotype_heritability_error"] = (
        pail_phenotype["heritability_standard_error"]
    )
    # Organize the summary collection table.
    table = organize_metabolites_heritabilities_correlations_table(
        table=table_merge,
        table_reference=table_reference,
        threshold_metabolite_heritability=threshold_metabolite_heritability,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(path_correlations)
        print(table)
    # Return information.
    return table


def read_collect_combine_phenotype_metabolites_studies(
    table_reference_shin_2014=None,
    table_reference_panyard_2021=None,
    threshold_metabolite_heritability=None,
    paths=None,
    report=None,
):
    """
    Reads, collects, combines, and organizes heritabilities and genetic
    correlations between phenotypes and metabolites.

    arguments:
        table_reference_shin_2014 (object): Pandas data frame of metabolites'
            identifiers and names from study
        table_reference_panyard_2021 (object): Pandas data frame of metabolites'
            identifiers and names from study
        threshold_metabolite_heritability (float): threshold for metabolite
            heritability
        paths (dict<str>): collection of paths to directories for procedure's
            files
        report (bool): whether to print reports

    raises:

    returns:
        (dict): heritability estimations for metabolites from multiple GWAS

    """

    # Collect metabolites' heritabilities from each GWAS.
    pail = dict()
    pail["table_yengo_2018_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=os.path.join(
            paths["dock"], "heritability", "30124842_yengo_2018",
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["30124842_yengo_2018"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )
    pail["table_pulit_2018_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=os.path.join(
            paths["dock"], "heritability", "30239722_pulit_2018",
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["30239722_pulit_2018"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )
    pail["table_stahl_2019_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=os.path.join(
            paths["dock"], "heritability", "31043756_stahl_2019",
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["31043756_stahl_2019"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )
    pail["table_howard_2019_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=os.path.join(
            paths["dock"], "heritability", "30718901_howard_2019",
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["30718901_howard_2019"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )
    pail["table_walters_2018_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=os.path.join(
            paths["dock"], "heritability", "30482948_walters_2018",
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["30482948_walters_2018"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )
    pail["table_pgc3_2021_all_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=(
            paths["heritability"]["00000000_pgc3_2021_all"]
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["00000000_pgc3_2021_all"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )
    pail["table_pgc3_2021_bd1_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=(
            paths["heritability"]["00000000_pgc3_2021_bd1"]
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["00000000_pgc3_2021_bd1"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )
    pail["table_pgc3_2021_bd2_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=(
            paths["heritability"]["00000000_pgc3_2021_bd2"]
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["00000000_pgc3_2021_bd2"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )

    pail["table_ruderfer_2018_scz_bpd_vs_ctl_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=(
            paths["heritability"]["29906448_ruderfer_2018_scz_bpd_vs_ctl"]
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["29906448_ruderfer_2018_scz_bpd_vs_ctl"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )
    pail["table_ruderfer_2018_scz_vs_bpd_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=(
            paths["heritability"]["29906448_ruderfer_2018_scz_vs_bpd"]
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["29906448_ruderfer_2018_scz_vs_bpd"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )

    pail["table_ripke_2021_shin_2014"] = read_collect_combine_study(
        table_reference=table_reference_shin_2014,
        file_phenotype_heritability="heritability_report.log",
        path_phenotype_heritability=(
            paths["heritability"]["00000000_ripke_2021"]
        ),
        path_metabolite_heritabilities=(
            paths["heritability"]["24816252_shin_2014"]
        ),
        path_correlations=(
            paths["correlation"]["00000000_ripke_2021"]["24816252_shin_2014"]
        ),
        threshold_metabolite_heritability=threshold_metabolite_heritability,
        report=report,
    )

    # Return information.
    return pail


def write_product_study_table(
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
        index=False,
    )
    pass


def write_product_studies(
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
        write_product_study_table(
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

    # Cohort tables in PLINK format.
    write_product_studies(
        information=information["studies"],
        path_parent=paths["genetic_correlation"],
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

    utility.print_terminal_partition(level=1)
    print(path_dock)
    print("version check: 6")

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
    # Read, collect, and combine estimations of heritability and genetic
    # correlations between phenotypes and metabolites.
    pail_studies = read_collect_combine_phenotype_metabolites_studies(
        table_reference_shin_2014=source["table_reference_shin_2014"],
        table_reference_panyard_2021=source["table_reference_panyard_2021"],
        threshold_metabolite_heritability=0.05, # metabolite heritability
        paths=paths,
        report=True,
    )

    # Collect information.
    information = dict()
    information["studies"] = pail_studies
    # Write product information to file.
    write_product(
        paths=paths,
        information=information
    )

    pass



if (__name__ == "__main__"):
    execute_procedure()
