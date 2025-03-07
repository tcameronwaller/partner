"""
Query MyGene.info for annotation information about genes.

Author:

    T. Cameron Waller, Ph.D.
    tcameronwaller@gmail.com
    Rochester, Minnesota 55902
    United States of America

License:

    This module file is part of the project package directory 'partner'
    (https://github.com/tcameronwaller/partner/).

    Project 'partner' supports data analysis in multiple other projects.
    Copyright (C) 2024 Thomas Cameron Waller

    The code within project 'partner' is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the GNU
    General Public License, or (at your option) any later version.

    The code within project 'partner' is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with project 'partner'. If not, see <http://www.gnu.org/licenses/>.
"""

################################################################################
# Author: T. Cameron Waller, Ph.D.
# Date, first execution: 31 October 2024
# Date, last execution: 1 November 2024
# Review: TCW; 1 November 2024
################################################################################
# Note

# tool: MyGene.info
# site, home: https://www.mygene.info/
# site, live API: https://www.mygene.info/v3/api
# site, documentation: https://docs.mygene.info/en/latest/index.html
# site, host PyPi: https://pypi.org/project/mygene/

# For simplicity, the current design of query only supports a single type of
# source identifier or name and a single type of product identifier or name.
# It is possible to use multiple types of identifiers or names in the query
# source and to deliver multiple annotations, attributes, or features from each
# gene's record.

# Simple process of converting delimiters in a list of gene identifiers.
# $ sed 's/,/\n/g' input.txt > output.txt

################################################################################
# Installation and importation

# Standard
import sys
import os
import copy

# Relevant
import pandas
import scipy
import numpy
import mygene

# Custom
import partner.utility as putly

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_source_extract_items(
    path_file=None,
    delimiter=None,
    report=None,
):
    """
    Reads and extracts from a source file the identifiers or names of items.

    arguments:
        path_file_source (str): path to source file in text format from which
            to read identifiers or names of items
        delimiter (str): delimiter between items in text representation of list
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers or names of genes for a query

    """

    # Read information from file.
    items = putly.read_file_text_list(
        delimiter=delimiter,
        path_file=path_file,
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: script_mygene_convert_gene_identifiers_names.py")
        function = "read_source_extract_items"
        print(str("function: " + function + "()"))
        putly.print_terminal_partition(level=4)
        count_items = len(items)
        print("count of items in set or list: " + str(count_items))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return items


def query_mygene_information(
    items_source=None,
    type_source=None,
    type_product=None,
    species=None,
    report=None,
):
    """
    Make query request to MyGene.info and return delivery in response.

    While the "querymany" method of MyGene.info does offer the option to
    deliver information in response to a query as a data-frame table, this
    table has a complex structure (list of dictionaries within cells) when a
    single identifier from Entrez matches multiple identifiers from Ensembl.

    arguments:
        items_source (list<str>): identifiers or names of genes for a query
        type_source (str): type of identifiers or names of genes in source
            query
        type_product (str): type of identifiers or names of genes in product
            delivery
        species (str): name of species, either "human", "mouse", or another
            relevant option
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): delivery of information in response to query

    """

    mygene_class = mygene.MyGeneInfo()
    delivery = mygene_class.querymany(
        items_source,
        scopes=type_source,
        fields=type_product,
        species=species,
        as_dataframe=False,
        df_index=False,
        #email="tcameronwaller@gmail.com",
    )
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("module: script_mygene_convert_gene_identifiers_names.py")
        function = "query_mygene_information"
        print(str("function: " + function + "()"))
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return delivery


def parse_delivery_mygene_information(
    delivery=None,
    items_source=None,
    type_source=None,
    type_product=None,
    species=None,
    report=None,
):
    """
    Parse and organize delivery of information in response to query at
    MyGene.info.

    arguments:
        delivery (list<dict>): delivery of information in response to query
        items_source (list<str>): identifiers or names of genes for a query
        type_source (str): type of identifiers or names of genes in source
            query
        type_product (str): type of identifiers or names of genes in product
            delivery
        species (str): name of species, either "human", "mouse", or another
            relevant option
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): identifiers of genes from a set

    """

    # Collect records of information, which will become rows in table.
    records = list()
    for record_delivery in delivery:
        if (
            ("notfound" not in record_delivery.keys())
        ):
            print(record_delivery)
            query_source = record_delivery["query"]
            if (str(type_product).strip() == "ensembl.gene"):
                if ("ensembl" in record_delivery.keys()):
                    if (isinstance(record_delivery["ensembl"], dict)):
                        # Collect information.
                        record = dict()
                        record[type_source] = record_delivery["query"]
                        record[type_product] = (
                            record_delivery["ensembl"]["gene"]
                        )
                        records.append(record)
                    elif (isinstance(record_delivery["ensembl"], list)):
                        for item in record_delivery["ensembl"]:
                            # Collect information.
                            record = dict()
                            record[type_source] = record_delivery["query"]
                            record[type_product] = item["gene"]
                            records.append(record)
                            pass
                        pass
                    pass
            elif (str(type_product).strip() == "entrezgene"):
                if ("entrezgene" in record_delivery.keys()):
                    identifier_product = record_delivery["entrezgene"]
                else:
                    identifier_product = record_delivery["_id"]
                    pass
                # Collect information.
                record = dict()
                record[type_source] = record_delivery["query"]
                record[type_product] = identifier_product
                records.append(record)
            else:
                # Collect information.
                record = dict()
                record[type_source] = record_delivery["query"]
                record[type_product] = record_delivery[type_product]
                records.append(record)
                pass
            pass
        pass
    # Extract identifiers or names.
    items_product = list()
    for record in records:
        items_product.append(record[type_product])
    # Collect unique items.
    items_product_unique = putly.collect_unique_elements(
        elements=items_product,
    )

    # Organize information in a table.
    table = pandas.DataFrame(data=records)

    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["items"] = items_product_unique
    # Report.
    if report:
        putly.print_terminal_partition(level=4)
        print("module: script_mygene_convert_gene_identifiers_names.py")
        function = "parse_delivery_mygene_information"
        print(str("function: " + function + "()"))
        putly.print_terminal_partition(level=5)
        print(table)
        putly.print_terminal_partition(level=5)
        pass
    # Return information.
    return pail


################################################################################
# Procedure


def execute_procedure(
    path_file_source=None,
    path_file_product=None,
    delimiter_source=None,
    delimiter_product=None,
    type_source=None,
    type_product=None,
    species=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_file_source (str): path to source file in text format from which
            to read identifiers or names of genes
        path_file_product (str): path to product file to which to write in text
            format the identifiers or names of genes
        delimiter_source (str): text delimiter between items in source file,
            not white space
        delimiter_product (str): text delimiter between items in product file,
            not white space
        type_source (str): type of identifiers or names of genes in source
            query
        type_product (str): type of identifiers or names of genes in product
            delivery
        species (str): name of species, either "human", "mouse", or another
            relevant option

    raises:

    returns:

    """

    ##########
    # Parameters.
    report = True

    ##########
    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: script_mygene_convert_gene_identifiers_names.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("path_file_source: " + str(path_file_source))
        print("path_file_product: " + str(path_file_product))
        print("delimiter_source: " + str(delimiter_source))
        print("delimiter_product: " + str(delimiter_product))
        print("type_source: " + str(type_source))
        print("type_product: " + str(type_product))
        print("species: " + str(species))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read and extract identifiers or names in source.
    delimiter_source = delimiter_source.replace("tab", "\t") # "\\t"
    delimiter_source = delimiter_source.replace("newline", "\n") # "\\n"
    items_source = read_source_extract_items(
        path_file=path_file_source,
        delimiter=delimiter_source,
        report=report,
    )

    ##########
    # Query records for genes and their features in MyGene.info.
    delivery = query_mygene_information(
        items_source=items_source,
        type_source=type_source,
        type_product=type_product,
        species=species,
        report=report,
    )
    print("!!!!!!!!!!!!!!!! mygene delivery")
    print(delivery)
    # Use temporary file for testing to avoid repetitive queries at
    # MyGene.info.
    if False:
        # Write product information to file.
        name_file = str(os.path.basename(path_file_source) + "_temp")
        path_directory = os.path.dirname(path_file_product)
        putly.write_object_to_file_pickle(
            object=delivery,
            name_file=name_file,
            path_directory=path_directory,
        )
    if False:
        # Read source information from file.
        delivery = putly.read_object_from_file_pickle(
            path_file=str(path_file_source + "_temp.pickle"),
        )

    ##########
    # Parse delivery from response to query at MyGene.info and organize
    # information in a table and in a list of identifiers or names.
    pail = parse_delivery_mygene_information(
        delivery=delivery,
        items_source=items_source,
        type_source=type_source,
        type_product=type_product,
        species=species,
        report=report,
    )

    ##########
    # Write product information to file.
    #name_file = os.path.basename(path_file_product).split()
    #path_directory = os.path.dirname(path_file_product)
    # Resolve misinterpretation of delimiter strings.
    delimiter_product = delimiter_product.replace("tab", "\t")
    delimiter_product = delimiter_product.replace("newline", "\n")
    putly.write_list_to_file_text(
        elements=pail["items"],
        delimiter=delimiter_product,
        path_file=path_file_product,
    )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_source = sys.argv[1]
    path_file_product = sys.argv[2]
    delimiter_source = sys.argv[3]
    delimiter_product = sys.argv[4]
    type_source = sys.argv[5]
    type_product = sys.argv[6]
    species = sys.argv[7]

    # Call function for procedure.
    execute_procedure(
        path_file_source=path_file_source,
        path_file_product=path_file_product,
        delimiter_source=delimiter_source,
        delimiter_product=delimiter_product,
        type_source=type_source,
        type_product=type_product,
        species=species,
    )

    pass



#
