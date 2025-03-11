"""
Query SQLite database of Molecular Signatures Database (MSigDB) for a
selection of reference gene sets.

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
# Date, first execution: 7 March 2025
# Date, last execution: 7 March 2025
# Review: TCW; 7 March 2025
################################################################################
# Note

# data format: SQLite
# site: https://www.sqlite.org
# site: https://sqlite.org/cli.html
#   - official command-line terminal interface for SQLite
# site: https://docs.python.org/3/library/sqlite3.html
# site: https://www.sqlalchemy.org/

# database: Molecular Signatures Database (MSigDB)
# site: https://www.gsea-msigdb.org/gsea/msigdb
# site: https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
# site: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Msigdb_browser

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
import sqlite3

# Custom
import partner.utility as putly

#dir()
#importlib.reload()

###############################################################################
# Functionality


def read_organize_source_parameter_sets(
    path_file_table_selection=None,
    report=None,
):
    """
    Reads and organizes source information about parameters for selection of
    sets from MSigDB.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        path_file_table_selection (str): path to source file of table with
            parameters for selection of a collection of gene sets from the
            information information from MSigDB
        report (bool): whether to print reports

    raises:

    returns:
        (dict): collection of source information about parameters for selection
            of sets from MSigDB

    """

    # Specify variable types of columns within table.
    types_columns = dict()
    types_columns["inclusion"] = "string" # "int32"
    types_columns["sequence"] = "string" # "int32"
    types_columns["groups"] = "string"
    types_columns["species"] = "string"
    types_columns["collection"] = "string"
    types_columns["set"] = "string"

    # Read information from file.

    # Table of parameters for parallel instances.
    table = pandas.read_csv(
        path_file_table_selection,
        sep="\t",
        header=0,
        dtype=types_columns,
        na_values=[
            "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>",
        ],
        encoding="utf-8",
    )
    # Organize information.
    table["inclusion"] = pandas.to_numeric(
        table["inclusion"],
        downcast="integer",
        errors="coerce",
    )
    table["sequence"] = pandas.to_numeric(
        table["sequence"],
        downcast="integer",
        errors="coerce",
    )

    # Collect information.
    records = list()
    for index, row in table.iterrows():
        if (int(row["inclusion"]) == 1):
            # Collect information and parameters for current row in table.
            pail = dict()
            pail["sequence"] = row["sequence"]
            pail["groups"] = str(row["groups"]).strip()
            pail["species"] = str(row["species"]).strip()
            pail["collection"] = str(row["collection"]).strip()
            pail["set"] = str(row["set"]).strip()
            # Collect information and parameters for current row in table.
            records.append(pail)
            pass
        pass
    # Organize information.
    #sets_names = list(map(lambda d: d["set"] if "set" in d, records))
    species = [d["species"] for d in records if "species" in d]
    species = putly.collect_unique_elements(
        elements=species
    )
    collections_names = [d["collection"] for d in records if "collection" in d]
    collections_names = putly.collect_unique_elements(
        elements=collections_names
    )
    sets_names = [d["set"] for d in records if "set" in d]
    sets_names = putly.collect_unique_elements(
        elements=sets_names
    )
    # Collect information.
    pail = dict()
    pail["table"] = table
    pail["records"] = records
    pail["species"] = species
    pail["collections_names"] = collections_names
    pail["sets_names"] = sets_names
    # Report.
    if report:
        count_records = len(records)
        putly.print_terminal_partition(level=3)
        print("module: script_sqlite_msigdb_query_organize_gene_sets.py")
        print("function: read_organize_source_parameter_sets()")
        putly.print_terminal_partition(level=5)
        print("collections:")
        print(collections_names)
        print("count: " + str(len(collections_names)))
        putly.print_terminal_partition(level=5)
        print("sets:")
        print(sets_names)
        print("count: " + str(len(sets_names)))
        pass
    # Return information.
    return pail


def report_description_msigdb(
    connection=None,
    cursor=None,
    report=None,
):
    """
    Reports basic description of the MSigDB database.

    arguments:
        connection (object): handle to an open and active connection to a
            database in SQLite format
        cursor (object): handle to an open and active cursor with a connection
            to a database in SQLite format
        report (bool): whether to print reports

    raises:

    returns:


    """

    # Report.
    if report:
        putly.print_terminal_partition(level=3)
        print("module: script_sqlite_msigdb_query_organize_gene_sets.py")
        print("function: report_description_msigdb()")
        # Extract names of tables within database.
        putly.print_terminal_partition(level=5)
        print("Names of tables within database:")
        reply = cursor.execute("SELECT name FROM sqlite_master")
        print(reply.fetchone())
        print(reply.fetchall())
        # Extract version and description of database.
        #reply = cursor.execute("SELECT name FROM MSigDB")
        putly.print_terminal_partition(level=5)
        print("names of columns from table 'MSigDB':")
        reply = cursor.execute("SELECT name FROM PRAGMA_TABLE_INFO('MSigDB')")
        print(reply.fetchall())
        putly.print_terminal_partition(level=5)
        print("names and details of columns from table 'MSigDB':")
        reply = cursor.execute("PRAGMA TABLE_INFO('MSigDB')")
        table_reply = pandas.DataFrame(
            reply.fetchall(),
            columns=[
                "cid",
                "name",
                "type",
                "notnull",
                "dflt_value",
                "pk",
            ]
        )
        print(table_reply)
        putly.print_terminal_partition(level=5)
        #print("everything from table 'MSigDB':")
        #reply = cursor.execute("SELECT * FROM MSigDB")
        #print(reply.fetchall())
    pass




############

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
    path_file_source_msigdb=None,
    path_file_table_selection=None,
    path_file_product=None,
    identifier_type=None,
    delimiter_name=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_file_source_msigdb (str): path to source file of information from
            MSigDB in SQLite format
        path_file_table_selection (str): path to source file of table with
            parameters for selection of a collection of gene sets from the
            information information from MSigDB
        path_file_product (str): path to product file to which to write in text
            format the identifiers or names of genes from a selection of sets
            in MSigDB
        identifier_type (str): type of identifiers or names of genes to extract
            from MSigDB and organize in the product file
        delimiter_name (str): text delimiter between items in product file,
            specified by name to avoid problems passing white space between
            scripts

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
        print("module: script_sqlite_msigdb_query_organize_gene_sets.py")
        print("function: execute_procedure()")
        putly.print_terminal_partition(level=5)
        print("system: local")
        print("path_file_source_msigdb: " + str(path_file_source_msigdb))
        print("path_file_table_selection: " + str(path_file_table_selection))
        print("path_file_product: " + str(path_file_product))
        print("identifier_type: " + str(identifier_type))
        print("delimiter_name: " + str(delimiter_name))
        putly.print_terminal_partition(level=5)
        pass

    ##########
    # Read, parse table of parameters for selection.
    pail_parameters = read_organize_source_parameter_sets(
        path_file_table_selection=path_file_table_selection,
        report=report,
    )
    #pail["species"] = species
    #pail["collections_names"] = collections_names
    #pail["sets_names"] = sets_names

    putly.print_terminal_partition(level=2)
    print("now testing SQLite3 queries...")

    # Connect to database.
    connection = sqlite3.connect(path_file_source_msigdb)
    cursor = connection.cursor()
    # Report basic description of the database.
    report_description_msigdb(
        connection=connection,
        cursor=cursor,
        report=report,
    )
    sets_names = pail_parameters["sets_names"]
    connection.set_trace_callback(print)
    # "WHERE gset.id IN (?) " # does not work
    # "WHERE gset.id IN ({', '.join([?]*len(sets_names))}) "
    #sets_names_expand = ", ".join(['?'] * len(sets_names))
    #reply = cursor.execute(
    #    f"SELECT standard_name 'na', group_concat(symbol, '	') "
    #    f"FROM gene_set gset "
    #    f"INNER JOIN gene_set_gene_symbol gsgs on gset.id = gene_set_id "
    #    f"INNER JOIN gene_symbol gsym on gsym.id = gene_symbol_id "
    #    f"WHERE gset.id IN ({sets_names_expand}) "
    #    f"GROUP BY standard_name ORDER BY standard_name ASC;",
    #    sets_names_expand,
    #) # variable parameters
    reply = cursor.execute(
        "SELECT standard_name 'na', group_concat(symbol, '	') "
        "FROM gene_set gset "
        "INNER JOIN gene_set_gene_symbol gsgs on gset.id = gene_set_id "
        "INNER JOIN gene_symbol gsym on gsym.id = gene_symbol_id "
        "WHERE standard_name IN ("
        "'HALLMARK_ADIPOGENESIS', 'HALLMARK_ANGIOGENESIS'"
        ") "
        "GROUP BY standard_name ORDER BY standard_name ASC;",
    ) # variable parameters
    connection.set_trace_callback(None)
    print(reply.fetchall())

    # TODO: simplify the query and figure out the way to insert the expanded
    # list for the "IN" statement.




    if False:

        ##########
        # Parse delimiter.
        delimiter = delimiter_name
        delimiter = delimiter.replace("newline", "\n") # "\\n"
        delimiter = delimiter.replace("tab", "\t") # "\\t"
        delimiter = delimiter.replace("space", " ") # " "
        delimiter = delimiter.replace("semicolon", ";")
        delimiter = delimiter.replace("colon", ":")
        delimiter = delimiter.replace("comma", ",")

        ##########
        # Query records for genes and their features in MyGene.info.
        delivery = query_mygene_information(
            items_source=items_source,
            type_source=type_source,
            type_product=type_product,
            species=species,
            report=report,
        )

    pass


if (__name__ == "__main__"):
    # Parse arguments from terminal.
    path_file_script = sys.argv[0] # always the first argument
    path_file_source_msigdb = sys.argv[1]
    path_file_table_selection = sys.argv[2]
    path_file_product = sys.argv[3]
    identifier_type = sys.argv[4]
    delimiter_name = sys.argv[5]

    # Call function for procedure.
    execute_procedure(
        path_file_source_msigdb=path_file_source_msigdb,
        path_file_table_selection=path_file_table_selection,
        path_file_product=path_file_product,
        identifier_type=identifier_type,
        delimiter_name=delimiter_name,
    )
    pass



#
