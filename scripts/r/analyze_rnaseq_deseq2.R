# Script for execution in R.

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 25 July 2024
# Date, last execution or modification: 25 July 2024
# Review: TCW; 25 July 2024
###############################################################################
# Note

# documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

###############################################################################
# Organize libraries and packages.

#library("optparse") # alternative for management of arguments

###############################################################################
# Organize arguments.

arguments = commandArgs(trailingOnly=TRUE)
if (length(arguments)==0) {
    # There are not any arguments.
    stop("This script requires arguments.n", call.=FALSE)
} else if (length(arguments)==3) {
  # There are a correct count of arguments.
  print("Thank you for providing the correct count of arguments. Good job!")
  path_file_source_table_sample <- arguments[1]
  path_file_source_table_signal <- arguments[2]
  path_file_product_table <- arguments[3]
} else {
  # There are an incorrect count of arguments.
  stop("There seem to be an incorrect count of arguments.n", call.=FALSE)
}

###############################################################################
# Read source information from file.

table_sample <- read.table(
    path_file_source_table_sample,
    header = TRUE,
    sep = "\t",
    quote = "\"'",
    dec = ".",
    na.strings = c(
        "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>"
    ),
    encoding = "UTF-8",
)
print(table_sample)



###############################################################################
# End
