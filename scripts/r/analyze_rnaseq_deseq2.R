# Script for execution in R.

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 24 July 2024
# Date, last execution or modification: 24 July 2024
# Review: TCW; 24 July 2024
###############################################################################
# Note

# documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

###############################################################################
# Organize libraries and packages.

#library("optparse") # alternative management of arguments

###############################################################################
# Organize arguments.

arguments = commandArgs(trailingOnly=TRUE)
if (length(arguments)==0) {
  stop("This script requires arguments.n", call.=FALSE)
} else if (length(arguments)==1) {
  # if there are too few arguments, then specify defaults
  arguments[2] = "test_default"
}

###############################################################################
# End
