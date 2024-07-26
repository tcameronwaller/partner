# Script for execution in R.

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 25 July 2024
# Date, last execution or modification: 25 July 2024
# Review: TCW; 25 July 2024
###############################################################################
# Note

# https://www.biostars.org/p/368158/
# It might be a method improvement to use HTSeq or featureCounts to allocate
# and quantify reads to transcripts and then to use "tximport" (R Bioconductor
# package) to combine transcripts to gene level.



# documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#se

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
# Organize libraries and packages.

#library("optparse") # alternative for management of arguments
library("DESeq2")

###############################################################################
# Organize parameters and settings.

#names(options()) # print names from the options menu
options(width=200) # control print wrapping behavior per width of the console

###############################################################################
# Read and organize source information from file.

table_sample <- read.table(
    path_file_source_table_sample,
    header = TRUE,
    row.names = "identifier",
    check.names = FALSE,
    sep = "\t",
    quote = "\"'",
    dec = ".",
    na.strings = c(
        "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>"
    ),
    encoding = "UTF-8",
)
# Create row names from integer series.
#row.names(table_sample) <- 1:nrow(table_sample)
# Create row names from existing column in table.
#row.names(table_sample) <- table_sample$identifier
#table_sample$identifier <- NULL
# Transfer row names to a new column in the table.
#table_sample$identifier <- row.names(table_sample)
# Report.
cat("\n----------\n----------\n----------\n\n")
print("Table of information and attributes about samples.")
print(table_sample[1:10, ])
print(paste("columns: ", ncol(table_sample)))
print(paste("rows: ", nrow(table_sample)))

table_signal <- read.table(
    path_file_source_table_signal,
    header = TRUE,
    row.names = "identifier_gene",
    check.names = FALSE,
    sep = "\t",
    quote = "\"'",
    dec = ".",
    na.strings = c(
        "nan", "na", "NAN", "NA", "<nan>", "<na>", "<NAN>", "<NA>"
    ),
    encoding = "UTF-8",
)
cat("\n----------\n----------\n----------\n\n")
print("Table of signals for genes across samples.")
print(table_signal[1:10, 1:10])
print(paste("columns: ", ncol(table_signal)))
print(paste("rows: ", nrow(table_signal)))

###############################################################################
# Create DESeq2 data set.

# Check coherence.
cat("\n----------\n----------\n----------\n\n")
print("Confirm that both tables have identical sequences of samples.")
samples_sample <-rownames(table_sample)
samples_signal <- colnames(table_signal)
print(samples_sample)
print(samples_signal)
all(rownames(table_sample) %in% colnames(table_signal))
all(rownames(table_sample) == colnames(table_signal))

# DESeq2 requires integer counts for the signals.
#table_signal_deseq <- as.matrix(round(table_signal, digits = 0))
table_signal_deseq <- round(table_signal, digits = 0)

data_set <- DESeqDataSetFromMatrix(
    countData = table_signal_deseq,
    colData = table_sample,
    design = ~ condition
)
cat("\n----------\n----------\n----------\n\n")
print("Data set in DESeq2.")
print(data_set)



###############################################################################
# End
