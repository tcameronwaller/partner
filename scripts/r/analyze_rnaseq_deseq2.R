# Script for execution in R.

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 25 July 2024
# Date, last execution or modification: 28 July 2024
# Review: TCW; 28 July 2024
###############################################################################
# Note

# https://www.biostars.org/p/368158/
# It might be a method improvement to use HTSeq or featureCounts to allocate
# and quantify reads to transcripts (without "fractional" allocation) and then
# to use "tximport" (R Bioconductor package) to combine transcripts to gene level.
# https://bioconductor.org/packages/release/bioc/html/tximport.html

# DESeq2
# documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#se
# Paired samples from same subject in different experimental groups
#    documentation: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#can-i-use-deseq2-to-analyze-paired-samples
#    documentation: https://support.bioconductor.org/p/9135209/
#    - attribution of variance to pairs of samples from the same experimental subject
#    - Apparently DESeq2 does employ a mixed effects model, according to
#      Michael Love (2021).



###############################################################################
# Organize arguments.

arguments = commandArgs(trailingOnly=TRUE)
if (length(arguments)==0) {
    # There are not any arguments.
    stop("This script requires arguments.n", call.=FALSE)
} else if (length(arguments)==5) {
  # There are a correct count of arguments.
  print("Thank you for providing the correct count of arguments. Good job!")
  path_file_source_table_sample <- arguments[1]
  path_file_source_table_signal <- arguments[2]
  path_file_product_table <- arguments[3]
  threads <- arguments[4]
  report <- arguments[5]
} else {
  # There are an incorrect count of arguments.
  stop("There seem to be an incorrect count of arguments.n", call.=FALSE)
}



###############################################################################
# Organize libraries and packages.

#library("optparse") # alternative for management of arguments
library("BiocParallel")
library("DESeq2")



###############################################################################
# Organize parameters and settings.

#names(options()) # print names from the options menu
options(width=200) # control print wrapping behavior per width of the console
register(MulticoreParam(strtoi(threads)))
# Report.
cat("\n----------\n----------\n----------\n\n")
print(paste("threads: ", strtoi(threads)))



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
inclusion <- all(rownames(table_sample) %in% colnames(table_signal))
equality <- all(rownames(table_sample) == colnames(table_signal))
print(paste("sample mutual inclusion: ", inclusion))
print(paste("sample mutual equality: ", equality))

# DESeq2 requires integer counts for the signals.
# Round float values to integer representations.
#table_signal_deseq <- as.matrix(round(table_signal, digits = 0))
table_signal_deseq <- round(table_signal, digits = 0)

# Initialize data set in DESeq2.
data_deseq <- DESeqDataSetFromMatrix(
    countData = table_signal_deseq,
    colData = table_sample,
    design = ~ subject + condition
)
# Define levels of factors explicitly.
#data_deseq$condition <- rlevel(data_deseq$condition, ref = "control")
data_deseq$condition <- factor(
    data_deseq$condition,
    levels = c("control", "intervention_1"),
    exclude = NA
)
data_deseq$condition <- droplevels(data_deseq$condition)
data_deseq$subject <- factor(data_deseq$subject)
# Report.
cat("\n----------\n----------\n----------\n\n")
print("Data set in DESeq2.")
print(data_deseq)

# Merge additional attributes of samples.
# See documentation for DESeq2.



###############################################################################
# Perform differential expression analysis.
# DESeq2 supports analyses on pairs of samples in different experimental
# groups.
# If registration of multiple parallel processing cores occurred previously,
# then the argument "BPPARAM=MulticoreParam(threads)" is unnecessary.

data_set <- DESeq(
    data_deseq,
    parallel=TRUE
)
table_result <- results(
    data_deseq,
    contrast=c("condition", "intervention_1", "control"),
    alpha=0.05,
    parallel=TRUE
)
table_result_sort <- table_result[order(table_result$pvalue),]
table_result_sort
summary(table_result_sort)
count_significant <- sum(
    table_result_sort$padj<0.05,
    na.rm=TRUE
)
print(paste("count of significant differences: ", count_significant))
table_result_sort_significant <- subset(
    table_result_sort,
    padj<0.1
)



###############################################################################
# Write product information to file.

write.table(
    as.data.frame(table_result_sort_significant),
    file=path_file_product_table,
    sep="\t",
    eol="\n",
    na="NA",
    dec=".",
    row.names=TRUE,
    col.names=TRUE
)



###############################################################################
# End
