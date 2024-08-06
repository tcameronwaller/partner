# Script for execution in R.

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 25 July 2024
# Date, last execution or modification: 6 August 2024
# Review: TCW; 6 August 2024
###############################################################################
# Note

# https://www.biostars.org/p/368158/
# It might be a method improvement to use HTSeq or featureCounts to allocate
# and quantify reads to transcripts (without "fractional" allocation) and then
# to use "tximport" (R Bioconductor package) to combine transcripts to gene level.
# https://bioconductor.org/packages/release/bioc/html/tximport.html

# DESeq2
# documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
#    - how to specify contrast conditions for differential expression analysis
# documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs

# define formulaic design of analysis with interactions between categorical factors
# - documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
# - documentation: https://www.biostars.org/p/472591/
#    - how to specify and interpret interaction terms in differential expression analysis

# define formulaic design of analysis for pairs of samples between experimental factor groups
# - documentation: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#can-i-use-deseq2-to-analyze-paired-samples
# - documentation: https://support.bioconductor.org/p/9135209/
# - attribution of variance to pairs of samples from the same experimental subject
# - Apparently DESeq2 does employ a mixed effects model, according to
#   Michael Love (2021).
# adjust for covariate factor with nested pairs of samples for each subject or individual
# - documentation: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#group-specific-condition-effects-individuals-nested-within-groups
#    - Careful organization of categorical factor variables enables this analysis design.
# - documentation: https://support.bioconductor.org/p/100828/
#    - When subject or individual for pairs of samples is a factor in the
#      model, it is not strictly necessary to adjust further for the
#      categorical covariate groups that nest those subjects or individuals.
#      Essentially, the subject or individual already adjusts for these and
#      other unknown covariates to focus instead on the pairwise changes.



###############################################################################
# Organize arguments.

# Explain arguments.
#  1. path_file_source_table_sample: full path to file for source table
#  2. path_file_source_table_gene: full path to file for source table
#  3. path_file_source_table_signal: full path to file for source table
#  4. path_file_product_table: full path to file for product table of results
#      from differential expression analysis in DESeq2
#  5. formula_text: formulaic design of analysis in DESeq2
#  6. condition: name of column within table for factor variable with
#      categorical values corresponding to experimental conditions of primary
#      interest
#  7. levels_condition: categorical values of condition factor variable
#  8. supplement: name of column within table for factor variable with
#      categorical values corresponding to groups of secondary interest as
#      covariates
#  9. levels_supplement: categorical values of supplement factor variable
# 10. subject: name of column within table for factor variable with categorical
#      values corresponding to pairs of samples or observations between the
#      experimental conditions of primary interest
# 11. threads: count of concurrent or parallel process threads on node cores
# 12. report: whether to print report information to terminal

# Parse arguments.
arguments = commandArgs(trailingOnly=TRUE)
print(paste("count of arguments: ", length(arguments)))
if (length(arguments)==0) {
    # There are not any arguments.
    stop("This script requires 12 arguments.n", call.=FALSE)
} else if (length(arguments)==12) {
  # There are a correct count of arguments.
  print("correct count of arguments: 12")
  path_file_source_table_sample <- arguments[1]
  path_file_source_table_gene <- arguments[2]
  path_file_source_table_signal <- arguments[3]
  path_file_product_table <- arguments[4]
  formula_text <- paste(as.vector(unlist(strsplit(arguments[5], ","))), collapse=" + ")
  condition <- arguments[6]
  levels_condition <- as.vector(unlist(strsplit(arguments[7], ",")))
  supplement <- arguments[8]
  levels_supplement <- as.vector(unlist(strsplit(arguments[9], ",")))
  subject <- arguments[10]
  threads <- arguments[11]
  report <- arguments[12]
} else {
  # There are an incorrect count of arguments.
  print("There seem to be an incorrect count of arguments.")
  stop("This script requires 12 arguments.", call.=FALSE)
}

# Report.
cat("\n--------------------------------------------------\n")
cat("--------------------------------------------------\n")
cat("--------------------------------------------------\n\n")
print("Arguments.")
cat("\n----------\n----------\n----------\n\n")
print(paste(
    "1. path_file_source_table_sample: ", path_file_source_table_sample
))
print(paste(
    "2. path_file_source_table_gene: ", path_file_source_table_gene
))
print(paste(
    "3. path_file_source_table_signal: ", path_file_source_table_signal
))
print(paste("4. path_file_product_table: ", path_file_product_table))
print(paste("5. formula_text: ", formula_text))
print(paste("6. condition: ", condition))
print(paste("7. levels_condition: ", paste(levels_condition, collapse=", ")))
print(paste("8. supplement: ", supplement))
print(paste(
    "9. levels_supplement: ", paste(levels_supplement, collapse=", ")
))
print(paste("10. subject: ", subject))
print(paste("11. threads: ", threads))
print(paste("12. report: ", report))
cat("\n----------\n----------\n----------\n\n")



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

# image <- "
#     \/
#    _--_
#   /    \
#  / ^  ^ \
# |  *  *  |
#  \      /
#  /      \
#  |      |
#  |      |
#   \ .. /\
#       |_/
# "
print(image)

###############################################################################
# Read and organize source information from file.

# Report.
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n\n")
print("Read source and organize source information from file.")
cat("\n----------\n----------\n----------\n\n")

# Table of samples.
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
cat("----------\n")
print(paste("path: ", path_file_source_table_sample))
cat("----------\n")
print(table_sample[1:10, ])
print(paste("columns: ", ncol(table_sample)))
print(paste("rows: ", nrow(table_sample)))
cat("----------\n")

# Table of genes.
table_gene <- read.table(
    path_file_source_table_gene,
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
# Report.
cat("\n----------\n----------\n----------\n\n")
print("Table of information and attributes about genes.")
cat("----------\n")
print(paste("path: ", path_file_source_table_gene))
cat("----------\n")
print(table_gene[1:10, ])
print(paste("columns: ", ncol(table_gene)))
print(paste("rows: ", nrow(table_gene)))

# Table of signals.
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
# Report.
cat("\n----------\n----------\n----------\n\n")
print("Table of signals for genes across samples.")
cat("----------\n")
print(paste("path: ", path_file_source_table_signal))
cat("----------\n")
print(table_signal[1:10, 1:10])
print(paste("columns: ", ncol(table_signal)))
print(paste("rows: ", nrow(table_signal)))



###############################################################################
# Create DESeq2 data set.

# Report.
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n\n")
print("Create DESeq2 data set.")
cat("\n----------\n----------\n----------\n\n")

# Report.
cat("\n----------\n----------\n----------\n\n")
print("Formulaic design and categorical factors for analysis in DESeq2.")
cat("----------\n")
print(paste("formulaic design: ~ ", formula_text))
cat("----------\n")
print(paste("factor for main experimental condition: ", condition))
print("values of experimental condition:")
print(levels_condition)
cat("----------\n")
print(paste("factor for covariate supplement: ", supplement))
print("values of covariate supplement:")
print(levels_supplement)
cat("----------\n")

##########
# Check coherence.
cat("\n----------\n----------\n----------\n\n")
print("Confirm that both tables have identical sequences of samples.")
cat("----------\n")
samples_sample <-rownames(table_sample)
samples_signal <- colnames(table_signal)
print("sample identifiers from table of sample attributes:")
print(samples_sample)
print("sample identifiers from table of signals:")
print(samples_signal)
cat("----------\n")
inclusion <- all(rownames(table_sample) %in% colnames(table_signal))
equality <- all(rownames(table_sample) == colnames(table_signal))
print(paste("sample lists mutual inclusion: ", inclusion))
print(paste("sample lists equality: ", equality))
cat("\n----------\n----------\n----------\n\n")

##########
# Simplify signals to counts.
# DESeq2 requires integer counts for the signals.
# Round float values to integer representations.
#table_signal_deseq <- as.matrix(round(table_signal, digits = 0))
table_signal_deseq <- round(table_signal, digits = 0)

##########
# Predefine categorical factor variables.
if (
    nchar(condition) > 0 &
    length(condition) > 1 &
    condition != "None"
) {
    table_sample[[condition]] <- factor(table_sample[[condition]])
}
if (
    nchar(supplement) > 0 &
    length(supplement) > 1 &
    supplement != "None"
) {
    table_sample[[supplement]] <- factor(table_sample[[supplement]])
}
if (
    nchar(subject) > 0 &
    length(subject) > 1 &
    subject != "None"
) {
    table_sample[[subject]] <- factor(table_sample[[subject]])
}

##########
# Initialize data set in DESeq2.
# Statistics (log2 fold change, p-values, etc) of priority in results will
# correspond to the last factor in the formula. The statistics for other
# factors (covariates, etc) are also accessible from the results
# ("results(data, name=...)").
# Use the function "formula()" to define a variable formulaic design.
data_deseq <- DESeqDataSetFromMatrix(
    countData = as.matrix(table_signal_deseq),
    colData = table_sample,
    design = formula(paste("~", formula_text))
)

##########
# Define levels of categorical factors explicitly.
# The "$" subset operator (example "data_deseq$condition") only works with
# literal arguments that do not need any prior evaluation.
if (
    nchar(condition) > 0 &
    length(condition) > 1 &
    condition != "None"
) {
    #data_deseq[[condition]] <- relevel(
    #    data_deseq[[condition]], ref = "control"
    #)
    data_deseq[[condition]] <- factor(
        data_deseq[[condition]],
        levels = levels_condition,
        exclude = NA
    )
    data_deseq[[condition]] <- droplevels(data_deseq[[condition]])
}
if (
    nchar(supplement) > 0 &
    length(supplement) > 1 &
    supplement != "None"
) {
    data_deseq[[supplement]] <- factor(
        data_deseq[[supplement]],
        levels = levels_supplement,
        exclude = NA
    )
    data_deseq[[supplement]] <- droplevels(data_deseq[[supplement]])
}
if (
    nchar(subject) > 0 &
    length(subject) > 1 &
    subject != "None"
) {
    data_deseq[[subject]] <- factor(data_deseq[[subject]])
}

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

# Report.
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n\n")
print("Perform differential expression analysis.")
cat("\n----------\n----------\n----------\n\n")

data_deseq <- DESeq(
    data_deseq,
    parallel=TRUE
)
table_result <- results(
    data_deseq,
    contrast=c(condition, levels_condition),
    alpha=0.05
)
table_result_sort <- table_result[order(table_result$pvalue),]
table_result_sort_significant <- subset(
    table_result_sort,
    padj<0.05
)
count_significant <- sum(
    table_result_sort$padj<0.05,
    na.rm=TRUE
)
# Report.
cat("\n----------\n----------\n----------\n\n")
print("Results of differential expression analysis in DESeq2.")
cat("----------\n")
print(table_result_sort)
cat("----------\n")
summary(table_result_sort)
cat("----------\n")
print(paste("count of significant differences: ", count_significant))
cat("\n----------\n----------\n----------\n\n")



###############################################################################
# Include information about genes in table of results.

# Report.
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n\n")
print("Prepare and organize report of results.")
cat("\n----------\n----------\n----------\n\n")

# Merge together tables.
table_merge <- merge(
    as.data.frame(table_result_sort_significant),
    table_gene,
    by="row.names"
)
colnames(
    table_merge
)[which(names(table_merge) == "Row.names")] <- "identifier_gene"

# Report.
cat("\n----------\n----------\n----------\n\n")
print("Table of results after merge with information about genes.")
cat("----------\n")
print(table_merge[1:10, 1:10])
print(paste("columns: ", ncol(table_merge)))
print(paste("rows: ", nrow(table_merge)))



###############################################################################
# Write product information to file.

# Report.
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n\n")
print("Write product information to file.")
cat("\n----------\n----------\n----------\n\n")

write.table(
    table_merge,
    file=path_file_product_table,
    sep="\t",
    eol="\n",
    na="NA",
    dec=".",
    row.names=TRUE,
    col.names=TRUE
)

# Report.
cat("\n----------\n----------\n----------\n\n")
print("Table of results of gene differential expression analysis in DESeq2.")
cat("----------\n")
print(paste("path: ", path_file_product_table))
cat("----------\n")


###############################################################################
# End
