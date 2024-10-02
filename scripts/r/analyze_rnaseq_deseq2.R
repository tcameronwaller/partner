# Script for execution in R.

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 25 July 2024
# Date, last execution or modification: 1 October 2024
# Review: TCW; 1 October 2024
###############################################################################
# Note

# General note about using R with files in Excel format.
# readxl()
# writexl()
# These functions make it possible to access separate sheets within the file.


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
# Organize libraries and packages.

#library("optparse") # alternative for management of arguments
library("BiocParallel")
library("DESeq2")



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
#  8. supplement_1: name of column within table for factor variable with
#      categorical values corresponding to groups of secondary interest as
#      covariates
#  9. levels_supplement_1: categorical values of supplement factor variable
# 10. supplement_2: name of column within table for factor variable with
#      categorical values corresponding to groups of secondary interest as
#      covariates
# 11. levels_supplement_2: categorical values of supplement factor variable
# 12. supplement_3: name of column within table for factor variable with
#      categorical values corresponding to groups of secondary interest as
#      covariates
# 13. levels_supplement_3: categorical values of supplement factor variable
# 14. subject: name of column within table for factor variable with categorical
#      values corresponding to pairs of samples or observations between the
#      experimental conditions of primary interest
# 15. results_contrast: designation of parameter for 'contrast' argument to
#      DESeq2 'results' function, for which value 'condition' designates to
#      specify the condition variable and its levels
# 16. results_name: parameter for 'name' argument to DESeq2 'results' function
# 17. threshold_significance: threshold alpha on p-value for determination of
#      significance
# 18. threads: count of concurrent or parallel process threads on node cores
# 19. report: whether to print report information to terminal

# Parse arguments.
arguments = commandArgs(trailingOnly=TRUE)
print(paste("count of arguments: ", length(arguments)))
if (length(arguments)==0) {
    # There are not any arguments.
    stop("This script requires 19 arguments.", call.=FALSE)
} else if (length(arguments)==19) {
  # There are a correct count of arguments.
  print("correct count of arguments: 19")
  path_file_source_table_sample <- arguments[1]
  path_file_source_table_gene <- arguments[2]
  path_file_source_table_signal <- arguments[3]
  path_file_product_table <- arguments[4]
  formula_text <- paste(as.vector(unlist(strsplit(arguments[5], ","))), collapse=" + ")
  condition <- arguments[6]
  levels_condition <- as.vector(unlist(strsplit(arguments[7], ",")))
  supplement_1 <- arguments[8]
  levels_supplement_1 <- as.vector(unlist(strsplit(arguments[9], ",")))
  supplement_2 <- arguments[10]
  levels_supplement_2 <- as.vector(unlist(strsplit(arguments[11], ",")))
  supplement_3 <- arguments[12]
  levels_supplement_3 <- as.vector(unlist(strsplit(arguments[13], ",")))
  subject <- arguments[14]
  results_contrast <- arguments[15]
  results_name <- arguments[16]
  threshold_significance <- as.double(arguments[17])
  threads <- arguments[18]
  report <- arguments[19]
} else {
  # There are an incorrect count of arguments.
  print("There seem to be an incorrect count of arguments.")
  stop("This script requires 17 arguments.", call.=FALSE)
}

# Report.
cat("\n--------------------------------------------------\n")
cat("--------------------------------------------------\n")
cat("--------------------------------------------------\n\n")
cat("\n--------------------------------------------------\n")
cat("--------------------------------------------------\n")
cat("--------------------------------------------------\n\n")
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
cat("----------\n")
print(paste("5. formula_text: ", formula_text))
cat("----------\n")
print(paste("6. condition: ", condition))
print(paste("7. levels_condition: ", paste(levels_condition, collapse=", ")))
print(paste("... effect is first level: ", levels_condition[1]))
print(paste(
    "... reference is last level: ",
    levels_condition[length(levels_condition)]
))
print("... fold change = (effect) / (reference)")
print(paste(
    "... see documentation and notes below about factor levels ",
    "and contrasts in DESeq2"
))
cat("----------\n")
print(paste("8. supplement_1: ", supplement_1))
print(paste(
    "9. levels_supplement_1: ", paste(levels_supplement_1, collapse=", ")
))
print(paste("10. supplement_2: ", supplement_2))
print(paste(
    "11. levels_supplement_2: ", paste(levels_supplement_2, collapse=", ")
))
print(paste("12. supplement_3: ", supplement_3))
print(paste(
    "13. levels_supplement_3: ", paste(levels_supplement_3, collapse=", ")
))
print(paste("14. subject: ", subject))
print(paste("15. results_contrast: ", results_contrast))
print(paste("16. results_name: ", results_name))


print(paste("17. threshold_significance: ", threshold_significance))
print(paste("18. threads: ", threads))
print(paste("19. report: ", report))
cat("\n----------\n----------\n----------\n\n")
print("Floating point precision on current computer system.")
print(paste("float minimum: ", .Machine$double.xmin))
print(paste("float maximum: ", .Machine$double.xmax))
cat("\n----------\n----------\n----------\n\n")

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
    row.names = "identifier_signal",
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
print(paste("categorical factor for covariate supplement_1: ", supplement_1))
print("values of covariate supplement_1:")
print(levels_supplement_1)
cat("----------\n")
print(paste("categorical factor for covariate supplement_2: ", supplement_2))
print("values of covariate supplement_2:")
print(levels_supplement_2)
cat("----------\n")
print(paste("categorical factor for covariate supplement_3: ", supplement_3))
print("values of covariate supplement_3:")
print(levels_supplement_3)
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
#nchar(condition) > 0 &
#length(condition) > 1 &

if (
    condition != "none" &
    length(levels_condition) > 1 &
    levels_condition[1] != "none"
) {
    table_sample[[condition]] <- factor(table_sample[[condition]])
}
if (
    supplement_1 != "none" &
    length(levels_supplement_1) > 1 &
    levels_supplement_1[1] != "none"
) {
    table_sample[[supplement_1]] <- factor(table_sample[[supplement_1]])
}
if (
    supplement_2 != "none" &
    length(levels_supplement_2) > 1 &
    levels_supplement_2[1] != "none"
) {
    table_sample[[supplement_2]] <- factor(table_sample[[supplement_2]])
}
if (
    supplement_3 != "none" &
    length(levels_supplement_3) > 1 &
    levels_supplement_3[1] != "none"
) {
    table_sample[[supplement_3]] <- factor(table_sample[[supplement_3]])
}
if (
    subject != "none"
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
#nchar(condition) > 0 &
#length(condition) > 1 &

if (
    condition != "none" &
    length(levels_condition) > 1 &
    levels_condition[1] != "none"
) {
    print("... setting explicit levels of 'condition' ...")
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
    supplement_1 != "none" &
    length(levels_supplement_1) > 1 &
    levels_supplement_1[1] != "none"
) {
    print("... setting explicit levels of 'supplement_1' ...")
    data_deseq[[supplement_1]] <- factor(
        data_deseq[[supplement_1]],
        levels = levels_supplement_1,
        exclude = NA
    )
    data_deseq[[supplement_1]] <- droplevels(data_deseq[[supplement_1]])
}
if (
    supplement_2 != "none" &
    length(levels_supplement_2) > 1 &
    levels_supplement_2[1] != "none"
) {
    print("... setting explicit levels of 'supplement_2' ...")
    data_deseq[[supplement_2]] <- factor(
        data_deseq[[supplement_2]],
        levels = levels_supplement_2,
        exclude = NA
    )
    data_deseq[[supplement_2]] <- droplevels(data_deseq[[supplement_2]])
}
if (
    supplement_3 != "none" &
    length(levels_supplement_3) > 1 &
    levels_supplement_3[1] != "none"
) {
    print("... setting explicit levels of 'supplement_3' ...")
    data_deseq[[supplement_3]] <- factor(
        data_deseq[[supplement_3]],
        levels = levels_supplement_3,
        exclude = NA
    )
    data_deseq[[supplement_3]] <- droplevels(data_deseq[[supplement_3]])
}
if (
    subject != "none"
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
# Contrast factors and levels for fold changes:
# Refer to the documentation of DESeq2 (vignettes at top of this file) for
# explanation of contrasts and comparisons between specific levels of the
# main factor variable. By default, DESeq2 reports the fold changes
# corresponding to the last variable in the formulaic design of the analysis.
# If this last variable is a categorical factor, the fold changes will
# correspond to the ratio of the first discrete level of that factor, which is
# the effect (numerator), to the last discrete level of that factor, which is
# the reference (denominator).
# When using interaction terms in the formulaic design of the analysis, there
# are special considerations.
# If registration of multiple parallel processing cores occurred previously,
# then the argument "BPPARAM=MulticoreParam(threads)" is unnecessary.

# Report.
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n\n")
print("Perform differential expression analysis.")
cat("\n----------\n----------\n----------\n\n")
print("fold changes and the formulaic design of analysis:")
print("effect: first level of the last factor in formulaic design")
print("reference: last level of the last factor in formulaic design")
print("(fold change) = (effect) / (reference)")
print("... which is synonymous to...")
print("(fold change) = (effect) vs (reference)")
cat("\n----------\n----------\n----------\n\n")

# Execute regression in DESeq2.
data_deseq <- DESeq(
    data_deseq,
    parallel=TRUE
)
cat("\n----------\n----------\n----------\n\n")
print("... names of individual effect coefficients from the analysis ...")
resultsNames(data_deseq)
cat("\n----------\n----------\n----------\n\n")

# Determine whether to extract results with the "contrast" argument or with the
# "name" argument.
if (
    results_contrast != "none" &
    results_contrast == "condition"
) {
    print("... using 'contrast' argument to DESeq2 'results' function ...")
    table_result <- results(
        data_deseq,
        contrast=c(condition, levels_condition),
        alpha=threshold_significance
    )
} else if (
    results_name != "none"
) {
    print("... using 'name' argument to DESeq2 'results' function ...")
    table_result <- results(
        data_deseq,
        name=results_name,
        alpha=threshold_significance
    )
}
#else {
#  table_result <- results(
#      data_deseq,
#      alpha=threshold_significance
#  )
#}

# Calculate negative base-ten logarithm of p-value.
table_result$neglog10pvalue <- (log(table_result$pvalue, base=10) * -1)
# Calculate rank metric.
table_result$rank_metric <- (
    table_result$log2FoldChange * table_result$neglog10pvalue
)
# Filter rows in table.
#table_result_filter <- subset(
#    table_result,
#        (pvalue < 0.1) &
#        (log2FoldChange < -0.3 | log2FoldChange > 0.3)
#)
# Sort rows in table.
table_result_sort <- table_result[order(table_result$pvalue),]
# Prepare information for summary.
count_significant <- sum(
    table_result_sort$padj<threshold_significance,
    na.rm=TRUE
)
table_result_sort_significant <- subset(
    table_result_sort,
    padj<threshold_significance
)
# Report.
cat("\n----------\n----------\n----------\n\n")
print("Results of differential expression analysis in DESeq2.")
cat("----------\n")
print(table_result)
cat("----------\n")
summary(table_result_sort)
cat("----------\n")
cat("----------\n")
cat("----------\n")
print(paste(
    "value of alpha for p-value significance: ", threshold_significance
))
cat("----------\n")
print(paste("count of significant differences: ", count_significant))
cat("----------\n")
print(table_result_sort_significant)
cat("----------\n")
summary(table_result_sort_significant)
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
    as.data.frame(table_result), # whichever version here determines output
    table_gene,
    by="row.names"
)
colnames(
    table_merge
)[which(names(table_merge) == "Row.names")] <- "identifier_gene"
row.names(table_merge) <- NULL

# Report.
cat("\n----------\n----------\n----------\n\n")
print("Table of results after filter and merge with information about genes.")
cat("----------\n")
print(table_merge[1:10, 1:10])
print(paste("columns: ", ncol(table_merge)))
print(paste("rows: ", nrow(table_merge)))

table_merge_sort <- table_merge[order(table_merge$pvalue),]

# Report.
cat("\n----------\n----------\n----------\n\n")
print("Table of results after filter and sort with information about genes.")
cat("----------\n")
print(table_merge_sort[1:10, 1:10])
print(paste("columns: ", ncol(table_merge_sort)))
print(paste("rows: ", nrow(table_merge_sort)))



###############################################################################
# Write product information to file.

# Report.
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n")
cat("\n--------------------------------------------------\n\n")
print("Write product information to file.")
cat("\n----------\n----------\n----------\n\n")

write.table(
    table_merge_sort,
    file=path_file_product_table,
    sep="\t",
    eol="\n",
    na="NA",
    dec=".",
    row.names=FALSE,
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
