#!/usr/bin/Rscript

################################################################################
# Organize arguments.

arguments = commandArgs(trailingOnly=TRUE)
path_file_source = arguments[1]
path_file_product = arguments[2]

#table_source = read.table(path_file_source, header=TRUE)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install("MungeSumstats")
library(MungeSumstats)

# Reference files for human genome assembly GRCh38.
#BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh38")
#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
#BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

# Reference files for human genome assembly GRCh37.
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

path_handle <- MungeSumstats::format_sumstats(
    path=path_file_source,
    ref_genome="GRCh37",
    dbSNP=155,
    impute_beta=FALSE,
    impute_se=FALSE,
    sort_coordinates=FALSE,
    rmv_chrPrefix=TRUE,
    rmv_chr=NULL,
    save_path=path_file_product,
    write_vcf=FALSE,
    nThread=4
)

# Print the variable for the save path.
write("hello world!", stdout())
write(path_handle, stdout())

# Clear garbage from memory.
#gc()

# Clear objects in R memory.
#rm(list = ls())

# Restart R session.
# This can be a strong strategy to clear the R memory.
#.rs.restartR()
