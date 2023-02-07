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

#path_handle <- MungeSumstats::format_sumstats(
#    path=path_file_source,
#    ref_genome="GRCh37",
#    dbSNP=155,
#    allele_flip_check=1,
#    allele_flip_drop=TRUE,
#    bi_allelic_filter=1,
#    remove_multi_rs_snp=1,
#    impute_beta=0,
#    impute_se=0,
#    sort_coordinates=FALSE,
#    rmv_chrPrefix=TRUE,
#    rmv_chr=NULL,
#    save_path=path_file_product,
#    write_vcf=FALSE,
#    nThread=4
#)

# Operation "rmv_chr" removes chromosomes "X", "Y", and "MT" by default.
# Operation "on_ref_genome" is computational intensive and slow.
# on_ref_genome=FALSE,
# Operation "snp_ids_are_rs_ids" determines SNP reference sequence identifiers
# (rsIDs) from chromosome and base position.

path_handle <- MungeSumstats::format_sumstats(
    path=path_file_source,
    ref_genome="GRCh37",
    sort_coordinates=FALSE,
    rmv_chrPrefix=TRUE,
    rmv_chr=NULL,
    snp_ids_are_rs_ids=FALSE,
    on_ref_genome=TRUE,
    save_path=path_file_product,
    write_vcf=FALSE,
    nThread=4
)

# Print the variable for the save path.
# R "print" and "cat" functions write to "stdout".
# R "message", "warning", and "stop" functions write to "stderr".
#write("hello world!", stdout())
#write(path_handle, stdout())
print("File path to GWAS summary statistics after Munge procedure:")
print(path_handle)

# Clear garbage from memory.
#gc()

# Clear objects in R memory.
#rm(list = ls())

# Restart R session.
# This can be a strong strategy to clear the R memory.
#.rs.restartR()
