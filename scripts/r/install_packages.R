# Script for execution in R.

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 24 July 2024
# Date, last execution or modification: 24 July 2024
# Review: TCW; 24 July 2024
###############################################################################
# Note

# Use "renv" to manage virtual environment libraries of specific packages and
# their versions.
# renv
#   description: https://posit.co/blog/renv-project-environments-for-r/
#   site: https://rstudio.github.io/renv/
#   installation from CRAN:
#      from R console: install.packages("renv")

#Rscript -e 'renv::run("/path/to/myscript.R")'
# renv::init(
#    project = /path/to/project/environment/directory,
#)

###############################################################################
# Bioconductor installation
# Bioconductor provides a specialized repository of packages for
# bioinformatics.
#   site: https://www.bioconductor.org/install/

##########
# Bioconductor v3.19
# - description: Bioconductor provides a specialized repository of packages for
#   bioinformatics.
# - site: https://www.bioconductor.org/install/
# - site: https://support.bioconductor.org/
# version: 1.44.0
# date, release: ___
# installation:
# - system: halyard
# - date, installation: 24 July 2024
# - version, installation: v3.19
# Execution from R console or R script.
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(version = "3.19")
BiocManager::version() # "3.19"; TCW; 24 July 2024
BiocManager::valid() # "TRUE"; TCW; 24 July 2024



##########
# DESeq2 v1.44.0
# - description: DESeq2 analyzes differential gene expression using models for
#   negative binomial distribution.
# - site: https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html
# installation:
# - system: halyard
# - date, installation: 24 July 2024
# - version, installation: v1.44.0
# Execution from R console or R script.
BiocManager::install("DESeq2")




##########
# edgeR v4.2.1
# - description: edgeR analyzes differential gene expression using models for
#   negative binomial distribution.
# - site: https://bioconductor.org/packages/release/bioc/html/edgeR.html
# installation:
# - system: halyard
# - date, installation: 24 July 2024
# - version, installation: v4.2.1
# Execution from R console or R script.
BiocManager::install("edgeR")



##########
# End
