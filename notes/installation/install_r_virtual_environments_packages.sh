

################################################################################
# Organize paths.

# Directories.
cd ~
path_directory_paths="./Downloads/paths_process_local"
path_directory_tool=$(<"$path_directory_paths/path_directory_tool.txt")
path_directory_repository_partner=$(<"$path_directory_paths/path_directory_repository_partner.txt")
path_directory_scripts="$path_directory_repository_partner/scripts"


###############################################################################
# Installation of R
# site: https://www.r-project.org/
# CRAN mirror host: Iowa State University in Ames, Iowa
# site: https://mirror.las.iastate.edu/CRAN/

##########
# R v4.4.1
# description: https://mirror.las.iastate.edu/CRAN/doc/manuals/r-release/NEWS.html
# version: v4.4.1
# date, release: 14 June 2024
# installation:
# - installation from source
# - system: halyard
# - date, installation: 23 July 2024
# - version, installation: v4.4.1
# Install dependencies.
# https://docs.posit.co/resources/install-r-source.html
sudo sed -i.bak "/^#.*deb-src.*universe$/s/^# //g" /etc/apt/sources.list
sudo apt update
sudo apt build-dep r-base
# Install R itself.
cd "${path_directory_tool}/r"
wget https://mirror.las.iastate.edu/CRAN/src/base/R-4/R-4.4.1.tar.gz
tar -xzvf R-4.4.1.tar.gz # Extract from tar ball with GZip compression.
rm R-4.4.1.tar.gz
mv R-4.4.1 r-4.4.1
cd r-4.4.1
pwd
# Argument "--prefix" needs a full, explicit file path.
# This argument keeps installation within the specific directory path.
./configure \
--prefix="/.../tool/r/r-4.4.1" \
--enable-R-shlib \
--enable-memory-profiling
make
make check # synonym for "make test"?
make install # halyard: 24 July 2024
# Rscript.
"${path_directory_tool}/r/r-4.4.1/bin/Rscript" --help
"${path_directory_tool}/r/r-4.4.1/bin/Rscript" --version # "Rscript (R) version 4.4.1 (2024-06-14)"; TCW; 24 July 2024
# Initialize R console.
# Execution from Bash terminal.
# user@machine:/path$ "${path_tool}/r/r-4.4.1/bin/R"
# Execution from R console.
# > help.start()
# > ?mean
# Print path to library in which R installs packages.
# > getOption("repos") # Print repositories from which R will install packages.
# > .libPaths() # "$path_tool/r/r-4.4.1/lib/R/library"
# > library()
# > library(lib.loc='~/r-packages')
# > lapply(.libPaths(), list.files) # Print packages available in each environmental library.
# > q()
# Remove installation.
#rm -rf $path_r_441



###############################################################################
# Installation of Bioconductor repository and other packages

# Execution from Bash terminal.
"${path_directory_tool}/r/r-4.4.1/bin/Rscript" "$path_directory_scripts/r/install_packages.R"

#Rscript -e 'renv::run("/path/to/myscript.R")'


# End
