#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 13 July 2024
# Date, last execution or modification: 13 July 2024
# Review: TCW; 13 July 2024
###############################################################################
# Note




###############################################################################
# Organize arguments.

paths_file_source=${1} # full paths to files for source genomic or transcriptomic sequence information in BAM format
path_directory_source=${1} # full path to parent directory of files for source genomic or transcriptomic sequence information in BAM format

###############################################################################
# Organize paths.


###############################################################################
# Organize parameters.

###############################################################################
# Execute procedure.

echo "...Test variable expansion in parameter to script..."
echo "Here are the expanded file paths from parameter 'paths_file_source'."
echo $paths_file_source
echo "Here are the contents of the directory 'path_directory_source'."
echo "File suffix is '.bam'."
echo $path_directory_source/*.bam

###############################################################################
# End.
