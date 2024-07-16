#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 13 July 2024
# Date, last execution or modification: 14 July 2024
# Review: TCW; 14 July 2024
###############################################################################
# Note


# Documentation for subread featureCounts
# https://rdrr.io/bioc/Rsubread/man/featureCounts.html
# https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html

# Publication
# https://academic.oup.com/bioinformatics/article/30/7/923/232889?login=false


###############################################################################
# Organize arguments.

path_directory_source=${1} # full path to parent directory of files for source genomic or transcriptomic sequence information in BAM format
path_file_product=${2} # full path to file for product quantification of reads allocatable to specific genomic features
path_file_annotation_gtf_gzip=${3} # full path to file for annotation of reference genome with GZip compression
threads=${4} # count of concurrent or parallel process threads on node cores
report=${5} # whether to print reports to terminal
path_environment_htseq=${6} # full path to Python virtual environment

###############################################################################
# Organize paths.

stamp_date=$(date +%Y-%m-%d)
name_base_file_product="$(basename $path_file_product .tsv)"
path_directory_product="$(dirname $path_file_product)"
path_directory_temporary="${path_directory_product}/temporary_${name_base_file_product}_${stamp_date}" # hopefully unique
path_file_temporary_1="${path_directory_temporary}/${name_base_file_product}_genome_annotation_temporary_1.gtf"

# Initialize directory.
mkdir -p $path_directory_product
rm -r $path_directory_temporary
mkdir -p $path_directory_temporary

# Remove any previous version of the product file.
rm $path_file_temporary_1
rm $path_file_product


###############################################################################
# Organize parameters.


###############################################################################
# Execute procedure.

##########
# Deompress file format.
gzip -dcvf $path_file_annotation_gtf_gzip > $path_file_temporary_1

# Quantify reads at genomic features.
if true; then
  # Filter.
  $path_execution_featurecounts \
  view \
  --bam \
  --with-header \
  --require-flags 0x1,0x2 \
  --excl-flags 0x4,0x8,0x200 \
  --threads $threads \
  --output $path_file_temporary_1 \
  $path_file_source
  # Sort coordinates of file in BAM format.
  $path_execution_samtools \
  sort \
  -n \
  -@ $threads \
  -o $path_file_product \
  $path_file_temporary_1
  # Create index for file in BAM format.
  $path_execution_samtools \
  index \
  --bai \
  --threads $threads \
  -o $path_file_product_index \
  $path_file_product
fi



if true; then

  ##########
  # Activate Python Virtual Environment.
  source "${path_environment_htseq}/bin/activate"
  echo "----------"
  echo "...Confirm Python Virtual Environment path..."
  which python3
  sleep 1s
  echo "----------"
  # Regulate concurrent or parallel process threads on node cores.
  # Force Python program (especially SciPy) not to use all available cores on a
  # cluster computation node.
  export MKL_NUM_THREADS=$threads
  export NUMEXPR_NUM_THREADS=$threads
  export OMP_NUM_THREADS=$threads

  #paths_file_source=()
  #while IFS= read -r -d $'\0'; do
  #  paths_file_source+=("$REPLY")
  #done < <(find $path_directory_source -maxdepth 1 -mindepth 1 -type f -name "*.bam" -print0)
  #paths_file_source_epansion="${paths_file_source[*]}"

  ##########
  # Quantify reads allocatable to specific genomic features.
  python3 -m HTSeq.scripts.count \
  --format "bam" \
  --minaqual 10 \
  --idattr "gene_id" \
  --type "gene" \
  --order "name" \
  --stranded "no" \
  --mode "union" \
  --nonunique "fraction" \
  --with-header \
  --additional-attr "gene_id" \
  --additional-attr "gene_name" \
  --additional-attr "exon_number" \
  --additional-attr "gene_type" \
  --add-chromosome-info \
  --delimiter "\t" \
  --counts_output $path_file_product \
  --quiet \
  --nprocesses $threads \
  $path_directory_source/*.bam \
  $path_file_temporary_1

  ##########
  # Deactivate Python Virtual Environment.
  deactivate
  echo "----------"
  echo "confirm deactivation of Python Virtual Environment..."
  which python3
fi

##########
# Remove temporary, intermediate files.
rm -r $path_directory_temporary

###############################################################################
# End.
