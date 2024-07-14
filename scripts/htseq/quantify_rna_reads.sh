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

# Notice that the default operation of the sort method in SamTools is to sort
# aligned sequence reads by position coordinates and not by name. For paired
# end sequence reads in sort sequence by position coordinates rather than by
# name, it is necessary to specify the "--order=pos" option of HTSeq's count
# method.
# If HTSeq's count method has difficulty (due to the buffer necessary to find
# paired end reads of either direction), then try sorting the reads instead by
# name, using the "-n" parameter option. Then use the "--order=name" parameter
# option of HTSeq's count method.

# Deompress file format of the genome annotation.
#gzip -dcvf $path_file_compressed > $path_file_decompressed


# TODO: TCW; 11 July 2024
# It might be necessary for the GTF annotation file to use UCSC gene identifiers
# since these will match the BAM aligned reads.

# TODO: TCW; 11 July 2024
# The HTSeq documentation specifies that it's necessary to use the "end_included"
# argument of the "HTSeq.GFF_Reader()" to specify whether or not the GTF file
# specifies position ranges that include the end position. If the length of the
# range of an exon is divisible by 3 (codon) then the "end_included" parameter
# must be "False".

# Documentation for htseq-count
# https://htseq.readthedocs.io/en/latest/htseqcount.html#htseqcount



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

if true; then

  ##########
  # Deompress file format.
  gzip -dcvf $path_file_annotation_gtf_gzip > $path_file_temporary_1

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
