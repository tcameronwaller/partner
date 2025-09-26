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

# Notice that the default operation of the sort method in SamTools is to sort
# aligned sequence reads by position coordinates and not by name. For paired
# end sequence reads in sort sequence by position coordinates rather than by
# name, it is necessary to specify the "--order=pos" option of HTSeq's count
# method. For paired end sequence reads in sort sequence by name, it is
# necessary to specify the "--order=name" option of HTSeq's count method. It
# tends to be more efficient to sort reads by name so that paired end reads for
# each fragment are more or less adjacent to each other in the BAM file.
# Otherwise, HTSeq's count method can have unnecessary computational burden to
# find paired end reads of either direction. Sort the records using the "-n"
# parameter in SamTools' "sort" method.

# For paired-end RNA sequence technologies, it is most accurate to quantify
# fragments, each of which are represented by two reads in opposite directions.
# Before quantification in HTSeq, use SamTools to filter to proper pairs of
# paired end sequence reads for each fragment. These fragments will offer
# greater confidence and accuracy in quantification.

# Notice that it is necessary for the GTF file annotating genomic features to
# use the same format of chromosome identifiers as were in the reference
# genome used in prior alignment.

# The HTSeq documentation specifies that it's necessary to use the
# "end_included" argument of the "HTSeq.GFF_Reader()" to specify whether or not
# the GTF file specifies position ranges that include the end position. It can
# be ambiguous and difficult to determine whether the GTF annotation file even
# uses a consistent criterion.

# Note: TCW; 27 August 2025
# The "stranded" argument can make a difference in the quantification
# algorithm. It is advisable to make sure that this argument is consistent with
# the sequencing technology and method. It is also advisable to check the
# results of HTSeq quantification to see how many reads it is missing from the
# quantification, I think. Notice that when I quantified the RNAseq data for
# muscle and adipose in 2024, I used the argument setting "--stranded 'no'".
# I do not know for sure whether this was the most accurate method for the
# sequencing technology and method that produced our sequencing data; however,
# I do think that "--stranded 'no'" might be the most versatile and forgiving
# option, especially in our situation in which I do not know all details about
# the origin of the data.

# Decompress file format of the genome annotation.
#gzip -dcvf $path_file_compressed > $path_file_decompressed

# Documentation for htseq-count
# https://htseq.readthedocs.io/en/latest/htseqcount.html#htseqcount
# https://cloud.genepattern.org/gp/module/doc/urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00412:3

# Tools to merge together multiple reports of results from HTSeq Count.
# https://www.genepattern.org/modules/docs/MergeHTSeqCounts/1#gsc.tab=0

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

# Initialize file.
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
