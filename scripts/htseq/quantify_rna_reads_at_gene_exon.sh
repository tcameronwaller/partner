#!/bin/bash

#chmod u+x script.sh
#chmod -R 777

###############################################################################
# Author: T. Cameron Waller
# Date, first execution: 11 July 2024
# Date, last execution or modification: 11 July 2024
# Review: TCW; 11 July 2024
###############################################################################
# Note


###############################################################################
# Organize arguments.

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
# Organize paths.



###############################################################################
# Organize parameters.


###############################################################################
# Execute procedure.

# TODO: need to activate appropriate Python environment
# TODO: need to export global variables, including any necessary paths

# --order=pos, # specify that sort sequence is by alignment position

python3 -m HTSeq.scripts.count \
--format "bam" \
--minaqual 10 \
--idattr "gene_id" \
--type "gene" \
--order "pos" \
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
$path_file_source_bam \
$path_file_annotation_gtf


###############################################################################
# End.
