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
# name, it is necessary to specify the "--order=pos" option.


# TODO: TCW; 11 July 2024
# The HTSeq documentation specifies that it's necessary to use the "end_included"
# argument of the "HTSeq.GFF_Reader()" to specify whether or not the GTF file
# specifies position ranges that include the end position. If the length of the
# range of an exon is divisible by 3 (codon) then the "end_included" parameter
# must be "False".

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
--order=pos

###############################################################################
# End.
