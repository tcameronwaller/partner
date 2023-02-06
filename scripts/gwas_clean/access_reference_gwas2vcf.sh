#!/bin/bash

################################################################################
# Notes:

################################################################################

################################################################################
# Organize arguments.

path_directory_parent=${1} # full path to parent directory within which to create child directories and save files
report=${2} # whether to print reports

################################################################################
# Organize paths.

path_directory_sequence="${path_directory_parent}/genome_sequence"
path_directory_dbsnp="${path_directory_parent}/dbsnp"

# Initialize directory.
#rm -r $path_directory_sequence
#rm -r $path_directory_dbsnp
mkdir -p $path_directory_sequence
mkdir -p $path_directory_dbsnp

###########################################################################
# Execute procedure.



################################################################################
# Option 1: Access references files organized specifically for GWAS2VCF.



##########
# Genomic sequence.
cd $path_directory_sequence
# GRCh36/hg18/b36
#wget http://fileserve.mrcieu.ac.uk/ref/2.8/b36/human_b36_both.fasta.gz; gzip -d human_b36_both.fasta.gz
#wget http://fileserve.mrcieu.ac.uk/ref/2.8/b36/human_b36_both.fasta.fai
#wget http://fileserve.mrcieu.ac.uk/ref/2.8/b36/human_b36_both.dict

# GRCh37/hg19/b37
#wget http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.gz; gzip -d human_g1k_v37.fasta.gz
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.gz # <-- need to access again... TCW; 6 February 2023
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.fai
wget http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.dict

# GRCh38/hg38/b38
#wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
#wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
#wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict

##########
# dbSNP.
cd $path_directory_dbsnp
# GRCh37/hg19/b37
wget http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.b37.vcf.gz .
wget http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.b37.vcf.gz.tbi .

# GRCh38/hg38/b38
#wget http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.gz .
#wget http://fileserve.mrcieu.ac.uk/dbsnp/dbsnp.v153.hg38.vcf.gz.tbi .



################################################################################
# Option 2: Access references files from original sources.


# dbSNP.
# Need to translate names of chromosomes in the dbSNP reference.
# Use the "annotate" funtion in BCFTools.

# See: https://mrcieu.github.io/gwas2vcf/install/

# See also script: ".../promiscuity/scripts/access_genetic_reference/call_access_translate_sequence_dbsnp_human_grch37_grch38_vcf"
