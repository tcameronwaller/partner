#!/bin/bash

################################################################################
# Notes:

# Download newest genetic reference files:
# genomic sequences in GRCh37 and GRCh38
# dbSNP in GRCh37 and GRCh38
# Use BCFTools to annotate the reference sequences with SNP rs IDs... is that right?
# UCSC and Ensembl chain files from GRCh37 to GRCh38 and vice versa


# update installation of BCFTools and HTSlib
# install GWAS2VCF
# install BCFTools plugin for GWAS-VCF

# Steps on GWAS sum stats before SBayesR or LDPred2
# 1. convert GWAS sum stats to team standard format
# 2. convert GWAS sum stats to GWAS-VCF format
# 3. run BCFTools +Munge plugin on GWAS in GWAS-VCF format
# 4. convert from GWAS-VCF format to GWAS catalog format
# 5. convert from GWAS catalog format to format for SBayesR and LDPred2

# Steps on output from SBayesR and LDPred2
# 1. convert output to format for GWAS2VCF
# 2. convert to GWAS-VCF format
# 3. run BCFTools +Liftover plugin from GRCh37 to GRCh38
# 4. convert from GWAS-VCF format to GWAS catalog format
# 5. convert from GWAS catalog format to format readable by PLINK2 score function


################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")
path_environment_gwas2vcf="${path_directory_tools}/python/environments/gwas2vcf"
path_gwas2vcf="${path_directory_tools}/gwas2vcf/gwas2vcf/main.py"
path_bcftools=$(<"./tools_bcftools.txt")

#path_plink2=$(<"./tools_plink2.txt")
#path_gctb=$(<"./tools_waller_gctb.txt")
#path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
#path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_product="${path_directory_dock}/test_gwas_clean"
path_directory_temporary="${path_directory_product}/temporary_tcw_2431687" # hopefully unique
path_directory_reference_gwas2vcf="${path_directory_product}/reference_gwas2vcf"

# Files.
#identifier_gwas="32042192_ruth_2020_testosterone_female"
identifier_gwas="30367059_teumer_2018_tsh_female"
#path_file_gwas_source=
path_file_gwas_standard_source="${path_directory_dock}/hormone_genetics/gwas_format_standard/${identifier_gwas}.txt.gz"
path_file_gwas_product="${path_directory_product}/${identifier_gwas}.txt.gz"
path_file_gwas2vcf_parameter="${path_directory_process}/promiscuity/scripts/gwas_clean/parameter_gwas_standard_to_gwas2vcf.json"
#path_file_reference_genome_sequence="${path_directory_reference_gwas2vcf}/genome_sequence/human_g1k_v37.fasta.gz"
path_file_reference_genome_sequence="${path_directory_reference_gwas2vcf}/genome_sequence/human_g1k_v37_test.fasta.gz"
path_file_reference_dbsnp="${path_directory_reference_gwas2vcf}/dbsnp/dbsnp.v153.b37.vcf.gz"


# Temporary files.
name_base_file_gwas_product="$(basename $path_file_gwas_product .txt.gz)"
path_file_temporary_gwas_decompress="${path_directory_temporary}/${name_base_file_gwas_product}.txt"
path_file_temporary_gwas_vcf="${path_directory_temporary}/${name_base_file_gwas_product}.vcf"
name_base_file_genome_sequence="$(basename $path_file_reference_genome_sequence .fasta.gz)"
path_file_temporary_genome_decompress="${path_directory_temporary}/${name_base_file_genome_sequence}.fasta"
path_file_temporary_gwas_nhgriebi_vcf="${path_directory_temporary}/gwas_nhgri_ebi_gwas_catalog_format.vcf.gz"
path_file_temporary_gwas_nhgriebi_tsv="${path_directory_temporary}/gwas_nhgri_ebi_gwas_catalog_format.tsv"

# Scripts.
#path_script_gwas_format_source_to_standard="${path_directory_process}/promiscuity/scripts/gwas_format/format_gwas_team/translate_gwas_30367059_teumer_2018.sh"
path_script_access_reference_gwas2vcf="${path_directory_process}/promiscuity/scripts/gwas_clean/access_reference_gwas2vcf.sh"
#path_script_gwas_format_standard_to_gwas2vcf="${path_directory_process}/promiscuity/scripts/utility/gwas_clean/translate_gwas_standard_to_gwas2vcf.sh" # <-- do not need
#path_script_gwas_clean="${path_directory_process}/promiscuity/scripts/utility/gwas_clean/clean_gwas.sh"

# Initialize directories.
#rm -r $path_directory_product
rm -r $path_directory_temporary
mkdir -p $path_directory_product
mkdir -p $path_directory_temporary
cd $path_directory_product



################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.



##########
# Translate GWAS summary statistics to standard format.
# Product Format (Team Standard)
# delimiter: white space
# columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT
if false; then
  /usr/bin/bash "${path_script_gwas_format}" \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $report
fi



##########
# Installation: GWAS2VCF
# PubMed: 33441155
# GWAS-VCF format specification: https://github.com/MRCIEU/gwas-vcf-specification
# Host of GWAS2VCF: https://github.com/MRCIEU/gwas2vcf
# Documentation for GWAS2VCF: https://mrcieu.github.io/gwas2vcf
# Host of GWASGlue: https://github.com/MRCIEU/gwasglue
# Examples of analyses: https://mrcieu.github.io/gwasglue/articles/
# Refer to notes in script "install_local_software.sh".
# Refer to notes in script "install_python_virtual_environments_packages.sh"
if false; then
  # Activate Virtual Environment.
  source "${path_environment_gwas2vcf}/bin/activate"
  echo "confirm Python Virtual Environment path..."
  which python3
  sleep 5s
  # Test installation of GWAS2VCF.
  python3 $path_gwas2vcf -h
  # Deactivate Virtual Environment.
  deactivate
  which python3
fi



##########
# Access genomic reference information for GWAS2VCF.
if false; then
  /usr/bin/bash $path_script_access_reference_gwas2vcf \
  $path_directory_reference_gwas2vcf \
  $report
fi



##########
# Translate GWAS summary statistics from standard format to GWAS-VCF format.
# Use tool GWAS2VCF.
# Functions.
# 1. Verify and introduce SNPs' rs identifiers from dbSNP reference.

if true; then
  # Decompress the GWAS summary statistics.
  gzip -dcvf $path_file_gwas_standard_source > $path_file_temporary_gwas_decompress
  # Decompress the reference genome sequence.
  gzip -dcvf $path_file_reference_genome_sequence > $path_file_temporary_genome_decompress
  # Define parameters for GWAS2VCF.
  # Index of columns in source GWAS summary statistics bases on zero.
  # Use "ncontrol_col" to designate the column for count of observations in
  # linear GWAS or for the column for count of controls in logistic GWAS.
  # Can also use the parameter "--cohort_cases" to designate a single value for
  # the count of cases in logistic GWAS.
  # Subsequent use of the parameter "--cohort_controls" to will rewrite any
  # SNP-specific counts of controls.
  #echo "{
  #  'chr_col': 1,
  #  'pos_col': 2,
  #  'snp_col': 0,
  #  'ea_col': 3,
  #  'oa_col': 4,
  #  'beta_col': 6,
  #  'se_col': 7,
  #  'ncontrol_col': 9,
  #  'pval_col': 8,
  #  'eaf_col': 5,
  #  'delimiter': '\t',
  #  'header': true,
  #  'build': 'GRCh37'
  #}" > $path_file_temporary_gwas2vcf_parameter
  # Activate Virtual Environment.
  source "${path_environment_gwas2vcf}/bin/activate"
  echo "confirm Python Virtual Environment path..."
  which python3
  sleep 5s
  # Call GWAS2VCF.
  python3 $path_gwas2vcf \
  --data $path_file_temporary_gwas_decompress \
  --json $path_file_gwas2vcf_parameter \
  --id $identifier_gwas \
  --ref $path_file_temporary_genome_decompress \
  --dbsnp $path_file_reference_dbsnp \
  --out $path_file_temporary_gwas_vcf \
  --log INFO
  # Deactivate Virtual Environment.
  deactivate
  which python3
fi



##########
# Translate GWAS summary statistics from GWAS-VCF format to standard format.
# Use tool GWAS2VCF.
# documentation: https://mrcieu.github.io/gwas2vcf/downstream/#convert

if true; then
  # Translate from GWAS-VCF format to NHGRI-EBI GWAS Catalog format.
  $path_bcftools query \
  -e 'ID == "."' \
  -f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t%AF\t[%ES\t%SE]\n' \
  $path_file_temporary_gwas_nhgriebi_vcf | \
  awk 'BEGIN {print "variant_id\tp_value\tchromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tbeta\tstandard_error"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' > $path_file_temporary_gwas_nhgriebi_tsv
fi


# TODO: next translate columns to standard format

# TODO: then remove the temporary directory



##########
# Munge GWAS summary statistics in Bioconductor package "MungeSumstats".
# Documentation: https://bioconductor.org/packages/release/bioc/manuals/MungeSumstats/man/MungeSumstats.pdf
# Note: MungeSumstats removes any SNPs on chromosomes X, Y, or mitochondrion.


##########
# TCW; 24 January 2023
# The Bioconductor packange "MungeSumstats" might be preferrable over the
# "+munge" plugin for BCFTools.
##########
# Munge GWAS summary statistics in BCFTools "+munge" plugin.
# documentation: https://github.com/freeseek/score#convert-summary-statistics


#
