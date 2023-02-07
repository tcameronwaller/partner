#!/bin/bash

################################################################################
# Notes:



################################################################################



################################################################################
# Organize arguments.

path_file_gwas_source=${1} # full path to file for source GWAS summary statistics with GZip compression
path_file_gwas_product=${2} # full path to file for product GWAS summary statistics in format with GZip compression
threads=${3} # count of processing threads to use
report=${4} # whether to print reports

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")
path_environment_gwas2vcf="${path_directory_tools}/python/environments/gwas2vcf"
path_gwas2vcf="${path_directory_tools}/gwas2vcf/gwas2vcf/main.py"
path_bcftools=$(<"./tools_bcftools.txt")
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_reference_gwas2vcf="${path_directory_reference}/gwas2vcf"
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_product="$(dirname $path_file_gwas_product)"
path_directory_temporary="${path_directory_product}/temporary_tcw_test_5271693" # hopefully unique

# Files.
name_base_file_gwas_product="$(basename $path_file_gwas_product .txt.gz)"
identifier_gwas=name_base_file_gwas_product
path_file_munge_report="${path_directory_product}/${name_base_file_gwas_product}_munge_report.log"
path_file_gwas2vcf_report="${path_directory_product}/${name_base_file_gwas_product}_gwas2vcf_report.log"
path_file_gwas2vcf_parameter="${path_directory_process}/promiscuity/scripts/gwas_clean/parameter_gwas_standard_to_gwas2vcf.json"
path_file_reference_genome_sequence="${path_directory_reference_gwas2vcf}/genome_sequence/human_g1k_v37.fasta.gz"
path_file_reference_dbsnp="${path_directory_reference_gwas2vcf}/dbsnp/dbsnp.v153.b37.vcf.gz"

# Temporary files.
#path_ftemp_gwas_premunge="${path_directory_temporary}/${identifier_gwas}_before_munge.txt"
#path_ftemp_gwas_premunge_gz="${path_directory_temporary}/${identifier_gwas}_before_munge.txt.gz"
#path_ftemp_gwas_postmunge_gz="${path_directory_temporary}/${identifier_gwas}_after_munge.txt.gz"
#path_ftemp_gwas_postmunge_standard="${path_directory_temporary}/${identifier_gwas}_munge_standard.txt"
path_ftemp_gwas_source_decomp="${path_directory_temporary}/${name_base_file_gwas_product}_source.txt"
name_base_file_genome_sequence="$(basename $path_file_reference_genome_sequence .fasta.gz)"
path_ftemp_genome_decomp="${path_directory_temporary}/${name_base_file_genome_sequence}.fasta"
path_ftemp_gwas_vcf="${path_directory_temporary}/${name_base_file_gwas_product}.vcf"
path_ftemp_gwas_vcf_gz="${path_directory_temporary}/${name_base_file_gwas_product}.vcf.gz"
path_ftemp_gwas_postvcf_tsv="${path_directory_temporary}/${name_base_file_gwas_product}_postvcf.tsv"
path_ftemp_gwas_postvcf_standard_text="${path_directory_temporary}/${name_base_file_gwas_product}_postvcf_standard_format.txt"

# Scripts.

# Initialize files.
#rm $path_file_munge_report
rm $path_file_gwas2vcf_report
rm $path_file_gwas_product

# Initialize directories.
#rm -r $path_directory_product
rm -r $path_directory_temporary
mkdir -p $path_directory_product
mkdir -p $path_directory_temporary
cd $path_directory_product



################################################################################
# Organize parameters.



################################################################################
# Execute procedure.


##########
# Apply any preparation to the GWAS summary statistics in the source file.
if true; then
  # Decompress the GWAS summary statistics.
  #gzip -dcvf $path_file_gwas_source > $path_ftemp_gwas_source_decomp
  # Keep same delimiters (field separators), but only keep first count of lines.
  zcat $path_file_gwas_source | awk 'NR < 100000 {
    print $0
  }' >> $path_ftemp_gwas_source_decomp
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "----------"
    echo "----------"
    echo "Count of lines in source GWAS summary statistics."
    echo "Count:"
    wc -l $path_ftemp_gwas_source_decomp
    echo "----------"
    echo "----------"
    echo "----------"
  fi
fi


##########
# Translate GWAS summary statistics from standard format to GWAS-VCF format.
# Use tool GWAS2VCF.
# PubMed: 33441155
# GWAS-VCF format specification: https://github.com/MRCIEU/gwas-vcf-specification
# Host of GWAS2VCF: https://github.com/MRCIEU/gwas2vcf
# Documentation for GWAS2VCF: https://mrcieu.github.io/gwas2vcf
# Host of GWASGlue: https://github.com/MRCIEU/gwasglue
# Examples of analyses: https://mrcieu.github.io/gwasglue/articles/
# Refer to notes in script "install_local_software.sh".
# Refer to notes in script "install_python_virtual_environments_packages.sh"

if true; then
  # Decompress the reference genome sequence.
  gzip -dcvf $path_file_reference_genome_sequence > $path_ftemp_genome_decomp
  # Define parameters for GWAS2VCF within a text file in "json" format.
  # Parameters do need to be within a separate text file.
  # Index of columns in source GWAS summary statistics bases on zero.
  # Use "ncontrol_col" to designate the column for count of observations in
  # linear GWAS or for the column for count of controls in logistic GWAS.
  # Can also use the parameter "--cohort_cases" to designate a single value for
  # the count of cases in logistic GWAS.
  # Subsequent use of the parameter "--cohort_controls" to will rewrite any
  # SNP-specific counts of controls.
  # Activate Virtual Environment.
  source "${path_environment_gwas2vcf}/bin/activate"
  echo "confirm Python Virtual Environment path..."
  which python3
  sleep 5s
  # Force Python program (especially SciPy) not to use all available cores on a
  # cluster computation node.
  export MKL_NUM_THREADS=$threads
  export NUMEXPR_NUM_THREADS=$threads
  export OMP_NUM_THREADS=$threads
  # Call GWAS2VCF.
  # GWAS2VCF automatically applies GZip compression to the file in VCF format
  # and calculates a Tabix index.
  python3 $path_gwas2vcf \
  --data $path_ftemp_gwas_source_decomp \
  --json $path_file_gwas2vcf_parameter \
  --id $identifier_gwas \
  --ref $path_ftemp_genome_decomp \
  --dbsnp $path_file_reference_dbsnp \
  --out $path_ftemp_gwas_vcf \
  --log INFO 2>&1 | tee $path_file_gwas2vcf_report
  # Deactivate Virtual Environment.
  deactivate
  which python3
  # Examine file.
  #$path_bcftools head $path_ftemp_gwas_vcf_gz
fi



##########
# Extract information from GWAS summary statistics in GWAS-VCF format.
# Convert format to a tabular text file.
# Use tool GWAS2VCF.
# documentation: https://mrcieu.github.io/gwas2vcf/downstream/#convert
# ID: "Study variant identifier"; reference sequence identifier (rsID)
# ES: "Effect size estimate relative to the alternative allele"
# SE: "Standard error of effect size estimate"
# LP: "-log10 p-value for effect estimate"
# AF: "Alternate allele frequency in the association study"
# SS: "Sample size used to estimate genetic effect"; count of observations per SNP
# EZ: "Z-score provided if it was used to derive the EFFECT and SE fields"; z-score
# NC: "Number of cases used to estimate genetic effect"; count of cases per SNP
# Note: TCW; 7 February 2023
# This extraction converts the logarithm of the p-value to the p-value itself.
# It seems to be a problem to request a field that does not exist in the
# specific GWAS-VCF file.
# It might be necessary to query GWAS-VCF files differently for those with or
# without counts of cases and controls.
if true; then
  $path_bcftools query \
  -e 'ID == "."' \
  -f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t%AF\t[%ES]\t[%SE]\t[%SS]\n' \
  $path_ftemp_gwas_vcf_gz | \
  awk 'BEGIN {print "variant_id\tp_value\tchromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tbeta\tstandard_error\tobservations"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' > $path_ftemp_gwas_postvcf_tsv
fi



##########
# Translate GWAS summary statistics after GWAS2VCF clean to standard format.
# Translate GWAS summary statistics from GWAS2VCF export format to
# standard format.
# Source Format: Export from GWAS2VCF GWAS-VCF
# Effect allele: "effect_allele"
# Delimiter: tab
# Columns: variant_id p_value chromosome base_pair_location effect_allele other_allele effect_allele_frequency beta standard_error observations
# Columns: 1          2       3          4                  5             6            7                       8    9              10
# Product Format: Team Standard
# Effect allele: "A1"
# Delimiter: white space
# Columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT
if true; then
  echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_ftemp_gwas_postvcf_standard_text
  cat $path_ftemp_gwas_postvcf_tsv | awk 'BEGIN {FS = "\t"; OFS = " "} NR > 1 {
    print $1, $3, $4, $5, $6, $7, $8, $9, $2, $10, "NA", (1.0), "NA", "NA"
  }' >> $path_ftemp_gwas_postvcf_standard_text
  # Compress file format.
  gzip -cvf $path_ftemp_gwas_postvcf_standard_text > $path_file_gwas_product
fi



##########
# Remove temporary directories and files.
# Suppress this block for debugging.
if true; then
  rm -r $path_directory_temporary
fi





#
