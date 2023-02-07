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
path_directory_temporary="${path_directory_product}/temporary_tcw_5271693" # hopefully unique

# Files.
name_base_file_gwas_product="$(basename $path_file_gwas_product .txt.gz)"
identifier_gwas=name_base_file_gwas_product
path_file_munge_report="${path_directory_product}/${name_base_file_gwas_product}_munge_report.log"
path_file_gwas2vcf_report="${path_directory_product}/${name_base_file_gwas_product}_gwas2vcf_report.log"
path_file_gwas2vcf_parameter="${path_directory_process}/promiscuity/scripts/gwas_clean/parameter_gwas_standard_to_gwas2vcf.json"
path_file_reference_genome_sequence="${path_directory_reference_gwas2vcf}/genome_sequence/human_g1k_v37.fasta.gz"
#path_file_reference_genome_sequence="${path_directory_reference_gwas2vcf}/genome_sequence/human_g1k_v37_test.fasta.gz"
path_file_reference_dbsnp="${path_directory_reference_gwas2vcf}/dbsnp/dbsnp.v153.b37.vcf.gz"

# Temporary files.
path_ftemp_gwas_premunge="${path_directory_temporary}/${identifier_gwas}_before_munge.txt"
path_ftemp_gwas_premunge_gz="${path_directory_temporary}/${identifier_gwas}_before_munge.txt.gz"
path_ftemp_gwas_postmunge_gz="${path_directory_temporary}/${identifier_gwas}_after_munge.txt.gz"
path_ftemp_gwas_postmunge_standard="${path_directory_temporary}/${identifier_gwas}_munge_standard.txt"


path_ftemp_gwas_vcf="${path_directory_temporary}/${name_base_file_gwas_product}.vcf"
path_ftemp_gwas_vcf_gz="${path_directory_temporary}/${name_base_file_gwas_product}.vcf.gz"


name_base_file_genome_sequence="$(basename $path_file_reference_genome_sequence .fasta.gz)"
path_file_temporary_genome_decompress="${path_directory_temporary}/${name_base_file_genome_sequence}.fasta"
path_file_temporary_gwas_nhgriebi_vcf="${path_directory_temporary}/gwas_nhgri_ebi_gwas_catalog_format.vcf.gz"
path_file_temporary_gwas_nhgriebi_tsv="${path_directory_temporary}/gwas_nhgri_ebi_gwas_catalog_format.tsv"

# Scripts.
#path_script_gwas_format_source_to_standard="${path_directory_process}/promiscuity/scripts/gwas_format/format_gwas_team/translate_gwas_30367059_teumer_2018.sh"
path_script_access_reference_gwas2vcf="${path_directory_process}/promiscuity/scripts/gwas_clean/access_reference_gwas2vcf.sh"
path_script_mungesumstats="${path_directory_process}/promiscuity/scripts/gwas_clean/execute_bioconductor_mungesumstats_format_sumstats.R"
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



################################################################################
# Execute procedure.



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
  # Decompress the GWAS summary statistics.
  gzip -dcvf $path_file_gwas_standard_source > $path_ftemp_gwas_postmunge_standard
  # Decompress the reference genome sequence.
  gzip -dcvf $path_file_reference_genome_sequence > $path_file_temporary_genome_decompress
  # Define parameters for GWAS2VCF within a text file in "json" format.
  # Parameters do need to be within a separate text file.
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
  # Force Python program (especially SciPy) not to use all available cores on a
  # cluster computation node.
  export MKL_NUM_THREADS=$threads
  export NUMEXPR_NUM_THREADS=$threads
  export OMP_NUM_THREADS=$threads
  # Call GWAS2VCF.
  # GWAS2VCF automatically applies GZip compression to the file in VCF format
  # and calculates a Tabix index.
  python3 $path_gwas2vcf \
  --data $path_ftemp_gwas_postmunge_standard \
  --json $path_file_gwas2vcf_parameter \
  --id $identifier_gwas \
  --ref $path_file_temporary_genome_decompress \
  --dbsnp $path_file_reference_dbsnp \
  --out $path_ftemp_gwas_vcf \
  --log INFO 2>&1 | tee $path_file_gwas2vcf_report
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
  $path_ftemp_gwas_vcf_gz | \
  awk 'BEGIN {print "variant_id\tp_value\tchromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tbeta\tstandard_error"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' > $path_file_temporary_gwas_nhgriebi_tsv
fi


# TODO: next translate columns to standard format




# TODO: then remove the temporary directory






#
