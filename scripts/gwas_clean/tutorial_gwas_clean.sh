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
path_directory_temporary="${path_directory_product}/temporary_tcw_test_output_format" # hopefully unique
path_directory_reference_gwas2vcf="${path_directory_product}/reference_gwas2vcf"

# Files.
#identifier_gwas="32042192_ruth_2020_testosterone_female"
identifier_gwas="30367059_teumer_2018_tsh_female"
#path_file_gwas_source=
path_file_gwas_standard_source="${path_directory_dock}/hormone_genetics/gwas_format_standard/${identifier_gwas}.txt.gz"
path_file_munge_report="${path_directory_product}/${identifier_gwas}_munge_report.log"
path_file_gwas_product="${path_directory_product}/${identifier_gwas}.txt"
path_file_gwas_product_gz="${path_directory_product}/${identifier_gwas}.txt.gz"
path_file_gwas2vcf_parameter="${path_directory_process}/promiscuity/scripts/gwas_clean/parameter_gwas_standard_to_gwas2vcf.json"
#path_file_reference_genome_sequence="${path_directory_reference_gwas2vcf}/genome_sequence/human_g1k_v37.fasta.gz"
path_file_reference_genome_sequence="${path_directory_reference_gwas2vcf}/genome_sequence/human_g1k_v37_test.fasta.gz"
path_file_reference_dbsnp="${path_directory_reference_gwas2vcf}/dbsnp/dbsnp.v153.b37.vcf.gz"


# Temporary files.
name_base_file_gwas_product="$(basename $path_file_gwas_product .txt.gz)"
path_file_temporary_gwas_munge_source="${path_directory_temporary}/${identifier_gwas}_before_munge.txt"
path_file_temporary_gwas_munge_source_compress="${path_directory_temporary}/${identifier_gwas}_before_munge.txt.gz"
path_file_temporary_gwas_munge_product="${path_directory_temporary}/${identifier_gwas}_after_munge.txt.gz"
path_file_temporary_gwas_munge_standard="${path_directory_temporary}/${identifier_gwas}_munge_standard.txt"
path_file_temporary_gwas_munge_standard_compress="${path_directory_temporary}/${identifier_gwas}_munge_standard.txt.gz"

#path_file_temporary_gwas_munge="${path_directory_temporary}/${name_base_file_gwas_product}_munge.txt"
path_file_temporary_gwas_decompress="${path_directory_temporary}/${name_base_file_gwas_product}.txt"
path_file_temporary_gwas_vcf="${path_directory_temporary}/${name_base_file_gwas_product}.vcf"
path_file_temporary_gwas_vcf_compress="${path_directory_temporary}/${name_base_file_gwas_product}.vcf.gz"
name_base_file_genome_sequence="$(basename $path_file_reference_genome_sequence .fasta.gz)"
path_file_temporary_genome_decompress="${path_directory_temporary}/${name_base_file_genome_sequence}.fasta"
path_file_temporary_gwas_nhgriebi_vcf="${path_directory_temporary}/gwas_nhgri_ebi_gwas_catalog_format.vcf.gz"
path_file_temporary_gwas_nhgriebi_tsv="${path_directory_temporary}/gwas_nhgri_ebi_gwas_catalog_format.tsv"

path_file_temporary_gwas_postvcf_standard="${path_directory_temporary}/gwas_postvcf_format.txt"
path_file_temporary_gwas_postvcf_standard_gz="${path_directory_temporary}/gwas_postvcf_format.txt.gz"

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

# TODO: TCW; 6 February 2023
# Plan:
# 1. install GZ-Sort
# 2. set up Python 3.8 environment for SumStatsRehab
# 3. install SumStatsRehab
# 4. test SumStatsRehab


# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04920-7#Sec1
# "SumStatsRehab" might be able to introdue rsIDs?
# PubMed: 36284273
# https://github.com/Kukuster/SumStatsRehab

# Dependency of SumStatsRehab: GZ-Sort
# http://kmkeen.com/gz-sort/
# https://github.com/keenerd/gz-sort

# TODO: TCW; 6 February 2023

# TODO: Need to use Bioconductor MungeSumstats first to introduce rsIDs

##########
# Munge GWAS summary statistics in Bioconductor package "MungeSumstats".
# PubMed: 34601555
# Host of MungeSumstats: https://github.com/neurogenomics/MungeSumstats
# Getting Started Guide: https://neurogenomics.github.io/MungeSumstats/articles/MungeSumstats.html
# Documentation: https://bioconductor.org/packages/release/bioc/html/MungeSumstats.html
# Documentation: https://bioconductor.org/packages/release/bioc/vignettes/MungeSumstats/inst/doc/MungeSumstats.html
# Documentation: https://bioconductor.org/packages/release/bioc/manuals/MungeSumstats/man/MungeSumstats.pdf
# Note: MungeSumstats by default removes any SNPs on chromosomes X, Y, or the
# mitochondrial chromosome, but it is possible to suppress this behavior.

if false; then
  # Translate GWAS summary statistics from standard format to a format
  # compatible with MungeSumstats.
  # Source Format (Team Standard)
  # Effect allele: "A1"
  # Delimiter: white space
  # Columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT
  # Product Format (MungeSumstats)
  # Effect allele: "A2"
  # Delimiter: white space
  # Columns: SNP CHR BP A1 A2 FRQ BETA SE P N Z INFO N_CAS N_CON
  echo "SNP CHR BP A1 A2 FRQ BETA SE P N Z INFO N_CAS N_CON" > $path_file_temporary_gwas_munge_source
  zcat $path_file_gwas_standard_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    print $1, $2, $3, toupper($5), toupper($4), $6, $7, $8, $9, $10, $11, $12, $13, $14
  }' >> $path_file_temporary_gwas_munge_source
  # Compress file format.
  gzip -cvf $path_file_temporary_gwas_munge_source > $path_file_temporary_gwas_munge_source_compress
  # Call MungeSumstats.
  Rscript $path_script_mungesumstats $path_file_temporary_gwas_munge_source_compress $path_file_temporary_gwas_munge_product > $path_file_munge_report
  # Translate GWAS summary statistics from MungeSumstats format to standard
  # format.
  # Source Format (MungeSumstats)
  # Effect allele: "A2"
  # Delimiter: white space
  # Columns: SNP CHR BP A1 A2 FRQ BETA SE P N Z INFO N_CAS N_CON
  # Product Format (Team Standard)
  # Effect allele: "A1"
  # Delimiter: white space
  # Columns: SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT
  echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_gwas_munge_standard
  zcat $path_file_temporary_gwas_munge_product | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    print $1, $2, $3, toupper($5), toupper($4), $6, $7, $8, $9, $10, $11, $12, $13, $14
  }' >> $path_file_temporary_gwas_munge_standard
  # Compress file format.
  gzip -cvf $path_file_temporary_gwas_munge_standard > $path_file_temporary_gwas_munge_standard_compress
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
  #gzip -dcvf $path_file_gwas_standard_source > $path_file_temporary_gwas_decompress
  # Switch to tab delimiters (field separators).
  #zcat $path_file_gwas_standard_source | awk 'BEGIN {FS = " "; OFS = "\t"} {
  #  print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14
  #}' >> $path_file_temporary_gwas_decompress
  # Keep same delimiters (field separators), but only keep first count of lines.
  zcat $path_file_gwas_standard_source | awk 'NR < 10000 {
    print $0
  }' >> $path_file_temporary_gwas_decompress
  wc -l $path_file_temporary_gwas_decompress
  # Keep same delimiters (field separators), but only keep chromosomes 1-22.
  #chromosomes=($(for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; do echo $i; done))
  #zcat $path_file_gwas_standard_source | awk 'NR < 10000 {
  #  print $0
  #}' >> $path_file_temporary_gwas_decompress
  #wc -l $path_file_temporary_gwas_decompress


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
  export MKL_NUM_THREADS=8
  export NUMEXPR_NUM_THREADS=8
  export OMP_NUM_THREADS=8
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
  # Examine file.
  #$path_bcftools head $path_file_temporary_gwas_vcf

fi



##########
# Translate GWAS summary statistics from GWAS-VCF format to standard format.
# Use tool GWAS2VCF.
# documentation: https://mrcieu.github.io/gwas2vcf/downstream/#convert

# Note: TCW; 7 February 2023
# It seems to be a problem to request a field that does not exist in the
# specific GWAS-VCF file.
# It might be necessary to query GWAS-VCF files differently for those with or
# without counts of cases and controls.

if false; then
  # This block uses the same extraction code from the GWAS2VCF documentation.
  # The documentation claims that this extraction uses the NHGRI-EBI GWAS
  # Catalog format.
  # This extraction omits the count of observations used to determine effect for
  # each SNP.
  # Translate from GWAS-VCF format to NHGRI-EBI GWAS Catalog format.
  # Notice that the extraction converts the logarithm of the p-value to the
  # p-value itself.
  $path_bcftools query \
  -e 'ID == "."' \
  -f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t%AF\t[%ES\t%SE]\n' \
  $path_file_temporary_gwas_vcf_compress | \
  awk 'BEGIN {print "variant_id\tp_value\tchromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tbeta\tstandard_error"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' > $path_file_temporary_gwas_nhgriebi_tsv
fi

if true; then
  # Attempt to keep sample size (count).
  # Translate from GWAS-VCF format to NHGRI-EBI GWAS Catalog format.
  # ID: "Study variant identifier"; reference sequence identifier (rsID)
  # ES: "Effect size estimate relative to the alternative allele"
  # SE: "Standard error of effect size estimate"
  # LP: "-log10 p-value for effect estimate"
  # AF: "Alternate allele frequency in the association study"
  # SS: "Sample size used to estimate genetic effect"; count of observations per SNP
  # EZ: "Z-score provided if it was used to derive the EFFECT and SE fields"; z-score
  # NC: "Number of cases used to estimate genetic effect"; count of cases per SNP
  $path_bcftools query \
  -e 'ID == "."' \
  -f '%ID\t[%LP]\t%CHROM\t%POS\t%ALT\t%REF\t%AF\t[%ES]\t[%SE]\t[%SS]\n' \
  $path_file_temporary_gwas_vcf_compress | \
  awk 'BEGIN {print "variant_id\tp_value\tchromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tbeta\tstandard_error\tobservations"}; {OFS="\t"; if ($2==0) $2=1; else if ($2==999) $2=0; else $2=10^-$2; print}' > $path_file_temporary_gwas_nhgriebi_tsv
fi

if true; then
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
  echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_temporary_gwas_postvcf_standard
  cat $path_file_temporary_gwas_nhgriebi_tsv | awk 'BEGIN {FS = "\t"; OFS = " "} NR > 1 {
    print $1, $3, $4, $5, $6, $7, $8, $9, $2, $10, "NA", (1.0), "NA", "NA"
  }' >> $path_file_temporary_gwas_postvcf_standard
  # Compress file format.
  gzip -cvf $path_file_temporary_gwas_postvcf_standard > $path_file_temporary_gwas_postvcf_standard_gz
fi






# TODO: next translate columns to standard format

# TODO: then remove the temporary directory




##########
# TCW; 24 January 2023
# The Bioconductor packange "MungeSumstats" might be preferrable over the
# "+munge" plugin for BCFTools.
##########
# Munge GWAS summary statistics in BCFTools "+munge" plugin.
# documentation: https://github.com/freeseek/score#convert-summary-statistics


#
