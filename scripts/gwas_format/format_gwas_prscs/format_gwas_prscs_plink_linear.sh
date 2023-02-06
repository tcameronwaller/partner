#!/bin/bash

###########################################################################
###########################################################################
###########################################################################

# "Organize GWAS summary statistics for analysis in PRS-CS."

# review:
# TCW, 10 May 2022

###########################################################################
###########################################################################
###########################################################################
# Notes

################################################################################
# Organize variables.
path_gwas_source=${1} # full path to source file with GWAS summary statistics, with gzip compression
path_gwas_constraint=${2} # full path to temporary file for constraint of GWAS summary statistics to bounds
path_gwas_format=${3} # full path to file for formatted GWAS summary statistics
path_gwas_format_compress=${4} # full path to file for formatted GWAS summary statistics after compression
report=${5} # whether to print reports

###########################################################################
# Execute procedure.

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "path to source file: " $path_gwas_source
  echo "path to target file: " $path_gwas_format_compress
fi

# Here are some useful regular expressions to evaluate values in "awk".
# ( $12 ~ /^[0-9]+$/ ); ( $12 ~ /^[[:alpha:]]+$/ ); ( $12 ~ /^[[:punct:]]+$/ )

# Remove any previous versions of temporary files.
rm $path_gwas_constraint
rm $path_gwas_format
rm $path_gwas_format_compress

##########
# Constrain GWAS probability values within bounds for LDSC.
# Probability values in column "P" from PLINK2 have values of "NA" for missing
# or floating point values between zero and one.
# LDSC uses 64-bit floating point precision (double precision) to represent
# values from +/- 2.23E-308 to +/- 1.80E308 (https://github.com/bulik/ldsc/issues/144).
# Constrain probability values from 1.0E-305 to 1.0.
zcat $path_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR == 1' > $path_gwas_constraint
zcat $path_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ( NF != 13)
    # Skip any rows with incorrect count of column fields.
    next
  else if ( ( $13 != "NA" ) && ( ($13 + 0) < 1.0E-305 ) )
    # Constrain probability value.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, ( 1.0E-305 )
  else if ( ( $13 != "NA" ) && ( ($13 + 0) > 1.0 ) )
    # Constrain probability value.
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, ( 1.0 )
  else
    print $0
  }' >> $path_gwas_constraint

##########
# Format of GWAS summary statistics for PRS-CS.
# https://github.com/getian107/PRScs

# Format of GWAS reports by PLINK2 for linear regression (".glm.linear").
# https://www.cog-genomics.org/plink/2.0/formats
# PLINK2 report format is similar for logistic regression (".glm.logistic").
# Logistic report has "OR" in place of "BETA", but positions are the same.

# PLINK2's format documentation describes column "REF" as "Reference allele" and
# column "ALT" as "All alternate alleles, comma-separated". Column "REF" and
# column "ALT" are always different (TCW, 7 July 2021).
# PLINK2's format documentation also describes column "A1" as "Counted allele in
# regression". Column "A1" always matches either column "REF" or column "ALT"
# (TCW, 7 July 2021).

# The designations of "reference" and "alternate" alleles are irrelevant to
# downstream analyses on GWAS summary statistics (such as LDSC).
# What matters is which allele PLINK2 counted as the "effect" allele in GWAS
# regression and which allele was the other allele.

# Column "A1" is the "effect" allele that PLINK2 counted in the regression,
# corresponding to the coefficient (beta).
# Determine the "non-effect" or "other" allele from either the "reference" or
# "alternate" allele, whichever differs from the "effect" allele.

# description: ............................ PRS-CS column ......... PLINK column .......... position
# variant identifier (RS ID): .............  "SNP" ................  "ID" .................. 3
# alternate allele (effect allele): .......  "A1" .................  "A1" .................. 6
# reference allele (non-effect allele): ...  "A2" .................  "REF" or "ALT" ........ 4 or 5
# effect (coefficient or odds ratio): .....  "BETA" or "OR" .......  "BETA" ................ 10
# probability (p-value): ..................  "P" ..................  "P" ................... 13

# Organize information from linear GWAS.
# Select relevant columns and place them in the correct order.
echo "SNP A1 A2 BETA P" > $path_gwas_format
cat $path_gwas_constraint | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
  if ($6 == $5 && $6 != $4)
    print $3, toupper($6), toupper($4), $10, $13
  else if ($6 == $4 && $6 != $5)
    print $3, toupper($6), toupper($5), $10, $13
  else
    next
  }' >> $path_gwas_format

# Compress file format.
gzip -cvf $path_gwas_format > $path_gwas_format_compress

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "after bound constraints:"
  head -10 $path_gwas_constraint
  echo "after format:"
  head -10 $path_gwas_format
  echo "----------"
fi

# Remove temporary files.
rm $path_gwas_constraint
rm $path_gwas_format
