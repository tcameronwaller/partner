#!/bin/bash

# Note:
# Each linear GWAS (30,000 - 200,000 records; 22 chromosomes) requires about
# 5-7 hours to run on the grid.

################################################################################
# Organize argument variables.

path_table_phenotypes_covariates=${1} # full path to file for table with phenotypes and covariates
path_report=${2} # full path to parent directory for GWAS summary statistics
analysis=${3} # unique name for association analysis
phenotypes=${4} # names of table's column or columns for single or multiple phenotypes, dependent variables
covariates=${5} # name of table's columns for covariates, independent variables
threads=${6} # count of processing threads to use
maf=${7} # minor allele frequency threshold filter

###########################################################################
# Organize variables.

# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_plink2=$(<"./tools_plink2.txt")
path_ukb_genotype=$(<"./ukbiobank_genotype.txt")

# Iterate on chromosomes.
#chromosomes=("x")
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "x")
for chromosome in "${chromosomes[@]}"; do
  #echo "chromosome: ${chromosome}"

  # Specify target directory.
  path_chromosome="$path_report/chromosome_${chromosome}"
  # Determine whether the temporary directory structure already exists.
  if [ ! -d $path_chromosome ]; then
      # Directory does not already exist.
      # Create directory.
      mkdir -p $path_chromosome
  fi
  cd $path_chromosome

  # Specify source directory and files.
  if [[ "$chromosome" == "x" ]]; then
    # On 25 October 2021, Gregory Jenkins told me:
    # "I believe you are right on 'ukb22828_cX_b0_v3_s486620.sample' and 'ukb22828_cX_b0_v3.bgen'"
    # "They were missing chrom X and I asked them to download that in March/April of this year I think?  So that's why it looks different."
    path_genotype="$path_ukb_genotype/Chromosome/ukb22828_cX_b0_v3.bgen"
    path_sample="$path_ukb_genotype/Chromosome/ukb22828_cX_b0_v3_s486620.sample"
  else
    path_genotype="$path_ukb_genotype/Chromosome/ukb_imp_chr${chromosome}_v3.bgen"
    path_sample="$path_ukb_genotype/Chromosome/ukb46237_imp_chr${chromosome}_v3_s487320.sample"
  fi

  if true; then
    # Call PLINK2.
    # This call to PLINK2 includes "--freq" and produces a report of allele
    # frequencies.
    # PLINK2 command "--glm" drives genotypic association analyses, either
    # linear or logistic regressions across Simple Nucleotide Polymorphisms
    # (SNPs).
    # 90,000 Mebibytes (MiB) is 94.372 Gigabytes (GB)
    # --pfilter 1 \
    # --pfilter drops SNPs with null p-values and any beyond threshold (such as 1)
    # But, maybe pfilter is actually a problem.
    # Argument "--reference-allele" allows to give an explicit list of alleles to designate as "A1".
    # "Warning: No --bgen REF/ALT mode specified ('ref-first', 'ref-last', or 'ref-unknown'). This will be required as of alpha 3."
    $path_plink2 \
    --memory 90000 \
    --threads $threads \
    --bgen $path_genotype ref-first \
    --sample $path_sample \
    --keep $path_table_phenotypes_covariates \
    --pheno $path_table_phenotypes_covariates \
    --pheno-name $phenotypes \
    --covar $path_table_phenotypes_covariates \
    --covar-name $covariates \
    --maf $maf \
    --xchr-model 1 \
    --freq \
    --glm omit-ref no-x-sex hide-covar \
    --out report
  fi
done



##########
