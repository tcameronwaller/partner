
################################################################################
# Note

# Specify desired version of dbSNP reference using cooridnates from appropriate
# assembly of the human genome.

# BCFTools can work with dbSNP reference vcf file in BGZF (bgzip) format with a
# Tabix index (.tbi) for more efficient performance.


# BGZIP documentation
# http://www.htslib.org/doc/bgzip.html

# BCFTools documentation
# https://samtools.github.io/bcftools/bcftools.html

# Example of simple script call format for using BCFTools Annotate to introduce rsIDs
# https://gist.github.com/obenshaindw/99440ac43de07548e453
# "/usr/bin/htslib/bcftools/bcftools annotate -a /reference/dbsnp_137.b37.vcf.gz -c ID vcf_to_add_id_to.vcf"


# "Cheat sheet" on BCFTools
# https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b

################################################################################
