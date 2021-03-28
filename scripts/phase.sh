
SCRIPT=`basename $0`
# Ensure that arguments are passed to script, if not display help
if [ "$#" -ne 5 ]; then
cat << EOF
Usage: sh ${SCRIPT} input mapfile outprefix chrom threads
Parse a plink file by chromosome, convert to vcf, phase vcf with shapeit4, and store result as BCF (and index BCF).

input: prefix of plink files

mapfile: genetic map positions for shapeit4

outprefix: prefix to use for output bcf.  Will be overwritten if they already exist.

chrom: which chromosome to parse

threads: number of processors to use for phasing.  Do not need more than 6.

EOF
  exit 1
fi

# check that input file exists
if ! [ -e "$1" ]; then
  echo "$1 not found" >&2
  exit 1
fi

# Check that mapfile exists
if ! [ -e "$2" ]; then
  echo "$2 not found" >&2
  exit 1
fi

#Load modules
module load plink/1.90b6.10
module load htslib/1.6
module load bcftools/1.9

#Input info
#Path to directory containing master PLINK files
INPUT="$1"

#Location of genetic map files
MAPFILE="$2"

OUTFILE="$3"

CHROM="$4"

#Number of threads used for phasing.  12 may be overkill.  Phasing completes in < 30 minutes.
THREADS="$5"

#Convert plink files to vcf and separate by chromosome
plink --bfile ${INPUT} --chr ${CHROM} --out ${OUTFILE} --recode vcf-iid

#compress and index VCF
bgzip ${OUTFILE}.vcf
tabix -p vcf ${OUTFILE}.vcf.gz

#Phase VCF
shapeit4 -I ${OUTFILE}.vcf.gz -M ${MAPFILE} --region ${CHROM} -O ${OUTFILE}.phz.bcf --log ${OUTFILE}.phz.log -T ${THREADS} 

#Index the BCF file
bcftools index ${OUTFILE}.phz.bcf
