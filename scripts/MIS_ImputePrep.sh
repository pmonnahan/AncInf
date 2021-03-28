#!/bin/bash

SCRIPT=`basename $0`
# Ensure that arguments are passed to script, if not display help
if [ "$#" -ne 3 ]; then
cat << EOF
Usage: sh ${SCRIPT} <plink_prefix> <working_directory> <t/f delete_ref>

plink_prefix: full path to PREFIX of plink bed file
working_directory: Where are work will be performed and stored
delete_ref: t for true f for false

EOF
  exit 1
fi

PLINK="$1"
WorkingDir="$2"
Delete_ref="$3"
BASE=$(basename "${PLINK}")

module load plink/1.90b6.10
module load samtools

mkdir -p ${WorkingDir}
cd ${WorkingDir}

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai

wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
wget http://qbrc.swmed.edu/zhanxw/software/checkVCF/checkVCF-20140116.tar.gz

unzip HRC-1000G-check-bim-v4.2.7.zip
rm HRC-1000G-check-bim-v4.2.7.zip
gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
tar -xzvf checkVCF-20140116.tar.gz
rm checkVCF-20140116.tar.gz
gunzip human_g1k_v37.fasta.gz
samtools faidx human_g1k_v37.fasta


plink --bfile ${PLINK} --freq --make-bed --out ${BASE}

perl HRC-1000G-check-bim.pl -b "${BASE}.bim" -f "${BASE}.frq" -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
sh Run-plink.sh

mkdir ./impute_rdy_vcfs
cd ./impute_rdy_vcfs

for file in ../*-chr*.bim; do
  prefix=$(basename "${file}" | rev | cut -d "." -f 2- | rev);
  plink --bfile "../${prefix}" --recode vcf --out ${prefix};
  bcftools sort "${prefix}.vcf" -Oz -o "${prefix}.vcf.gz";
  python2.7 ../checkVCF.py -r ../human_g1k_v37.fasta -o "${prefix}.vcf.gz.out" "${prefix}.vcf.gz";
  rm "${prefix}.vcf"
  done

cd ..
mkdir plink_files
mkdir RefDat
mv *bim plink_files; mv *fam plink_files; mv *log plink_files; mv *bed plink_files; mv *.txt plink_files
mv *.tab RefDat; mv *.fasta RefDat; mv *.fai RefDat; mv *.fa RefDat
rm example.vcf.gz

if [ ${Delete_ref} = 't' ]; then
  rm -r -f RefDat
fi