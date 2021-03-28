#!/bin/bash
#PBS -l mem=64gb,nodes=1:ppn=12,walltime=24:00:00
#PBS -A spectorl
#PBS -m a
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi

#Load modules
module load plink/1.90b6.10
module load htslib/1.6
module load bcftools/1.9

#Input info
#Path to directory containing master PLINK files
INDIR="/home/spectorl/pmonnaha/data/QC_allRaces"
#Prefix of master PLINK files containing genotype data for all individuals and all chromosomes
BFILE="all_aric_merged"
#Location of genetic map files
MAPDIR="/home/pankrat2/public/resources/Genome/hg19/chr/"
MAPBASE="genetic_map_hg19_chr${PBS_ARRAYID}.txt.gz"
#Number of threads used for phasing.  12 may be overkill.  Phasing completes in < 30 minutes.
THREADS="12"

#Output info
OUTDIR="/panfs/roc/scratch/pmonnaha/ALL/"
OUTFILE=${BFILE}.chr${PBS_ARRAYID}

MAPFILE=${MAPDIR}/${MAPBASE}

#Convert plink files to vcf and separate by chromosome
plink --bfile ${INDIR}/${BFILE} --chr ${PBS_ARRAYID} --out ${OUTDIR}/${OUTFILE} --recode vcf-iid

#compress and index VCF
bgzip ${OUTDIR}/${OUTFILE}.vcf
tabix -p vcf ${OUTDIR}/${OUTFILE}.vcf.gz

#Phase VCF
shapeit4 -I ${OUTDIR}/${OUTFILE}.vcf.gz -M ${MAPFILE} --region ${PBS_ARRAYID} -O ${OUTDIR}/${OUTFILE}.phz.bcf --log ${OUTDIR}/${OUTFILE}.phz.log -T ${THREADS} 

#Index the BCF file
bcftools index ${OUTDIR}/${OUTFILE}.phz.bcf
