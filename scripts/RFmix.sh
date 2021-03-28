#!/bin/bash
#PBS -l mem=32gb,nodes=1:ppn=6,walltime=12:00:00
#PBS -A spectorl
#PBS -m abe
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi

#Load modules
module load bcftools/1.9

#Input info
#Path to directory containing phased BCF files
INDIR="/panfs/roc/scratch/pmonnaha/ALL/"
#Prefix of phased BCF files.  Ensure that task array index is appropriately referencing chromosome number
PREFIX="all_aric_merged.chr${PBS_ARRAYID}.phz.bcf"
#Location VCF file containing the genotype information for the reference populations
REFDIR="/home/pankrat2/public/bin/ref/1000G_v5/"
#Again ensure that task array index is appropriately referencing chromosome number
REFBASE="ALL.chr${PBS_ARRAYID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
#This file links the sample IDs in the reference VCF to populations.
POPFILE="/home/pankrat2/public/bin/ref/1000G_v5/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.superpopulations.txt"

#Location of genetic map files
MAPFILE="/home/spectorl/pmonnaha/accessory/genetic_map_hg19.txt"

#Number of threads used for ancestry inference.  Set equal to ppn requested above.
THREADS="6"

#Output info
OUTDIR="/panfs/roc/scratch/pmonnaha/ALL/rfmix/"

# make output directory if it doesn't already exist
if ! [ -e "$OUTDIR" ]; then
  mkdir -p $OUTDIR
fi

OUTPRE=${PREFIX%.phz.bcf}

INFILE=${INDIR}/${PREFIX}
REFFILE=${REFDIR}/${REFBASE}
OUTFILE=${OUTDIR}/${OUTPRE}

#Run RFMix
rfmix -f ${INFILE} -r ${REFFILE} -m ${POPFILE} -g ${MAPFILE} -o ${OUTFILE} --chromosome=${PBS_ARRAYID} --n-threads=${THREADS} -G 12 --reanalyze-reference -e 10
