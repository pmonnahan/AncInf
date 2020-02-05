#!/bin/bash

SCRIPT=`basename $0`
# Ensure that arguments are passed to script, if not display help
if [ "$#" -ne 7 ]; then
cat << EOF
Usage: sh ${SCRIPT} <PlinkPrefix> <OutDir> <ChromKey> <DownloadRef> <DataPrepStep1> <DataPrepStep2> <DataPrepStep3> <SaveDataPrepIntermeds>

ChromKey is text file specifying how to rename chromosomes
Subsequent arguments are all t/f and should likely stay set to t

EOF
  exit 1
fi

RawData="$1"
WorkingDir="$2"
ChromKey="$3"
DownloadRef="$4"
DataPrepStep1="$4"
DataPrepStep2="$5"
DataPrepStep3="$6"
SaveDataPrepIntermeds="$7"

BASE=$(basename RawData)

# =================
## IMPORTANT NOTES:
# =================

# The bed/bim/fam trio must have the proper variant ID's as identified by NCBI (otherwise fixing the data to the reference will likely not work)
# You also need to make sure that you have the proper reference build in relation to your genetic data you are trying to fix (don't try and fix GRCh37 data to a GRCh38 reference)

# =================
## DEPENDENCIES:
# =================

# BCFtools v1.8 or later and the BCFtools plugin +fixref
# htslib v1.8 or later -- which is a BCFTools dependency
module load bcftools/1.9
module load htslib/1.9
module load plink/1.90b6.10

bcftools_Exec="bcftools"
Plink_Exec="plink"
bgzip_Exec="bgzip"
gunzip_Exec="gunzip"

# Set Working Directory
# -------------------------------------------------
echo
echo Changing to Working Directory
echo ----------------------------------------------
mkdir -p ${WorkingDir}
echo ${WorkingDir}

	cd ${WorkingDir}


# Controls whether BCFTools +Fixref is performed on the dataset
echo "Performing BCFTools +Fixref on dataset in ${RawData}"
  echo ----------------------------------------------

#Make Temp Directory in which all Temp files will populate
  mkdir -p ./TEMP
  mkdir -p ./RefAnnotationData/

# Download all the Reference Data to Reformat the files
# ----------------------------------------------------------------------------

if [ "${DownloadRef}" == "t" ]; then

  echo
  echo Downloading Reference Data and index files from 1K Genomes and NCBI
  echo ----------------------------------------------

  echo Downloading: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
  echo
  wget --directory-prefix=./RefAnnotationData/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
  wget --directory-prefix=./RefAnnotationData/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai

  # Unzip Fasta file
  ${gunzip_Exec} -d ./RefAnnotationData/human_g1k_v37.fasta.gz

  # Rezip Fasta File in bgzip
  ${bgzip_Exec} ./RefAnnotationData/human_g1k_v37.fasta
  rm ./RefAnnotationData/human_g1k_v37.fasta


# Download the annotation files (make sure the the build version is correct) to flip/fix the alleles
  echo
  echo Downloading ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
  echo
  wget --directory-prefix=./RefAnnotationData/ ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
  wget --directory-prefix=./RefAnnotationData/ ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz.tbi


fi



# STEP 1: Convert the Plink BED/BIM/FAM into a VCF into a BCF so that it may be fixed with BCFtools
# --------------------------------------------------------------------------------------------------

if [ "${DataPrepStep1}" == "t" ]; then

  # Convert Plink file into a VCF

  printf "\n\nConverting $RawData Plink files into VCF format \n"
  echo ----------------------------------------------
  echo
  echo

  ${Plink_Exec} --bfile $RawData --recode vcf --out ./TEMP/DataFixStep1_${BASE}

# Convert from a VCF into a BCF and also rename the chromosomes to match the reference fasta (where [chr]23 is X, 24 is Y, etc.)

  printf "\n\nConverting VCF into a BCF with chromosome names that match the reference .fasta annotation \n\nNOTE: You may need to manually adjust ./Odyssey/0_DataPrepModule/RefAnnotationData/PlinkChrRename.txt depending on the fasta reference you use in order to match the chromosome names \n"
  echo ----------------------------------------------
  echo
  echo

  ${bcftools_Exec} annotate -Ob --rename-chrs ${ChromKey} ./TEMP/DataFixStep1_${BASE}.vcf > ./TEMP/DataFixStep1_${BASE}.bcf

fi

# STEP 2: Align Input File to the Reference Annotation (Fix with BCFtools)
# --------------------------------------------------------------------------------------------------

if [ "${DataPrepStep2}" == "t" ]; then

# Run bcftools +fixref to see the number of wrong SNPs
  printf "\n\nRun bcftools +fixref to first view the number of correctly annotated/aligned variants to the Reference annotation \n"
  echo ----------------------------------------------
  echo
  echo

  ${bcftools_Exec} +fixref ./TEMP/DataFixStep1_${BASE}.bcf -- -f ./0_DataPrepModule/RefAnnotationData/human_g1k_v37.fasta.gz

# Run bcftools to fix/swap the allels based on the downloaded annotation file
  printf "\n\nRun bcftools +fixref to fix allels based on the downloaded annotation file \n"
  echo ----------------------------------------------
  echo
  echo

  ${bcftools_Exec} +fixref ./TEMP/DataFixStep1_${BASE}.bcf -Ob -o ./TEMP/DataFixStep2_${BASE}-RefFixed.bcf -- -d -f ./RefAnnotationData/human_g1k_v37.fasta.gz -i ./RefAnnotationData/All_20170710.vcf.gz

# Rerun the bcftool +fixref check to see if the file has been fixed and all unmatched alleles have been dropped
  printf "\n\nRun bcftools +fixref to see if the file has been fixed - all alleles are matched and all unmatched alleles have been dropped \n"
  echo ----------------------------------------------
  echo
  echo

  ${bcftools_Exec} +fixref ./TEMP/DataFixStep2_${BASE}-RefFixed.bcf -- -f ./RefAnnotationData/human_g1k_v37.fasta.gz

fi


# STEP 3: Sort the Ref-Aligned BCF output and convert back into Plink format for Odyssey Pipeline
# --------------------------------------------------------------------------------------------------

if [ "${DataPrepStep3}" == "t" ]; then


# Sort the BCF output
  printf "\n\nSorting the BCF output since fixing it may have made it unsorted \n"
  echo ----------------------------------------------

  (${bcftools_Exec} view -h ./TEMP/DataFixStep2_${BASE}-RefFixed.bcf; ${bcftools_Exec} view -H ./TEMP/DataFixStep2_${BASE}-RefFixed.bcf | sort -k1,1d -k2,2n;) | ${bcftools_Exec} view -Ob -o ./TEMP/DataFixStep3_${BASE}-RefFixedSorted.bcf

  printf "Done \n\n\n"

# Convert BCF back into Plink .bed/.bim/.fam for Shapeit2 Phasing
  printf "\n\nConverting Fixed and Sorted BCF back into Plink .bed/.bim/.fam \n"
  echo ----------------------------------------------
  echo
  echo

  ${Plink_Exec} --bcf ./TEMP/DataFixStep3_${BASE}-RefFixedSorted.bcf --make-bed --out ./TEMP/DataFixStep3_${BASE}-RefFixSorted


# Finally Remove any positional duplicates
  # i.e. same position and alleles, but differently named variants since Shapeit will not tolerate these


  printf "\n\nFinding Positional and Allelic Duplicates \n"
  echo ----------------------------------------------
  echo
  echo

  ${Plink_Exec} --bfile ./TEMP/DataFixStep3_${BASE}-RefFixSorted --list-duplicate-vars ids-only suppress-first --out ./TEMP/Dups2Remove

  # Report Number of duplicates:
  DuplicateNumber="$(wc ./TEMP/Dups2Remove.dupvar | awk '{print $1}')"

  printf "\n\nRemoving Positional and Allelic Duplicates if they exist\nFound ${DuplicateNumber} Duplicate Variant/s\n"
  echo ----------------------------------------------
  echo
  echo

  ${Plink_Exec} --bfile ./TEMP/DataFixStep3_${BASE}-RefFixSorted --exclude ./TEMP/Dups2Remove.dupvar --make-bed --out ./TEMP/DataFixStep4_${BASE}-RefFixSortedNoDups




# Add back in the sex information
  printf "\n\nRestoring Sample Sex Information \n"
  echo ----------------------------------------------
  echo
  echo

  ${Plink_Exec} --bfile ./TEMP/DataFixStep4_${BASE}-RefFixSortedNoDups --update-sex ${RawData}.fam 3 --make-bed --out ${BASE}-FixRef


  echo
  echo
  echo ----------------------------------------------
  printf "Analysis Ready Data -- ${BASE}-FixRef \n"
  echo ----------------------------------------------



fi

# After Step: Cleanup File Intermediates
# --------------------------------------------------------------------------------------------------

if [ "${SaveDataPrepIntermeds}" == "f" ]; then

  echo
  echo ----------------------------------------------
  echo Tidying Up -- Cleanup Intermediate Files
  echo ----------------------------------------------

  if [ -d "./TEMP" ]; then rm -r ./TEMP; fi


fi




