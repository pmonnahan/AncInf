#!/bin/bash

SCRIPT=`basename $0`
# Ensure that arguments are passed to script, if not display help
if [ "$#" -ne 1 ]; then
cat << EOF
Usage: sh ${SCRIPT} <WorkingDir>

EOF
  exit 1
fi

WorkingDir="$1"

mkdir ./Reference

# Retrieves the (default) Reference Genome from the IMPUTE Website
				# ----------------------------------------------------------------------------------
				# Collects the 1000Genome Reference Build from the Impute Site
					#(https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html)
					# Reference Build Specs: 1,000 Genomes haplotypes -- Phase 3 integrated variant set release in NCBI build 37 (hg19) coordinates
					# Ref Build Updated Aug 3 2015

					printf "\n\nRetrieving 1K Genome Phase 3 Ref Panel and hg19 Genetic Map from Impute2 Website \n-------------------------------------------------------------------------------\n\n\n"
						wget --directory-prefix=${WorkingDir}/Reference/ https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
						wget --directory-prefix=${WorkingDir}/Reference/ https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz


				#Unzip the packaged ref panel
					printf "\n\nUnpackaging Ref Panel \n--------------------------\n\n"
						tar -xzf ${WorkingDir}/Reference/1000GP_Phase3.tgz -C ${WorkingDir}Reference/
						tar -xzf ${WorkingDir}/Reference/1000GP_Phase3_chrX.tgz -C ${WorkingDir}Reference/

				# Since untar makes an additional directory, move all the files from the 1000GP_Phase3 folder and move it into the Ref Directory
					printf "\n\nCleaning Up \n-------------------\n\n"
						mv ${WorkingDir}/Reference/1000GP_Phase3/* ${WorkingDir}Reference/

				# Delete the now empty directory and the tgz zipped Ref panel
					rmdir ${WorkingDir}/Reference/1000GP_Phase3/
					rm ${WorkingDir}/Reference/*.tgz

        # Download BCF files
          wget --directory-prefix=${WorkingDir}/Reference/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/bcf_files/ALL.chr\*


