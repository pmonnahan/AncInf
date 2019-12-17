

import allel
import subprocess
import argparse
import os
import pdb



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script prepares the data contained in multiple SV VCFs (corresponding to the same samples called against different reference genomes) for subsequent analysis with SVMap.py (which links the SVs across reference genomes).\nRequirements: SUVIVOR_ant, vcftools, and bgzip/tabix')
    parser.add_argument('-f', type=str, metavar='vcf_info_file', required=True, help='tab delimited file with each line containing 1.) the reference genotype ID (e.g. B73; used for naming), 2.) a bed file with gene locations ONLY and 3.) the vcf file.')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-s', type=str, metavar='output_suffix', required=True, help='Output Suffix')
    parser.add_argument('-sp', type=str, metavar='SURVIVOR_ant_path', required=False, default='/Users/pmonnahan/Documents/Research/Maize/MaizeSV/software/SURVIVOR_ant/bin/survivor_ant-core-0.1.0/survivor_ant')
    parser.add_argument('-vp', type=str, metavar='vcftools_perl_folder', required=False, default='/usr/local/bin/')
    parser.add_argument('-ad', type=str, metavar='annotation_distance', required=False, default='0', help = "distance from SV to buffer for looking for gene overlap when annotating merged vcf with gene info; this is accomodated by -b, which provides more explicit info regarding the location where the overlap is found, so this can be left at 0 (default)")
    parser.add_argument('-b', type=int, metavar='buffer', required=False, default=2000, help = "Adds this amount to either side of gene boundary for annotation of variants")
    args = parser.parse_args()
