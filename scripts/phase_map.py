"""
Author: Patrick Monnahan
Purpose: Convert the genetic map format for RFMix into a format usable by Shapeit4.  Also, parse map by chromosome.
Requirements:
 Takes the following arguments:
    -i: input map
    -o: output prefix.
Date:
"""

import argparse

def phase_map(input_file, output_prefix):
	#The genetic map file for shapeit4 has different column ordering than for rfmix. 
	old_chrom = "-9"
	with open(input_file, 'r') as mapfile: #Loop over lines in input genetic map
		for line in mapfile:
			line = line.split()
			curr_chrom = line[0]
			if curr_chrom != old_chrom: #Check if we have moved on to next chromosome.  If so, close current file and open next one.
				if old_chrom != "-9": newmap.close()
				newmap = open(output_prefix + curr_chrom, 'w')
				old_chrom = curr_chrom
			if len(line) == 3:
				newmap.writelines(f"{line[1]}\t{line[0]}\t{line[2]}\n")
	newmap.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-i', type=str, metavar='input_map', required=True, help='')
	parser.add_argument('-o', type=str, metavar='output_prefix', default="accessory/Shapeit4_genetic_map")
	args = parser.parse_args()

	phase_map(args.i, args.o)