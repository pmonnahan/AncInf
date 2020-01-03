import argparse

def phase_map(input_file, output_prefix):
	#The genetic map file for shapeit4 has different column ordering than for rfmix. 
	old_chrom = "-9"
	with open(input_file, 'r') as mapfile:
		for line in mapfile:
			line = line.split()
			curr_chrom = line[0]
			if curr_chrom != old_chrom:
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