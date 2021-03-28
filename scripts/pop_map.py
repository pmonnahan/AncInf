"""
Author: Patrick Monnahan
Purpose: Link population IDs to superpopulations and create a file that lists each individual along with their (super) population name for RFMix.
Requirements:
 Takes the following arguments:
    -i: input map.  File containing list of individuals and their population IDs.  One sample per line with first column being sample name and second column being pouplation ID
    -p: pop_ids. Comma-separated list (no-spaces) of the population name/IDs as specified in input_map.
    -n: pop_names.  Comma-separated string containing the names to be used for each pop_id in output.  These are typically the super-population names that group together the various populations in the pop_ids.
    -o: output filename.
Date:
"""

import argparse
from shutil import copyfile

def pop_map(input_file, pop_ids, pop_names, output):
	POPS = dict(zip(pop_ids, pop_names)) #Create dictionary where each key is a population name and each value is the superpopulation name.
	with open(output, 'w') as new_file:
		with open(input_file, mode='rU') as pop_file: #Open and read the input_map.
			for line in pop_file:
				if len(line.split()) == 2: #Ensure that each line has the expected 2 elements
					ind = line.split()[0]
					pop = line.split()[1]
					if pop in POPS.keys(): #Check if population listed for this individual in the input_map is present in the dictionary created above.
						new_pop = POPS[pop] #Retrieve superpopulation name for the population ID for this individual.
						new_file.writelines(f"{ind}\t{new_pop}\n")
	return()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-i', type=str, metavar='input_map', required=True, help='File containing list of individuals and their pop_id')
	parser.add_argument('-p', type=str, metavar='pop_ids', required=True, help='Comma-separated list (no-spaces) of the population name/IDs as specified in input_map')
	parser.add_argument('-n', type=str, metavar='pop_names', required=True, help='comma separated string containing the names to be used for each pop_id in output')
	parser.add_argument('-o', type=str, metavar='outfile', default="Population_Map_File.txt", help='Output file name')
	args = parser.parse_args()

	pop_ids = args.p.split(",")
	pop_names = args.n.split(",")

	assert len(pop_ids) == len(pop_names), print("A name should be provided for every population ID requested")

	if len(pop_ids) == 1:
		if pop_ids[0] == "all": copyfile(args.i, output)
		else: print("Only one reference population provided")
	else:
		pop_map(args.i, pop_ids, pop_names, args.o)