import argparse
from shutil import copyfile

def pop_map(input_file, pop_ids, pop_names, output):
	POPS = dict(zip(pop_ids, pop_names))
	with open(output, 'w') as new_file:
		with open(input_file, mode='rU') as pop_file:
			for line in pop_file:
				if len(line.split()) == 2:
					ind = line.split()[0]
					pop = line.split()[1]
					if pop in POPS.keys():
						new_pop = POPS[pop]
						new_file.writelines(f"{ind}\t{new_pop}\n")
	return()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-i', type=str, metavar='input_map', required=True, help='File containing list of individuals and their pop_id')
	parser.add_argument('-p', type=str, metavar='pop_ids', required=True, help='Population name/IDs as specified in input_map')
	parser.add_argument('-n', type=str, metavar='pop_names', default=-9, help='comma separated string containing the names to be used for each pop_id in output')
	parser.add_argument('-o', type=str, metavar='output', default="Population_Map_File.txt", help='Output file name')
	args = parser.parse_args()

	pop_ids = args.p.split(",")
	if args.n == -9:
		pop_names = ["POP" for x in pop_ids]
	else:
		pop_names = args.n.split(",")

	assert len(pop_ids) == len(pop_names), print("A name should be provided for every population ID requested")

	if len(pop_ids) == 1:
		if pop_ids[0] == "all": copyfile(args.i, output)
		else: print("Only one reference population provided")
	else:
		pop_map(args.i, pop_ids, pop_names, args.o)