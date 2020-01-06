#takes in rfmix output and writes ancestry plink tped files for admixture mapping, etc
import argparse
from math import ceil
import re
import pdb
import subprocess

def getMSPinfo(msp_file):
    # reading order of individuals
    ind_list = []
    pop_labels = []
    with open(msp_file) as ind_info:
        for line in ind_info:
            if line.startswith("#Sub"):
                line = line.strip().split()
                for pop in line[2:]:
                    pop_labels.append(pop.split("=")[0])
            elif line.startswith("#chm"):
                line = line.strip().split()
                for ind in line[7:]:
                    if ind.endswith(".0"):
                        ind = ind[:-2]
                        ind_list.append(ind)
            else:
                break
    return ind_list, pop_labels


def main(args):
    chrs = map(str, range(22,23))
    inds = open(args.fam)
    in_fam = []
    for ind in inds:
        in_fam.append(ind.strip().split())
    
    #go through every chromosome and population, printing tped along the way. counting number of tracts of each ancestry: G = absent, A = present
    for ii, i in enumerate(chrs):
        ind_list, pop_labels = getMSPinfo(args.rfmix + ".chr" + str(i) + ".msp.tsv")

        if ii == 0:
            for pop in pop_labels:
                with open(args.out + '_' + pop + '.tfam', 'w') as out_fam:
                    # Sort information in in_fam based on order in ind_list
                    idxs = [ind_list.index(x[0]) for x in in_fam]
                    new_list = [x for _, x in sorted(zip(idxs, in_fam), key=lambda pair: pair[0])]
                    for ind in new_list:
                        line = ' '.join(ind) + '\n'
                        out_fam.writelines(line)

            out_tped = []
            for pop in pop_labels:
                out_tped.append(open(args.out + '_' + pop + '.tped', 'w'))


        if args.snps:
            # TODO: need to elaborate rfmix2plink to use fb files instead of msp files (i.e. output SNPs instead of window midpoints)
            pass
        else:
            with open(args.rfmix + ".chr" + str(i) + ".msp.tsv") as rfmix:
                for k, line in enumerate(rfmix):
                    if line.startswith("#") or line.startswith("chrom"):
                        continue
                    myLine = line.strip().split()  # start index is hardcoded at 6
                    start, end, gstart, gend, num_snps = map(float, myLine[1:6])
                    myLine = myLine[6:]
                    for out in out_tped:
                        out.write(' '.join(map(str, [i, k-1, (gstart + gend) / 2, int(ceil((start + end) / 2))])) + ' ')
                    for j in range(len(myLine)/2):
                        current_anc = [myLine[2*j], myLine[2*j+1]]
                        for pop in range(len(pop_labels)):
                            pop_count = current_anc.count(str(pop))
                            if pop_count == 0:
                                out_tped[pop].write('G G ')
                            elif pop_count == 1:
                                out_tped[pop].write('G A ')
                            else:
                                out_tped[pop].write('A A ')
                    for pop in range(len(pop_labels)):
                        out_tped[pop].write('\n')

    for pop in range(len(pop_labels)):
        out_tped[pop].close()
        if args.transpose:
            cmd = "plink --tfile " + out_tped[pop] + " --recode --out " + out_tped[pop].replace('.tped', '')
            pp = subprocess.call(cmd)


    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rfmix')
    parser.add_argument('--fam')
    parser.add_argument('--out')
    parser.add_argument('--snps', help='if this flag is set, the fb file (containing SNPs) will be used instead of '
                                       'the msp file, the "position" of which is taken as the midpoint of the '
                                       'window', default=False, action="store_true")
    parser.add_argument('--transpose', help='if this flag is set, the tped files will be transposed to ped format', default=False, action="store_true")
    args = parser.parse_args()
    main(args)