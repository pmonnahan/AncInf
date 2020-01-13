#takes in rfmix output and writes ancestry plink tped files for admixture mapping, etc
import argparse
from math import ceil
import re
import pdb
import subprocess
import os
from itertools import izip_longest

def grouper(n, iterable, fillvalue=None, start_idx = 4):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(map(float, iterable[start_idx:]))] * n
    return izip_longest(fillvalue=fillvalue, *args)

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
    chrs = map(str, range(1, 23))
    inds = open(args.fam)
    in_fam = []
    ped_key = open(args.outdir + "ped_key.txt", 'w')
    for ind in inds:
        in_fam.append(ind.strip().split())
    
    #go through every chromosome and population, printing tped along the way. counting number of tracts of each ancestry: G = absent, A = present
    for ii, i in enumerate(chrs):
        ind_list, pop_labels = getMSPinfo(args.rfmix + ".chr" + str(i) + ".msp.tsv")
        npop = len(pop_labels)

        if ii == 0:
            for pop in pop_labels:
                with open(args.outdir + args.out + '_' + pop + '.tfam', 'w') as out_fam:
                    # Sort information in in_fam based on order in ind_list
                    idxs = [ind_list.index(x[0]) for x in in_fam]
                    new_list = [x for _, x in sorted(zip(idxs, in_fam), key=lambda pair: pair[0])]
                    for ind in new_list:
                        line = ' '.join(ind) + '\n'
                        out_fam.writelines(line)

            out_tped = []
            for pop in pop_labels:
                out_tped.append(open(args.outdir + args.out + '_' + pop + '.tped', 'w'))
                key_line = pop + "\t" + args.outdir + args.out + '_' + pop + '.ped\n'
                ped_key.write(key_line)
            ped_key.close()

        if args.snps:
            with open(args.rfmix + ".chr" + str(i) + ".fb.tsv") as fbk:
                for k, line in enumerate(fbk):
                    if line.startswith("#") or line.startswith("chrom"):
                        continue
                    myLine = line.strip().split()
                    pos = myLine[1]
                    gpos = myLine[2]
                    fbk_max = []
                    for group in grouper(npop, myLine):
                        # pdb.set_trace()
                        idx = group.index(max(group))
                        if max(group) > args.fbk_threshold:
                            fbk_max.append(idx)
                        else:
                            fbk_max.append(-9)
                    for out in out_tped:
                        snp_id = str(i) + "-" + str(k-1)
                        out.write(' '.join(map(str, [i, snp_id, gpos, pos])) + ' ')
                    for j in range(len(fbk_max) / 2):
                        for pop in range(len(pop_labels)):
                            current_anc = [fbk_max[2 * j], fbk_max[2 * j + 1]]
                            if -9 in current_anc:
                                out_tped[pop].write('0 0 ')
                            else:
                                pop_count = current_anc.count(pop)
                                if pop_count == 0:
                                    out_tped[pop].write('G G ')
                                elif pop_count == 1:
                                    out_tped[pop].write('G A ')
                                else:
                                    out_tped[pop].write('A A ')
                    for pop in range(len(pop_labels)):
                        out_tped[pop].write('\n')

                        # for j in fbk_max:
                        #     if j == -9:
                        #         out_line += "0 "
                        #     elif j == pop:
                        #         out_line += "A "
                        #     else:
                        #         out_line += "G "

        else:
            with open(args.rfmix + ".chr" + str(i) + ".msp.tsv") as rfmix:
                for k, line in enumerate(rfmix):
                    if line.startswith("#") or line.startswith("chrom"):
                        continue
                    myLine = line.strip().split()  # start index is hardcoded at 6
                    start, end, gstart, gend, num_snps = map(float, myLine[1:6])
                    myLine = myLine[6:]

                    for out in out_tped:
                        snp_id = str(i) + "-" + str(k-1)
                        out.write(' '.join(map(str, [i, snp_id, (gstart + gend) / 2, int(ceil((start + end) / 2))])) + ' ')
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
            Name = args.outdir + args.out + '_' + pop_labels[pop]
            cmd = "plink --tfile " + Name + " --recode --out " + Name
            print(cmd)
            pp = subprocess.call(cmd, shell=True)


    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rfmix')
    parser.add_argument('--fam')
    parser.add_argument('--out')
    parser.add_argument('--outdir', default = os.getcwd() + '/plink/')
    parser.add_argument('--snps', help='if this flag is set, the fb file (containing SNPs) will be used instead of '
                                       'the msp file, the "position" of which is taken as the midpoint of the '
                                       'window', default=False, action="store_true")
    parser.add_argument('--transpose', help='if this flag is set, the tped files will be transposed to ped format', default=False, action="store_true")
    parser.add_argument('--fbk_threshold', type=float, default=0.9)
    args = parser.parse_args()

    main(args)

