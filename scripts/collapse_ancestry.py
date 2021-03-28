__author__ = 'armartin'

# Taken from https://github.com/armartin/ancestry_pipeline/blob/master/collapse_ancestry.py
# takes in RFMix output and creates collapsed 2 haploid bed files per individual

from datetime import datetime
import time
import argparse
from itertools import izip_longest
import gzip
import re
import pdb
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rfmix', help='Prefix for collection of rfmix output files, chr expected in filename',
                        required=True)
    parser.add_argument('--fbk_threshold', type=float, default=0.9)
    parser.add_argument('--ind', help='Individual ID, must match an individual in the rfmix Viterbi output',
                        required=True)
    parser.add_argument('--chrX', help='include chrX?', default=False, action="store_true")
    parser.add_argument('--out', help='prefix to bed file, _A.bed and _B.bed will be appended', required=True)

    args = parser.parse_args()
    return (args)


def grouper(n, iterable, fillvalue=None, start_idx = 4):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(map(float, iterable[start_idx:]))] * n
    return izip_longest(fillvalue=fillvalue, *args)


def check_gt_posterior(fbk_max, fbk_threshold, index, add, hap_anc, line, current_anc):
    if fbk_max[index * 2 + add] >= fbk_threshold:
        hap_anc = myLine[index * 2 + add]
    else:
        hap_anc = -9
        current_anc[add] = -9
    return current_anc, hap_anc


def checkFiles(rfmix_prefix, chrom):
    found = 0
    try:
        rfmix = open(rfmix_prefix + ".chr" + str(chrom) + ".msp.tsv")
        found += 1
        rfmix.close()
    except IOError:
        print("Did not find Viterbi output (.msp file) from RFmix for chromosome " + str(chrom))
    try:
        rfmix = open(rfmix_prefix + ".chr" + str(chrom) + ".sis.tsv")
        found += 1
        rfmix.close()
    except IOError:
        print("Did not find StayInState file (.sis file) from RFmix for chromosome " + str(chrom))
    try:
        rfmix = open(rfmix_prefix + ".chr" + str(chrom) + ".fb.tsv")
        found += 1
        rfmix.close()
    except IOError:
        print("Did not find ForwardBackward file (.fb file) from RFmix for chromosome " + str(chrom))

    return found


def find_haplotype_bounds(index, add, chrom, pop_order, hap):
    print("Starting chromosome " + str(chrom))
    npop = len(pop_order)
    msp = open(args.rfmix + ".chr" + str(chrom) + ".msp.tsv")
    fbk = open(args.rfmix + ".chr" + str(chrom) + ".fb.tsv")
    # snp_map = open(snp_file.replace('snp_locations', 'map')) #map of physical position -> genetic position

    last_anc_pos_cm = [None, None, 0]

    # Skip header lines in sis and fb files
    my_reg = msp.readline().strip().split()
    my_reg = msp.readline().strip().split()
    assert my_reg[0] == "#chm", "Error parsing header of msp file"

    # Get first non-header line of msp file
    my_reg = msp.readline().strip().split()
    start, end, gstart, gend = my_reg[1:5]
    my_reg = my_reg[6:]
    pop_label = my_reg[index * 2 + add]

    counter = 0
    curr_pos = 0
    for line in fbk:
        if line.startswith("#") or line.startswith("chrom"):
            continue

        counter += 1
        myLine = line.strip().split()  # start index is hardcoded at 6

        # pdb.set_trace()
        pos, gpos = myLine[1:3]

        # Check that SNP from fbk file resides within the current region given by msp file
        if int(start) <= int(pos) <= int(end):
            pass
        else:  # If not, increment msp file and HOPE that next window captures the SNP.
            my_reg = msp.readline().strip().split()
            start, end, gstart, gend = my_reg[1:5]
            assert int(start) < int(pos) < int(end), "window increment did not capture current SNP"
            my_reg = my_reg[6:]
            pop_label = my_reg[index * 2 + add]

        fbk_max = []
        for i in grouper(npop, myLine):
            fbk_max.append(max(i))
        # pdb.set_trace()
        if fbk_max[index * 2 + add] < args.fbk_threshold:
            # pdb.set_trace()
            pop_label = -9
        else:
            pop_label = my_reg[index * 2 + add]
        # fencepost for start of the chromosome
        if counter == 1:
            last_anc_pos_cm = [pop_label, pos, gpos]
            post_anc_pos_cm = [pop_label, pos, gpos]
            continue

        # start regular iterations
        current_anc_pos_cm = [pop_label, pos, gpos]
        if current_anc_pos_cm[0] == last_anc_pos_cm[0]:
            last_anc_pos_cm = current_anc_pos_cm
            continue
        else:
            # we've reached the end of a region. Need to print.
            if last_anc_pos_cm[0] == -9:
                hap.write(str(chrom) + '\t' + post_anc_pos_cm[1] + '\t' + last_anc_pos_cm[1] +
                          '\tUNK\t' + post_anc_pos_cm[2] + '\t' + last_anc_pos_cm[2] + '\n')
            else:
                hap.write(str(chrom) + '\t' + post_anc_pos_cm[1] + '\t' + last_anc_pos_cm[1] + '\t' +
                          pop_order[int(last_anc_pos_cm[0])] + '\t' +
                          post_anc_pos_cm[2] + '\t' + last_anc_pos_cm[2] + '\n')
            post_anc_pos_cm = current_anc_pos_cm

        last_anc_pos_cm = current_anc_pos_cm

    # last iteration, still need to print
    if last_anc_pos_cm[0] == -9:
        hap.write(str(chrom) + '\t' + post_anc_pos_cm[1] + '\t' + current_anc_pos_cm[1] +
                  '\tUNK\t' + post_anc_pos_cm[2] + '\t' + current_anc_pos_cm[2] + '\n')
    else:
        hap.write(str(chrom) + '\t' + post_anc_pos_cm[1] + '\t' + current_anc_pos_cm[1] +
                  '\t' + pop_order[int(current_anc_pos_cm[0])] + '\t' +
                  post_anc_pos_cm[2] + '\t' + current_anc_pos_cm[2] + '\n')
    print("Finished chromosome " + str(chrom))


def main(individual, chrs):
    # open bed files (2 haplotypes per individual)
    hap_a = open(args.out + '_A.bed', 'w')
    hap_b = open(args.out + '_B.bed', 'w')
    for chrom in chrs:
        found_files = checkFiles(args.rfmix, chrom)
        if found_files != 3:
            print(
                "Did not find all of the necessary input files using the provided rfmix output prefix for chromosome " + str(
                    chrom))
        else:
            index, pop_labels = getMSPinfo(args.rfmix + ".chr" + str(chrom) + ".msp.tsv", individual)
            if not os.path.exists(os.path.dirname(args.out) + "/pop_order.txt"):
                with open(os.path.dirname(args.out) + "pop_order.txt", 'w') as pop_order_file:
                    pop_order_file.write(",".join(pop_labels))
            find_haplotype_bounds(index, 0, chrom, pop_labels, hap_a)
            find_haplotype_bounds(index, 1, chrom, pop_labels, hap_b)

    hap_a.close()
    hap_b.close()


def getMSPinfo(msp_file, current_ind):
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
                        ind = ind.strip(".0")
                        ind_list.append(ind)
            else:
                break
    try:
        ind_index = ind_list.index(current_ind)
    except ValueError:
        raise Exception('Individual is not in the list provided')

    return ind_index, pop_labels


if __name__ == '__main__':
    # load parameters and files
    args = parse_args()

    if args.ind is None:
        raise Exception('individual not set')
    current_ind = args.ind

    # set up chromosome variables
    chrs = range(22, 23)
    if args.chrX:
        chrs.append('X')

    print('Starting [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']')
    main(current_ind, chrs)
    print('Finished [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']')
