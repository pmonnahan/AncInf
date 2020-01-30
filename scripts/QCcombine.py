#!/usr/bin/env python
"""
Author: Patrick Monnahan
Purpose: 
Requirements: 
 Takes the following arguments:
    -i: 
    -o:
Date:
"""

# Import Modules
import subprocess
import argparse
import os
import pdb
import pandas as pd
import shutil
from functools import reduce

# Define Functions
def var_miss(bfile, threshold, output, plink):
    cmd = f"{plink} --bfile {bfile} --geno {threshold} --make-bed --out {output}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(out1, err1)

def geno_miss(bfile, threshold, output, plink):
    cmd = f"{plink} --bfile {bfile} --mind {threshold} --make-bed --out {output}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(out1, err1)

def geno_flt(bfile, keptSamples, output, plink):
    cmd = f"{plink} --bfile {bfile} --keep {keptSamples} --make-bed --out {output}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(out1, err1)

def calc_stats(bfile, plink):
    '''calculate frequency, missingness, HWE prob, and missingness by sex'''

    cmd = f"{plink} --bfile {bfile} --freq --missing --hardy --mpheno 3 --pheno {bfile}.fam --test-missing --out {bfile}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish

    return(out1, err1)

def wrapQC(input, output, tvm1, tgm, tvm2, maf, hwe, mbs, plink):

    # TODO: STILL NEED TO DO REFALT SWITCHING...maybe this is best left to user to do prior to running pipeline.  Each array needs its own key to update.  Hard to generalize without requireing each input file to have associated key
    # TODO: parallelize using pool
    out1, err1 = var_miss(input, tvm1, f"{output}.var_miss{tvm1}", plink)
    out2, err2 = geno_miss(f"{output}.var_miss{tvm1}", tgm, f"{output}.geno_miss{tgm}", plink)
    out3, err3 = geno_flt(input, f"{output}.geno_miss{tgm}.fam",
                         f"{output}.geno_flt{tgm}", plink)

    out2, err2 = calc_stats(f"{output}.geno_flt{tgm}", plink)
    DF_list = []
    for stat in ['frq', 'hwe', 'lmiss', 'missing']:
        DF_list.append(pd.read_csv(f"{output}.geno_flt{tgm}.{stat}", sep = r"\s+"))
    df = reduce(lambda df1, df2: pd.merge(df1, df2, on='SNP'), DF_list)
    df2 = df.query(f"MAF < {maf} | P_x < {hwe} | F_MISS > {tvm2} | P_y < {mbs}")
    df2['SNP'].to_csv(f"{output}.fltdSNPs.txt", header=False)
    cmd = f"{plink} --bfile {output}.geno_flt{tgm} --exclude {output}.fltdSNPs.txt --make-bed --out {output}.fin"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish

    return(out2, err2)

# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='comma-separated string (no space) of the PLINK files to be QCed and combined')
    parser.add_argument('-o', type=str, metavar='output', required=True, help='output name for combined PLINK file')
    parser.add_argument('-d', type=str, metavar='output_directory', required=True, help='directory to store all output')
    parser.add_argument('-t', type=str, metavar='threads', default=1, help='Number of threads to use')
    parser.add_argument('-p', type=str, metavar='plink_path', default="plink", help='')
    parser.add_argument('-tvm1', type=str, metavar='treshold_variant_missingness1', default=0.2, help='Used for generated filtered dataset for calculating genotype missingness')
    parser.add_argument('-tgm', type=str, metavar='treshold_genotype_missingness', default=0.1, help='')
    parser.add_argument('-tvm2', type=str, metavar='treshold_variant_missingness2', default=0.02, help='Missingness allowed for variants post genotype filtration')
    parser.add_argument('-maf', type=str, metavar='minor_allele_frequency', default=0.01,
                        help='')
    parser.add_argument('-mbs', type=str, metavar='missingness_by_sex', default=0.00001,
                        help='')
    parser.add_argument('-hwe', type=str, metavar='p_value_cutoff_HWEtest', default=0.00001,
                        help='')
    args = parser.parse_args()

    files = args.i.split(",")

    if not os.path.exists(args.d):
        os.mkdir(args.d)

    with open(f"{args.d}/files2merge.txt", 'w') as outfiles:
        for file in files:
            base = os.path.basename(file)
            dirname = args.d + "/" + base
            if not os.path.exists(dirname): os.mkdir(dirname)
            out, err = wrapQC(file, f"{dirname}/{base}", args.tvm1, args.tgm, args.tvm2, args.maf, args.hwe, args.mbs, args.p)
            outfiles.write(f"{dirname}/{base}.fin\n")

    if len(files) != 1:
        DF_list = []
        with open(f"{args.d}/files2merge.txt", 'r') as files:
            for file in files:
                DF_list.append(pd.read_csv(f"{file.strip()}.bim", sep=r"\s+", header=None))
        pdb.set_trace()
        df = reduce(lambda df1, df2: pd.merge(df1, df2, on=1), DF_list)
        df['1'].to_csv(f"{args.d}/mergeSNPs.txt", header=False)
        cmd = f"{args.p} --merge-list {args.d}/files2merge.txt --extract {args.d}/mergeSNPs.txt --make-bed --out {args.d}/{args.o}"
        pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
        out1, err1 = pp1.communicate()  # Wait for it to finish
    else:
        cmd = f"{args.p} --bfile {files[0]} --make-bed --out {args.d}/{args.o}"
        pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
        out1, err1 = pp1.communicate()  # Wait for it to finish

