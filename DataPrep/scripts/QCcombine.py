
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
import itertools
import shutil
from functools import reduce
import multiprocessing

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

def var_flt(bfile, keptVariants, output, plink):
    cmd = f"{plink} --bfile {bfile} --extract {keptVariants} --make-bed --out {output}"
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

    # TODO: Add first round of QC, filtering for related samples.  Perhaps this should be separate script
    # TODO: Add last round of QC, filtering for related samples in the combined dataset
    out1, err1 = var_miss(input, tvm1, f"{output}.var_miss{tvm1}", plink)
    out2, err2 = geno_miss(f"{output}.var_miss{tvm1}", tgm, f"{output}.geno_miss{tgm}", plink)
    out3, err3 = geno_flt(input, f"{output}.geno_miss{tgm}.fam",
                         f"{output}.geno_flt{tgm}", plink)

    out2, err2 = calc_stats(f"{output}.geno_flt{tgm}", plink)
    DF_list = []
    for stat in ['frq', 'hwe', 'lmiss', 'missing']:
        tdf = pd.read_csv(f"{output}.geno_flt{tgm}.{stat}", sep = r"\s+")
        if stat=="hwe":
            tdf = tdf[tdf['TEST']=="UNAFF"]
        DF_list.append(tdf)
    df = reduce(lambda df1, df2: pd.merge(df1, df2, on='SNP'), DF_list)
    # ARIC is failing HWE massively.c
    df2 = df.query(f"MAF < {maf} | P_x < {hwe} | F_MISS > {tvm2} | P_y < {mbs}")
    df3 = df.query(f"MAF < {maf}")
    df4 = df.query(f"P_x < {hwe}")
    df5 = df.query(f"F_MISS > {tvm2}")
    df6 = df.query(f"P_y < {mbs}")
    df2['SNP'].to_csv(f"{output}.fltdSNPs.txt", header=False, index=False)
    df6.to_csv(f"{output}.fltdSNPs.mbs.txt", index=False)
    df3.to_csv(f"{output}.fltdSNPs.maf.txt", index=False)
    df4.to_csv(f"{output}.fltdSNPs.hwe.txt", index=False)
    df5.to_csv(f"{output}.fltdSNPs.miss.txt", index=False)
    cmd = f"{plink} --bfile {output}.geno_flt{tgm} --list-duplicate-vars suppress-first --out {output}; tail -n +2 {output}.dupvar | cut -f 4 >> {output}.fltdSNPs.txt;"
    cmd += f"cat {output}.geno_flt{tgm}.bim | " + "awk '{if($5 == \"N\" || $6 == \"N\") print $2}' >> " + output + ".fltdSNPs.txt; "
    cmd += f"cut -f 2 {output}.geno_flt{tgm}.bim | sort | uniq -d >> {output}.fltdSNPs.txt; {plink} --bfile {output}.geno_flt{tgm} --exclude {output}.fltdSNPs.txt --make-bed --out {output}.QC; sed 's/^24/23/' {output}.QC.bim | sed 's/^X/23/' > {output}.QC.bim.tmp; mv {output}.QC.bim.tmp {output}.QC.bim"
    print(cmd)
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish

    return(output + ".QC")


# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='input plink file')
    parser.add_argument('-o', type=str, metavar='output', required=True, help='output name for QCed PLINK file')
    parser.add_argument('-d', type=str, metavar='output_directory', required=True, help='directory to store all output')
    parser.add_argument('-p', type=str, metavar='plink_path', default="plink", help='')
    parser.add_argument('-tvm1', type=str, metavar='treshold_variant_missingness1', default=0.2, help='Used for generated filtered dataset for calculating genotype missingness')
    parser.add_argument('-tgm', type=str, metavar='treshold_genotype_missingness', default=0.1, help='')
    parser.add_argument('-tvm2', type=str, metavar='treshold_variant_missingness2', default=0.05, help='Missingness allowed for variants post genotype filtration')
    parser.add_argument('-maf', type=str, metavar='minor_allele_frequency', default=0.01,
                        help='')
    parser.add_argument('-mbs', type=str, metavar='missingness_by_sex', default=0.0000001,
                        help='')
    parser.add_argument('-hwe', type=str, metavar='p_value_cutoff_HWEtest', default=-1,
                        help='')
    parser.add_argument('-s', type=str, metavar='sample_file', default='all',
                        help='File containing sample IDs to retain from merged plink files. One sample per line')
    args = parser.parse_args()

    if not os.path.exists(args.d):
        os.mkdir(args.d)


    base = os.path.basename(args.i).strip(".bed")
    dirname = args.d + "/" + base
    if not os.path.exists(dirname): os.mkdir(dirname)
    wrapQC(args.i, f"{dirname}/{base}", args.tvm1, args.tgm, args.tvm2, args.maf, args.hwe, args.mbs, args.p)




