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

    # TODO: STILL NEED TO DO REFALT SWITCHING...maybe this is best left to user to do prior to running pipeline.  Each array needs its own key to update.  Hard to generalize without requireing each input file to have associated key
    # TODO: Add filtering for related samples
    # TODO: parallelize using pool
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


def flipstrand_and_updateID(input, names_file, snps2flip, output, plink):
    cmd = f"{plink} --bfile {input} --update-name {names_file} --flip {snps2flip} --make-bed --out {output}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()
    return(output)


def var_flt(bfile, keptVariants, output, plink):
    cmd = f"{plink} --bfile {bfile} --extract {keptVariants} --make-bed --out {output}"
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(output)


def parallelQC(arg_list, threads, function="wrapQC"):
    pool = multiprocessing.Pool(threads)
    if function == "wrapQC":
        outfiles = pool.starmap(wrapQC, arg_list)
    if function == "flipstrand_and_updateID":
        outfiles = pool.starmap(flipstrand_and_updateID, arg_list)
    if function == "var_flt":
        outfiles = pool.starmap(var_flt, arg_list)
    pool.close()
    pool.join()
    return(outfiles)


def merge_files(file, file_number, output, samples='all', plink='plink'):
    pdb.set_trace()
    if file_number == 1:
        if samples == 'all': cmd = f"{plink} --bfile {file} --make-bed --allow-no-sex --out {output}"
        else: cmd = f"{plink} --bfile {file} --keep {samples} --make-bed --allow-no-sex --out {output}"
    else:
        if samples == 'all': cmd = f"{plink} --merge-list {file} --make-bed --allow-no-sex --out {output}"
        else: cmd = f"{plink} --merge-list {file} --keep {samples} --make-bed --allow-no-sex --out {output}"
    print(cmd)
    pp1 = subprocess.Popen(cmd, shell=True)  # Run cmd1
    out1, err1 = pp1.communicate()  # Wait for it to finish
    return(out1, err1)


def read_bims(file_key):
    DF_dict = {}
    with open(file_key, 'r') as files:
        for file in files:
            df = pd.read_csv(f"{file.strip()}.bim", sep=r"\s+", header=None, dtype=str)
            DF_dict[os.path.basename(file.strip())] = df
    return(DF_dict)


# Set up command line execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', type=str, metavar='input', required=True, help='comma-separated string (no space) of the PLINK files to be QCed and combined')
    parser.add_argument('-o', type=str, metavar='output', required=True, help='output name for combined PLINK file')
    parser.add_argument('-d', type=str, metavar='output_directory', required=True, help='directory to store all output')
    parser.add_argument('-t', type=int, metavar='threads', default=1, help='Number of threads to use')
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

    files = args.i.split(",")

    if not os.path.exists(args.d):
        os.mkdir(args.d)

    arg_list = []
    with open(f"{args.d}/files2merge.txt", 'w') as outfiles:
        for file in files:
            base = os.path.basename(file)
            dirname = args.d + "/" + base
            if not os.path.exists(dirname): os.mkdir(dirname)
            arg_list.append((file, f"{dirname}/{base}", args.tvm1, args.tgm, args.tvm2, args.maf, args.hwe, args.mbs, args.p))
            outfiles.write(f"{dirname}/{base}.QC\n")

    # Run QC on each input file
    outs = parallelQC(arg_list, args.t)

    if len(files) != 1:
        DF_dict = read_bims(f"{args.d}/files2merge.txt")

        # Identify merge issues before they occur.  SNPs with matching positions but different names.  Matching SNPs where strand is flipped. e.g. A-G in one is T-C in the other
        Names_dict = {}
        Flips_dict = {}
        for i, comb in enumerate(itertools.combinations(DF_dict.keys(), r=2)):
            if i==0:  # Hold first data frame as the reference and flip needed sites in subsequent data frames with respect to constant reference
                reference = comb[0]
            df1 = DF_dict[comb[0]]
            df2 = DF_dict[comb[1]]
            df = pd.merge(df1, df2, on=(0,3))
            # pdb.set_trace()
            # TODO: need to catch when names match, alleles match, but positions do not.
            bad_names = df[df['1_x'] != df['1_y']]
            bad_names_flips = bad_names[(bad_names['4_x'] != bad_names['4_y']) & (bad_names['4_x'] != bad_names['5_y']) & (bad_names['5_x'] != bad_names['4_y']) & (bad_names['5_x'] != bad_names['5_y'])]

            df_b = pd.merge(df1, df2, on=1)
            good_names_flips = df_b[(df_b['4_x'] != df_b['4_y']) & (df_b['4_x'] != df_b['5_y']) & (df_b['5_x'] != df_b['4_y']) & (df_b['5_x'] != df_b['5_y'])]

            # try:
            #     bad_names = bad_names[~bad_names['1_y'].isin(Names_dict[comb[1]]['1_y'])]
            #     # Names_dict[comb[1]] = pd.concat([Names_dict[comb[1]], bad_names[['1_y', '1_x']]])
            # except KeyError:
            #
            #     Names_dict[comb[1]] = bad_names[['1_y', '1_x']]
            if comb[0] == reference:
                Names_dict[comb[1]] = bad_names[['1_y', '1_x']]
                flips = list(good_names_flips[1]) + list(bad_names_flips['1_x'])  # SNP ids to be flipped should be specified in terms of reference ids bc
                Flips_dict[comb[1]] = flips

        # harmonize SNPids across data sets
        with open(f"{args.d}/files2premerge.txt", 'w') as mergefile:
            mergefile.write(f"{args.d}/{reference.strip('.QC')}/{reference}\n")
            arg_list4 = []
            for name, data in Names_dict.items():
                data[['1_y', '1_x']].to_csv(f"{args.d}{name.strip('.QC')}/{name}_id_updates", sep="\t", index=False, header=False)
                with open(f"{args.d}/{name.strip('.QC')}/{name}_flips", 'w') as flip_file:
                    flip_file.write('\n'.join(Flips_dict[name]))
                arg_list4.append((f"{args.d}{name.strip('.QC')}/{name}", f"{args.d}/{name.strip('.QC')}/{name}_id_updates", f"{args.d}/{name.strip('.QC')}/{name}_flips", f"{args.d}{name.strip('.QC')}/{name}.ids", args.p))
                mergefile.write(f"{args.d}{name.strip('.QC')}/{name}.ids\n")

        outs = parallelQC(arg_list4, args.t, function="flipstrand_and_updateID")
        # pdb.set_trace()
        DF_dict = read_bims(f"{args.d}/files2premerge.txt")

        # for i, comb in enumerate(itertools.combinations(DF_dict.keys(), r=2)):
        #     if i==0:  # Hold first data frame as the reference and flip needed sites in subsequent data frames with respect to constant reference
        #         reference = comb[0]
        #     df1 = DF_dict[comb[0]]
        #     df2 = DF_dict[comb[1]]
        #     tdf = pd.merge(df1, df2, on=1)

        #Losing LOTS of SNPs here.
        # Poor overlap between aric and stjude/cog9906
        df = reduce(lambda df1, df2: pd.merge(df1, df2, on=1), DF_dict.values())
        df[1].to_csv(f"{args.d}/mergeSNPs.txt", header=False, index=False)

        arg_list5 = []
        with open(f"{args.d}/files2merge.txt", 'w') as mergefile:
            with open(f"{args.d}/files2premerge.txt", 'r') as premerge:
                for input in premerge:
                    input = input.strip()
                    arg_list5.append((input, f"{args.d}/mergeSNPs.txt", f"{input}.flt", args.p))
                    mergefile.write(f"{input}.flt\n")

        outs = parallelQC(arg_list5, args.t, function="var_flt")
        # # df = reduce(lambda df1, df2: pd.merge(df1, df2, on=1, how='outer', indicator=True), DF_list)
        # # same_ids = df[df._merge=="both"][1]
        # # Need to find and fix all conflicts by filtering df.  very difficult to accomplish with plink
        # # df = reduce(lambda df1, df2: pd.merge(df1, df2, on=(0,3)), DF_list)  # Merges on chromosome and position

        # First merge attempt is likely to have strand flip errors
        out1, err1 = merge_files(f"{args.d}/files2merge.txt", len(files), f"{args.d}/{args.o}", samples=args.s, plink=args.p)

        # TODO: perform flip-scan via plink to look for incorrect strand assignment in a subset of sample after finding best-flipped combon
    else:
        out1, err1 = merge_files(files[0], len(files), f"{args.d}/{args.o}",
                                 snp_file='all', samples=args.s, plink=args.p)
        # if os.path.exists(f"{args.d}/{args.o}-merge.missnp"):
        #     with open(f"{args.d}/files2merge.txt", 'r') as mergefiles:
        #         arg_list2 = []
        #         for file in mergefiles:
        #             arg_list2.append((file.strip(), f"{args.d}/{args.o}-merge.missnp", args.p))
        #     # How about instead just filtering the .missnp from all files 2 be merged
        #     pool = multiprocessing.Pool(args.t)
        #     outs,errs = pool.starmap(flipStrand, arg_list2)
        #     pool.close()
        #     pool.join()
        #     mergefiles = [x[0] for x in arg_list2]
        #     flip_files = [f"{x}.flip" for x in mergefiles]
        #     zip_tup = tuple([(x,y) for x,y in zip(mergefiles,flip_files)])
        #     combos = [x for x in itertools.product(*zip_tup)]
        #
        #     arg_list3 = []
        #
        #     for i, combo in enumerate(combos):
        #         with open(f"{args.d}/files2merge_{i}.txt", 'w') as tmp_file:
        #             for entry in combo:
        #                 tmp_file.write(entry + "\n")
        #         arg_list3.append((f"{args.d}/files2merge_{i}.txt", len(files), f"{args.d}/{args.o}_{i}",
        #                           f"{args.d}/mergeSNPs.txt", args.s, args.p))
        #     pool = multiprocessing.Pool(args.t)
        #     pdb.set_trace()
        #     outs, errs = pool.starmap(merge_files, arg_list3)
        #     pool.close()
        #     pool.join()
        #





