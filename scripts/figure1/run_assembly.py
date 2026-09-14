import os
import argparse
import re
from tree_class import *
from preprocess import *
from store_junction import *
from plot_assembly import *
import pickle
import pprint 

parser = argparse.ArgumentParser()

parser.add_argument("-b", "--bam", required=True)
parser.add_argument("-s", "--sample", required=True)
parser.add_argument(
    "-c",
    "--chr",
    required=True,
    help="chromosome(s), e.g. 9 or 9,10,11 or all"
)
parser.add_argument("-o", "--outdir", required=True)

args = parser.parse_args()

bam = args.bam
sample_name = args.sample
target_chr_input = args.chr
if target_chr_input.lower() == "all":
    target_chrs = [str(i) for i in range(1,23)] + ["X","Y"]
else:
    target_chrs = target_chr_input.split(",")
res_path = args.outdir

for target_chr in target_chrs:

    print(f"Processing chr{target_chr}")

    junctionRes, dataFrameList = extract_junction(
        bam,
        target_chr,
        "gi|333031|lcl|HPV16REF.1|"
    )

    junctionList, clusteredJunctionRes = merge_junction(junctionRes)

    write_junction(
        junctionList,
        sample=f'{sample_name}.{target_chr}.merge',
        outputPath=res_path
    )

    definedJunctionList = define_junction(junctionList)

    write_junction(
        definedJunctionList,
        sample=f'{sample_name}.{target_chr}',
        outputPath=res_path
    )

    finalJunction, summary = \
        flip_reads_by_weighted_host_strand_grouped(
            definedJunctionList
        )

    summary.to_csv(
        f"{res_path}/{sample_name}.{target_chr}.read_flip_summary.csv",
        index=False
    )

    write_junction(
        finalJunction,
        sample=f'{sample_name}.{target_chr}.final',
        outputPath=res_path
    )

    mirrors_summary, mirrors_detail = \
        detect_mirrors_after_flip(finalJunction)

    mirrors_summary.to_csv(
        f"{res_path}/{sample_name}.{target_chr}.mirror_summary.csv",
        index=False
    )

    mirrors_detail.to_csv(
        f"{res_path}/{sample_name}.{target_chr}.mirror_detail.csv",
        index=False
    )

    allAssemblies = build_tree(finalJunction)

    write_assembly(
        allAssemblies,
        f'{res_path}/{sample_name}.{target_chr}.filtered.assembile.pickle'
    )

    outdir = f'{res_path}/{sample_name}/{target_chr}'
    os.system(f'mkdir -p {outdir}')

    plot_backward(allAssemblies, outdir, sample_name)
    plot_forward(allAssemblies, outdir, sample_name)