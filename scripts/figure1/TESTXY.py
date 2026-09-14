import os
import re
from tree_class import *
from preprocess import *
from store_junction import *
from plot_assembly import *
import pickle
import pprint 
sr_samples = ['2655-CB-90','2655-CB-96','2655-CB-93','81279','2655-CB-90','2655-CB-96','2655-CB-93','81279','2655-CB-90','2655-CB-96','2655-CB-93','81279','2655-CB-90','2655-CB-96','2655-CB-93','81279','2655-CB-90','2655-CB-96','2655-CB-93','81279','2655-CB-90','2655-CB-96','2655-CB-93','81279','2655-CB-90','2655-CB-96','2655-CB-93','81279','2655-CB-90','2655-CB-96','2655-CB-93','81279']
# lr_samples = ['5219-CB-1.bam','4102-CB_UPCI-152.bam', '4745-CB-1_UDSCC-2.bam','UM-SCC-47.bam','2610-CB-3_Tumor-A.bam','3409-CB-1_UM104.bam','3816-CB_PDX-294R.bam','4745-CB-2_VUSCC-147.bam','5219-CB-10.bam','5219-CB-11.bam','5219-CB-12.bam','5219-CB-14.bam','5219-CB-15.bam','5219-CB-16.bam','5219-CB-17.bam','5219-CB-18.bam','5219-CB-2.bam','5219-CB-3.bam']
lr_samples = ['90_chr9_reads.allAligns.sorted.bam', '152_chr9_reads.allAligns.sorted.bam' ]
sample_name = ['UPCI:90','UPCI:152','UD:2','UM:47','TumorA','UM104','PDX-294R','VU147','2453','1410','2352','3749','3720','3536', '3173','3744','UPCI:154','UMCV:6']

sr_path = '/home/wenjingu/scratch/HPV_fusion/targeted_exom/our_res/allRes/'
# lr_path = '/nfs/turbo/oto-brenner-lab/Xinyi/long-DNA/Long_read_DNA_alignment/'
lr_path = '/nfs/turbo/oto-brenner-lab/Xinyi/UPCI90152/'
res_path = "chr910_24_2025"
os.system(f'mkdir -p {res_path}')
# depth_path = '/home/wenjingu/remillsscr/HPV_fusion/HPV_nanopore_assembler/'
#'2453','1410','2352','3749','3720','3536', '3173','3744','UPCI:154','UMCV:6'
targets = {'UPCI:90','UPCI:152'}
target_chrs = ['9']


for i in range(min(len(sample_name), len(lr_samples), len(sr_samples))):
    name = sample_name[i]
    if name not in targets:
        continue

    bam = f'{lr_path}/{lr_samples[i]}'
    lr_sample_name = lr_samples[i].split('.')[0]
    # depthFile = f'{depth_path}/{lr_sample_name}.txt'  

    for chro in target_chrs:
        junctionRes, dataFrameList = extract_junction(bam, chro, "gi|333031|lcl|HPV16REF.1|")
        junctionList, clusteredJunctionRes = merge_junction(junctionRes)
        write_junction(junctionList, sample=f'{name}.{chro}.merge', outputPath=res_path)
        definedJunctionList = define_junction(junctionList)
        write_junction(definedJunctionList, sample=f'{name}.{chro}', outputPath=res_path)
        # finalJunction = flip_reads_by_nonhpv_minus(definedJunctionList)
        finalJunction, summary = flip_reads_by_weighted_host_strand_grouped(definedJunctionList)
        summary.to_csv(f"{res_path}/{name}.{chro}.read_flip_summary.csv", index=False)
        write_junction(finalJunction, sample=f'{name}.{chro}.final', outputPath=res_path)
        pickle_path = f"{res_path}/{name}.{chro}.final.junctionListDic.pickle"
        # rev_df = find_reverse_junctions_from_pickle(pickle_path, save_csv=True) 
        mirrors_summary, mirrors_detail = detect_mirrors_after_flip(finalJunction)
        mirrors_summary.to_csv(f"{res_path}/{name}.{chro}.mirror_summary.csv", index=False)
        mirrors_detail.to_csv(f"{res_path}/{name}.{chro}.mirror_detail.csv", index=False)
        allAssemblies = build_tree(finalJunction)
        # depthDic = depth_dic(definedJunctionList)
        # cov = read_cov(chro, depthFile)
        # uniqfilteredAllAssemblies = filter_read(allAssemblies, cov, depthDic)
        # write_assembly(uniqfilteredAllAssemblies, f'{res_path}/{name}.{chro}.filtered.assembile.pickle')
        write_assembly(allAssemblies, f'{res_path}/{name}.{chro}.filtered.assembile.pickle')

# plot
for file in os.listdir(res_path):
    for name in targets:
        if 'filtered.assembile' in file and name in file:
            with open(f'{res_path}/{file}', 'rb') as inputFile:
                uniqFilteredAllAssembiles = pickle.load(inputFile)
                chro = file.split('.')[1]
            sampleName = '.'.join(file.split('.')[0:2])
            outdir = f'{res_path}/{sampleName}/{chro}'
            os.system(f'mkdir -p {outdir}')
            plot_backward(uniqFilteredAllAssembiles, outdir, sampleName)
            plot_forward(uniqFilteredAllAssembiles, outdir, sampleName)

