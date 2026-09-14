#!/bin/bash
#SBATCH --job-name=T
#SBATCH --mail-user=ppxinyi@umich.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100gb
#SBATCH --time=100:00:00
#SBATCH --account=chadbren99
#SBATCH --partition=standard
#SBATCH --output=1.log
#SBATCH --error=1.err   

# source ~/miniconda3/etc/profile.d/conda.sh
conda activate assembly

cd /nfs/turbo/oto-brenner-lab/Xinyi/HPV_nanopore_assembler

/home/ppxinyi/miniconda3/envs/assembly/bin/python -c "import sys, pysam, pandas, matplotlib; print(sys.executable); print('OK')"

/home/ppxinyi/miniconda3/envs/assembly/bin/python TESTXY.py
