#!/bin/bash
#SBATCH --job-name=hpvAsm
#SBATCH --account=xxxxx
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=1.out
#SBATCH --error=1.err

source activate assembly_env

while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--bam)
            bam="$2"
            shift 2
            ;;
        -s|--sample)
            sample="$2"
            shift 2
            ;;
        -c|--chr)
            chr="$2"
            shift 2
            ;;
        -o|--outdir)
            outdir="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

echo "BAM: $bam"
echo "Sample: $sample"
echo "Chr: $chr"
echo "Output: $outdir"

python run_assembly.py \
    -b "$bam" \
    -s "$sample" \
    -c "$chr" \
    -o "$outdir"
