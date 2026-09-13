#!/bin/bash
#SBATCH --job-name=orthofinder_thiotrichales
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=120:00:00

date
echo "===== Starting OrthoFinder Thiotrichales ====="

eval "$(conda shell.bash hook)"
conda activate orthofinder

INPUT_DIR="/mnt/data/sur/users/mreyes/data/thiotrichales/proteinas"
OUTPUT_DIR="/mnt/data/sur/users/mreyes/exp/thiotrichales/results/orthofinder_${SLURM_JOB_ID}"

orthofinder \
    -f "$INPUT_DIR" \
    -o "$OUTPUT_DIR" \
    -t $SLURM_CPUS_PER_TASK \
    -a $SLURM_CPUS_PER_TASK

echo "===== OrthoFinder Done ====="
date
