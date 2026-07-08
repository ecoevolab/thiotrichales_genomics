#!/bin/env bash
#$ -N beggiatoas2
#$ -o $JOB_NAME.log
#$ -e $JOB_NAME.error
#$ -cwd
#$ -S /bin/bash
#$ -l h_rss=16G
#$ -l h_rt=24:00:00
#$ -pe openmp 10

date
echo "===== Starting Prokka Batch ====="

eval "$(conda shell.bash hook)"
conda activate prokka

INPUT_DIR="/mnt/atgc-d3/sur/shared_data/beggiatoas/fna_links"
OUTPUT_DIR="/mnt/atgc-d3/sur/users/mreyesr/exp/beggiatoas_outputs/outputs_prokka"

for INPUT_GENOME in "$INPUT_DIR"/*.{fna,fa,fasta}; do
        [[ -e "$INPUT_GENOME" ]] || continue

        SAMPLE_NAME=$(basename "$INPUT_GENOME" | sed -E 's/\.(fna|fa|fasta)$//')

        echo ">>> Ejecutando Prokka para $SAMPLE_NAME"

        prokka --outdir "$OUTPUT_DIR/$SAMPLE_NAME" \
        --prefix "$SAMPLE_NAME" \
        --cpus 10 \
        --kingdom Bacteria \
        --centre X \
        --compliant \
        "$INPUT_GENOME"
done

echo "===== Prokka Batch Done ====="
date
