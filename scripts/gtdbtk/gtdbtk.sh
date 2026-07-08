#!/bin/bash
#SBATCH --job-name=gtdbtk_all
#SBATCH --output=gtdbtk_%j.log
#SBATCH --error=gtdbtk_%j.error
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

eval "$(conda shell.bash hook)"
conda activate gtdbtk-2.3.2

export GTDBTK_DATA_PATH="/mnt/data/sur/shared_data/release220"

GENOME_DIR="/mnt/data/sur/shared_data/beggiatoas/fna_links"
OUT_DIR="/mnt/data/sur/users/mreyes/data/beggiatoas/results/gtdbtk"

for GENOME in "$GENOME_DIR"/*.fna; do
    BASENAME=$(basename "$GENOME" .fna)
    OUT_FOLDER="$OUT_DIR/$BASENAME"
    TEMP_DIR=$(mktemp -d)
    
    cp "$GENOME" "$TEMP_DIR/"
    
    gtdbtk classify_wf \
        --genome_dir "$TEMP_DIR" \
        --out_dir "$OUT_FOLDER" \
        --extension fna \
        --cpus 8 \
        --prefix "$BASENAME" \
        --skip_ani_screen
    
    rm -rf "$TEMP_DIR"
done

date
