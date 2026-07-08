#!/bin/bash
#SBATCH --job-name=antismash_missing
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

echo "===== Starting antiSMASH batch ====="
date

eval "$(conda shell.bash hook)"
conda activate antismash5

PROKKA_DIR="/mnt/data/sur/users/mreyes/data/beggiatoas/results/prokka"
OUT_DIR="/mnt/data/sur/users/mreyes/data/beggiatoas/results/antismash"
mkdir -p "$OUT_DIR"

TARGETS=(
3001836915
640963011
640963012
8002131595
8002135289
8056884285
8066494198
8071263300
8074213743
8125999560
)

for ID in "${TARGETS[@]}"; do
    INPUT_FOLDER="${PROKKA_DIR}/${ID}"
    OUTPUT_FOLDER="${OUT_DIR}/${ID}"
    
    rm -rf "$OUTPUT_FOLDER"
    mkdir -p "$OUTPUT_FOLDER"
    
    GBK=$(find "$INPUT_FOLDER" -maxdepth 1 -name "*.gbk" -o -name "*.gbff")
    
    if [ -z "$GBK" ]; then
        echo "ERROR: No se encontró archivo .gbk para $ID"
        continue
    fi
    
    echo "Procesando $ID..."
    
    antismash \
        --output-dir "$OUTPUT_FOLDER" \
        --cpus 4 \
        --genefinding-tool prodigal\
        "$GBK"
done

echo "===== antiSMASH batch DONE ====="
date
