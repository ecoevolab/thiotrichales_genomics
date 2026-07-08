#!/bin/env bash
#$ -N antismashbeggiatoas2
#$ -o $JOB_NAME.log
#$ -e $JOB_NAME.error
#$ -cwd
#$ -S /bin/bash
#$ -l h_rss=32G
#$ -l h_rt=72:00:00
#$ -pe openmp 10

date
echo "===== Starting antiSMASH batch ====="

eval "$(conda shell.bash hook)"
conda activate antismash

INPUT_DIR="/mnt/atgc-d3/sur/users/mreyesr/exp/beggiatoas_outputs/gbk_files"
OUTPUT_DIR="/mnt/atgc-d3/sur/users/mreyesr/exp/beggiatoas_outputs/outputs_antismash"

mkdir -p "$OUTPUT_DIR"

for GBK_FILE in "$INPUT_DIR"/*.gbk; do
	[[ -e "$GBK_FILE" ]] || continue  
	BASENAME=$(basename "$GBK_FILE" .gbk)
    
	echo ">>> Corriendo antiSMASH para $BASENAME"
    
	antismash --cpus 10 \
          	--genefinding-tool prodigal \
          	--output-dir "$OUTPUT_DIR/$BASENAME" \
          	"$GBK_FILE"
done

echo "===== antiSMASH batch DONE ====="
date


