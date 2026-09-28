#!/bin/bash

DIR="/mnt/data/sur/users/mreyes/exp/thiotrichales/results/antismash"

echo -e "genoma\ttotal_bp_bgcs"

for genoma_dir in "$DIR"/*/; do
    genoma_id=$(basename "$genoma_dir")
    total_genoma=0

    while IFS= read -r -d '' gbk; do
        bp=$(awk '/^LOCUS/{print $3}' "$gbk")
        total_genoma=$((total_genoma + bp))
    done < <(find "$genoma_dir" -name "*.region*.gbk" -print0)

    echo -e "${genoma_id}\t${total_genoma}"
done
