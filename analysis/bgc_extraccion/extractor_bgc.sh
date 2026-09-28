#!/bin/bash

# Objetivo:
# Toma un directorio donde se encuentren subcarpetas de salida de antismash y calcula para cada una de las subcarpetas (genomas) 
# el total de pares de bases que ocupan los BGC en el genoma
#
# Salida:
# Un archivo .csv, donde la primera columna corresponder al nombre de la subcarpeta que proceso y la segunda es el total de pares de bases





# Directorio donde se encuentra la salida de antismash
DIR="/mnt/data/sur/users/mreyes/exp/thiotrichales/results/antismash"


echo -e "genoma\ttotal_bp_bgcs"

# iteracion sobre cada uno de los genomas 
for genoma_dir in "$DIR"/*/; do
    genoma_id=$(basename "$genoma_dir")
    total_genoma=0

# Se busca la primera linea de cada archivo gbk y se extrae el numero de pares de bases que esta en la variable LOCUS
    while IFS= read -r -d '' gbk; do
        bp=$(awk '/^LOCUS/{print $3}' "$gbk")
        total_genoma=$((total_genoma + bp))
    done < <(find "$genoma_dir" -name "*.region*.gbk" -print0)

    echo -e "${genoma_id}\t${total_genoma}"
done
