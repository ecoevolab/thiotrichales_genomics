#!/bin/bash

# uso: ./top_hits.sh resultados_ddl_tblout.txt "D-alanine-D-alanine ligase" salida_ddl.tsv

INPUT="$1"
PRODUCTO_ESPERADO="$2"
OUTPUT="$3"

if [ -z "$INPUT" ] || [ -z "$PRODUCTO_ESPERADO" ] || [ -z "$OUTPUT" ]; then
    echo "uso: $0 archivo_tblout.txt \"producto esperado\" salida.tsv"
    exit 1
fi

echo -e "genome_locustag\tgenome_id\tevalue\texp_dom\tdescription\tproduct_match" > "$OUTPUT"

grep -v "^#" "$INPUT" | awk -v producto="$PRODUCTO_ESPERADO" '
{
    target = $1
    evalue = $5
    exp_dom = $11

    desc = ""
    for (i = 19; i <= NF; i++) {
        desc = desc $i " "
    }
    sub(/ $/, "", desc)

    split(target, partes, "_")
    genome_id = partes[1]

    if (tolower(desc) ~ tolower(producto)) {
        match_flag = "TRUE"
    } else {
        match_flag = "FALSE"
    }

    print target"\t"genome_id"\t"evalue"\t"exp_dom"\t"desc"\t"match_flag
}' | sort -t$'\t' -k2,2 -k3,3g | awk -F'\t' '!seen[$2]++' >> "$OUTPUT"

echo "listo, tabla guardada en $OUTPUT"
