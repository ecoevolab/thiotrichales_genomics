PROKKA_DIR="$1"

if [ -z "$PROKKA_DIR" ]; then
    echo "uso: $0 carpeta_prokka"
    exit 1
fi

OUT_DIR="${PROKKA_DIR}/Gff"
mkdir -p "$OUT_DIR"

n_procesados=0
n_faltantes=0

for carpeta in "$PROKKA_DIR"/*/; do
    genome_id=$(basename "$carpeta")
    input_gff="${carpeta}${genome_id}.gff"
    output_gff="${OUT_DIR}/${genome_id}_limpio.gff"

    if [ ! -f "$input_gff" ]; then
        echo "aviso: no encontre gff para $genome_id, se omite"
        n_faltantes=$((n_faltantes + 1))
        continue
    fi

    awk '/^##FASTA/{exit} !/^#/' "$input_gff" | sort -k1,1V -k4,4n > "$output_gff"

    n_procesados=$((n_procesados + 1))
done

echo "listo: $n_procesados genomas procesados, $n_faltantes sin gff encontrado"
