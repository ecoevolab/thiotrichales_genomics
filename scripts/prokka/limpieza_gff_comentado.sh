#!/bin/bash # =
------------------------------------------------------------------------------
# limpiar_ordenar_gff.sh # 
------------------------------------------------------------------------------
# Objetivo
# Recorrer las subcarpetas de genomas dentro de la carpeta prokka, extraer el archivo .gff de cada una
# quitar las lineas de comentarios, la secuencia ##FASTA que prokka agrega por defecto, luego
# ordenar por contig y por posición genomica (START). Finalmente guardarlos en un nuevo directorio
#
#
# Entrada: carpeta de prokka donde se tienen todos los subdirectiorios de cada genoma
# cada una con un archivo gff
#
# Salida: un subdirectorio Gff, dentro de la carpeta prokka, que dentro tendra los archivos
# .gff limpios y ordenados.
#
# Ejemplo de uso: /limpiar_ordenar_gff.sh <direccion_carpeta_prokka>

PROKKA_DIR="$1"  #ruta de la carpeta de prokka

# Si no se pasa un directorio se termina el programa con error
if [ -z "$PROKKA_DIR" ]; then 
    echo "uso: $0 carpeta_prokka"
    exit 1
fi

# directorio de salida
OUT_DIR="${PROKKA_DIR}/Gff"
mkdir -p "$OUT_DIR" 

# se inicializan los contadores
n_procesados=0
n_faltantes=0

# se itera sobre cada carpeta del directorio
for carpeta in "$PROKKA_DIR"/*/; do
# Se toma el ID de la carpeta del genoma que se esta iterando
    genome_id=$(basename "$carpeta")
# Se construye la ruta completa al gff del genoma que estamos iterando
    input_gff="${carpeta}${genome_id}.gff"
# se construye la ruta dende se va guardar
    output_gff="${OUT_DIR}/${genome_id}_limpio.gff"

# Caso donde no se encuentre dentro de la n subcarpeta un archivo .gff y continua con la siguiente subcarpeta
    if [ ! -f "$input_gff" ]; then
        echo "aviso: no encontre gff para $genome_id, se omite"
        n_faltantes=$((n_faltantes + 1))
        continue
    fi

# Cuando se encuentra un gff se remueven las lineas de comentarios y los ##FASTA, se remuerven todas features excepto CDS
# luego  se ordena la columna 1  que es la de los contigs y luego la 4 por posicion de incio finalmente se guarda todo en la salida
awk -F'\t' '/^##FASTA/{exit} !/^#/ && $3=="CDS"' "$input_gff" | sort -k1,1V -k4,4n > "$output_gff"
    n_procesados=$((n_procesados + 1))
done
 
# Mensaje para decir que ya se ejecuto todo y en cuales no se encontro los gff
echo "listo: $n_procesados genomas procesados, $n_faltantes sin gff encontrado"
