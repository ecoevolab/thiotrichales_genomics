#!/bin/bash

# uso: ./top_hits.sh resultados_hmmer.txt "producto_esperado" salida.tsv

INPUT="$1"   #resultados provenientes de hmmer
PRODUCTO_ESPERADO="$2" # el producto proteico que se espera obtener ejemplo "D-alanine--D-alanine ligase" 
OUTPUT="$3" # nombre del archivo donde se va a guardar el resultado en formato tsv


# En general lo hace es verificar si se le pasaron los 3 argumentos solicitados al momento de ejecutar el script
# -z significa: zero lenght
# Primero se verifica si el argumento del input esta vacio o si el argumento de producto esperado esta vacio o si el argumento de ouput esta vacio.
if [ -z "$INPUT" ] || [ -z "$PRODUCTO_ESPERADO" ] || [ -z "$OUTPUT" ]; then

# si alguno de los argumentos esta vacio se imprime en pantalla el formato de uso del script
    echo "uso: $0 archivo_de_salida_hmmer.txt \"producto esperado\" salida.tsv"
    exit 1
# Finalizacion con error
fi # se cierra el if en bash

# Se imprime los encabezados en el archivo de salida, el -e es para decirle a echo que separe 
# encabezados por un tabulador

echo -e "genome_locustag\tgenome_id\tevalue\texp_dom\tdescription\tproduct_match" > "$OUTPUT"


# Se lee el archivo de entrada y con -v se excluyen todas las lineas que empiecen con #, luego 
# se crea una variable llamada producto que dentro tiene el producto que se espera y la ultima # comilla es para indicarle a awk que todo lo que sigue son intrucciones que se ejecutan una 
# vez por cada linea

grep -v "^#" "$INPUT" | awk -v producto="$PRODUCTO_ESPERADO" '
{
    target = $1. # almacena la columna 1 del archivo de salida de Hmmer que corresponde al ID
    evalue = $5. # almacena el e-value de "full sequence"
    exp_dom = $11. # almacena el número de dominios esperados

# La descripción del producto puede tener varias palabras (ej. "D-alanine--D-alanine ligase # # B"), y awk separa por espacios, así que cada palabra cae en una columna distinta
# Este ciclo junta todas las columnas desde la 19 hasta la última (NF)  para reconstruir la 
# descripción completa en una sola variable. 
    desc = "" 
    for (i = 19; i <= NF; i++) {
        desc = desc $i " "
    } # fin del ciclo for

# busca un patro dentro de la variable 
    sub(/ $/, "", desc)

# Se corta en pedazos el ID ej. 2835558653_GGNGOHCK_00698) guarda los pedazos en el arreglo
# partes (partes[1] = "2835558653", partes[2] = "GGNGOHCK", partes[3] = "00698")

    split(target, partes, "_")
    genome_id = partes[1]

# Convierte todo el texto a minuscula de la descrición de hmmer y la compara con el producto 
# que se puso en consola que se espera 
    if (tolower(desc) ~ tolower(producto)) {
         match_flag = "TRUE". # la condicion es verdadera (coincide)
    } else {
        match_flag = "FALSE" # la descripcion no coincide
    }

# Se imprimen todas las variables con la informacion guardad separados por un tab y se forma  # guarda en el archivo de salida
    print target"\t"genome_id"\t"evalue"\t"exp_dom"\t"desc"\t"match_flag
}' | sort -t$'\t' -k2,2 -k3,3g | awk -F'\t' '!seen[$2]++' >> "$OUTPUT"

# mensaje en pantalla para decir donde se guardo la salida/resultados
echo "listo, tabla guardada en $OUTPUT"
