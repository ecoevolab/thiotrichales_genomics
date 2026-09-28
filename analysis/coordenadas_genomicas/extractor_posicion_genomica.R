library(tidyverse)

# Objetivo:
# Toma los archivos .gff y los archivos de salida del alineamiento con hmmer
# por cada gen, busca las coordenadas de los locus_tag para obtener la posición genomica
# (contig, start, end, strand) de cada gen confirmado
#
# Entrada:
# - Carpeta con subcarpetas por cada gen, cada uno con su el formato
# salida_<gen>_filtrado.tsv
# - Carpeya con los GFF individuales ordenados y procesados (solo líneas CDS
# ordenados por contig y posición)
#
# Salida: 
# - Archivo coordenadas_genes.csv, una fila por cada gen confirmado 
# NOTA: confirmado -> gen del cual se encontro el locus_tag y coordenadas en el gff.
# Columnas: gen, genome_id, locustag, evalue, contig, start, end, strand, 
# coordenada_encontrada (TRUE/FALSE)
# 


# Directorio donde estan las subcarpetas por gen - de los alineamientos/hmmer
base_dir <- "/Users/monicareyes/Desktop/aligments"


# Carpeta donde estan los GFFs individuales ordenados y limpios (por posicion, sin comentarios y ordenados por locus)
gff_dir <- "/Users/monicareyes/Desktop/gff_procesados"

# Lista de los 14 genes (de interes - que queremos buscar, tienen que coincidir con el nombre de la carpeta de alineamiento)
genes <- c("murG", "murC", "murB", "ddl", "ftsQ", "ftsA",
           "ftsZ", "lpxC", "mreB", "mreC", "mreD", "mrdA", "rodA", "mltB", "rlpA", "rodZ")


leer_filtrado <- function(gen) {
  archivo <- file.path(base_dir, gen, paste0("salida_", gen, "_filtrado.tsv"))
  read_tsv(archivo, show_col_types = FALSE,
           col_types = cols(genome_id = col_character())) %>%
    mutate(
      gen = gen,
      # locus_tag = lo que sigue despues del genome_id y el guion bajo
      locus_tag = str_remove(genome_locustag, paste0("^", genome_id, "_"))
    ) %>%
    select(gen, genome_id, locus_tag, evalue)
}

hits_focales <- map_dfr(genes, leer_filtrado)

# leer todos los GFF (una sola vez) y extraer contig/posicion/strand/locus_tag ----

leer_gff <- function(archivo_gff) {
  genome_id <- basename(archivo_gff) %>% str_remove("_limpio\\.gff$")
  
  read_tsv(archivo_gff, show_col_types = FALSE, col_names = FALSE,
           col_types = cols(.default = "c")) %>%
    rename(contig = X1, source = X2, type = X3, start = X4,
           end = X5, score = X6, strand = X7, phase = X8, attributes = X9) %>%
    filter(type == "CDS") %>%
    mutate(
      genome_id = genome_id,
      start = as.integer(start),
      end = as.integer(end),
      locus_tag = str_extract(attributes, "locus_tag=[^;]+") %>% str_remove("locus_tag=")
    ) %>%
    select(genome_id, contig, start, end, strand, locus_tag)
}

archivos_gff <- list.files(gff_dir, pattern = "_limpio\\.gff$", full.names = TRUE)
todos_los_cds <- map_dfr(archivos_gff, leer_gff)

#cruzar hits confirmados con sus coordenadas del GFF 

coordenadas_focales <- hits_focales %>%
  left_join(todos_los_cds, by = c("genome_id", "locus_tag")) %>%
  mutate(coordenada_encontrada = !is.na(contig))

#verificacion: deberia haber coordenadas para todos los hits

sin_coordenadas <- coordenadas_focales %>% filter(!coordenada_encontrada)
if (nrow(sin_coordenadas) > 0) {
  cat("hits sin coordenadas encontradas en el gff :( \n")
  print(sin_coordenadas)
} else {
  cat("todos los hits confirmados tienen coordenadas, alegria")
}

write_csv(coordenadas_focales, "coordenadas_genes_focales.csv")

cat("ya termine :)")

