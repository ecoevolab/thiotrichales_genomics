library(tidyverse)

# Directorio donde estan las subcarpetas por gen - de los alineamientos
base_dir <- "/Users/monicareyes/Desktop/aligments"

# Lista de los 14 genes 
genes <- c("murG", "murC", "murB", "ddl", "ftsQ","ftsA", 
           "ftsZ", "lpxC", "mreB", "mreC", "mreD", "mrdA", "rodA", "mltB", "rlpA", "rodZ")

# Opcional, se puede pasar los ID de los genes que se quieren directamente en la variable genome_ID
# En este caso se tomaron los GFFs por que se intentaban cruzar las coordenadas, pero se descarto esa ópcion 
# Y se transfirio esa idea a un nuevo script
# Carpeta donde estan los GFFs individuales
gff_dir <- "/Users/monicareyes/Desktop/gff_procesados"

# Lista de ID de los genomas, puede ser reemplazado por un archivo con los ID directamente
genoma_ids <- list.files(gff_dir, pattern = "_limpio\\.gff$") %>%
  str_remove("_limpio\\.gff$") # remueve el .gff de los gff individuales

# Primero se lee y  cuentan las copias de cada gen por cada genoma y se almacenan 
leer_conteo <- function(gen) {
  archivo <- file.path(base_dir, gen, paste0("salida_", gen, "_filtrado.tsv"))
  read_tsv(archivo, show_col_types = FALSE,
           col_types = cols(genome_id = col_character())) %>%
    count(genome_id, name = "n_copias") %>%
    mutate(gen = gen)
}

conteos <- map_dfr(genes, leer_conteo)

# los genomas sin ningun hit para ese gen reciben el valor de 0
tabla_completa <- conteos %>%
  complete(genome_id = genoma_ids, gen = genes, fill = list(n_copias = 0))

# resumen en tabla de los genes y su relacion con los genomas - para verificar a mano
# o para un uso posterior
tabla_resumen <- tabla_completa %>%
  pivot_wider(names_from = gen, values_from = n_copias)

write_csv(tabla_resumen, "resumen_copias_por_gen.csv")

# Heatmap
tabla_plot <- tabla_completa %>%
  mutate(
    categoria = case_when(
      n_copias == 0 ~ "absent",
      n_copias == 1 ~ "1 copy",
      n_copias == 2 ~ "2 copies",
      n_copias >= 3 ~ "3+ copies"
    ),
    categoria = factor(categoria, levels = c("absent", "1 copy", "2 copies", "3+ copies")),
    gen = factor(gen, levels = genes)
  )

p <- ggplot(tabla_plot, aes(x = gen, y = genome_id, fill = categoria)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(
    "absent"   = "grey90",
    "1 copy"   = "steelblue3",
    "2 copies"  = "sienna1",
    "3+ copies" = "red2"
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7)
  ) +
  labs(
    x = "Gene", y = "Genome ID", fill = "number of copies",
    title = "Preliminary - Number of gene copies"
  )
p

