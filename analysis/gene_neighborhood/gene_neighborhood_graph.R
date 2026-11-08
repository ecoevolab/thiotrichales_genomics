rm(list = ls())          
library(tidyverse)       
library(gggenes)         # is a ggplot2 extension for drawing gene arrow maps.



gff_dir    <- "/Users/monicareyes/Desktop/gff_procesados"  # carpeta con los gff de los genomas
csv_coords <- "coordenadas_genes_focales.csv"              # coordenadas genomicas de los genes de interes 
orden_txt  <- "/Users/monicareyes/Desktop/ID.txt"          # orden de las filas

ventana    <- 10000     # pares de bases a los lados de cada ancla
anclas     <- c("ddl", "mreB", "rodZ")   # anclas 
listones   <- TRUE      # bandas que unen al mismo gen en las diferentes filas
escala_pb  <- 3000      # largo de la barra de escala 

# genes que queremos colorear/resaltar, todos los demás que no esten en la lista van a ser grises
genes <- c("murG", "murC", "murB", "ddl", "ftsZ", "ftsQ", "ftsA", "lpxC", "mreB", "mreC",
           "mreD", "mrdA", "rodA", "mltB", "rodZ", "rlpA")

# leer el archivo txt para definir el orden de las filas
orden <- read_tsv(orden_txt, col_names = FALSE, show_col_types = FALSE,
                  col_types = cols(.default = "c"))
orden <- tibble(
  genome_id = str_trim(orden$X1),                               
  etiqueta  = if (ncol(orden) >= 2) coalesce(orden$X2, orden$X1) # caso donde no hay ID 
  else orden$X1
) %>%
  mutate(y = n():1)   


leer_gff <- function(archivo_gff) {
  genome_id <- basename(archivo_gff) %>% str_remove("_limpio\\.gff$")  
  read_tsv(archivo_gff, col_names = FALSE, show_col_types = FALSE,
           col_types = cols(.default = "c")) %>%                       
    rename(contig = X1, source = X2, type = X3, start = X4,            
           end = X5, score = X6, strand = X7, phase = X8, attributes = X9) %>%
    filter(type == "CDS") %>%                                          
    mutate(genome_id = genome_id,
           start = as.integer(start),                                 
           end   = as.integer(end),
           locus_tag = str_extract(attributes, "locus_tag=[^;]+") %>%
             str_remove("locus_tag=")) %>%
    select(genome_id, contig, start, end, strand, locus_tag)
}


archivos_gff <- list.files(gff_dir, pattern = "_limpio\\.gff$", full.names = TRUE)
cds <- map_dfr(archivos_gff, leer_gff)   


focales <- read_csv(csv_coords, show_col_types = FALSE,
                    col_types = cols(genome_id = col_character())) %>%
  filter(coordenada_encontrada, gen %in% genes) %>%  
  group_by(genome_id, locus_tag) %>%                 
  slice_min(evalue, n = 1, with_ties = FALSE) %>%    
  ungroup()


vecindario <- function(gen_ancla) {
  
  anc <- focales %>%
    filter(gen == gen_ancla) %>%
    group_by(genome_id) %>%
    slice_min(evalue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(genome_id, contig, a_start = start, a_end = end, a_strand = strand)  
  
  datos <- cds %>%
    inner_join(anc, by = c("genome_id", "contig")) %>%
    filter(end >= a_start - ventana, start <= a_end + ventana) %>%
    mutate(
      xmin = if_else(a_strand == "+", start - a_start, a_end - end),
      xmax = if_else(a_strand == "+", end - a_start,   a_end - start),
      forward = (strand == "+") == (a_strand == "+")
    ) %>%
    left_join(select(focales, genome_id, locus_tag, gen),
              by = c("genome_id", "locus_tag")) %>%
    mutate(gen = replace_na(gen, "otro")) %>%       
    group_by(genome_id, gen) %>%
    arrange(xmin, .by_group = TRUE) %>%
    mutate(copia = if_else(gen != "otro", row_number(), NA_integer_),  
           n_copias = n()) %>%                                         
    ungroup() %>%
    mutate(clave = if_else(gen == "otro", NA_character_,
                           if_else(n_copias > 1, paste0(gen, copia), gen))) %>%
    right_join(orden, by = "genome_id")
  
  datos
}

armar_listones <- function(datos, alto = 0.15) {
  d <- datos %>% filter(!is.na(clave)) %>% select(y, clave, xmin, xmax)  # solo genes de interés
  
  d_ab <- d                       
  d_ar <- mutate(d, y = y - 1)  
  pares <- inner_join(d_ab, d_ar, by = c("y", "clave"), suffix = c("_ab", "_ar")) %>%
    mutate(id = row_number())     
  bind_rows(
    pares %>% transmute(id, x = xmin_ab, y = y + alto),
    pares %>% transmute(id, x = xmax_ab, y = y + alto),
    pares %>% transmute(id, x = xmax_ar, y = y + 1 - alto),
    pares %>% transmute(id, x = xmin_ar, y = y + 1 - alto)
  )
}

# colores para los genes 
paleta <- c(
  murG = "#668B8B",
  murC = "#8B7765",
  murB = "#CD96CD",
  ddl  = "black",
  ftsZ = "#9BCD9B",
  ftsQ = "#009ACD",
  ftsA = "red",
  lpxC = "#8968CD",
  mreB = "black",
  mreC = "#9BCD9B",
  mreD = "red",
  mrdA = "#009ACD",
  rodA = "#4F6D6D",
  mltB = "#C77DCA",
  rodZ = "#E63946",
  rlpA = "#5C4D9A",
  otro = "grey85"
)

# dibuja los paneles
dibujar <- function(gen_ancla) {
  datos <- vecindario(gen_ancla)
  
  p <- ggplot()
  if (listones) {
    lis <- armar_listones(datos)
    if (nrow(lis) > 0) {
      p <- p + geom_polygon(data = lis, aes(x = x, y = y, group = id),
                            fill = "grey88", colour = NA)
    }
  }
  
  p <- p +
    geom_hline(data = orden, aes(yintercept = y),
               colour = "grey60", linetype = "dotdash", linewidth = 0.3) +
    # flechas de genes
    geom_gene_arrow(data = filter(datos, !is.na(xmin)),
                    aes(xmin = xmin, xmax = xmax, y = y, fill = gen, forward = forward),
                    arrowhead_height = unit(3, "mm"), arrow_body_height = unit(2, "mm")) +
    scale_fill_manual(values = paleta, breaks = genes, name = NULL) +  
    scale_y_continuous(breaks = orden$y, labels = orden$etiqueta,
                       expand = expansion(add = c(1.2, 0.6))) +        
    labs(title = paste0("Vecindario de ", gen_ancla), x = "pb respecto al ancla", y = NULL) +
    theme_classic() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  
  p +
    annotate("segment", x = -ventana, xend = -ventana + escala_pb, y = 0, yend = 0) +
    annotate("text", x = -ventana + escala_pb / 2, y = 0.45,
             label = paste0(escala_pb, " bp"), size = 3)
}

# generar los paneles
paneles <- map(anclas, dibujar) %>% set_names(anclas)   

# Guardamos un pdf de los paneles
iwalk(paneles, ~ ggsave(paste0("vecindario_", .y, ".pdf"), .x,
                        width = 12, height = max(5, 0.25 * nrow(orden) + 2)))

paneles$ddl    
paneles$mreB
paneles$rodZ

# Tabla para verificar si el orden en el que se graficaron los genes es correcta
tabla_vecindario <- function(gen_ancla) {
  vecindario(gen_ancla) %>%
    filter(!is.na(xmin), gen != "otro") %>%     # quita filas vacías y genes que no nos interesan 
    arrange(desc(y), xmin) %>%                  # orden del txt
    transmute(panel = gen_ancla, genome_id,
              volteado = a_strand == "-",       
              contig, locus_tag, gen = clave,
              start_genomico = start, end_genomico = end, hebra_genoma = strand,
              xmin, xmax)
}

tablas <- map_dfr(anclas, tabla_vecindario)
write_csv(tablas, "vecindarios_detalle.csv")
