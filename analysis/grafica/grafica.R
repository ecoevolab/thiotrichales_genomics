rm(list = ls())
library(ggplot2)
library(cowplot)

setwd("/Users/monicareyes/Desktop/input_grafica/grafica/")

archivo_ids  <- "ID.txt"
archivo_completeness <- "input_Completeness.csv"
archivo_assembly <- "input_genome_size.csv"
archivo_cds  <- "input_CDS_count.csv"
archivo_bgcs <- "input_bgc.csv"
archivo_genes <- "vecindarios_detalle.csv"  
archivo_pares_bgc <- "pb_d_bgc.csv"

ruta_salida <- "grafica_completa.png"

genome_ids <- trimws(readLines(archivo_ids))
genome_ids <- genome_ids[nzchar(genome_ids)]
n_muestras <- length(genome_ids)

cat("Filas:", n_muestras, "\n")

archivo_valido <- function(archivo) {
  !is.na(archivo) && nzchar(archivo) && file.exists(archivo) && file.size(archivo) > 0
}

cargar_columna <- function(archivo, n_esperado, nombre) {
  if (!archivo_valido(archivo)) return(NULL)
  valores <- read.csv(archivo, header = FALSE)$V1
  if (length(valores) != n_esperado) {
    warning(nombre, ": tiene ", length(valores), " filas, se esperaban ", n_esperado, ". Se omite.")
    return(NULL)
  }
  valores
}

completeness <- cargar_columna(archivo_completeness, n_muestras, "completeness")
assembly <- cargar_columna(archivo_assembly, n_muestras, "assembly size")
cds   <- cargar_columna(archivo_cds, n_muestras, "CDS count")
bgcs     <- cargar_columna(archivo_bgcs, n_muestras, "BGCs")
pb_bgcs <- NULL
if (archivo_valido(archivo_pares_bgc)) {
  pb_raw <- read.delim(archivo_pares_bgc, stringsAsFactors = FALSE,
                       colClasses = c(genoma = "character"))
  pb_raw$genoma <- trimws(pb_raw$genoma)
  sin_bgc <- setdiff(genome_ids, pb_raw$genoma)
  if (length(sin_bgc) > 0) {
    warning("estos genomas no aparecen en el archivo de pares de bases de BGCs: ",
            paste(sin_bgc, collapse = ", "))
  }
  pb_bgcs <- pb_raw$total_bp_bgcs[match(genome_ids, pb_raw$genoma)]
}


genes_pattern1 <- c("murG", "murC", "murB", "ddl", "ftsZ", "ftsQ", "lpxC", "ftsA")  
genes_pattern2 <- c("mreD", "mrdA", "rodZ")                                         

if (!archivo_valido(archivo_genes)) {
  stop("no encuentro el csv de genes '", archivo_genes, "' en ", getwd())
}
pattern <- NULL
if (archivo_valido(archivo_genes)) {
  genes_long <- read.csv(archivo_genes, stringsAsFactors = FALSE,
                         colClasses = c(genome_id = "character"))
  genes_long$genome_id <- trimws(genes_long$genome_id)
  
  sin_datos <- setdiff(genome_ids, genes_long$genome_id)
  if (length(sin_datos) > 0) {
    warning("estos genomas de ID.txt no aparecen en el csv de genes: ",
            paste(sin_datos, collapse = ", "))
  }
  
  # el script de vecindarios numera las copias (mreD1, mreD2): se quita el numero final
  genes_long$gen_base <- sub("[0-9]+$", "", genes_long$gen)
  
  faltan_genes <- setdiff(c(genes_pattern1, genes_pattern2), unique(genes_long$gen_base))
  if (length(faltan_genes) > 0) {
    warning("estos genes no aparecen en el csv y se cuentan como ausentes (revisa la lista `genes` del script de vecindarios): ",
            paste(faltan_genes, collapse = ", "), call. = FALSE)
  }
  
  # un mismo gen puede salir en varios paneles: se cuenta por locus_tag unico
  genes_unicos <- unique(genes_long[, c("genome_id", "locus_tag", "gen_base")])
  
  copias <- as.matrix(table(
    factor(genes_unicos$genome_id, levels = genome_ids),
    factor(genes_unicos$gen_base, levels = c(genes_pattern1, genes_pattern2))
  ))
  
  pattern1 <- ifelse(rowSums(copias[, genes_pattern1, drop = FALSE] >= 1) == length(genes_pattern1), "+", "-")
  pattern2 <- ifelse(rowSums(copias[, genes_pattern2, drop = FALSE] >= 2) == length(genes_pattern2), "+", "-")
  
  pattern <- data.frame(V1 = unname(pattern1), V2 = unname(pattern2))
  stopifnot(nrow(pattern) == n_muestras)
}

colores <- list(assembly = '#5f8a8b', cds = '#8b7355', bgcs = '#c96a9e')

titulos <- list(
  id_completeness = 'genome ID   completeness [%]',
  assembly = 'assembly size\n[Mbp]',
  cds = 'CDS count\n(x1000)',
  bgcs = 'BGCs\n[%]',
  pb_bgcs = 'BGCs\n[bp]',
  pattern = 'pattern'
)

plots  <- list()
anchos <- c()

if (!is.null(completeness)) {
  df_plot <- data.frame(y = 1:n_muestras, id = genome_ids, comp = round(completeness, 1))
  p1 <- ggplot(df_plot, aes(y = y)) +
    geom_text(aes(x = 0, label = id), hjust = 1, vjust = 0.5, size = 3.3) +
    geom_text(aes(x = 0.15, label = sprintf("%.1f", comp)),
              hjust = 0, vjust = 0.5, size = 4, fontface = "bold") +
    scale_y_reverse(limits = c(n_muestras + 0.5, 0.5)) +
    xlim(-1, 1) +
    labs(title = titulos$id_completeness) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
          plot.margin = unit(c(0.5, 0.2, 1.5, 0.2), "cm"))
} else {
  df_plot <- data.frame(y = 1:n_muestras, id = genome_ids)
  p1 <- ggplot(df_plot, aes(y = y)) +
    geom_text(aes(x = 0, label = id), hjust = 0.5, vjust = 0.5, size = 3.3) +
    scale_y_reverse(limits = c(n_muestras + 0.5, 0.5)) +
    xlim(-1, 1) +
    labs(title = "genome ID") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
          plot.margin = unit(c(0.5, 0.2, 1.5, 0.2), "cm"))
}
plots <- append(plots, list(p1)); anchos <- c(anchos, 1.5)

agregar_barras <- function(valores, titulo, color, x_max, x_by) {
  if (is.null(valores)) return(NULL)
  df_plot <- data.frame(y = 1:n_muestras, valor = valores)
  ggplot(df_plot, aes(x = valor, y = y)) +
    geom_col(fill = color, color = "white", width = 0.8, linewidth = 0.5, orientation = "y") +
    scale_y_reverse(limits = c(n_muestras + 0.5, 0.5)) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, by = x_by)) +
    labs(title = titulo, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linetype = "dashed", color = "gray", linewidth = 0.3),
      plot.margin = unit(c(0.5, 0.2, 1.5, 0.2), "cm")
    )
}

p2 <- agregar_barras(if (!is.null(assembly)) assembly / 1e6 else NULL, titulos$assembly, colores$assembly, 15, 5)
if (!is.null(p2)) { plots <- append(plots, list(p2)); anchos <- c(anchos, 1.2) }

cds_k <- if (!is.null(cds)) { if (max(cds, na.rm = TRUE) > 1000) cds / 1000 else cds } else NULL
p3 <- agregar_barras(cds_k, titulos$cds, colores$cds, 15, 5)
if (!is.null(p3)) { plots <- append(plots, list(p3)); anchos <- c(anchos, 1.2) }

p4 <- agregar_barras(bgcs, titulos$bgcs, colores$bgcs, 35, 10)
if (!is.null(p4)) { plots <- append(plots, list(p4)); anchos <- c(anchos, 1.2) }

if (!is.null(pb_bgcs)) {
  df_plot <- data.frame(
    y = 1:n_muestras,
    valor = format(round(pb_bgcs), big.mark = ",", scientific = FALSE, trim = TRUE)
  )
  p4b <- ggplot(df_plot, aes(x = 0, y = y)) +
    geom_text(aes(label = valor), hjust = 0.5, vjust = 0.5, size = 4, fontface = "bold") +
    scale_y_reverse(limits = c(n_muestras + 0.5, 0.5)) +
    xlim(-1, 1) +
    labs(title = titulos$pb_bgcs) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.margin = unit(c(0.5, 0.2, 1.5, 0.2), "cm"))
  plots <- append(plots, list(p4b)); anchos <- c(anchos, 1.0)
}

if (!is.null(pattern)) {
  df_plot <- data.frame(
    y = 1:nrow(pattern),
    pattern1 = as.character(pattern$V1),
    pattern2 = as.character(pattern$V2)
  )
  p5 <- ggplot(df_plot, aes(x = 0, y = y)) +
    geom_text(aes(x = 0.3, label = pattern1), hjust = 0.5, vjust = 0.5,
              size = 6, family = "mono", fontface = "bold") +
    geom_text(aes(x = 0.7, label = pattern2), hjust = 0.5, vjust = 0.5,
              size = 6, family = "mono", fontface = "bold") +
    scale_y_reverse(limits = c(n_muestras + 0.5, 0.5)) +
    xlim(0, 1) +
    labs(title = titulos$pattern) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.margin = unit(c(0.5, 0.2, 1.5, 0.2), "cm"))
  plots <- append(plots, list(p5)); anchos <- c(anchos, 0.8)
}

combined_plot <- plot_grid(
  plotlist = plots,
  ncol = length(plots),
  align = "h",
  axis = "tb",
  rel_widths = anchos
)

print(combined_plot)

ggsave(ruta_salida, combined_plot,
       width = sum(anchos) * 2.2, height = max(6, n_muestras * 0.35),
       units = "in", dpi = 300, bg = "white")
