if (!require(cowplot)) install.packages("cowplot", repos='http://cran.rstudio.com/')
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)

archivo_completeness <- "/mnt/data/sur/users/mreyes/data/beggiatoas/results/checkm/storage/Completeness_input.csv"
archivo_assembly     <- "/mnt/data/sur/users/mreyes/data/beggiatoas/results/checkm/storage/Genomesize_input.csv"
archivo_cds          <- "/mnt/data/sur/users/mreyes/data/beggiatoas/results/checkm/storage/CDS_input.csv"
archivo_bgcs         <- "/mnt/data/sur/users/mreyes/data/beggiatoas/results/checkm/storage/BGC_percentage.csv"
archivo_pattern      <- "/mnt/data/sur/users/mreyes/data/beggiatoas/results/checkm/storage/patterns_input.csv"

ruta_salida <- "/mnt/data/sur/users/mreyes/data/beggiatoas/data_processed/test_with_selected_genomes/csv_extraccion/grafica_completa.png"

completeness <- read.csv(archivo_completeness, header = FALSE)$V1
assembly     <- read.csv(archivo_assembly, header = FALSE)$V1
cds          <- read.csv(archivo_cds, header = FALSE)$V1
bgcs         <- read.csv(archivo_bgcs, header = FALSE)$V1 
pattern      <- read.csv(archivo_pattern, header = FALSE)

n_muestras <- length(completeness)

cat("Datos cargados:\n")
cat("  Muestras:", n_muestras, "\n")
cat("  Assembly: min =", round(min(assembly)/1e6, 2), "max =", round(max(assembly)/1e6, 2), "Mbp\n")
cat("  CDS: min =", round(min(cds)/1000, 2), "max =", round(max(cds)/1000, 2), "k\n")
cat("  BGCs: min =", round(min(bgcs), 2), "max =", round(max(bgcs), 2), "%\n\n")

colores <- list(
  assembly = '#5f8a8b',
  cds = '#8b7355',
  bgcs = '#c96a9e'
)

titulos <- list(
  completeness = 'completeness\n[%]',
  assembly = 'assembly size\n[Mbp]',
  cds = 'CDS count\n(x1000)',
  bgcs = 'BGCs\n[%]',
  pattern = 'pattern'
)

plots <- list()

df_plot <- data.frame(y = 1:n_muestras, valor = completeness)

p1 <- ggplot(df_plot, aes(x = 0, y = y)) +
  geom_text(aes(label = sprintf("%.1f", valor)), 
            hjust = 0.5, vjust = 0.5, size = 4, fontface = "bold") +
  scale_y_reverse(limits = c(n_muestras + 0.5, 0.5)) +
  xlim(-1, 1) +
  labs(title = titulos$completeness) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.margin = unit(c(0.5, 0.2, 1.5, 0.2), "cm")
  )

plots <- append(plots, list(p1))

assembly_mbp <- assembly / 1e6
df_plot <- data.frame(y = 1:n_muestras, valor = assembly_mbp)

p2 <- ggplot(df_plot, aes(x = valor, y = y)) +
  geom_col(fill = colores$assembly, color = "white", width = 0.8, linewidth = 0.5, orientation = "y") +
  scale_y_reverse(limits = c(n_muestras + 0.5, 0.5)) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5)) +
  labs(title = titulos$assembly, x = NULL, y = NULL) +
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

plots <- append(plots, list(p2))

cds_dividido <- cds / 1000
df_plot <- data.frame(y = 1:n_muestras, valor = cds_dividido)

p3 <- ggplot(df_plot, aes(x = valor, y = y)) +
  geom_col(fill = colores$cds, color = "white", width = 0.8, linewidth = 0.5, orientation = "y") +
  scale_y_reverse(limits = c(n_muestras + 0.5, 0.5)) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5)) +
  labs(title = titulos$cds, x = NULL, y = NULL) +
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

plots <- append(plots, list(p3))

df_plot <- data.frame(y = 1:n_muestras, valor = bgcs)

p4 <- ggplot(df_plot, aes(x = valor, y = y)) +
  geom_col(fill = colores$bgcs, color = "white", width = 0.8, linewidth = 0.5, orientation = "y") +
  scale_y_reverse(limits = c(n_muestras + 0.5, 0.5)) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5)) +
  labs(title = titulos$bgcs, x = NULL, y = NULL) +
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

plots <- append(plots, list(p4))

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
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.margin = unit(c(0.5, 0.2, 1.5, 0.2), "cm")
  )

plots <- append(plots, list(p5))

combined_plot <- plot_grid(
  plotlist = plots, 
  ncol = 5, 
  align = "v", 
  axis = "lr", 
  rel_widths = c(0.8, 1.2, 1.2, 1.2, 0.8)
)

ggsave(ruta_salida, combined_plot, 
       width = 20, height = 12, units = "in", dpi = 300,
       bg = "white")

cat("\nGrafica guardada en:", ruta_salida, "\n")
