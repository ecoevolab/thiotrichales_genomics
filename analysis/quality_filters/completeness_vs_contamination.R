#Limpiar memoria
rm(list = ls())

# Cargar librerias necesarias
library(ggplot2)

# Leer el csv combinado 
checkm_df <- read.table("/Users/monicareyes/Desktop/checkm_input_all.csv", header = FALSE, sep = "",
                        col.names = c("Completeness", "Contamination", "Genome_size"))

# Asegurar que las columnas sean numéricas
checkm_df$Completeness  <- as.numeric(checkm_df$Completeness)
checkm_df$Contamination <- as.numeric(checkm_df$Contamination)
checkm_df$Genome_size   <- as.numeric(checkm_df$Genome_size)

# Score de calidad estándar (CheckM/GTDB-Tk)
checkm_df$quality_score <- checkm_df$Completeness - 5 * checkm_df$Contamination

# Categorías de calidad basadas en ese score (puedes ajustar los cortes después de ver el gráfico)
checkm_df$categoria <- cut(checkm_df$quality_score,
                           breaks = c(-Inf, 25, 50, 75, Inf),
                           labels = c("low", "medium", "good", "high"))

# Scatter plot: Completeness vs Contamination - complet
ggplot(checkm_df, aes(x = Completeness, y = Contamination,
                      color = categoria, size = Genome_size)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(50, 90), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = c(5, 10), linetype = "dashed", color = "gray40") +
  labs(title = "Completeness vs Contaminación (CheckM) - without pre-filter",
       x = "Completeness (%)", y = "Contaminación (%)",
       color = "Category", size = "Genome size (pb)") +
  theme_minimal()

# Opcional guardar en la memoria
ggsave("/Users/monicareyes/Desktop/completeness_vs_contamination_completo.png", width = 8, height = 6, dpi = 300)



# Scatter plot: Completeness vs Contamination zoom
ggplot(checkm_df, aes(x = Completeness, y = Contamination,
                      color = categoria, size = Genome_size)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(50, 90), linetype = "dashed", color = "gray40") + # lineas punteadas en el eje x
  geom_hline(yintercept = c(5, 10), linetype = "dashed", color = "gray40") + # lineas punteadas en el eje y
  coord_cartesian(ylim = c(0, 30)) + # Aqui puedes ajustar el zoom en este momento esta de 0 30 en el y
  coord_cartesian(ylim = c(0, 30)) + # zoom en el x 
  labs(title = "Completeness vs Contaminación (zoom 0-30%)",
       x = "Completeness (%)", y = "Contaminación (%)",
       color = "Category", size = "Genome size (pb)") +
  theme_minimal()

# Opcional guardar en la memoria
ggsave("/Users/monicareyes/Desktop/completeness_vs_contamination_zoom.png", width = 8, height = 6, dpi = 300)

