# limpiar memoria
rm(list = ls())
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

#ruta
OG <- read_tsv("/Users/monicareyes/Desktop/Orthogroups/Orthogroups.GeneCount.tsv")


# Hacemos un archivo copia donde se guarden solamente las columnas de los genomas que nos interesa
conservar <- OG %>% 
  select("Orthogroup","2894041767", "2890541225", "2890630416", "2894040609", "2890631485", "2890616370", "2772190905","2711768587", "8074213743", "2772190906",
         "2721755989", "2772190729", "2806310698", "2711768588", "2710724217", "2767802434","2835562463", "2835554630", "8002131595","2684622566",
         "2556921650","2828901765","2619618925","2710724216","2514752027","2835551090","2767802425","2599185268","2835558653","3001836915", "2767802433",
         "2585428147","2828749676","2508501047","2788500402","2654588010","8002135289","8066494198","8056884285", "8125999560", "2802428810", "2506520049",
         "2882495133","2502790011","2619619276","2889727904","2517093012","2236661048","8071263300","640963011","2600255314","2788499801","2772190875",
         "2955952073","2956002012","2955963837","2955976453", "2955989171", "ralstonia_solanacearum")

# guardamos el archivo copia
write_tsv(conservar, "OG_HC_genomes.tsv")


OG_HC <- read_tsv("OG_HC_genomes.tsv")
# Comenzamos a seleccionar las filas que tengan 1 o menos copias en cada genoma
cols <- setdiff(names(OG_HC), c("Orthogroup", "Total"))
filtrado <- OG_HC %>%
  filter(if_all(all_of(cols), ~ .x <= 1)) %>%
  mutate(Total = rowSums(across(all_of(cols)))) #agrega una columna extra llamada total y se cuenta el número de genomas en el cual esta presente cada copia. 

# guardamos la copia
write_tsv(filtrado, "OG_single_copy.tsv")
nrow(OG_HC); nrow(filtrado) #comparamos con cuantos OG comenzamos y con cuantos nos quedamos

presencia <- read_tsv("OG_single_copy.tsv")
# El numero de copias puede cambiar dependiendo el numero de porcentaje que se ponga 
p_copias <- ceiling(0.8 * length(cols))
filtrado_80_per <- presencia %>% 
  filter(Total >= p_copias)
# guardamos los datos
write_tsv(filtrado_80_per, "OG_80_percent_in_genomes.tsv")
# Verificamos cuantos OG quedaron 
nrow(filtrado_80_per)  

# Histograma para visualizar en cuantos genomas esta presente cada OG 
# ruta del archivo de arriba
OG_80 <- read_tsv("OG_80_percent_in_genomes.tsv")

OG_freq <- OG_80 %>% 
  mutate(n_genomes = rowSums(across(all_of(cols)) > 0)) %>% 
  select(Orthogroup, n_genomes) %>% 
  arrange(desc(n_genomes)) %>% 
  mutate(Orthogroup = factor(Orthogroup, levels = Orthogroup))

ggplot(OG_freq, aes(x = Orthogroup, y = n_genomes, fill = n_genomes)) +
  geom_col() +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Genomas") + #degradado de colores
  scale_y_continuous(limits = c(0, 59), breaks = seq(0, 59, 5)) +
  labs(x = "Orthogroup", y = "Número de genomas") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4))

# histograma para visualizar cuantos OG tiene cada genoma
genoma_freq <- OG_80 %>% 
  summarise(across(all_of(cols), ~ sum(.x > 0))) %>% 
  pivot_longer(everything(), names_to = "genoma", values_to = "n_OG") %>% 
  arrange(desc(n_OG)) %>% 
  mutate(genoma = factor(genoma, levels = genoma))

ggplot(genoma_freq, aes(x = genoma, y = n_OG, fill = n_OG)) +
  geom_col() +
  geom_hline(yintercept = nrow(OG_80), linetype = "dashed") +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "OG") + #degradado de colores
  scale_y_continuous(limits = c(0, nrow(OG_80)), breaks = seq(0, nrow(OG_80), 25)) +
  labs(x = "Genoma", y = paste0("Número de OG conservados (Max en ", nrow(OG_80), ")")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))



