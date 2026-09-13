# Librerias
library(readr)
library(here) #define la ruta en la que estas trabajando - se toma el here muy literal
library(stringr)
library(dplyr)


# El objetivo de este script es convertir una tabla de meta datos en formato Rmd a un archivo csv 

#-------------------------------------------------------------------------------------------------

# Input:
# Una tabla en formado Rmd con el siguiente formatos:
# |Genome ID | Organism name (putative) | Ubication | Ecosystem type | Sequencing Strategy| Study name
# |:----|:----|:----|:----|:----|:----|
# |2894041767 | Thiothrix sp. AG-917-F21 | Atlantic Ocean  | Marine water | Whole Genome Sequencing | Charting the complexity 

# basicamente una tabla separada por "|", si el segundo reglon no se tiene puede omitirse la linea de la variable table_lines


# Output:
# un archivo csv (separado por ",") del archivo input.

#-------------------------------------------------------------------------------------------------

lines <- readLines("summary_table.Rmd") # Solamente se necesita ajustar la ruta del archivo Rmd

# extrae las lineas que empiezan con "|"
table_lines <- lines[str_detect(lines, "^\\s*\\|")]

# elimina la fila separadora tipo |:----|:----|
table_lines <- table_lines[!str_detect(table_lines, "^[|:\\-\\s]+$")]

# Ajuste del formato
meta <- read_delim(paste(table_lines, collapse = "\n"), delim = "|", trim_ws = TRUE)

# Eliminacion de las lineas vacias del inviio y el final 
meta <- meta %>% select(where(~ !all(is.na(.)) && !all(. == "", na.rm = TRUE)))

names(meta) <- trimws(names(meta))

# Guardado el archivo en la ruta
write_csv(meta, "metadata.csv")
