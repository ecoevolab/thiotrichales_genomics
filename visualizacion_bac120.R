# Librerias
library(ggtree)
library(treeio)
library(ape)
library(here)
library(dplyr)

metadata <- read.csv("metadata.csv")
metadata <- metadata %>% rename(label = Genome.ID)
tree <- read.tree("/Users/monicareyes/Desktop/Thiotrichales_bac120_IQtree/gtdbtk.treefile")

setdiff(metadata$label, tree$tip.label)
setdiff(tree$tip.label, metadata$label)

p1 <- ggtree(tree) %<+% metadata 

# Agregar los nombres de las bacterias
# Colorear las puntas por el tipo de ecosistema que viven las bacterias
p1 +
  geom_tiplab(size = 2, aes (label = Organism.name..putative.)) +
  geom_tippoint(aes(color = Ecosystem.type), size = 2) +
  theme(legend.position = "right")
