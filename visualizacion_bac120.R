# Librerias
library(ggtree)
library(treeio)

tree <- read.tree("/Users/monicareyes/Desktop/Thiotrichales_bac120_IQtree/gtdbtk.treefile")
p1 <- ggtree(tree) + 
  geom_tiplab(size=2) + 
  geom_nodelab(aes(label=label), size=1.5)


ggsave("test.png", width = 5, height = 8)
