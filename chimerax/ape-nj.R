# ape - Neighbor-Joining
library(ape)
nboot <- 1000
dist_matrix <- data$`Alignment score`
bootstrap_trees <- boot.phylo(x = dist_matrix, B = nboot)