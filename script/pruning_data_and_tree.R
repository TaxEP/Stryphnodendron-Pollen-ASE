library(ape)
library(phytools)

data <- read.csv("data/pollen_data.csv", na.strings = "", 
                 stringsAsFactors = F)
tree <- read.nexus("data/pip_group.tre")

#======#
# TREE #
#======#

tree <- drop.tip(tree, setdiff(tree$tip.label, data$name_phylogeny))

# Setting all branch lengths equal to one

tree <- compute.brlen(tree, 1)

tree <- force.ultrametric(tree, method = "extend")

write.nexus(tree, file = "output/data/pruned_tree.nex")

# Subtree

tree_subset <- extract.clade(tree, getMRCA(tree, c("Lachesiodendron_viridiflorum", "Microlobius_foetidus")))
write.nexus(tree_subset, file = "output/data/tree_subset.nex")

#=====#
# CAT #
#=====#

data_pruned <- data[data$name_phylogeny %in% tree$tip.label, ]

data_cat <- data_pruned[, c(4, 9:11)]

write.csv(data_cat, "output/data/data_cat.csv", row.names = F)

#======#
# CONT #
#======#

data_cont <- data_pruned[, c(4, 12:22)]

write.csv(data_cont, "output/data/data_cont.csv", row.names = F)

