library(ape)
library(phytools)

data <- read.csv("pollen_data.csv", na.strings = "", 
                 stringsAsFactors = F)
tree <- read.nexus("pip_group.tre")

#======#
# TREE #
#======#

tab.not.tree <- data$name_phylogeny[!data$name_phylogeny %in% tree$tip.label]

tree <- drop.tip(tree, setdiff(tree$tip.label, data$name_phylogeny))

node_descendants <- tree$tip.label[getDescendants(tree, getMRCA(tree, c("Mimosa_oedoclada", "Mimosa_dichroa")))]

tree <- drop.tip(tree, node_descendants)

# Setting all branch lengths equal to one

tree <- compute.brlen(tree, 1)

write.nexus(tree, file = "pruned_tree.nex")

# Subtree

tree_subset <- extract.clade(tree, getMRCA(tree, c("Lachesiodendron_viridiflorum", "Microlobius_foetidus")))
write.nexus(tree_subset, file = "tree_subset.nex")

#=====#
# CAT #
#=====#

data_pruned <- data[data$name_phylogeny %in% tree$tip.label, ]

data_cat <- data_pruned[, c(6, 12:17)]

write.csv(data_cat, "data_cat.csv", row.names = F)

#======#
# CONT #
#======#

data_cont <- data_pruned[, c(6, 18:30)]

write.csv(data_cont, "data_cont.csv", row.names = F)

