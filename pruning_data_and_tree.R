library(ape)

#=====#
# CAT #
#=====#

data <- read.csv("piptadenia_pollen - cat.csv", na.strings = "", 
                 stringsAsFactors = F)
tree <- read.tree("stryphnod_clean_updated.tre")

species.names.tab <- as.vector(data$name_phylogeny)
species.names.tree <- as.vector(tree$tip.label)
tab.not.tree <- data$name_phylogeny[!data$name_phylogeny %in% tree$tip.label]

data2 <- data[data$name_phylogeny %in% species.names.tree, ]

write.csv(data2, "data_cat.csv", row.names = F)

#======#
# CONT #
#======#

data <- read.csv("piptadenia_pollen - cont.csv", na.strings = "", 
                 stringsAsFactors = F)
tree <- read.tree("stryphnod_clean_updated.tre")

species.names.tab <- as.vector(data$name_phylogeny)
species.names.tree <- as.vector(tree$tip.label)
tab.not.tree <- data$name_phylogeny[!data$name_phylogeny %in% tree$tip.label]

data2 <- data[data$name_phylogeny %in% species.names.tree, ]

tree2 <- drop.tip(tree, setdiff(species.names.tree, species.names.tab))

write.csv(data2, "data_cont.csv", row.names = F)

# Setting all branch lengths equal to one

tree2 <- compute.brlen(tree2, 1)

write.nexus(tree2, file = "pruned_tree.nex")

