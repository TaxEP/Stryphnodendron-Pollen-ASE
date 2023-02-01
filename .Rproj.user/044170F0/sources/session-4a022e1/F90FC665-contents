#==============================================================================#

# Loading phytools and dependencies
require(phytools)

#==============================================================================#

#====================#
# CATEGORICAL TRAITS #
#====================#

# Loading tree
stryphnod.tree <- read.nexus("pruned_tree.nex")

# Loading character states (categorical traits must be formatted as a 
# data.frame, with tip labels as row names and trait labels as column names)
traits <- read.csv("data_cat.csv", header = TRUE, row.names = 1,
                   na.strings = "")

# The data.frame needs to be transformed into a named vector
# containing tip labels as names and character states as factors.
plyr::count(traits$SEM_Ornamentation)
cat.trait <- factor(traits$SEM_Ornamentation, levels = c("areolate", "psilate", "rugulate", "verrucate", "verrucate-scabrate", "fossulate", "reticulate",
                                                         "areolate_&_fossulate", "areolate_&_psilate", "areolate_&_reticulate", "areolate_&_rugulate", "areolate_&_verrucate", "psilate_&_verrucate",
                                                         "missing"))
cat.trait[is.na(cat.trait)] <- "missing"
names(cat.trait) <- rownames(traits)

# The named vector, then, needs to be transformed into a binary matrix
# with tip labels as row names and character states as columns
bin.matrix <- to.matrix(cat.trait, levels(cat.trait))

# Dealing with polymorphisms
## When taxa with polymorphisms are present, we must assign a value
# equal to 1/(number of states) for each possible state and then remove 
# the column that refers to the polymorphic state. 

# Polymorphism areolate_&_fossulate
polymorph8 <- rownames(bin.matrix)[bin.matrix[ , c(8)] == 1]
bin.matrix[polymorph8, c(1, 6)] <- 1/2

# Polymorphism areolate_&_psilate
polymorph9 <- rownames(bin.matrix)[bin.matrix[ , c(9)] == 1]
bin.matrix[polymorph9, c(1, 2)] <- 1/2

# Polymorphism areolate_&_reticulate
polymorph10 <- rownames(bin.matrix)[bin.matrix[ , c(10)] == 1]
bin.matrix[polymorph10, c(1, 7)] <- 1/2

# Polymorphism areolate_&_rugulate
polymorph11 <- rownames(bin.matrix)[bin.matrix[ , c(11)] == 1]
bin.matrix[polymorph11, c(1, 3)] <- 1/2

# Polymorphism areolate_&_verrucate
polymorph12 <- rownames(bin.matrix)[bin.matrix[ , c(12)] == 1]
bin.matrix[polymorph12, c(1, 4)] <- 1/2

# Polymorphism psilate_&_verrucate
polymorph13 <- rownames(bin.matrix)[bin.matrix[ , c(13)] == 1]
bin.matrix[polymorph13, c(2, 4)] <- 1/2

# Removing columns indicating polymorphism
bin.matrix <- bin.matrix[ , -c(8, 9, 10, 11, 12, 13)]

# Which taxa show missing data? 
missing.data <- names(which(bin.matrix[ , "missing"] == 1))

# For missing data, we must assign a prior probability distribution on the tips 
# that is flat across all possible states.
## Removing the column 'missing'
bin.matrix <- bin.matrix[ , -c(ncol(bin.matrix))]
## Assigning 1/(number of states) for all possible states
bin.matrix[row.names(bin.matrix) %in% missing.data, ] <- 1/ncol(bin.matrix)

# putting the columns in the same order of the tip.labels
bin.matrix[match(stryphnod.tree$tip.label, rownames(bin.matrix)),] -> bin.matrix

# Plotting
cols <- setNames(palette()[1:length(colnames(bin.matrix))], colnames(bin.matrix))
plotTree(stryphnod.tree, fsize = 0.7, lwd = 1, ftype="i", offset = 0.5, type = "phylogram")
tiplabels(pie=bin.matrix,piecol=cols,
          cex=0.4)
legend(x="topleft", legend=colnames(bin.matrix),pt.cex=2,cex=0.9,pch=21,
       pt.bg=cols)


# Saving as pdf
pdf("all_analyses/cat/SEM_Ornamentation/plot_polymorph.pdf")
plotTree(stryphnod.tree, fsize = 0.7, lwd = 1, ftype="i", offset = 0.5, type = "phylogram")
tiplabels(pie=bin.matrix,piecol=cols,
          cex=0.4)
dev.off()
