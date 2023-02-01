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
plyr::count(traits$Outline)
cat.trait <- factor(traits$Outline, levels = c("circular", "elliptical", "oval", 
                                               "circular_&_elliptical", "elliptical_&_oval", 
                                               "oval_&_circular_&_elliptical",
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

# Polymorphism 1_&_2
polymorph4 <- rownames(bin.matrix)[bin.matrix[ , c(4)] == 1]
bin.matrix[polymorph4, c(1, 2)] <- 1/2

# Polymorphism 2_&_3
polymorph5 <- rownames(bin.matrix)[bin.matrix[ , c(5)] == 1]
bin.matrix[polymorph5, c(2, 3)] <- 1/2

# Polymorphism 1_&_2_&_3
polymorph6 <- rownames(bin.matrix)[bin.matrix[ , c(6)] == 1]
bin.matrix[polymorph6, c(1, 2, 3)] <- 1/3

# Removing columns indicating polymorphism
bin.matrix <- bin.matrix[ , -c(4:6)]

# Which taxa show missing data? 
missing.data <- names(which(bin.matrix[ , "missing"] == 1))

# For missing data, we must assign a prior probability distribution on the tips 
# that is flat across all possible states.
## Removing the column 'missing'
bin.matrix <- bin.matrix[ , -c(ncol(bin.matrix))]
## Assigning 1/(number of states) for all possible states
bin.matrix[row.names(bin.matrix) %in% missing.data, ] <- 1/ncol(bin.matrix) 

# Running stochastic mapping (this can take a while)
trees <- make.simmap(stryphnod.tree, bin.matrix, model = "ER", nsim = 100)
obj <- summary(trees, plot = FALSE)

# Plotting
cols <- setNames(palette()[1:length(colnames(bin.matrix))], colnames(bin.matrix))
plot(obj, colors = cols, fsize = 0.7, cex = c(0.35,0.4), lwd = 1, ftype="i", offset = 0.5, type = "phylogram")
add.simmap.legend(colors = cols, x = 0.1, y = 37, prompt = FALSE, fsize=0.9)

# Saving as pdf
pdf("all_analyses/cat/Outline/plot.pdf")
plot(obj, colors = cols, fsize = 0.7, cex = c(0.35,0.4), lwd = 1, ftype="i", offset = 0.5, type = "phylogram")
add.simmap.legend(colors = cols, x = 0.02, y = 37, prompt = FALSE, fsize=0.9)
dev.off()

