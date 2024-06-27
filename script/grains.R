#==============================================================================#

# Loading phytools and dependencies
require(phytools)

#==============================================================================#

#====================#
# CATEGORICAL TRAITS #
#====================#

# Loading tree
stryphnod.tree <- read.nexus("output/data/pruned_tree.nex")

# Loading character states (categorical traits must be formatted as a 
# data.frame, with tip labels as row names and trait labels as column names)
traits <- read.csv("output/data/data_cat.csv", header = TRUE, row.names = 1,
                   na.strings = "")

## Testing evolutionary models ------------------------------------------------

# Seting seed for replicability
set.seed(7)

# adjusting data for fitting
traits_test <- as.data.frame(lapply(traits, function(x) gsub("&", "+", x)))
states <- traits_test$grains
names(states) <- rownames(traits)

# Fitting models

# equal rates
unordered.er <- fitpolyMk(stryphnod.tree,states,model="ER")
# symmetrical rates
unordered.sym <- fitpolyMk(stryphnod.tree,states,model="SYM")
# all rates different
unordered.ard <- fitpolyMk(stryphnod.tree,states,model="ARD")

# Evaluating results

# equal rates
round(AIC(unordered.er),2)
# symmetrical rates
round(AIC(unordered.sym),2)
# all rates different
round(AIC(unordered.ard),2)

# ER - 181.27 ; SYM - 187.87 ; ARD - 218.73

## Preparing data -------------------------------------------------------------

# The data.frame needs to be transformed into a named vector
# containing tip labels as names and character states as factors.
plyr::count(traits$grains)
cat.trait <- factor(traits$grains, levels = c("8", "12", "16", "32", "8&12",
                                              "8&16", "12&16", "missing"))
cat.trait[is.na(cat.trait)] <- "missing"
names(cat.trait) <- rownames(traits)

# The named vector, then, needs to be transformed into a binary matrix
# with tip labels as row names and character states as columns
bin.matrix <- to.matrix(cat.trait, levels(cat.trait))

# Dealing with polymorphisms
## When taxa with polymorphisms are present, we must assign a value
# equal to 1/(number of states) for each possible state and then remove 
# the column that refers to the polymorphic state. 

# Polymorphism 8&12
polymorph5 <- rownames(bin.matrix)[bin.matrix[ , c(5)] == 1]
bin.matrix[polymorph5, c(1, 2)] <- 1/2

# Polymorphism 8&16
polymorph6 <- rownames(bin.matrix)[bin.matrix[ , c(6)] == 1]
bin.matrix[polymorph6, c(1, 3)] <- 1/2

# Polymorphism 12&16
polymorph7 <- rownames(bin.matrix)[bin.matrix[ , c(7)] == 1]
bin.matrix[polymorph7, c(2, 3)] <- 1/2

# Removing columns indicating polymorphism
bin.matrix <- bin.matrix[ , -c(5,6,7)]

# Which taxa show missing data? 
missing.data <- names(which(bin.matrix[ , "missing"] == 1))

# For missing data, we must assign a prior probability distribution on the tips 
# that is flat across all possible states.
## Removing the column 'missing'
bin.matrix <- bin.matrix[ , -c(ncol(bin.matrix))]
## Assigning 1/(number of states) for all possible states
bin.matrix[row.names(bin.matrix) %in% missing.data, ] <- 1/ncol(bin.matrix) 

## SIMMAP --------------------------------------------------------------------

# Running stochastic mapping (this can take a while)
trees <- make.simmap(stryphnod.tree, bin.matrix, model = "ER", nsim = 100)
obj <- summary(trees, plot = FALSE)

# putting the columns in the same order of the tip.labels, so we can plot polymophisms
bin.matrix <- bin.matrix[match(stryphnod.tree$tip.label, rownames(bin.matrix)),]

# putting the columns of obj$ace in the same order of bin.matrix
obj$ace <- obj$ace[, colnames(bin.matrix)]

# plotting and saving
cols <- setNames(palette()[1:length(colnames(bin.matrix))], colnames(bin.matrix))

pdf("output/plots/grains.pdf")
par(lwd = 0.1)
plotTree(stryphnod.tree,
         type = "phylogram",
         fsize = 0.7, 
         offset = 0.5,
         cex = 0.2, 
         lwd = 1, 
         ftype = "i")
nodelabels(pie = obj$ace,
           piecol = cols, 
           cex=0.5)
tiplabels(pie = bin.matrix,
          piecol = cols,
          cex=0.45)
add.simmap.legend(colors = cols, x = 0, y = 5, prompt = FALSE, fsize=0.5)
dev.off()

