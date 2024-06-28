#==============================================================================#

# Loading phytools and dependencies
require(phytools)

#==============================================================================#

#====================#
# CATEGORICAL TRAITS #
#====================#

# Loading tree
phy <- read.nexus("output/data/pruned_tree.nex")

# Loading character states (categorical traits must be formatted as a 
# data.frame, with tip labels as row names and trait labels as column names)
traits <- read.csv("output/data/data_cat.csv", header = TRUE, row.names = 1)

## Testing evolutionary models ------------------------------------------------

# Seting seed for replicability
set.seed(7)

# adjusting data for fitting
traits_test <- as.data.frame(lapply(traits, function(x) gsub("&", "+", x)))
states <- traits_test$outline
names(states) <- rownames(traits)

# Standardize branch lengths to fit models
phy$edge.length <- phy$edge.length / max(phy$edge.length)

# Fitting models (this part can take a while...)

# equal rates
unordered.er <- fitpolyMk(phy,states,model="ER")
# symmetrical rates
unordered.sym <- fitpolyMk(phy,states,model="SYM")
# all rates different
unordered.ard <- fitpolyMk(phy,states,model="ARD")

# Evaluating results

aic.w(c(AIC(unordered.er), AIC(unordered.sym), AIC(unordered.ard)))

## Preparing data -------------------------------------------------------------

# The data.frame needs to be transformed into a named vector
# containing tip labels as names and character states as factors.
plyr::count(traits$outline)
cat.trait <- factor(traits$outline, levels = c("circular", "elliptical", "oval", 
                                               "circular&elliptical", "elliptical&oval", 
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

# Polymorphism 1&2
polymorph4 <- rownames(bin.matrix)[bin.matrix[ , c(4)] == 1]
bin.matrix[polymorph4, c(1, 2)] <- 1/2

# Polymorphism 2&3
polymorph5 <- rownames(bin.matrix)[bin.matrix[ , c(5)] == 1]
bin.matrix[polymorph5, c(2, 3)] <- 1/2

# Removing columns indicating polymorphism
bin.matrix <- bin.matrix[ , -c(4:5)]

# Which taxa show missing data? 
missing.data <- names(which(bin.matrix[ , "missing"] == 1))

# For missing data, we must assign a prior probability distribution on the tips 
# that is flat across all possible states.
## Removing the column 'missing'
bin.matrix <- bin.matrix[ , -c(ncol(bin.matrix))]
## Assigning 1/(number of states) for all possible states
bin.matrix[row.names(bin.matrix) %in% missing.data, ] <- 1/ncol(bin.matrix) 

# Seting seed for replicability
set.seed(7)

## SIMMAP ---------------------------------------------------------------------

# Running stochastic mapping (this can take a while)
trees <- make.simmap(phy, bin.matrix, model = "SYM", nsim = 100)
obj <- summary(trees, plot = FALSE)

# putting the columns in the same order of the tip.labels, so we can plot polymophisms
bin.matrix <- bin.matrix[match(phy$tip.label, rownames(bin.matrix)),]

# putting the columns of obj$ace in the same order of bin.matrix
obj$ace <- obj$ace[, colnames(bin.matrix)]

# plotting and saving
cols <- setNames(palette()[1:length(colnames(bin.matrix))], colnames(bin.matrix))

pdf("output/plots/outline_SYM.pdf")
par(lwd = 0.1)
plotTree(phy,
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

