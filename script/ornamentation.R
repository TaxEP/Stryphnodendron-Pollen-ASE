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
traits <- read.csv("output/data/data_cat.csv", header = TRUE, row.names = 1)

## Testing evolutionary models ------------------------------------------------

# Seting seed for replicability
set.seed(7)

# adjusting data for fitting
traits_test <- as.data.frame(lapply(traits, function(x) gsub("&", "+", x)))
states <- traits_test$ornamentation
names(states) <- rownames(traits)

# Fitting models (this part can take a while...)

# equal rates
unordered.er <- fitpolyMk(stryphnod.tree,states,model="ER")
# symmetrical rates
unordered.sym <- fitpolyMk(stryphnod.tree,states,model="SYM")
# all rates different
unordered.ard <- fitpolyMk(stryphnod.tree,states,model="ARD")

# Evaluating results

aic.w(c(AIC(unordered.er), AIC(unordered.sym), AIC(unordered.ard)))

## Preparing data -------------------------------------------------------------

# The data.frame needs to be transformed into a named vector
# containing tip labels as names and character states as factors.
plyr::count(traits$ornamentation)
cat.trait <- factor(traits$ornamentation, levels = c("areolate", "psilate", "rugulate", "verrucate","verrucate-scabrate", "fossulate", "reticulate",
                                                         "areolate&fossulate", "areolate&psilate", "areolate&reticulate", "areolate&rugulate", "areolate&verrucate", "microverrucate&rugulate", "psilate&verrucate",
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

# Polymorphism areolate&fossulate
polymorph8 <- rownames(bin.matrix)[bin.matrix[ , c(8)] == 1]
bin.matrix[polymorph8, c(1, 6)] <- 1/2

# Polymorphism areolate&psilate
polymorph9 <- rownames(bin.matrix)[bin.matrix[ , c(9)] == 1]
bin.matrix[polymorph9, c(1, 2)] <- 1/2

# Polymorphism areolate&reticulate
polymorph10 <- rownames(bin.matrix)[bin.matrix[ , c(10)] == 1]
bin.matrix[polymorph10, c(1, 7)] <- 1/2

# Polymorphism areolate&rugulate
polymorph11 <- rownames(bin.matrix)[bin.matrix[ , c(11)] == 1]
bin.matrix[polymorph11, c(1, 3)] <- 1/2

# Polymorphism areolate&verrucate
polymorph12 <- rownames(bin.matrix)[bin.matrix[ , c(12)] == 1]
bin.matrix[polymorph12, c(1, 4)] <- 1/2

# Polymorphism microverrucate&rugulate
# putting microverrucate as verrucate
polymorph13 <- rownames(bin.matrix)[bin.matrix[ , c(13)] == 1]
bin.matrix[polymorph13, c(1, 4)] <- 1/2

# Polymorphism psilate&verrucate
polymorph14 <- rownames(bin.matrix)[bin.matrix[ , c(14)] == 1]
bin.matrix[polymorph14, c(2, 4)] <- 1/2

# Removing columns indicating polymorphism
bin.matrix <- bin.matrix[ , -c(8:14)]

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
trees <- make.simmap(stryphnod.tree, bin.matrix, model = "ER", nsim = 100)
obj <- summary(trees, plot = FALSE)

# putting the columns in the same order of the tip.labels, so we can plot polymophisms
bin.matrix <- bin.matrix[match(stryphnod.tree$tip.label, rownames(bin.matrix)),]

# putting the columns of obj$ace in the same order of bin.matrix
obj$ace <- obj$ace[, colnames(bin.matrix)]

# plotting and saving
cols <- setNames(palette()[1:length(colnames(bin.matrix))], colnames(bin.matrix))

pdf("output/plots/ornamentation.pdf")
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

