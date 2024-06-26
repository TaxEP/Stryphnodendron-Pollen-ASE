#==============================================================================#

# Loading phytools and dependencies
require(phytools)

#==============================================================================#

#===================#
# CONTINUOUS TRAITS #
#===================#

# Loading tree
stryphnod.tree <- read.nexus("output/data/pruned_tree.nex")

# Loading trait values (continuous traits must be formatted as a data.frame,
# with tip labels as row names and trait labels as column names)
traits <- read.csv("output/data/data_cont.csv", header = TRUE, row.names = 1)

## Preparing data -------------------------------------------------------------

# contMap requires as input a named vector containing the character values and
# respective tip labels
cont.trait <- traits[ , "shorter_diameter_mean"]
names(cont.trait) <- rownames(traits)

# Removing NA (the analysis will recognize discrepancies between the tree
# and the matrix as missing data)
cont.trait <- cont.trait[!is.na(cont.trait)]

# Checking if the tree contain all taxa
missing.names <- names(cont.trait)[!names(cont.trait) %in% stryphnod.tree$tip.label]

# Seting seed for replicability
set.seed(7)

## contMap --------------------------------------------------------------------

# Mapping continuous character by estimating states at internal nodes using
# the method anc.ML, which estimates trait values for tips with missing data
obj <- contMap(stryphnod.tree, cont.trait, method = "anc.ML", plot = FALSE)

# Inverting colours 
obj <- setMap(obj, invert = TRUE)

# Plotting
pdf("output/plots/sd_mean.pdf")

plot(obj, fsize = c(0.7, 0.7), 
     outline = FALSE, lwd = c(3,7), 
     leg.txt = "shorter diameter (mean)")

dev.off()

## Log contMap ----------------------------------------------------------------

cont.trait.log <- log(cont.trait)

# Mapping continuous character by estimating states at internal nodes using
# the method anc.ML, which estimates trait values for tips with missing data
obj.log <- contMap(stryphnod.tree, cont.trait.log, method = "anc.ML", plot = FALSE)

# Inverting colours 
obj.log <- setMap(obj.log, invert = TRUE)

# Plotting
pdf("output/plots/SUPP-sd_mean_log.pdf")

plot(obj.log, fsize = c(0.7, 0.7), 
     outline = FALSE, lwd = c(3,7), 
     leg.txt = "shorter diameter (mean-log)")

dev.off()
