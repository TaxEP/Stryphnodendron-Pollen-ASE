#==============================================================================#

# Loading phytools and dependencies
require(phytools)

#==============================================================================#

#===================#
# CONTINUOUS TRAITS #
#===================#

# Loading tree
stryphnod.tree <- read.nexus("pruned_tree.nex")

# Loading trait values (continuous traits must be formatted as a data.frame,
# with tip labels as row names and trait labels as column names)
traits <- read.csv("data_cont.csv", header = TRUE, row.names = 1)

# contMap requires as input a named vector containing the character values and
# respective tip labels
cont.trait <- traits[ , "et_mean"]
names(cont.trait) <- rownames(traits)

# Removing NA (the analysis will recognize discrepancies between the tree
# and the matrix as missing data)
cont.trait <- cont.trait[!is.na(cont.trait)]

# Checking if the tree contain all taxa
missing.names <- names(cont.trait)[!names(cont.trait) %in% stryphnod.tree$tip.label]

# Mapping continuous character by estimating states at internal nodes using
# the method anc.ML, which estimates trait values for tips with missing data
obj <- contMap(stryphnod.tree, cont.trait, method = "anc.ML", plot = FALSE)

# Inverting colours 
obj <- setMap(obj, invert = TRUE)

# Plotting
plot(obj, fsize = c(0.4,1), 
     outline = FALSE, lwd = c(3,7), 
     leg.txt = "et_mean")

# Saving as pdf
pdf("all_analyses/cont/et_mean/plot.pdf"); plot(obj, fsize = c(0.4,1), 
                                                outline = FALSE, lwd = c(3,7), 
                                                leg.txt = "et_mean"); dev.off()

# What about using log transformed data?
cont.trait.log <- log(cont.trait)

# Mapping continuous character by estimating states at internal nodes using
# the method anc.ML, which estimates trait values for tips with missing data
obj.log <- contMap(stryphnod.tree, cont.trait.log, method = "anc.ML", plot = FALSE)

# Inverting colours 
obj.log <- setMap(obj.log, invert = TRUE)

# Plotting
plot(obj.log, fsize = c(0.4, 0.8), 
     outline = FALSE, lwd = c(3,7), 
     leg.txt = "et_mean (log)")

# Saving as pdf
pdf("all_analyses/cont/et_mean/plot_log.pdf"); plot(obj.log, fsize = c(0.4, 0.8), 
                                                    outline = FALSE, lwd = c(3, 7), 
                                                    leg.txt = "et_mean (log)"); dev.off()
