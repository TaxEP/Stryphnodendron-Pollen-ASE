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

## Testing evolutionary models ------------------------------------------------

# creating species column
traits_species <- rownames_to_column(traits, var = "species")
# selecting species and trait columns
traits_species <- traits_species[ , c(1,7)]

# checking trait distribution
plot(density(na.omit(traits_species$longer_diameter_mean)))

# transforming data
traits_species$longer_diameter_mean <- log(traits_species$longer_diameter_mean)

p_BM <- phylopars(traits_species, stryphnod.tree, model = "BM")
p_OU <- phylopars(traits_species, stryphnod.tree, model = "OU")
p_EB <- phylopars(traits_species, stryphnod.tree, model = "EB")
aic.w(c(AIC(p_BM), AIC(p_OU), AIC(p_EB)))

## Preparing data -------------------------------------------------------------

# contMap requires as input a named vector containing the character values and
# respective tip labels
anc_recon <- as.data.frame(p_BM$anc_recon)
cont.trait <- anc_recon$longer_diameter_mean[1:44]
names(cont.trait) <- rownames(anc_recon)[1:44]

# Checking if the tree contain all taxa
missing.names <- names(cont.trait)[!names(cont.trait) %in% stryphnod.tree$tip.label]

## contMap --------------------------------------------------------------------

# Seting seed for replicability
set.seed(7)

# Mapping continuous character by estimating states at internal nodes using
# the method anc.ML, which estimates trait values for tips with missing data

anc <- c(p_BM$anc_recon[45:84], p_BM$logLik)
names(anc) <- c(45:84, "logLik")

obj <- contMap(stryphnod.tree, 
               cont.trait, 
               method = "user",
               anc.states = anc,
               plot = FALSE)

# Inverting colours 
obj <- setMap(obj, invert = TRUE)

# Plotting and saving
pdf("output/plots/ld_mean_BM.pdf")

plot(obj, fsize = c(0.7, 0.7), 
     outline = FALSE, lwd = c(3,7), 
     leg.txt = "longer diameter (mean)")

dev.off()

