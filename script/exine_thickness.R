#==============================================================================#

# Loading phytools and dependencies
require(phytools)
library(Rphylopars)

#==============================================================================#

#===================#
# CONTINUOUS TRAITS #
#===================#

# Loading tree
phy <- read.nexus("output/data/pruned_tree.nex")

# Loading trait values (continuous traits must be formatted as a data.frame,
# with tip labels as row names and trait labels as column names)
traits <- read.csv("output/data/data_cont.csv", header = TRUE, row.names = 1)

## Testing evolutionary models ------------------------------------------------

# creating species column
traits <- cbind(species = rownames(traits), traits)
# selecting exine thickness
trait <- traits[ , c("species","exine_thickness_mean")]

# checking trait distribution
plot(density(na.omit(trait$exine_thickness_mean)))

# transforming data
trait$exine_thickness_mean <- log(trait$exine_thickness_mean)

# Standardize branch lengths to fit models
phy$edge.length <- phy$edge.length / max(phy$edge.length)

p_BM <- phylopars(trait, phy, model = "BM")
p_OU <- phylopars(trait, phy, model = "OU")
p_EB <- phylopars(trait, phy, model = "EB")
aic.w(c(AIC(p_BM), AIC(p_OU), AIC(p_EB)))

## Preparing data -------------------------------------------------------------

# contMap requires as input a named vector containing the character values and
# respective tip labels
anc_recon <- as.data.frame(p_OU$anc_recon)
cont.trait <- anc_recon$exine_thickness_mean[1:44]
names(cont.trait) <- rownames(anc_recon)[1:44]

# Checking if the tree contain all taxa
missing.names <- names(cont.trait)[!names(cont.trait) %in% phy$tip.label]

## contMap --------------------------------------------------------------------

# Seting seed for replicability
set.seed(7)

# Mapping continuous character by estimating states at internal nodes using
# the method anc.ML, which estimates trait values for tips with missing data

anc <- c(p_OU$anc_recon[45:84], p_OU$logLik)
names(anc) <- c(45:84, "logLik")

obj <- contMap(phy, 
               cont.trait, 
               method = "user",
               anc.states = anc,
               plot = FALSE)

# Inverting colours 
obj <- setMap(obj, invert = TRUE)

# Plotting and saving
pdf("output/plots/exine_OU.pdf")

plot(obj, fsize = c(0.7, 0.7), 
     outline = FALSE, lwd = c(3,7), 
     leg.txt = "exine thickness")

dev.off()

