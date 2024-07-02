#
# Ancestral state reconstruction - Continuous pollen traits
#
# Yago B. Souza and Rafael F. Barduzzi
#
# 2023-2024
#

# 1. Libraries ================================================================

require(phytools)
library(Rphylopars)

# 2. Get data =================================================================

# read tree
phy <- read.nexus("output/data/pruned_tree.nex")

# standardize branch lengths to fit models
phy$edge.length <- phy$edge.length / max(phy$edge.length)

# read continuous traits
traits_cont <- read.csv("output/data/data_cont.csv", 
                        header = TRUE, row.names = 1)

# read categorical traits
traits_cat <- read.csv("output/data/data_cat.csv", 
                   header = TRUE, row.names = 1)

# 3.1. Categorical traits evolutionary models =================================

# Set seed for replicability
set.seed(7)

## Grains number --------------------------------------------------------------

# adjust data for fitting
traits_test <- as.data.frame(lapply(traits_cat, function(x) gsub("&", "+", x)))
states <- traits_test$grains
names(states) <- rownames(traits_cat)

# Fit models (this part can take a while...)

# equal rates
unordered.er <- fitpolyMk(phy,states,model="ER")
# symmetrical rates
unordered.sym <- fitpolyMk(phy,states,model="SYM")
# all rates different
unordered.ard <- fitpolyMk(phy,states,model="ARD")

# evaluating AIC
aic.w(c(AIC(unordered.er), AIC(unordered.sym), AIC(unordered.ard)))

# RESULT: ER

## Outline --------------------------------------------------------------------

# adjust data for fitting
traits_test <- as.data.frame(lapply(traits_cat, function(x) gsub("&", "+", x)))
states <- traits_test$outline
names(states) <- rownames(traits_cat)

# deal with NAs (all states)
states[is.na(states)] <- "circular&elliptical&oval"

# Fit models (this part can take a while...)

# equal rates
unordered.er <- fitpolyMk(phy,states,model="ER")
# symmetrical rates
unordered.sym <- fitpolyMk(phy,states,model="SYM")
# all rates different
unordered.ard <- fitpolyMk(phy,states,model="ARD")

# evaluating AIC
aic.w(c(AIC(unordered.er), AIC(unordered.sym), AIC(unordered.ard)))

# RESULT: SYM

## Ornamentation --------------------------------------------------------------

# adjust data for fitting
traits_test <- as.data.frame(lapply(traits_cat, function(x) gsub("&", "+", x)))
states <- traits_test$ornamentation
names(states) <- rownames(traits_cat)

# deal with NAs (all states)
states[is.na(states)] <- 
  "areolate&psilate&rugulate&verrucate&fossulate%verrucate-fossulate"

# Fit models (this part can take a while...)

# equal rates
unordered.er <- fitpolyMk(phy,states,model="ER")
# symmetrical rates
unordered.sym <- fitpolyMk(phy,states,model="SYM")
# all rates different
unordered.ard <- fitpolyMk(phy,states,model="ARD")

# evaluating AIC
aic.w(c(AIC(unordered.er), AIC(unordered.sym), AIC(unordered.ard)))

# RESULT: ER

# 3.2. SIMMAP =================================================================

## Grains number (D.U.) -------------------------------------------------------

# The data.frame needs to be transformed into a named vector
# containing tip labels as names and character states as factors.
plyr::count(traits_cat$grains)
cat.trait <- factor(traits_cat$grains, levels = c("8", "12", "16", "32", "8&12",
                                              "8&16", "12&16", "missing"))
cat.trait[is.na(cat.trait)] <- "missing"
names(cat.trait) <- rownames(traits_cat)

# The named vector, then, needs to be transformed into a binary matrix
# with tip labels as row names and character states as columns
bin.matrix <- to.matrix(cat.trait, levels(cat.trait))

# Deal with polymorphisms
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

# remove columns indicating polymorphism
bin.matrix <- bin.matrix[ , -c(5,6,7)]

# which taxa show missing data? 
missing.data <- names(which(bin.matrix[ , "missing"] == 1))

# For missing data, we must assign a prior probability distribution on the tips 
# that is flat across all possible states.
## remove the column 'missing'
bin.matrix <- bin.matrix[ , -c(ncol(bin.matrix))]
## assign 1/(number of states) for all possible states
bin.matrix[row.names(bin.matrix) %in% missing.data, ] <- 1/ncol(bin.matrix) 

# run stochastic mapping (this can take a while)
trees <- make.simmap(phy, bin.matrix, model = "ER", nsim = 100)
obj_gn <- summary(trees, plot = FALSE)

# put columns in the same order of the tip.labels, so we can plot polymophisms
bin.matrix_gn <- bin.matrix[match(phy$tip.label, rownames(bin.matrix)),]

# put the columns of obj$ace in the same order of bin.matrix
obj_gn$ace <- obj_gn$ace[, colnames(bin.matrix_gn)]

# cols info
cols_gn <- setNames(palette()[1:length(colnames(bin.matrix_gn))], 
                    colnames(bin.matrix_gn))

## Outline (D.U.) -------------------------------------------------------------

# The data.frame needs to be transformed into a named vector
# containing tip labels as names and character states as factors.
plyr::count(traits_cat$outline)
cat.trait <- factor(traits_cat$outline, levels = c("circular", "elliptical", "oval", 
                                               "circular&elliptical", "elliptical&oval", 
                                               "missing"))
cat.trait[is.na(cat.trait)] <- "missing"
names(cat.trait) <- rownames(traits_cat)

# The named vector, then, needs to be transformed into a binary matrix
# with tip labels as row names and character states as columns
bin.matrix <- to.matrix(cat.trait, levels(cat.trait))

# Deal with polymorphisms
## When taxa with polymorphisms are present, we must assign a value
# equal to 1/(number of states) for each possible state and then remove 
# the column that refers to the polymorphic state. 

# Polymorphism 1&2
polymorph4 <- rownames(bin.matrix)[bin.matrix[ , c(4)] == 1]
bin.matrix[polymorph4, c(1, 2)] <- 1/2

# Polymorphism 2&3
polymorph5 <- rownames(bin.matrix)[bin.matrix[ , c(5)] == 1]
bin.matrix[polymorph5, c(2, 3)] <- 1/2

# remove columns indicating polymorphism
bin.matrix <- bin.matrix[ , -c(4:5)]

# which taxa show missing data? 
missing.data <- names(which(bin.matrix[ , "missing"] == 1))

# For missing data, we must assign a prior probability distribution on the tips 
# that is flat across all possible states.
## remove the column 'missing'
bin.matrix <- bin.matrix[ , -c(ncol(bin.matrix))]
## assign 1/(number of states) for all possible states
bin.matrix[row.names(bin.matrix) %in% missing.data, ] <- 1/ncol(bin.matrix) 

# run stochastic mapping (this can take a while)
trees <- make.simmap(phy, bin.matrix, model = "SYM", nsim = 100)
obj_out <- summary(trees, plot = FALSE)

# put columns in the same order of the tip.labels, so we can plot polymophisms
bin.matrix_out <- bin.matrix[match(phy$tip.label, rownames(bin.matrix)),]

# put the columns of obj$ace in the same order of bin.matrix
obj_out$ace <- obj_out$ace[, colnames(bin.matrix_out)]

# cols info
cols_out <- setNames(palette()[1:length(colnames(bin.matrix_out))], 
                     colnames(bin.matrix_out))

## Ornamentation (pollen) -----------------------------------------------------

# The data.frame needs to be transformed into a named vector
# containing tip labels as names and character states as factors.
plyr::count(traits_cat$ornamentation)
cat.trait <- factor(
  traits_cat$ornamentation, 
  levels = c("areolate", "psilate", "rugulate", "verrucate",
             "verrucate-scabrate", "fossulate", 
             "areolate&fossulate", "areolate&psilate", 
             "areolate&rugulate", "areolate&verrucate", "psilate&verrucate",
             "rugulate&verrucate", "missing"))

cat.trait[is.na(cat.trait)] <- "missing"
names(cat.trait) <- rownames(traits_cat)

# The named vector, then, needs to be transformed into a binary matrix
# with tip labels as row names and character states as columns
bin.matrix <- to.matrix(cat.trait, levels(cat.trait))

# Deal with polymorphisms
## When taxa with polymorphisms are present, we must assign a value
# equal to 1/(number of states) for each possible state and then remove 
# the column that refers to the polymorphic state. 

# Polymorphism areolate&fossulate
polymorph7 <- rownames(bin.matrix)[bin.matrix[ , c(7)] == 1]
bin.matrix[polymorph7, c(1, 6)] <- 1/2

# Polymorphism areolate&psilate
polymorph8 <- rownames(bin.matrix)[bin.matrix[ , c(8)] == 1]
bin.matrix[polymorph8, c(1, 2)] <- 1/2

# Polymorphism areolate&rugulate
polymorph9 <- rownames(bin.matrix)[bin.matrix[ , c(9)] == 1]
bin.matrix[polymorph9, c(1, 3)] <- 1/2

# Polymorphism areolate&verrucate
polymorph10 <- rownames(bin.matrix)[bin.matrix[ , c(10)] == 1]
bin.matrix[polymorph10, c(1, 4)] <- 1/2

# Polymorphism psilate&verrucate
polymorph11 <- rownames(bin.matrix)[bin.matrix[ , c(11)] == 1]
bin.matrix[polymorph11, c(2, 4)] <- 1/2

# Polymorphism rugulate&verrucate
polymorph12 <- rownames(bin.matrix)[bin.matrix[ , c(12)] == 1]
bin.matrix[polymorph12, c(3, 4)] <- 1/2

# remove columns indicating polymorphism
bin.matrix <- bin.matrix[ , -c(7:12)]

# which taxa show missing data? 
missing.data <- names(which(bin.matrix[ , "missing"] == 1))

# For missing data, we must assign a prior probability distribution on the tips 
# that is flat across all possible states.
## remove the column 'missing'
bin.matrix <- bin.matrix[ , -c(ncol(bin.matrix))]
## assign 1/(number of states) for all possible states
bin.matrix[row.names(bin.matrix) %in% missing.data, ] <- 1/ncol(bin.matrix) 

# run stochastic mapping (this can take a while)
trees <- make.simmap(phy, bin.matrix, model = "ER", nsim = 100)
obj_orn <- summary(trees, plot = FALSE)

# put columns in the same order of the tip.labels, so we can plot polymophisms
bin.matrix_orn <- bin.matrix[match(phy$tip.label, rownames(bin.matrix)),]

# put the columns of obj$ace in the same order of bin.matrix
obj_orn$ace <- obj_orn$ace[, colnames(bin.matrix_orn)]

# plot and save
cols_orn <- setNames(palette()[1:length(colnames(bin.matrix_orn))], 
                     colnames(bin.matrix_orn))

# 4.1. Continuous traits evolutionary models ==================================

# creating species column
traits_cont <- cbind(species = rownames(traits_cont), traits_cont)

## exine thickness ------------------------------------------------------------
trait_et <- traits_cont[ , c("species","exine_thickness_mean")]

# checking trait distribution
plot(density(na.omit(trait_et$exine_thickness_mean)))

# transforming data
trait_et$exine_thickness_mean <- log(trait_et$exine_thickness_mean)

# fit models
et_BM <- phylopars(trait_et, phy, model = "BM")
et_OU <- phylopars(trait_et, phy, model = "OU")
et_EB <- phylopars(trait_et, phy, model = "EB")
aic.w(c(AIC(et_BM), AIC(et_OU), AIC(et_EB)))

# RESULT: OU

## longer diameter ------------------------------------------------------------

trait_ld <- traits_cont[ , c("species","longer_diameter_mean")]

# checking trait distribution
plot(density(na.omit(trait_ld$longer_diameter_mean)))

# transforming data
trait_ld$longer_diameter_mean <- log(trait_ld$longer_diameter_mean)

# fit models
ld_BM <- phylopars(trait_ld, phy, model = "BM")
ld_OU <- phylopars(trait_ld, phy, model = "OU")
ld_EB <- phylopars(trait_ld, phy, model = "EB")
aic.w(c(AIC(ld_BM), AIC(ld_OU), AIC(ld_EB)))

# RESULTS: BM

## shorter diameter -----------------------------------------------------------

trait_sd <- traits_cont[ , c("species","shorter_diameter_mean")]

# checking trait distribution
plot(density(na.omit(trait_sd$shorter_diameter_mean)))

# transforming data
trait_sd$shorter_diameter_mean <- log(trait_sd$shorter_diameter_mean)

# fit models
sd_BM <- phylopars(trait_sd, phy, model = "BM")
sd_OU <- phylopars(trait_sd, phy, model = "OU")
sd_EB <- phylopars(trait_sd, phy, model = "EB")
aic.w(c(AIC(sd_BM), AIC(sd_OU), AIC(sd_EB)))

# RESULTS: EB

# 4.2. countMap ===============================================================

## Exine thickness (grains) ---------------------------------------------------

# contMap requires as input a named vector containing the character values and
# respective tip labels
anc_recon <- as.data.frame(et_OU$anc_recon)
cont.trait <- anc_recon$exine_thickness_mean[1:44]
names(cont.trait) <- rownames(anc_recon)[1:44]

# Checking if the tree contain all taxa
missing.names <- names(cont.trait)[!names(cont.trait) %in% phy$tip.label]

# Mapping continuous character by estimating states at internal nodes using
# the method anc.ML, which estimates trait values for tips with missing data

anc <- c(et_OU$anc_recon[45:84], et_OU$logLik)
names(anc) <- c(45:84, "logLik")

obj_exine <- contMap(phy, 
               cont.trait, 
               method = "user",
               anc.states = anc,
               plot = FALSE)

# Inverting colours 
obj_exine <- setMap(obj_exine, invert = TRUE)

## Longer diameter (D.U.) -----------------------------------------------------

# contMap requires as input a named vector containing the character values and
# respective tip labels
anc_recon <- as.data.frame(ld_BM$anc_recon)
cont.trait <- anc_recon$longer_diameter_mean[1:44]
names(cont.trait) <- rownames(anc_recon)[1:44]

# Checking if the tree contain all taxa
missing.names <- names(cont.trait)[!names(cont.trait) %in% phy$tip.label]

# Mapping continuous character by estimating states at internal nodes using
# the method anc.ML, which estimates trait values for tips with missing data

anc <- c(ld_BM$anc_recon[45:84], ld_BM$logLik)
names(anc) <- c(45:84, "logLik")

obj_longerd <- contMap(phy, 
               cont.trait, 
               method = "user",
               anc.states = anc,
               plot = FALSE)

# Inverting colours 
obj_longerd <- setMap(obj_longerd, invert = TRUE)

## Shorter diameter (D.U.) ----------------------------------------------------

# contMap requires as input a named vector containing the character values and
# respective tip labels
anc_recon <- as.data.frame(sd_EB$anc_recon)
cont.trait <- anc_recon$shorter_diameter_mean[1:44]
names(cont.trait) <- rownames(anc_recon)[1:44]

# Checking if the tree contain all taxa
missing.names <- names(cont.trait)[!names(cont.trait) %in% phy$tip.label]

# Mapping continuous character by estimating states at internal nodes using
# the method anc.ML, which estimates trait values for tips with missing data

anc <- c(sd_EB$anc_recon[45:84], sd_EB$logLik)
names(anc) <- c(45:84, "logLik")

obj_shorterd <- contMap(phy, 
               cont.trait, 
               method = "user",
               anc.states = anc,
               plot = FALSE)

# Inverting colours 
obj_shorterd <- setMap(obj_shorterd, invert = TRUE)

# 5. Plot and save ============================================================

# grains number and outline
pdf("output/plots/grains_ER-outline_SYM.pdf")

par(lwd = 0.1)

layout(matrix(c(1, 2, 3), 1, 3), widths = c(1, 0.6, 1))

plotTree(phy,
         type = "phylogram",
         ftype = "off",
         cex = 0.2, 
         lwd = 1, 
         ftype = "i")
nodelabels(pie = obj_gn$ace,
           piecol = cols_gn, 
           cex=1.2)
tiplabels(pie = bin.matrix_gn,
          piecol = cols_gn,
          cex=1)
legend("bottomleft", legend = names(cols_gn), fill = cols_gn, bty = "n", 
       cex = 1, 
       y.intersp = 1)
title("a. Dispersal unit number of grains", line = -1)

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), 
     ylim = get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)
text(0.5, seq_along(phy$tip.label), labels = gsub("_", " ", phy$tip.label), 
     cex = 1)

plotTree(phy,
         type = "phylogram",
         direction = "leftwards",
         ftype = "off",
         fsize = 0.5, 
         offset = 0.5,
         cex = 0.2, 
         lwd = 1, 
         ftype = "i")
nodelabels(pie = obj_out$ace,
           piecol = cols_out, 
           cex=1.2)
tiplabels(pie = bin.matrix_out,
          piecol = cols_out,
          cex=1)
legend("bottomright", legend = names(cols_out), fill = cols_out, bty = "n", 
       cex = 1, 
       y.intersp = 1)
title("b. Dispersal unit outline", line = -1)

dev.off()

# ornamentation and exine thickness
pdf("output/plots/ornamentation_ER-exine_OU.pdf")

par(lwd = 0.1)

layout(matrix(c(1, 2, 3), 1, 3), widths = c(1, 0.6, 1))

plotTree(phy,
         type = "phylogram",
         ftype = "off",
         cex = 0.2, 
         lwd = 1, 
         ftype = "i")
nodelabels(pie = obj_orn$ace,
           piecol = cols_orn, 
           cex=1.2)
tiplabels(pie = bin.matrix_orn,
          piecol = cols_orn,
          cex=1)
legend("bottomleft", legend = names(cols_orn), fill = cols_orn, bty = "n", 
       cex = 1, 
       y.intersp = 1)
title("a. Pollen ornamentation", line = -1)

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), 
     ylim = get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)
text(0.5, seq_along(phy$tip.label), labels = gsub("_", " ", phy$tip.label), 
     cex = 1)

plot(obj_exine, fsize = c(0.7, 0.7), 
     outline = FALSE, lwd = c(3,7), 
     leg.txt = "exine thickness (log)",
     direction="leftwards",
     ftype="off")
title("b. Exine thickness", line = -1)

dev.off()

# longer and shorter diameter

pdf("output/plots/logerd_BM-shorterd_EB.pdf")

par(lwd = 0.1)

layout(matrix(c(1, 2, 3), 1, 3), widths = c(1, 0.6, 1))

plot(obj_longerd, fsize = c(0.7, 0.7), 
     outline = FALSE, lwd = c(3,7), 
     leg.txt = "longer diameter (log)",
     ftype="off")
title("a. Dispersal unit longer diameter", line = -1)

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), 
     ylim = get("last_plot.phylo",envir=.PlotPhyloEnv)$y.lim)
text(0.5, seq_along(phy$tip.label), labels = gsub("_", " ", phy$tip.label), 
     cex = 1)

plot(obj_shorterd, fsize = c(0.7, 0.7), 
     outline = FALSE, lwd = c(3,7), 
     leg.txt = "shorter diameter (log)",
     direction="leftwards",
     ftype="off")
title("b. Dispersal unit shorter diameter", line = -1)

dev.off()

