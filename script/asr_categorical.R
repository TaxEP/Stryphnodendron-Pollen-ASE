#
# Ancestral state reconstruction - Categorical pollen traits
#
# Yago B. Souza and Rafael F. Barduzzi
#
# 2023-2024
#

# 1. Libraries ================================================================

require(phytools)

# 2. Get data =================================================================

# read tree
phy <- read.nexus("output/data/pruned_tree.nex")

# standardize branch lengths to fit models
phy$edge.length <- phy$edge.length / max(phy$edge.length)

# read traits
traits <- read.csv("output/data/data_cat.csv", header = TRUE, row.names = 1)

# 3. Grains number (D.U.) =====================================================

## Test evolutionary models ---------------------------------------------------

# set seed for replicability
set.seed(7)

# adjust data for fitting
traits_test <- as.data.frame(lapply(traits, function(x) gsub("&", "+", x)))
states <- traits_test$grains
names(states) <- rownames(traits)

# Fit models (this part can take a while...)

# equal rates
unordered.er <- fitpolyMk(phy,states,model="ER")
# symmetrical rates
unordered.sym <- fitpolyMk(phy,states,model="SYM")
# all rates different
unordered.ard <- fitpolyMk(phy,states,model="ARD")

# evaluating results
aic.w(c(AIC(unordered.er), AIC(unordered.sym), AIC(unordered.ard)))

## Prepare data ---------------------------------------------------------------

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

## SIMMAP --------------------------------------------------------------------

# run stochastic mapping (this can take a while)
trees <- make.simmap(phy, bin.matrix, model = "ER", nsim = 100)
obj <- summary(trees, plot = FALSE)

# put columns in the same order of the tip.labels, so we can plot polymophisms
bin.matrix <- bin.matrix[match(phy$tip.label, rownames(bin.matrix)),]

# put the columns of obj$ace in the same order of bin.matrix
obj$ace <- obj$ace[, colnames(bin.matrix)]

# plot and save
cols <- setNames(palette()[1:length(colnames(bin.matrix))], colnames(bin.matrix))

pdf("output/plots/grains_ER.pdf")
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

# 4. Outline (D.U.) ===========================================================

## Test evolutionary models ---------------------------------------------------

# set seed for replicability
set.seed(7)

# adjust data for fitting
traits_test <- as.data.frame(lapply(traits, function(x) gsub("&", "+", x)))
states <- traits_test$outline
names(states) <- rownames(traits)

# deal with NAs (all states)
states[is.na(states)] <- "circular&elliptical&oval"

# Fit models (this part can take a while...)

# equal rates
unordered.er <- fitpolyMk(phy,states,model="ER")
# symmetrical rates
unordered.sym <- fitpolyMk(phy,states,model="SYM")
# all rates different
unordered.ard <- fitpolyMk(phy,states,model="ARD")

# Evaluating results

aic.w(c(AIC(unordered.er), AIC(unordered.sym), AIC(unordered.ard)))

## Prepare data ---------------------------------------------------------------

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

# set seed for replicability
set.seed(7)

## SIMMAP ---------------------------------------------------------------------

# run stochastic mapping (this can take a while)
trees <- make.simmap(phy, bin.matrix, model = "SYM", nsim = 100)
obj <- summary(trees, plot = FALSE)

# put columns in the same order of the tip.labels, so we can plot polymophisms
bin.matrix <- bin.matrix[match(phy$tip.label, rownames(bin.matrix)),]

# put the columns of obj$ace in the same order of bin.matrix
obj$ace <- obj$ace[, colnames(bin.matrix)]

# plot and save
cols <- setNames(palette()[1:length(colnames(bin.matrix))], 
                 colnames(bin.matrix))

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

# 5. Ornamentation (pollen) ===================================================

## Test evolutionary models ---------------------------------------------------

# set seed for replicability
set.seed(7)

# adjust data for fitting
traits_test <- as.data.frame(lapply(traits, function(x) gsub("&", "+", x)))
states <- traits_test$ornamentation
names(states) <- rownames(traits)

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

# evaluating results
aic.w(c(AIC(unordered.er), AIC(unordered.sym), AIC(unordered.ard)))

## Prepare data ---------------------------------------------------------------

# The data.frame needs to be transformed into a named vector
# containing tip labels as names and character states as factors.
plyr::count(traits$ornamentation)
cat.trait <- factor(
  traits$ornamentation, 
  levels = c("areolate", "psilate", "rugulate", "verrucate",
             "verrucate-scabrate", "fossulate", 
             "areolate&fossulate", "areolate&psilate", 
             "areolate&rugulate", "areolate&verrucate", "psilate&verrucate",
             "rugulate&verrucate", "missing"))

cat.trait[is.na(cat.trait)] <- "missing"
names(cat.trait) <- rownames(traits)

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

## SIMMAP ---------------------------------------------------------------------

# run stochastic mapping (this can take a while)
trees <- make.simmap(phy, bin.matrix, model = "ER", nsim = 100)
obj <- summary(trees, plot = FALSE)

# put columns in the same order of the tip.labels, so we can plot polymophisms
bin.matrix <- bin.matrix[match(phy$tip.label, rownames(bin.matrix)),]

# put the columns of obj$ace in the same order of bin.matrix
obj$ace <- obj$ace[, colnames(bin.matrix)]

# plot and save
cols <- setNames(palette()[1:length(colnames(bin.matrix))], colnames(bin.matrix))

pdf("output/plots/ornamentation_ER.pdf")
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

