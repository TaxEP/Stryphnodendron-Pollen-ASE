library(ape)
library(tidyverse)

#---------------#
# Cleaning tree #
#---------------#

# getting tree and tip labels

read.tree("pseudopip_review_tree.nex") -> tree

# OPTIONAL - list of species to be corrected

as.data.frame(tree$tip.label) %>%
  filter(grepl("[0-9]", tree$tip.label) | grepl("gb", tree$tip.label))

# adjusting names manually (remove the "#" below if you need to run this part)

# edit(tree$tip.label) -> tips_adjusted_v1

# write.csv(tips_adjusted_v1, file = "tips_adjusted_v1.txt")

# finding duplicates

read.csv("tips_adjusted_v1.txt") -> tips_adjusted_v1

as.vector(tips_adjusted_v1$x) -> tips_adjusted_v1

tips_adjusted_v1[duplicated(tips_adjusted_v1)]

# classifying duplicates manually (I recommend do this with a text editor in the directory...)

read.csv("tips_adjusted_v2.txt") -> tips_adjusted_v2 # or edit(tips_adjusted_v1) -> tips_adjusted_v2

as.vector(tips_adjusted_v2$x) -> tips_adjusted_v2

# applying new tip labels

tree$tip.label <- tips_adjusted_v2

tree$tip.label[duplicated(tree$tip.label)] # duplicates must be 0 now

# removing duplicates randomly

as.data.frame(tree$tip.label) %>%
  filter(grepl("[2-9]", tree$tip.label))

new_tree <- drop.tip(tree, c(sample(c("Microlobius_foetidus1", "Microlobius_foetidus2"), 1), 
                             sample(c("Parapiptadenia_blanchetii1", "Parapiptadenia_blanchetii2"), 1),
                             sample(c("Parapiptadenia_pterosperma1", "Parapiptadenia_pterosperma2"), 1),
                             sample(c("Parapiptadenia_zehntneri1", "Parapiptadenia_zehntneri2", "Parapiptadenia_zehntneri3", "Parapiptadenia_zehntneri4"), 3),
                             sample(c("Pseudopiptadenia_bahiana1", "Pseudopiptadenia_bahiana2", "Pseudopiptadenia_bahiana3"), 2),
                             sample(c("Pseudopiptadenia_contorta1", "Pseudopiptadenia_contorta2"), 1),
                             sample(c("Pseudopiptadenia_sp1", "Pseudopiptadenia_sp2"), 1),
                             sample(c("Pityrocarpa_moniliformis1", "Pityrocarpa_moniliformis2", "Pityrocarpa_moniliformis3"), 2), 
                             sample(c("Pseudopiptadenia_brenanii1", "Pseudopiptadenia_brenanii2", "Pseudopiptadenia_brenanii3", "Pseudopiptadenia_brenanii4"), 3),
                             sample(c("Pseudopiptadenia_leptostachya1", "Pseudopiptadenia_leptostachya2"), 1)))

new_tree$tip.label

# removing the numbers...

new_tree$tip.label <- sub("[0-9]", "", new_tree$tip.label)

new_tree$tip.label

# all done. now, writing the cleared tree

write.tree(new_tree, file = "stryphnod_clean_new.tre")

# OBS: don't forget to update the taxon names if necessary:

write.csv(new_tree$tip.label, file = "new_tips_noupdate_new.txt")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------ #

library(ape)
library(tidyverse)

#---------------#
# Updating tips # - ! RUN THIS HAVING A UPDATED LIST OF TAXON NAMES !
#---------------#

# getting the updated taxon names

read.csv("new_tips_updated.txt") -> tips_updated

# if necessary, edit the "new names" column with: edit(tips_updated) -> tips_updated

as.vector(tips_updated$accepted.name) -> tips_updated

read.tree("stryphnod_clean.tre") -> stryphnod_clean

stryphnod_clean$tip.label <- tips_updated

# all done. now, writing the cleared AND updated tree

write.tree(stryphnod_clean, file = "stryphnod_clean_updated.tre")





