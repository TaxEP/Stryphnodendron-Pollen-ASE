library(ape)
library(tidyverse)

    #---------------#
    # Cleaning tree #
    #---------------#

# getting tree and tip labels

read.tree("stryphnod-mapping.tre") -> tree

# OPTIONAL - list of species to be corrected

as.data.frame(tree$tip.label) %>%
  filter(grepl("[0-9]", tree$tip.label) | grepl("gb", tree$tip.label))

# adjusting names manually (remove the "#" above if you need to run this part)

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

new_tree <- drop.tip(tree, c(sample(c("Anadenanthera_colubrina1", "Anadenanthera_colubrina2"), 1), 
                 sample(c("Pseudopiptadenia_brenanii1", "Pseudopiptadenia_brenanii2"), 1),
                 sample(c("Stryphnodendron_fissuratum1", "Stryphnodendron_fissuratum2"), 1),
                 sample(c("Stryphnodendron_coriaceum1", "Stryphnodendron_coriaceum2", "Stryphnodendron_coriaceum3"), 2),
                 sample(c("Stryphnodendron_paniculatum1", "Stryphnodendron_paniculatum2"), 1),
                 sample(c("Stryphnodendron_duckeanum1", "Stryphnodendron_duckeanum2", "Stryphnodendron_duckeanum3"), 2),
                 sample(c("Stryphnodendron_foreroi1", "Stryphnodendron_foreroi2"), 1),
                 sample(c("Stryphnodendron_sp1", "Stryphnodendron_sp2", "Stryphnodendron_sp3", "Stryphnodendron_sp4"), 3), 
                 sample(c("Stryphnodendron_pulcherrimum1", "Stryphnodendron_pulcherrimum2", "Stryphnodendron_pulcherrimum3", "Stryphnodendron_pulcherrimum4"), 3),
                 sample(c("Stryphnodendron_polyphyllum1", "Stryphnodendron_polyphyllum2"), 1),
                 sample(c("Stryphnodendron_velutinum1", "Stryphnodendron_velutinum2"), 1),
                 sample(c("Stryphnodendron_roseiflorum1", "Stryphnodendron_roseiflorum2"), 1),
                 sample(c("Stryphnodendron_obovatum1", "Stryphnodendron_obovatum2", "Stryphnodendron_obovatum3"), 2),
                 sample(c("Stryphnodendron_heringeri1", "Stryphnodendron_heringeri2"), 1),
                 sample(c("Stryphnodendron_rotundifolium1", "Stryphnodendron_rotundifolium2"), 1),
                 sample(c("Stryphnodendron_rotundifolium_villosum1", "Stryphnodendron_rotundifolium_villosum2"), 1),
                 sample(c("Stryphnodendron_adstringens1", "Stryphnodendron_adstringens2", "Stryphnodendron_adstringens3", "Stryphnodendron_adstringens4"), 3),
                 sample(c("Pseudopiptadenia_contorta1", "Pseudopiptadenia_contorta2"), 1),
                 sample(c("Microlobius_foetidus1", "Microlobius_foetidus2"), 1)))

new_tree$tip.label

# removing the numbers...

new_tree$tip.label <- sub("[0-9]", "", new_tree$tip.label)

new_tree$tip.label

# all done. now, writing the cleared tree

write.tree(new_tree, file = "stryphnod_clean.tre")

# OBS: don't forget to update the taxon names if necessary:

write.csv(new_tree$tip.label, file = "new_tips_noupdate.txt")

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





