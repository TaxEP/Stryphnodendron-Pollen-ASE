# Ancestral state estimation of pollen grain morphology in the Stryphnodendron clade

## Reference publication

For details on methodology and results, see the published paper:  
https://doi.org/10.1007/s40415-024-01059-y

## Data folder

- **pip_group.tre**: Phylogenetic hypothesis of the Stryphnodendron and Mimosa clades, synthesized from Vasconcelos et al. (2020), Borges et al. (2022), and Lima et al. (2022).
- **pollen_data.csv**: Literature data and new pollen morphology descriptions of Stryphnodendron clade taxa.
- **pseudopip_review_tree.nex**: Phylogeny of the Stryphnodendron clade from Borges et al. (2022).

## Output folder

- **data**: Intermediate files (morphological data subsets and pruned trees).
- **plots**: Figures showing the results of the ancestral state estimation analyses.

## Script folder

Scripts that perform ancestral state estimation analyses of categorical data:
- **grains.R**
- **ornamentation.R**
- **outline.R**

Scripts that perform ancestral state estimation analyses of continuous data:
- **exine_thickness_mean.R**
- **longer_diameter_mean.R**
- **shorter_diameter_mean.R**

Script used to subset pollen morphological data and prune the phylogeny:
- **pruning_data_and_tree.R**

### Contact

For questions, please feel free to contact rfbarduzzi@gmail.com.
