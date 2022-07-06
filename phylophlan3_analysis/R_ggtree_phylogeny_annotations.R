
setwd("D:/OneDrive - Baylor College of Medicine/Documents/Kellermayer_FMT/@Experiment/Phylogeny & ST/ggtree_phylogenetic tree")
getwd()

#install.packages("devtools")
#devtools::install_github("YuLab-SMU/ggtree")
#devtools::install_github("YuLab-SMU/ggtreeExtra")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")
#BiocManager::install("ggtreeExtra")

library(ggtree)
library(treeio)


tree <- read.tree(file = "RAxML_bestTree.input_genomes_refined.tre")

map <- read.csv("map.csv")
df_map <- as.data.frame(map)

tree_map <- full_join(tree, df_map, by='label')

ggtree(tree, layout = 'circular', branch.length='none') 

ggtree(tree_map, layout = 'circular', branch.length='none') + 
  geom_tippoint(aes(color=tcdA_structure))

ggtree(tree_map, layout = 'circular', branch.length='none') + 
  +   geom_tippoint(aes(color=tcdA_structure), size=1) + geom_tiplab(size=1.5) #showing strain IDs

ggtree(tree_map, layout = 'circular', branch.length='none', aes(color=isolate_source)) + 
  geom_tippoint(aes(shape=tcdA_structure))

library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)

first_custom <- c("gold1", "brown", "green3", "blue", "dimgray", "deeppink1", "burlywood", 
                   "lightslateblue")

second_custom <- c("black", "red1")


p1 <-
ggtree(tree_map, layout = 'fan', open.angle = 10, branch.length='none') + 
  geom_tippoint(aes(color=isolate_source), size = 1) + #increase the size of tippoint (sample point)
  scale_colour_manual(values=first_custom) + #use manual color
  guides(color = guide_legend(override.aes = list(size = 5))) + #adjust the symbol size in legends
  geom_fruit(geom=geom_tile, mapping=aes(fill=tcdA_structure), width=1.5, offset=0.05) + 
  scale_fill_manual(values=second_custom) + #use manual color
  new_scale_fill()

p1


third_custom <- c("gray", "gold1", "black", "cyan", "red1", "blue", "green")
  

#"blueviolet", "darkslateblue", "wheat2", "lightslateblue",  
#"coral2", "mediumvioletred", "greenyellow", "deeppink1", 
#"seagreen2", "mediumpurple1", "deepskyblue1", 
#"yellow", "dimgray", "cyan", "darkgreen", "slateblue", "gray"



p1+ geom_fruit_list(
    geom_fruit(geom=geom_tile, mapping=aes(fill=tcdA_structure), width=1.5, offset=0.05),
    scale_fill_manual(values=second_custom), #use manual color
    new_scale_fill(), # To initialize fill scale
    geom_fruit(geom=geom_tile, mapping=aes(fill=mlst_edited_highlight), width=1.5, offset=0.05),
    scale_fill_manual(values=third_custom)
    ) 




