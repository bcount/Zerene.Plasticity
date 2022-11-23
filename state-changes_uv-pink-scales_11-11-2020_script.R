## Vince Ficarrotta - For Jen and Brian

####################################################################
# Import data, libraries, and generate lists
####################################################################
library('ape'); library('geiger'); library('phytools'); library('castor'); library('gtools'); library('ggplot2')

### Full dataset from Fenner-Counterman, 10-22-2020, directors commentary cut
full_data <- read.csv('C:/Users/vfica/Documents/Counterman_uv/coliadinae color patterns_10-22-2020_Copy of vince blind.csv')
full_data$UV.num_wing <- as.factor(full_data$UV.num_wing)
full_data$Pig.num_wing <- as.factor(full_data$Pig.num_wing)

### split each species group into its own page in a list
### works
subset_full_data <- function(x){full_data[which(full_data$Species.Group == x), ]}
full_list <- lapply(unique(full_data$Species.Group), FUN = subset_full_data)

### Generate setNames list for UV parsimony - a very specific object structure for the asr
### works
set_data <- function(x){setNames(as.factor(full_list[[x]]$UV.num_wing), full_list[[x]]$GS)}
set_names_data <- lapply(1:length(full_list), FUN = set_data)

### Generate setNames list for pigment parsimony - a very specific object structure for the asr
### works
set_data_pig <- function(x){setNames(full_list[[x]]$Pig.num_wing, full_list[[x]]$GS)}
set_names_data_pig <- lapply(1:length(full_list), FUN = set_data_pig)

### generate vector of state changes
### works
grab_changes_rm <- function(x, spgr){spgr[[x]]$total_cost}

########################################################
# Phylogenetic trees
########################################################
### library
library('ape'); library('geiger'); library('phytools'); library('castor'); library('gtools'); library('ggplot2')

### read in data
# full_data <- read.csv('C:/Users/vfica/Documents/Counterman_uv/coliadinae color patterns - vince blind.csv')
tree_file_list <- mixedsort(list.files('C:/Users/vfica/Documents/Counterman_uv/Phylo-Tree-files'))
setwd('C:/Users/vfica/Documents/Counterman_uv/Phylo-Tree-files')
all_trees <- lapply(tree_file_list, FUN = read.nexus)
names(all_trees) <- tree_file_list

#check spelling of species against this list, ask ifyou need pierid species list file
taxa_list <- read.csv('C:/Users/vfica/Documents/Pierid_Phylogenetics/Lamas_Pieridae_04ii08_taxa-list.csv')
tl <- paste(taxa_list$ï..Genus, taxa_list$Sp, sep = '_')
checklist <- for(i in full_data$GS){i %in% tl} #s false returns are false negatives

# perform parsimonly reconstruction of ancestral states, branch lengths are unused
  ### UV apply version  
asr_f <- function(x){asr_max_parsimony(all_trees[[x]], transition_costs = 'all_equal', Nstates = NULL, tip_states = as.integer(set_names_data[[x]]))}
all_asr <- lapply(1:length(all_trees), FUN = asr_f)
names(all_asr) <- tree_file_list

### Pigment apply version  
asr_pig <- function(x){asr_max_parsimony(all_trees[[x]], transition_costs = 'all_equal', Nstates = NULL, tip_states = as.integer(set_names_data_pig[[x]]))}
all_asr_pig <- lapply(1:length(all_trees), FUN = asr_pig)
names(all_asr_pig) <- tree_file_list

### plot UV asr onto the trees
setwd('C:/Users/vfica/Documents/Counterman_uv/R objects/Phylo-tree_asr_plots')
colors <- c('black', 'cyan3', 'red2', 'green3', 'orange1')
state_names <- c('Missing Data', 'Absent', 'Hindwing', 'Forewing', 'Both Wings')

### UV plots
#E_1
pdf(file = 'E_1_phylo-tree.pdf')
plot.phylo(all_trees[[1]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[1]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[1]][all_trees[[1]]$tip.label], levels(set_names_data[[1]])), cex = 0.5, piecol = colors) 
legend(x = 'bottomleft', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 1: ', all_asr[[1]]$total_cost ,' State Changes'))
dev.off()

#E_2
pdf(file = 'E_2_phylo-tree.pdf') 
plot.phylo(all_trees[[2]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[2]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[2]][all_trees[[2]]$tip.label], levels(set_names_data[[2]])), cex = 0.5, piecol = colors) 
legend(x = 'center', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 2: ', all_asr[[2]]$total_cost ,' State Changes'))
dev.off()

#E_3
pdf(file = 'E_3_phylo-tree.pdf')
plot.phylo(all_trees[[3]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[3]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[3]][all_trees[[3]]$tip.label], levels(set_names_data[[3]])), cex = 0.5, piecol = colors) 
legend(x = 'center', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 3: ', all_asr[[3]]$total_cost ,' State Changes'))
dev.off()

#E_4
pdf(file = 'E_4_phylo-tree.pdf')
plot.phylo(all_trees[[4]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[4]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[4]][all_trees[[4]]$tip.label], levels(set_names_data[[4]])), cex = 0.5, piecol = colors) 
legend(x = 'bottomleft', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 4: ', all_asr[[4]]$total_cost ,' State Changes'))
dev.off()

#E_5
pdf(file = 'E_5_phylo-tree.pdf') 
plot.phylo(all_trees[[5]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[5]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[5]][all_trees[[5]]$tip.label], levels(set_names_data[[5]])), cex = 0.5, piecol = colors) 
legend(x = 'center', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 5: ', all_asr[[5]]$total_cost ,' State Changes'))
dev.off()

#E_6
pdf(file = 'E_6_phylo-tree.pdf') 
plot.phylo(all_trees[[6]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[6]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[6]][all_trees[[6]]$tip.label], levels(set_names_data[[6]])), cex = 0.5, piecol = colors) 
legend(x = 'bottomleft', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 6: ', all_asr[[6]]$total_cost ,' State Changes'))
dev.off()

#E_7
pdf(file = 'E_7_phylo-tree.pdf')
plot.phylo(all_trees[[7]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[7]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[7]][all_trees[[7]]$tip.label], levels(set_names_data[[7]])), cex = 0.5, piecol = colors) 
legend(x = 'center', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 7: ', all_asr[[7]]$total_cost ,' State Changes'))
dev.off()

#G_1
pdf(file = 'G_1_phylo-tree.pdf')
plot.phylo(all_trees[[8]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[8]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[8]][all_trees[[8]]$tip.label], levels(set_names_data[[8]])), cex = 0.5, piecol = colors) 
legend(x = 'bottomleft', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 1: ', all_asr[[8]]$total_cost ,' State Changes'))
dev.off()

#G_2
pdf(file = 'G_2_phylo-tree.pdf')
plot.phylo(all_trees[[9]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[9]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[9]][all_trees[[9]]$tip.label], levels(set_names_data[[9]])), cex = 0.5, piecol = colors) 
legend(x = 'bottomleft', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 2: ', all_asr[[9]]$total_cost ,' State Changes'))
dev.off()

#G_3
pdf(file = 'G_3_phylo-tree.pdf')
plot.phylo(all_trees[[10]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[10]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[10]][all_trees[[10]]$tip.label], levels(set_names_data[[10]])), cex = 0.5, piecol = colors) 
legend(x = 'center', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 3: ', all_asr[[10]]$total_cost ,' State Changes'))
dev.off()

#G_4
pdf(file = 'G_4_phylo-tree.pdf')
plot.phylo(all_trees[[11]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[11]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[11]][all_trees[[11]]$tip.label], levels(set_names_data[[11]])), cex = 0.5, piecol = colors) 
legend(x = 'center', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 4: ', all_asr[[11]]$total_cost ,' State Changes'))
dev.off()

#G_5
pdf(file = 'G_5_phylo-tree.pdf')
plot.phylo(all_trees[[12]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[12]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[12]][all_trees[[12]]$tip.label], levels(set_names_data[[12]])), cex = 0.5, piecol = colors) 
legend(x = 'center', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 5: ', all_asr[[12]]$total_cost ,' State Changes'))
dev.off()

#G_6
pdf(file = 'G_6_phylo-tree.pdf')
plot.phylo(all_trees[[13]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[13]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[13]][all_trees[[13]]$tip.label], levels(set_names_data[[13]])), cex = 0.5, piecol = colors) 
legend(x = 'center', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 6: ', all_asr[[13]]$total_cost ,' State Changes'))
dev.off()

#G_7
pdf(file = 'G_7_phylo-tree.pdf')
plot.phylo(all_trees[[14]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[14]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[14]][all_trees[[14]]$tip.label], levels(set_names_data[[14]])), cex = 0.5, piecol = colors) 
legend(x = 'center', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 7: ', all_asr[[14]]$total_cost ,' State Changes'))
dev.off()

#G_8
pdf(file = 'G_8_phylo-tree.pdf')
plot.phylo(all_trees[[15]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[15]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[15]][all_trees[[15]]$tip.label], levels(set_names_data[[15]])), cex = 0.5, piecol = colors) 
legend(x = 'bottomleft', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 8: ', all_asr[[15]]$total_cost ,' State Changes'))
dev.off()

#G_9
pdf(file = 'G_9_phylo-tree.pdf')
plot.phylo(all_trees[[16]], show.tip.label = T, label.offset = 0.09)
nodelabels(pie = all_asr[[16]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[16]][all_trees[[16]]$tip.label], levels(set_names_data[[16]])), cex = 0.5, piecol = colors) 
legend(x = 'bottomleft', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 9: ', all_asr[[16]]$total_cost ,' State Changes'))
dev.off()

#G_10
pdf(file = 'G_10_phylo-tree.pdf')
plot.phylo(all_trees[[17]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[17]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[17]][all_trees[[17]]$tip.label], levels(set_names_data[[17]])), cex = 0.5, piecol = colors) 
legend(x = 'bottomleft', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 10: ', all_asr[[17]]$total_cost ,' State Changes'))
dev.off()

#G_11
pdf(file = 'G_11_phylo-tree.pdf')
plot.phylo(all_trees[[18]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[18]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[18]][all_trees[[18]]$tip.label], levels(set_names_data[[18]])), cex = 0.5, piecol = colors) 
legend(x = 'center', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 11: ', all_asr[[18]]$total_cost ,' State Changes'))
dev.off()

#G_12
pdf(file = 'G_12_phylo-tree.pdf')
plot.phylo(all_trees[[19]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[19]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[19]][all_trees[[19]]$tip.label], levels(set_names_data[[19]])), cex = 0.5, piecol = colors) 
legend(x = 'bottomleft', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 12: ', all_asr[[19]]$total_cost ,' State Changes'))
dev.off()

#G_13
pdf(file = 'G_13_phylo-tree.pdf')
plot.phylo(all_trees[[20]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[20]]$ancestral_likelihoods, piecol = colors)
tiplabels(pie = to.matrix(set_names_data[[20]][all_trees[[20]]$tip.label], levels(set_names_data[[20]])), cex = 0.5, piecol = colors) 
legend(x = 'bottomleft', legend = state_names, pch = 22, pt.bg = colors, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 13: ', all_asr[[20]]$total_cost ,' State Changes'))
dev.off()

#######################
### Pigment plots
#######################
### plot UV asr onto the trees
setwd('C:/Users/vfica/Documents/Counterman_uv/R objects/Phylo-tree_asr_plots/asr_pig')
colors_pig <- c('black', 'cyan3', 'red2')
state_names_pig <- c('Absent', 'Forewing', 'Both Wings')

#E_1
pdf(file = 'E_1_phylo-tree_pig.pdf')
plot.phylo(all_trees[[1]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[1]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[1]][all_trees[[1]]$tip.label], levels(set_names_data_pig[[1]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'bottomleft', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 1: ', all_asr_pig[[1]]$total_cost ,' State Changes'))
dev.off()

#E_2
pdf(file = 'E_2_phylo-tree_pig.pdf') 
plot.phylo(all_trees[[2]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[2]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[2]][all_trees[[2]]$tip.label], levels(set_names_data_pig[[2]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'center', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 2: ', all_asr_pig[[2]]$total_cost ,' State Changes'))
dev.off()

#E_3
pdf(file = 'E_3_phylo-tree_pig.pdf')
plot.phylo(all_trees[[3]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[3]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[3]][all_trees[[3]]$tip.label], levels(set_names_data_pig[[3]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'center', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 3: ', all_asr_pig[[3]]$total_cost ,' State Changes'))
dev.off()

#E_4
pdf(file = 'E_4_phylo-tree_pig.pdf')
plot.phylo(all_trees[[4]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[4]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[4]][all_trees[[4]]$tip.label], levels(set_names_data_pig[[4]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'bottomleft', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 4: ', all_asr_pig[[4]]$total_cost ,' State Changes'))
dev.off()

#E_5
pdf(file = 'E_5_phylo-tree_pig.pdf') 
plot.phylo(all_trees[[5]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[5]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[5]][all_trees[[5]]$tip.label], levels(set_names_data_pig[[5]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'center', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 5: ', all_asr_pig[[5]]$total_cost ,' State Changes'))
dev.off()

#E_6
pdf(file = 'E_6_phylo-tree_pig.pdf') 
plot.phylo(all_trees[[6]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[6]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[6]][all_trees[[6]]$tip.label], levels(set_names_data_pig[[6]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'bottomleft', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 6: ', all_asr_pig[[6]]$total_cost ,' State Changes'))
dev.off()

#E_7
pdf(file = 'E_7_phylo-tree_pig.pdf')
plot.phylo(all_trees[[7]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[7]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[7]][all_trees[[7]]$tip.label], levels(set_names_data_pig[[7]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'center', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Euro Species Group 7: ', all_asr_pig[[7]]$total_cost ,' State Changes'))
dev.off()

#G_1
pdf(file = 'G_1_phylo-tree_pig.pdf')
plot.phylo(all_trees[[8]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[8]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[8]][all_trees[[8]]$tip.label], levels(set_names_data_pig[[8]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'bottomleft', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 1: ', all_asr_pig[[8]]$total_cost ,' State Changes'))
dev.off()

#G_2
pdf(file = 'G_2_phylo-tree_pig.pdf')
plot.phylo(all_trees[[9]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[9]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[9]][all_trees[[9]]$tip.label], levels(set_names_data_pig[[9]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'bottomleft', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 2: ', all_asr_pig[[9]]$total_cost ,' State Changes'))
dev.off()

#G_3
pdf(file = 'G_3_phylo-tree_pig.pdf')
plot.phylo(all_trees[[10]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[10]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[10]][all_trees[[10]]$tip.label], levels(set_names_data_pig[[10]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'center', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 3: ', all_asr_pig[[10]]$total_cost ,' State Changes'))
dev.off()

#G_4
pdf(file = 'G_4_phylo-tree_pig.pdf')
plot.phylo(all_trees[[11]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[11]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[11]][all_trees[[11]]$tip.label], levels(set_names_data_pig[[11]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'center', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 4: ', all_asr_pig[[11]]$total_cost ,' State Changes'))
dev.off()

#G_5
pdf(file = 'G_5_phylo-tree_pig.pdf')
plot.phylo(all_trees[[12]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[12]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[12]][all_trees[[12]]$tip.label], levels(set_names_data_pig[[12]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'center', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 5: ', all_asr_pig[[12]]$total_cost ,' State Changes'))
dev.off()

#G_6
pdf(file = 'G_6_phylo-tree_pig.pdf')
plot.phylo(all_trees[[13]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[13]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[13]][all_trees[[13]]$tip.label], levels(set_names_data_pig[[13]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'center', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 6: ', all_asr_pig[[13]]$total_cost ,' State Changes'))
dev.off()

#G_7
pdf(file = 'G_7_phylo-tree_pig.pdf')
plot.phylo(all_trees[[14]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[14]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[14]][all_trees[[14]]$tip.label], levels(set_names_data_pig[[14]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'center', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 7: ', all_asr_pig[[14]]$total_cost ,' State Changes'))
dev.off()

#G_8
pdf(file = 'G_8_phylo-tree_pig.pdf')
plot.phylo(all_trees[[15]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr[[15]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[15]][all_trees[[15]]$tip.label], levels(set_names_data_pig[[15]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'bottomleft', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 8: ', all_asr[[15]]$total_cost ,' State Changes'))
dev.off()

#G_9
pdf(file = 'G_9_phylo-tree_pig.pdf')
plot.phylo(all_trees[[16]], show.tip.label = T, label.offset = 0.09)
nodelabels(pie = all_asr_pig[[16]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[16]][all_trees[[16]]$tip.label], levels(set_names_data_pig[[16]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'bottomleft', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 9: ', all_asr_pig[[16]]$total_cost ,' State Changes'))
dev.off()

#G_10
pdf(file = 'G_10_phylo-tree_pig.pdf')
plot.phylo(all_trees[[17]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[17]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[17]][all_trees[[17]]$tip.label], levels(set_names_data_pig[[17]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'bottomleft', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 10: ', all_asr_pig[[17]]$total_cost ,' State Changes'))
dev.off()

#G_11
pdf(file = 'G_11_phylo-tree_pig.pdf')
plot.phylo(all_trees[[18]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[18]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[18]][all_trees[[18]]$tip.label], levels(set_names_data_pig[[18]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'center', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 11: ', all_asr_pig[[18]]$total_cost ,' State Changes'))
dev.off()

#G_12
pdf(file = 'G_12_phylo-tree_pig.pdf')
plot.phylo(all_trees[[19]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[19]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[19]][all_trees[[19]]$tip.label], levels(set_names_data_pig[[19]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'bottomleft', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 12: ', all_asr_pig[[19]]$total_cost ,' State Changes'))
dev.off()

#G_13
pdf(file = 'G_13_phylo-tree_pig.pdf')
plot.phylo(all_trees[[20]], show.tip.label = T, label.offset = 0.05)
nodelabels(pie = all_asr_pig[[20]]$ancestral_likelihoods, piecol = colors_pig)
tiplabels(pie = to.matrix(set_names_data_pig[[20]][all_trees[[20]]$tip.label], levels(set_names_data_pig[[20]])), cex = 0.5, piecol = colors_pig) 
legend(x = 'bottomleft', legend = state_names_pig, pch = 22, pt.bg = colors_pig, pt.cex = 1.5) # legend, ####### letters replaced with state #######
title(paste0('Grishin Species Group 13: ', all_asr_pig[[20]]$total_cost ,' State Changes'))
dev.off()

####### Bar plot the number of transitions of all asr by species group
#### collect state changes into dataframe for ggplot digesting
grab_changes_uv <- function(x){all_asr[[x]]$total_cost}
grab_changes_pig <- function(x){all_asr_pig[[x]]$total_cost}

state_changes_uv <- unlist(lapply(1:length(all_asr), FUN = grab_changes))
state_changes_uv <- as.data.frame(cbind(state_changes, tree_file_list))
state_changes_pig <- unlist(lapply(1:length(all_asr_pig), FUN = grab_changes_pig))
state_changes_pig <- as.data.frame(cbind(state_changes_pig, tree_file_list))

# UV fit x and y coordinates in ggplot aes and set stat as identity in geom_bar to produce bar graph
library('RColorBrewer')
ggplot(state_changes, aes(x = tree_file_list, y = state_changes_uv, fill = state_changes_uv)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = .07)) +
  xlab('Species Groups') +
  ylab('Number of Changes of UV') +
  ggtitle('State Changes across Species Groups')

# pig
ggplot(state_changes_pig, aes(x = tree_file_list, y = state_changes_pig, fill = state_changes_pig)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = .07)) +
  xlab('Species Groups') +
  ylab('Number of Changes of Pigment') +
  ggtitle('State Changes across Species Groups')

###########################
# Begin PCA https://www.datanovia.com/en/lessons/one-way-manova-in-r/ <- tips for data exploration
###########################
###################################################################
library('ape'); library('geiger'); library('phytools'); library('castor'); library('GGally'); library('car'); library('gtools'); library('dplyr'); library('ggplot2');library('pracma'); library('mvnormtest'); library('rstatix')

### Full dataset from Fenner-Counterman, 10-22-2020, directors commentary cut
full_data <- read.csv('C:/Users/vfica/Documents/Counterman_uv/coliadinae color patterns_10-22-2020_Copy of vince blind.csv')
pca_data <- read.csv('C:/Users/vfica/Documents/Counterman_uv/coliadinae color patterns_pca_10-28-2020_Copy of vince blind.csv')
pcadata <- as.data.frame(cbind.data.frame(pca_data$GS, pca_data$Species.Group, pca_data$Num_SG, pca_data$Forewing_UV, pca_data$Hindwing_UV, pca_data$Forewing_pig, pca_data$Hindwing_pig))
 
outpca <- prcomp(x = pcadata[3:7], center = T)
raw <- outpca$rotation %*% diag(outpca$sdev,length(outpca$sdev),length(outpca$sdev))
rot_raw_vmax <- varimax(raw)$loadings 
i_vmax <- t(pinv(rot_raw_vmax))
scores <- scale(pcadata[3:7]) %*% i_vmax #this is the verimax rotated scores analogous to outpca$x

  #pca plot
plot(outpca$x, col = pcadata$`pca_data$Species.Group`, pch = 20) #pch is the symbol on the plot
legend(x = 'topright', legend =  unique(pcadata$`pca_data$Species.Group`), col = unique(pcadata$`pca_data$Species.Group`), cex = .3, pch = 20)

  #verimax rotated loadings plot
plot(scores, col = pcadata$`pca_data$Species.Group`, pch = 20) #pch is the symbol on the plot
legend(x = 'topright', legend =  unique(pcadata$`pca_data$Species.Group`), col = unique(pcadata$`pca_data$Species.Group`), cex = .3, pch = 20)

### MANOVA
V1 <- scores$V1
V2 <- scores$V2
V3 <- scores$V3
V4 <- scores$V4
V5 <- scores$V5
model <- lm(cbind(V2, V3, V4, V5) ~ V1, scores)

Manova(model, test.statistic = 'Wilks')
  ### Results with Wilks
# Type II MANOVA Tests: Wilks test statistic
# Df test stat approx F num Df den Df    Pr(>F)    
# V1  1    0.2398   27.738      4     35 2.008e-10 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Manova(model, test.statistic = 'Pillai')
  ### Results with Pillai
# Type II MANOVA Tests: Pillai test statistic
# Df test stat approx F num Df den Df    Pr(>F)    
# V1  1    0.7602   27.738      4     35 2.008e-10 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# in each test statistic, significant p-values found between the species groups and uv and pigment scorings.

