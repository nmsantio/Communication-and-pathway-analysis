#########################################################################
# CELL CHAT Inference and analysis of cell-cell communication 
#########################################################################

# link https://github.com/sqjin/CellChat

# installing
install.packages("NMF")

# Python UMAP-learn package needs this
install.packages("reticulate")
devtools::install_github("jokergoo/circlize", force = TRUE)
devtools::install_github("jokergoo/ComplexHeatmap", force = TRUE)

# install cellchat
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocNeighbors")
devtools::install_github("sqjin/CellChat", force = TRUE)

# libraries
library(Seurat)
library(BiocNeighbors)
library(dplyr)
library(NMF)
library(CellChat)
library(patchwork)
library(ggalluvial)
library(ComplexHeatmap)
library(magrittr)
library(ggplot2)
library(igraph)
library(reticulate)
options(stringsAsFactors = FALSE)

#########################################################################
# PYTHON TESTING
# Should only be needed at first run to setup environment
# conda_create("r-reticulate")
# conda_install("r-reticulate", "umap-learn")

# This needs to be called before import
# use_condaenv("r-reticulate")
# reticulate::py_install(packages = "umap-learn")
# import does not work
# needed to use umap.method = "uwot" instead of "umap-learn"!!!

#########################################################################
# PREPROCESSING
#########################################################################

# SELECT DATA

#set database for ligand-receptor interactions
CellChatDB = CellChatDB.mouse # change to human if human data
showDatabaseCategory(CellChatDB)
CellChatDB.use = CellChatDB

#get the UMAP file with cluster names
seurat = readRDS(file = "Filepath\\Subset1_renamedC.rds")

# split into two sample sets
splitted = SplitObject(seurat, split.by = "sampleID")
ctrl = splitted$ctrl
cancer = splitted$cancer

# seurat to CellChat format
data.input = GetAssayData(ctrl, assay = "RNA", slot = "data") # normalized data matrix
labels = Idents(ctrl)
meta = data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat = createCellChat(object = data.input, meta = meta, group.by = "group")

# identify interactions based on the database
cellchat@DB = CellChatDB
cellchat = subsetData(cellchat) 
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)

# CELL-CELL COMMUNICATION NETWORK

# set network
cellchat = projectData(cellchat, PPI.mouse) #change to human if human data

# Filter 
cellchat = computeCommunProb(cellchat)
cellchat = filterCommunication(cellchat, min.cells = 10)

# identify pathways
cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)

# run, save and open each sample separately (here ctrl and cancer)
saveRDS(cellchat, file = "Filepath\\ctrlcellchat.rds")
ctrlcellchat = readRDS(file = "Filepath\\ctrlcellchat.rds")

saveRDS(cellchat, file = "Filepath\\cancercellchat.rds")
cancercellchat = readRDS(file = "Filepath\\cancercellchat.rds")

# rename if needed for further analysis (analysis to one at a time or both simultaneously)
cancercellchat = cellchat 
ctrlcellchat = cellchat 
cellchat = cancercellchat 
cellchat = ctrlcellchat

####################################################################
# ANALYSIS OF ONE SET AT A TIME
####################################################################

# VISUALIZATION OF CELL-CELL COMMUNICATION NETWORK

# general networks
groupSize = as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# single networks
pathways.show = c("NOTCH") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle") # circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord") # chord plot
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds") # heatmap

# ligand-receptor pair distance
pathways.show = c("NOTCH") 
netAnalysis_contribution(cellchat, signaling = pathways.show)

# ligand-receptor plots
pairLR.NOTCH = extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show = pairLR.NOTCH[1,] # show one ligand-receptor pair
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord") #change to circle if wanted

# ligand-receptor heatmaps (select interacting cluster nro by sources and targets commands)
netVisual_bubble(cellchat, remove.isolate = FALSE, sources.use = 4, targets.use = c(1:8)) 
netVisual_bubble(cellchat, remove.isolate = FALSE, sources.use = 5, targets.use = c(1:8)) 
netVisual_bubble(cellchat, signaling = c("NOTCH"), remove.isolate = FALSE) # selected pathways in all clusters

# ligand-receptor chord plot (select interacting cluster nro by sources and targets commands)
netVisual_chord_gene(cellchat, sources.use = c(7), targets.use = 5, legend.pos.x = 15) 
netVisual_chord_gene(cellchat, sources.use = c(7), targets.use = 3, legend.pos.x = 15) 
netVisual_chord_gene(cellchat, sources.use = c(7), targets.use = 8, legend.pos.x = 15) 

# gene expression blot
plotGeneExpression(cellchat, signaling = "CXCL") # only the ones with expression
plotGeneExpression(cellchat, signaling = "NOTCH", enriched.only = FALSE) # all

# SYSTEMS ANALYSIS

# sender and receiver clusters on one pathway
pathways.show = c("NOTCH") 
cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# all out- and incoming signals in all clusters
ht1 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# outgoing signals
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat = identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing") # colorful graph
netAnalysis_dot(cellchat, pattern = "outgoing")# dot plot table

# incoming signals
selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat = identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

# functional similarity
cellchat = computeNetSimilarity(cellchat, type = "functional") #calculate
cellchat = netEmbedding(cellchat, type = "functional", umap.method = "uwot") # umap
cellchat = netClustering(cellchat, type = "functional") # 2D classification
netVisual_embedding(cellchat, type = "functional", label.size = 3.5) # plot
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2) # zoom

# structural similarity
cellchat = computeNetSimilarity(cellchat, type = "structural") #calculate
cellchat = netEmbedding(cellchat, type = "structural", umap.method = "uwot") # umap
cellchat = netClustering(cellchat, type = "structural") # 2D classification
netVisual_embedding(cellchat, type = "structural", label.size = 3.5) # plot
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2) # zoom

##########################################################################
# ANALYSIS OF TWO SETS
##########################################################################

# COMBINE TWO SETS

# list of pathways
cellchat@netP$pathways

#combined analysis
object.list = list(ctrl = ctrlcellchat, cancer = cancercellchat)
cellchat = mergeCellChat(object.list, add.names = names(object.list))

# PREDICT CELL-CELL COMMUNICATION

# bar plot (overall view)
gg1 = compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 = compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# network (overall view) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# heatmap (overall view)
gg1 = netVisual_heatmap(cellchat)
gg2 = netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

# weighted cirle plot
weight.max = getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interactions - ", names(object.list)[i]))
}


# CONTEXT -DEPENDENT PATHWAYS

# functional similarity
cellchat = computeNetSimilarityPairwise(cellchat, type = "functional") # calculate
cellchat = netEmbedding(cellchat, type = "functional", umap.method = "uwot") # umap
cellchat = netClustering(cellchat, type = "functional") # map
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5) # plot

# structural similarity
cellchat = computeNetSimilarityPairwise(cellchat, type = "structural") # calculate
cellchat = netEmbedding(cellchat, type = "structural", umap.method = "uwot") # umap
cellchat = netClustering(cellchat, type = "structural") # map
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)  # plot
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2) 

#save and read
saveRDS(cellchat, file = "Filepath\\cellchat_twosamples.rds")
cellchat = readRDS(file = "Filepath\\cellchat_twosamples.rds")

# distances between datasets 
rankSimilarity(cellchat, type = "functional") 
rankSimilarity(cellchat, type = "structural") 

# context-dependent pathways
gg1 = rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 = rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

# outgoing signals 
i = 1
pathway.union = union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
ht1 = netAnalysis_signalingRole_heatmap(ht1, pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 14)
ht2 = netAnalysis_computeCentrality(object.list[[i + 1]], slot.name = "netP")
ht2 = netAnalysis_signalingRole_heatmap(ht2, pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# incoming signals 
i = 1
pathway.union = union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
ht1 = netAnalysis_signalingRole_heatmap(ht1, pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 14)
ht2 = netAnalysis_computeCentrality(object.list[[i + 1]], slot.name = "netP")
ht2 = netAnalysis_signalingRole_heatmap(ht2, pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 14)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# UP- AND DOWNREGULATED SIGNALLING PAIRS
# select cluster to compare (source) and clusters to which compare (targets)

# increased signals (ctrl vs cancer, 1 vs 2)
gg1 = netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 4 ,  comparison = c(1, 2), max.dataset = 2, title.name = "Upregulated in cancer", angle.x = 45, remove.isolate = T)
tiff("up_sourceall_targetcrec.tiff", units="in", width=5, height=6, res=300)
gg1 
dev.off()

gg1 = netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:8),  comparison = c(1, 2), max.dataset = 2, title.name = "Upregulated in cancer", angle.x = 45, remove.isolate = T)
gg1 


# decreased signals (ctrl vs cancer, 1 vs 2)
gg2 = netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 4,  comparison = c(1, 2), max.dataset = 1, title.name = "Downregulated in cancer", angle.x = 45, remove.isolate = T)
tiff("down_sourceall_targetcrec.tiff", units="in", width=5, height=7, res=300)
gg2
dev.off()

gg2 = netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:8),  comparison = c(1, 2), max.dataset = 1, title.name = "Downregulated in cancer", angle.x = 45, remove.isolate = T)
gg2

# individual pathways - circle plot (perform only functional, not structural, analysis before this!)
pathways.show = c("HIF") 
weight.max = getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # edge weight control
par(mfrow = c(1,2), xpd=TRUE)


for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

tiff("circle_NOTCH.tiff", units="in", width=10, height=7, res=300)
for (i in 1:length(object.list)) {netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))}
dev.off()


#Violinplot, one target pathway - RNA levels (figure labels: left set 1 ctrl, right set 2 cancer)
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("ctrl", "cancer")) # set factor level
plotGeneExpression(cellchat, signaling = "IFN", split.by = "datasets", colors.ggplot = T)

#save
saveRDS(cellchat, file = "cellchat.rds") 


