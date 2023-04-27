#######################################################################################################
# GSEA: GO AND KEGG

# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/ 

library(clusterProfiler)
library(pathview)
library(enrichplot)
library(ggplot2)
library(pathview)

# open file or import manually and rename
# cluster markers from integrated data identified by FindAllMarkers function
# cluster markers from comparison of sample 1 to sample 2 identified by FindMarkers function based on sampleID 
DEgenes = read.delim("Filepath\\Markerlist_cluster1", header = TRUE)
DEgenes = Markerlist_cluster1

# remove pval over 0.05 and sort by avg_log2FC
genelist1 = DEgenes[DEgenes$p_val < 0.05,]
genelist2 = DEgenes[order(-DEgenes$avg_log2FC),]

# select avg_log2FC values and name list 
gene_list = genelist2$avg_log2FC
names(gene_list) = rownames(genelist2)

# sort again just in case
gene_list = sort(gene_list, decreasing = TRUE)

##########################
#          GO
##########################

# select ontology, keyType, size, cutoff and organism according to your needs
gse = gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 10, 
             maxGSSize = 200, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Mm.eg.db", 
             pAdjustMethod = "none")

require(DOSE)

# select number of pathways by showCategory
# save as tiff or svg

# DOTPLOT
tiff("GOdot_cluster1.tiff", units="in", width=6, height=5, res=300)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

# CNETPLOTs

tiff("GOcnet_cluster1.tiff", units="in", width=10, height=15, res=300)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 8)
dev.off()

svg(filename="GOcnet_cluster1.svg",     width=12,     height=13,     pointsize=8)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 8)
dev.off()

##########################
#         KEGG
##########################

# prepare gene list as mentioned above
# convert gene IDs for gseKEGG function
# note: you might lose some genes if all IDs are not converted!
ids = bitr(names(genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

# remove duplicate IDS 
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe, which has only the genes which were successfully mapped using the bitr function above
# for this purpose, import manually and rename, choose pval, order FC, then sort by succesfull renaming

# import manually and rename
DEgenes = Cluster1_forKEGG #use numbers as rownames
genelist1 = DEgenes[DEgenes$p_val < 0.05,]
genelist2 = genelist1[order(-DEgenes$avg_log2FC),]
entrezlist = genelist1[genelist1$row.names %in% ids$SYMBOL,]

# Create a new column with the corresponding ENTREZ IDs
entrezlist$ENTREZID = dedup_ids$ENTREZID

# Create a vector 
kegg_gene_list = entrezlist$avg_log2FC

# Name vector with ENTREZ ids
names(kegg_gene_list) = entrezlist$ENTREZID

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "mmu"                            #mmu mouse, hsa human
kegg = gseKEGG(geneList     = kegg_gene_list,    # select organism, size, cutoff
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 5,
               maxGSSize    = 200,
               pvalueCutoff = 0.05, 
               pAdjustMethod = "none")

# select number of pathways by showCategory
# save as tiff or svg

dotplot(kegg, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
cnetplot(kegg, categorySize="pvalue", foldChange=gene_list, showCategory = 10)

# DOTPLOT
tiff("KEGGdot_cluster1.tiff", units="in", width=6, height=5, res=300)
dotplot(kegg, showCategory = 15, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
dev.off()

# CNETPLOT

tiff("KEGGcnet_cluster1.tiff", units="in", width=11, height=8, res=300)
cnetplot(kegg, categorySize="pvalue", foldChange=gene_list, showCategory = 10)
dev.off()

svg(filename="KEGGcnet_cluster1.svg",     width=11,     height=8,     pointsize=8)
cnetplot(kegg, categorySize="pvalue", foldChange=gene_list, showCategory = 10)
dev.off()

# PATHVIEW 
dme = pathview(gene.data=kegg_gene_list, pathway.id="mmu04218", species = kegg_organism)   
knitr::include_graphics("mmu04218.cellularsenescence.png") 

#######################################################################################################
# REACTOME PATHWAY ANALYSIS

library(ReactomePA)

# prepare gene list as mentioned above
# Convert gene IDs as for KEGG
# Create a vector with the entrez id list
reactomegenelist = ids$ENTREZID

# select number of pathways by showCategory, organism and size
x = enrichPathway(reactomegenelist, organism = "mouse", minGSSize = 10, maxGSSize =200, readable = FALSE)

barplot(x, showCategory=8)

#save as tiff
tiff("reactomePA_cluster1.tiff", units="in", width=6, height=6, res=300)
barplot(x, showCategory=10)
dev.off()


