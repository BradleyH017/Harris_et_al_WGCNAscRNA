# Bradley September 2021
# DE testing within clusters divided by high/low expressors across the turquoise hub cluster grouping (defined by WGCNA)

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, sep = "/")

# Load the object
load(paste(pathOut, "objects/Processed_seur.integrated.clusteredRds", sep = "/"))

# Add the turquoise hub cluster grouping into the Seurat meta data
temp <- seur.integrated@meta.data
low <- c("N20.EpiB", "N20.EpiA", "N8.EpiB", "N13.EpiA", "N13.EpiB")
high <- setdiff(levels(factor(seur.integrated@meta.data$Sample)), low)
temp$turq_hub_group <- rep("", nrow(temp))
for(r in 1:nrow(temp)){
	temp$turq_hub_group[r] <- ifelse(temp$Sample[r] %in% low, "Low", "High")
}
group_vec <- temp$turq_hub_group
names(group_vec) <- rownames(temp)
seur.integrated <- AddMetaData(seur.integrated, group_vec, col.name = "turq_hub_group")

# Now seperate Seurat object into clusters, to loop thorugh and perform DE within clusters, but across grouping
Idents(seur.integrated) <- "turq_hub_group"
seur_list <- SplitObject(seur.integrated, split.by = "seurat_clusters")
res.list <- vector("list", length = length(list))
for(cluster in seq_along(seur_list)){
	print(paste("Doing cluster", names(seur_list)[cluster]));
	res.list[[cluster]] <- FindMarkers(seur_list[[cluster]], ident.1 = "High", ident.2 = "Low" , test.use = "MAST");
	head(res.list[[cluster]]);
	write.csv(res.list[[cluster]], paste(pathOut, "/tables/DE_within_clusters/Turq_hub_High_Low_cluster_", names(seur_list)[cluster], ".csv", sep = ""));
	print(paste("Done cluster", names(seur_list)[cluster]))
}
res <- do.call(rbind, res.list)
write.csv(res, paste(pathOut, "tables/DE_within_clusters/All_clusters_Turq_hub_High_Low.csv", sep = "/"))

# Now do some enrichment of the results
# KEGG over representation test
# Load packages
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(biomaRt)
library(ggplot2)

# Set up
geneList <- vector("list", length = length(seur_list))
kegg.res <- vector("list", length = length(seur_list))
names(kegg.res) <- names(seur_list)

# Loop through -- FINISH
for(cluster in seq_along(seur_list)){
	# Extract
	temp <- data.frame(gene = rownames(res.list[[cluster]]), fc = res.list[[cluster]]$avg_log2FC)
	# Seperate
	geneList[[cluster]] <- vector("list", length = 2)
	names(geneList[[cluster]]) <- c("Up", "Down")
	geneList[[cluster]][['Up']] <- temp[temp$fc>0.5,]
	geneList[[cluster]][['Down']] <- temp[temp$fc<(-0.5),]
	kegg.res[[cluster]] <- vector("list", length = 2)
	names(kegg.res[[cluster]]) <- c("Up", "Down")
	# Convert to entrez ID
	for(set in c(1,2)){
		HgeneList2entrez_pLI <- AnnotationDbi::mapIds(keys = as.character(geneList[[cluster]][[set]]$gene),
                                          x = org.Hs.eg.db, keytype = "SYMBOL", column = "ENTREZID",
                                          multiVals = "first")
		HgeneList2entrez_pLI <- HgeneList2entrez_pLI[!is.na(names(HgeneList2entrez_pLI))]
		DF <- data.frame(gene = names(HgeneList2entrez_pLI), h_entrez = HgeneList2entrez_pLI)
		DF <- merge(geneList[[cluster]][[set]], DF, by = "gene")
		DF <- DF[!is.na(DF$h_entrez),] 
		geneList[[cluster]][[set]] <- DF$fc
		names(geneList[[cluster]][[set]]) <- DF$h_entrez
		# Now do the enrichment
		kegg.res[[cluster]][[set]] <- enrichKEGG(gene = names(geneList[[cluster]][[set]]), organism = 'hsa',
                                pvalueCutoff = 0.01)
		kegg.res[[cluster]][[set]] <- as.data.frame(kegg.res[[cluster]][[set]]@result);
		kegg.res[[cluster]][[set]] <- kegg.res[[cluster]][[set]][kegg.res[[cluster]][[set]]$p.adjust < 0.05,]
		write.csv(kegg.res[[cluster]][[set]], paste(pathOut, "/tables/DE_within_clusters/KEGG/cluster_", names(seur_list)[cluster], ".", names(kegg.res[[cluster]])[set], ".csv", sep = ""))
		}	
}













