# Bradley June 2021
# Summarising the results of WGCNA within scRNASeq clusters

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/")
source('analysis.r')
library('Seurat')
library(ggsci)
library(rstatix)
library(optparse)
library(patchwork)
library(fastcluster)
library(dynamicTreeCut)
library(WGCNA)

# Set up
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, "WGCNA_within_clusters", sep = "/")

# generate the cluster list
clusters <- seq(1:12)
clusters <- clusters-1
clusters <- paste("cluster_", clusters, sep = "")


# Load the correlations
# 04/02/22 edit: For the paper, maybe not interested in sample grouping until later on, so can remove that from the p value and correlation matrices
remove_grouping <- F
cor_list <- vector("list", length = length(clusters))
names(cor_list) <- clusters
p_list <- vector("list", length = length(clusters)) 
names(p_list) <- clusters
for(cluster in 1:length(clusters)){
	tempP <- paste(pathOut, clusters[cluster], sep = "/")
	cor_list[[cluster]] <- read.csv(paste(tempP, "tables/moduleTraitCor.csv", sep = "/"))
	p_list[[cluster]]  <- read.csv(paste(tempP, "tables/moduleTraitPvalue.csv", sep = "/"))
	if(remove_grouping == T){
		cor_list[[cluster]] <- cor_list[[cluster]][,-which(colnames(cor_list[[cluster]]) == "turq_hub_cluster")]
		p_list[[cluster]] <- p_list[[cluster]][,-which(colnames(p_list[[cluster]]) == "turq_hub_cluster")]
	}
	rm(tempP)
}
 
##### Experimental - BH correction of pvals within clusters BEFORE JOINING
for(cluster in seq_along(p_list)){
	rownames(p_list[[cluster]]) <- p_list[[cluster]]$X
	p_list[[cluster]] <- p_list[[cluster]][,-1]
	BH_pvals <- p.adjust(as.matrix(p_list[[cluster]]), method = "BH")
	BH_pvals_df <- matrix(BH_pvals, ncol = ncol(p_list[[cluster]]), nrow = nrow(p_list[[cluster]]))
        rownames(BH_pvals_df) <- rownames(p_list[[cluster]])
        colnames(BH_pvals_df) <- colnames(p_list[[cluster]])
	p_list[[cluster]] <- BH_pvals_df
}

# Plot cluster 11
c11cor <- cor_list[[12]]
rownames(c11cor) <- cor_list[[12]][,1]
c11cor <- c11cor[,-1]
c11BHp <- p_list[[12]]
library('WGCNA',)
textMatrix = paste(signif(c11cor, 2), "\n(",
                   signif(c11BHp, 1), ")", sep = "");
dim(textMatrix) = dim(c11cor)
par(mar = c(6, 10, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = textMatrix,
               xLabels = colnames(c11cor),
               yLabels = rownames(c11cor),
               ySymbols = colnames(c11cor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



# First, have a look at the clusters which contain the majority of my FDR significant trans-eQTLs. Have to first load the geneInfo
gIs <- vector("list", length = length(clusters))
names(gIs) <- clusters
for(cluster in 1:length(clusters)){
	tempP <- paste(pathOut, clusters[cluster], sep = "/")
	gIs[[cluster]] <- read.csv(paste(tempP, "tables/geneInfo.csv", sep = "/"), row.names = 1);
}
trans <- c("LRMP", "SH2D6", "PSTPIP2", "HTR3E", "TRPM5", "HTR3C", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "AZGP1", "GNG13", "CAMP", "ANKHD1", "GIN1", "SPAG6")
# Find the proportion in each module from each cluster
prop_sum <- vector("list", length = length(clusters))
names(prop_sum) <- clusters
for(cluster in 1:length(clusters)){
	prop_sum[[cluster]] <- vector("list", length = length(levels(factor(gIs[[cluster]]$moduleColor))));
	names(prop_sum[[cluster]]) <- levels(factor(gIs[[cluster]]$moduleColor));
	for(module in 1:length(prop_sum[[cluster]])){
		prop_sum[[cluster]][[module]] <- gIs[[cluster]][gIs[[cluster]]$moduleColor == levels(factor(gIs[[cluster]]$moduleColor))[module],];
		prop_sum[[cluster]][[module]] <- data.frame(cluster.color = paste(clusters[cluster], levels(factor(gIs[[cluster]]$moduleColor))[module], sep = "."), 
			prop_RNASeq_trans_eQTLs = length(intersect(trans, rownames(prop_sum[[cluster]][[module]])))/length(trans), count = length(intersect(trans, rownames(prop_sum[[cluster]][[module]]))), intersect = paste(intersect(trans, rownames(prop_sum[[cluster]][[module]])), sep = "", collapse = ","))
	};
	prop_sum[[cluster]] <- do.call(rbind, prop_sum[[cluster]])
}


# Now subset for those that are >3
for(cluster in 1:length(clusters)){
	 prop_sum[[cluster]] <- prop_sum[[cluster]][prop_sum[[cluster]]$count > 3,]
}

# Put it together
prop_sum <- do.call(rbind, prop_sum)
rownames(prop_sum) <- prop_sum$cluster.color

# Now do the module Trait matrix, but subsetting for those that have at least one trans-eQTL present
for(cluster in 1:length(cor_list)){
	cor_list[[cluster]]$X <- gsub("ME", "", cor_list[[cluster]]$X)
	rownames(cor_list[[cluster]]) <- cor_list[[cluster]]$X
	cor_list[[cluster]] <- cor_list[[cluster]][,-which(colnames(cor_list[[cluster]]) == "X")]
}
cor <- do.call(rbind, cor_list)
cor <- cor[rownames(cor) %in% rownames(prop_sum),]


for(cluster in 1:length(p_list)){
	rownames(p_list[[cluster]]) <- gsub("ME", "", rownames(p_list[[cluster]]))
	rownames(p_list[[cluster]]) <- paste(names(p_list)[cluster], ".", rownames(p_list[[cluster]]), sep = "")
}

p <- do.call(rbind, p_list)
p <- p[rownames(p) %in% rownames(prop_sum),]

# Order together
cor <- cor[order(rownames(cor)),]
p <- p[order(rownames(p)),]
all(rownames(p) == rownames(cor))

cor <- as.matrix(cor)
p <- as.matrix(p)

# Now plot 
# Will display correlations and their p-values
textMatrix = paste(signif(cor, 2), "\n(",
                   signif(p, 2), ")", sep = "");
dim(textMatrix) = dim(cor)
par(mar = c(6, 6, 3, 3));

colnames(cor)[which(colnames(cor) == "turq_hub_cluster")] <- "blue_hub_cluster"

# Put ylabs on 2 lines
ylab <- rownames(cor)
ylab <- gsub("\\.", "\n", ylab)
par(mar = c(8, 6, 3, 3));
# Display the correlation values within a heatmap plot
pathOut <- "BH_analysis/Seurat/August_2021/Epithelial/WGCNA_within_clusters_Peter"
pdf(file = paste(pathOut, "ModuleTraitMatrix_BH_pvals_over1_transeQTL.pdf", sep = "/"))
labeledHeatmap(Matrix = cor,
               xLabels = colnames(cor),
               yLabels = ylab,
               ySymbols = rownames(cor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               cex.lab = 1.3,
               cex.lab.x = 1.3,
               cex.lab.y = 1.3,
               cex.main = 1.4,
               zlim = c(-1,1),
               main = paste("Module-trait matrix. Modules >3 trans-eQTLs"))
dev.off()


# Now calculate the intramodule connectivity for this cluster and all modules
# Identifying hub genes
# calculating the intramodular connectivity (the degree of co-expression of a given gene with respect to the genes of a particular module)
# Generate colours, a vector of length nGenes which the colour label for each gene
# 1. Load the expression data from each object (Likely to be long and slow
kIM_list <- vector("list", length = length(clusters))
geneList <- vector("list", length = length(clusters))
names(geneList) <- clusters
#Function to remove variables if they exist
ifrm <- function(obj, env = globalenv()) {
    obj <- deparse(substitute(obj))
    if(exists(obj, envir = env)) {
        rm(list = obj, envir = env)
    }
}
# Generate all of the geneLists
for(cluster in seq_along(gIs)){
	rm(pathOut)
	pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, "WGCNA_within_clusters", sep = "/")
	tempP <- paste(pathOut, clusters[cluster], sep = "/")
	ifrm(t_exprdata)
	ifrm(modcols)
	ifrm(sft)
	load(paste(tempP, "objects/UpTo_geneTraitSignificance.Rds", sep ="/"))
	kIM_list[[cluster]] <- intramodularConnectivity.fromExpr(t_exprdata, colors = moduleColors, power=sft$powerEstimate)
	rownames(kIM_list[[cluster]]) <- colnames(t_exprdata)
	# Attach the colour of the genes
	want <- c("colnames.t_exprdata.", "moduleColor")
	col_df <- gIs[[cluster]][, colnames(gIs[[cluster]]) %in% want]
	col_df <- col_df[match(rownames(kIM_list[[cluster]]), rownames(col_df)),]
	all(rownames(col_df) == rownames(kIM_list[[cluster]]))
	kIM_list[[cluster]] <- cbind(kIM_list[[cluster]], col_df)
	# make a geneList for every module within every cluster
	geneList[[cluster]] <- vector("list", length = length(levels(factor(kIM_list[[cluster]]$moduleColor))))
	names(geneList[[cluster]]) <- levels(factor(kIM_list[[cluster]]$moduleColor))
	for(module in seq_along(levels(factor(kIM_list[[cluster]]$moduleColor)))){
		geneList[[cluster]][[module]] <- kIM_list[[cluster]][kIM_list[[cluster]]$moduleColor == levels(factor(kIM_list[[cluster]]$moduleColor))[module],]$kWithin
		names(geneList[[cluster]][[module]]) <- kIM_list[[cluster]][kIM_list[[cluster]]$moduleColor == levels(factor(kIM_list[[cluster]]$moduleColor))[module],]$colnames.t_exprdata.
		geneList[[cluster]][[module]] <- geneList[[cluster]][[module]][!is.na(geneList[[cluster]][[module]])]
	}
	rm(tempP)
}


# 2. Generate the sets
# Now do an enrichment of Peter's trans-eQTLs against the genes of each module, ranked by their kWithin
trans_list <- vector("list", length = 1)
names(trans_list) <- c("Vaughan_Shaw_11q23.1_eQTLs")
trans_list[[1]] <- trans

# Now do the enrichment of modules
library(gage, lib.loc = "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
library(fgsea, lib.loc = "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
library(ggplot2, lib.loc = "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")
library(ggridges, lib.loc = "/exports/igmm/eddie/dunlop-lab/BradleyH/R/x86_64-pc-linux-gnu-library/4.0.2")


# Want to first test the individual gene significance  to C11orf53 expression
test <- gIs
for(cluster in seq_along(test)){
	# Subset for those from the interesting modules
	test[[cluster]]$cluster.mod <- paste(names(test)[cluster], test[[cluster]]$moduleColor, sep = ".")
	test[[cluster]] <- test[[cluster]][test[[cluster]]$cluster.mod %in% rownames(cor),]
	# Now subset those that are significant for correlations with C11orf53 expression
	test[[cluster]] <- test[[cluster]][test[[cluster]]$p.GS.C11orf53 < 0.05 & test[[cluster]]$GS.C11orf53 > 0.5,]
	test[[cluster]] <- test[[cluster]][rownames(test[[cluster]]) %in% trans,]
}
# PRint these. Colour code by these genes to show which genes are actually
sapply(test, rownames)





fgseaRes <- vector("list", length=length(clusters))
for(cluster in seq_along(clusters)){
	fgseaRes[[cluster]] <- vector("list", length = length(geneList[[cluster]]))
	for(module in seq_along(fgseaRes[[cluster]])){
		fgseaRes[[cluster]][[module]] <- fgsea(trans_list, geneList[[cluster]][[module]], minSize = 7, maxSize = 500, scoreType = "pos", eps=0);
		fgseaRes[[cluster]][[module]] <- as.data.frame(fgseaRes[[cluster]][[module]]);
		fgseaRes[[cluster]][[module]]$cluster.color <- rep(paste(clusters[cluster], names(geneList[[cluster]])[module], sep = "."), nrow(fgseaRes[[cluster]][[module]]))
	}
	fgseaRes[[cluster]] <- do.call(rbind, fgseaRes[[cluster]])
}

# Put all together: Note the pvals are calculated from within comparisons, so are nt corrected for the total number of tests
fgseaRes <- do.call(rbind, fgseaRes)
fgseaRes <- fgseaRes[,-which(colnames(fgseaRes) == "leadingEdge")]
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, "WGCNA_within_clusters", sep = "/")
write.csv(fgseaRes, paste(pathOut, "fgsea.trans_in_all_modules.kIM.csv", sep = "/"))

	


##### Now do the same but for gIs, ranking by GS.C11orf53
geneList <- vector("list", length = length(clusters))
names(geneList) <- clusters
for(cluster in seq_along(clusters)){
	geneList[[cluster]] <- vector("list", length = length(levels(factor(gIs[[cluster]]$moduleColor))))
        names(geneList[[cluster]]) <- levels(factor(gIs[[cluster]]$moduleColor))
        for(module in seq_along(levels(factor(gIs[[cluster]]$moduleColor)))){
                geneList[[cluster]][[module]] <- gIs[[cluster]][gIs[[cluster]]$moduleColor == levels(factor(gIs[[cluster]]$moduleColor))[module], 
			which(colnames(gIs[[cluster]]) == "GS.C11orf53")]
                names(geneList[[cluster]][[module]]) <- gIs[[cluster]][gIs[[cluster]]$moduleColor == levels(factor(gIs[[cluster]]$moduleColor))[module],]$colnames.t_exprdata.
                geneList[[cluster]][[module]] <- geneList[[cluster]][[module]][!is.na(geneList[[cluster]][[module]])]
        }
}


fgseaRes <- vector("list", length=length(clusters))
for(cluster in seq_along(clusters)){
        fgseaRes[[cluster]] <- vector("list", length = length(geneList[[cluster]]))
        for(module in seq_along(fgseaRes[[cluster]])){
                fgseaRes[[cluster]][[module]] <- fgsea(trans_list, geneList[[cluster]][[module]], minSize = 7, maxSize = 500, scoreType = "std", eps=0);
                fgseaRes[[cluster]][[module]] <- as.data.frame(fgseaRes[[cluster]][[module]]);
                fgseaRes[[cluster]][[module]]$cluster.color <- rep(paste(clusters[cluster], names(geneList[[cluster]])[module], sep = "."), nrow(fgseaRes[[cluster]][[module]]))
        }
        fgseaRes[[cluster]] <- do.call(rbind, fgseaRes[[cluster]])
}

# Put all together: Recalculate the pvals for the number of tests where at least one eQTL was found
fgseaRes <- do.call(rbind, fgseaRes)
fgseaRes$BH.pval <- p.adjust(fgseaRes$pval, method = "BH")
fgseaRes <- fgseaRes[,-which(colnames(fgseaRes) == "leadingEdge")]
for(c in 1:ncol(fgseaRes)){
	if(class(fgseaRes[,c]) == "numeric"){
		fgseaRes[,c] <- signif(fgseaRes[,c], 3)
	}
}

pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, "WGCNA_within_clusters", sep = "/")
write.csv(fgseaRes, paste(pathOut, "fgsea.trans_in_all_modules.GS.C11orf53.csv", sep = "/"))



# Combining the table of GS and p.GS of GeneInfos for the interesting modules across clusters
trans_gis <- gIs
for(cluster in 1:length(trans_gis)){
	trans_gis[[cluster]]$cluster.colour <- paste(names(trans_gis)[cluster], trans_gis[[cluster]]$moduleColor, sep = ".")
	trans_gis[[cluster]] <- trans_gis[[cluster]][rownames(trans_gis[[cluster]]) %in% trans_list[[2]] & trans_gis[[cluster]]$cluster.colour %in% rownames(cor), 
		 colnames(trans_gis[[cluster]]) %in% c("cluster.colour", "GS.C11orf53", "p.GS.C11orf53")]
	trans_gis[[cluster]]$gene <- rownames(trans_gis[[cluster]])
	rownames(trans_gis) <- NULL
}
trans_gis <- do.call(rbind, trans_gis)
rownames(trans_gis) <- NULL
write.csv(trans_gis, paste(pathOut, "FDR_sig_trans_eQTL_GI_C11orf53_int_modules.csv", sep = "/"))


#### Exploratory - GSEA of individual gene correlation with C11orf53 and trans-eQTLs within clusters, not module
geneList <- vector("list", length = length(clusters))
names(geneList) <- clusters
for(cluster in seq_along(clusters)){
	geneList[[cluster]] <- gIs[[cluster]]$GS.C11orf53
	names(geneList[[cluster]]) <- rownames(gIs[[cluster]])
#	geneList <- geneList[order(-geneList)]
	}
	
# Perform enrichment
fgseaRes <- vector("list", length=length(clusters))
for(cluster in seq_along(clusters)){
	fgseaRes[[cluster]] <- fgsea(trans_list, geneList[[cluster]], minSize = 1, maxSize = 500, scoreType = "std", eps=0)
	fgseaRes[[cluster]] <- as.data.frame(fgseaRes[[cluster]])
	rownames(fgseaRes[[cluster]]) <- rep(clusters[cluster], nrow(fgseaRes[[cluster]]))
}
fgseaRes_all <- do.call(rbind, fgseaRes)

# manually p.adjust
fgseaRes_all$BH.padj <- p.adjust(fgseaRes_all$pval, method = "BH")

# Save
write.csv(fgseaRes_all[,-8], "BH_analysis/Seurat/August_2021/Epithelial/WGCNA_across_clusters_Peter/tables/gsea_trans_eQTL_cor_C11orf53.csv")

# plot
plotEnrichment(trans_list[[1]], geneList[[3]])


# Plot
library(enrichplot)
library(clusterProfiler)



