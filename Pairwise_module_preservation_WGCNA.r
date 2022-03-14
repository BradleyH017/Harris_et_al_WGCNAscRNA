# Bradley Feb 2022
# Module preservation of cluster 11 black module in all other clusters by differential WGCNA
# Following tutorial V from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/ModulePreservation/Tutorials/

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

# Load the data for cluster 11. To see how this has been generated, look at WGCNA_within_clusters.r script
# Note: For this, we had utilised the 5000 most variable genes from the pseudo-bulked expression 
c=11
pathOut <- paste(pathOut, "/cluster_", c, sep = "")
load(paste(pathOut, "objects/UpTo_geneTraitSignificance.Rds", sep = "/"))
# Extract expression data and transpose
datExpr11 <- exprdata
datExpr11 <- as.data.frame(t(datExpr11))
# Extract the module colours
colors11 <- moduleColors
names(colors11) <- colnames(datExpr11)

# optional - exploring preservation of 11q23.1 trans-eQTL targets in cluster 11 specifically. 
# Edit 10/02/22 - Have exluded CAMP, PSTPIP2 and PLCG2 as these are not significantly correlated (p<0.05, cor>0.5 with C11orf53 within the cluster 11 black module)
trans_only <- T
if(trans_only == T){
	# Put trans genes in a module called orange
	trans <- c("LRMP", "SH2D6", "HTR3E", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PTGS1", "IL17RB", "AZGP1", "GNG13")
	indx <- grep("black", colors11)
	colors11[names(colors11) %in% trans] <- "orange"
}

# Extract dissTOM
dissTOM11_old <- dissTOM
# Extract softPower
softPower11 <- softPower

# Now load the results for another cluster. Doing this as an array, so want to load the cluster list number  using optparse
option_list = list(make_option(c("-c", "--cluster"), action = "store", default = NA, type ="character", help="cluster number"))
opt = parse_args(OptionParser(option_list=option_list))
c = opt$c;
print(paste("The option for cluster to study is", c, sep = " "))

# Load data
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, "WGCNA_within_clusters", sep = "/")
pathOut <- paste(pathOut, "/cluster_", c, sep = "")
load(paste(pathOut, "objects/UpTo_geneTraitSignificance.Rds", sep = "/"))
# Extract expression data and transpose
datExpr_test <- exprdata
datExpr_test <- as.data.frame(t(datExpr_test))
# Extract the module colours
colors_test <- moduleColors
names(colors_test) <- colnames(datExpr_test)
# Extract softPower old
softPower_test_old <- softPower


# Check the common genes - This is very low, potentially too low for the analysis
print(paste("The percentage of common genes is ", (length(intersect(colnames(datExpr_test), colnames(datExpr11)))/ncol(datExpr11)*100), "%", sep = ""))

# Reduce the expression and colours by the common genes 
common <- intersect(colnames(datExpr_test), colnames(datExpr11))
datExpr11 <- datExpr11[,colnames(datExpr11) %in% common]
datExpr_test <- datExpr_test[,colnames(datExpr_test) %in% common]
colors11 <- colors11[names(colors11) %in% common]
colors_test <- colors_test[names(colors_test) %in% common]

# Cluster the genes from both modules. Use the same sft as with cluster 11
# Need to recalculate the dissTOM for cluster 11 as am now using reduced number of genes
dissTOM11 <- 1-TOMsimilarityFromExpr(datExpr11, power = softPower11, TOMType = "signed");
tree11 <- hclust(as.dist(dissTOM11), method = "average");
dissTOM_test = 1-TOMsimilarityFromExpr(datExpr_test, power = softPower11, TOMType = "signed");
tree_test <- hclust(as.dist(dissTOM_test), method = "average")

# For plots and results from here, regenerate the pathOut
pathOut <- paste("BH_analysis/Seurat/August_2021/Epithelial/WGCNA_within_clusters/MP_cluster11_in_others/", "cluster_", c, sep = "")

pdf(file=paste(pathOut, "Dendrogram_comparison.pdf", sep="/"), w=11, h=7)
# Set up appropriate screen sectioning 
layout(matrix(c(1:4), 4, 1), heights = rep(c(0.8, 0.2), 2)); 
# Plot the cluster 11 dendrogram 
plotDendroAndColors(tree11, colors11, "Cluster_11 modules", main = "Cluster_11 gene dendrogram and Cluster_11 module colors", dendroLabels = FALSE, setLayout = FALSE, marAll = c(2,10,3,0), cex.colorLabels = 1.4, cex.main = 2, cex.lab = 1.4, cex.axis = 1.4, addGuide = TRUE); 
# Plot the cluster 3 dendrogram with cluster 11 module colors 
plotDendroAndColors(tree_test, colors11, paste("Cluster_", c, " modules", sep = ""), main = paste("Cluster_", c, " gene dendrogram and Cluster_11 module colors", sep = ""), dendroLabels = FALSE, setLayout = FALSE, marAll = c(2,10,3,0), cex.colorLabels = 1.4, cex.main = 2, cex.lab = 1.4, cex.axis = 1.4, addGuide = TRUE);
dev.off()


# Set up module preservation stats
setLabels = c("Cluster_11", paste("Cluster", c, sep = "_"))
multiExpr = list(Cluster_11 = list(data = datExpr11), test_cluster = list(data = datExpr_test))
multiColor = list(Cluster_11 = colors11, test_cluster = colors_test)
nSets = 2

# Calculate module preservation
mp = modulePreservation(multiExpr, multiColor, referenceNetworks = c(1:2), nPermutations = 200, randomSeed = 1, verbose = 3)
save(mp, file = paste(pathOut, "/mp_cluster11_in_cluster_", c, ".Rds", sep = ""))

# Analyse preservation
ref = 1 
test = 2 
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1]) 
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")], signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
write.csv(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")], signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)), paste(pathOut, "statsObs_statsZ.csv", sep = "/"))

# Module labels and module sizes are also contained in the results 
modColors = rownames(mp$preservation$observed[[ref]][[test]]) 
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1]; 
# leave grey and gold modules out 
plotMods = !(modColors %in% c("grey", "gold")); 
# Text labels for points 
text = modColors[plotMods];

# Auxiliary convenience variable 
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2]) 
# Main titles for the plot 
mains = c("Preservation Median rank", "Preservation Zsummary"); 

# Start the plot 
pdf(file=paste(pathOut, "Pres_median_rank_z_sum_plot.pdf", sep = "/"), w=10, h=5)
par(mfrow = c(1,2)) 
par(mar = c(4.5,4.5,2.5,1)) 
for (p in 1:2) { 
	min = min(plotData[, p], na.rm = TRUE); 
	max = max(plotData[, p], na.rm = TRUE); 
	# Adjust ploting ranges appropriately 
	if (p==2) { 
		if (min > -max/10) {
			min = -max/10 
			ylim = c(min -0.1 * (max-min), max + 0.1 * (max-min)) 
			} else 
				ylim = c(max + 0.1 * (max-min), min -0.1 * (max-min)) 
				plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21, main = mains[p], cex = 2.4, ylab = mains[p], xlab = "Module size", log = "x", ylim = ylim, xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4) 
				labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08); 
			}	
			# For Zsummary, add threshold lines 
			if (p==2) { 
				abline(h=0) 
				abline(h=2, col = "blue", lty = 2) 
				abline(h=10, col = "darkgreen", lty = 2) 
			} 
		}
dev.off()



# plot all the density and connectivity plots in a single plot
# Re-initialize module color labels and sizes 
modColors = rownames(statsZ) 
moduleSizes = mp$quality$Z[[ref]][[test]][, 1]; 

# Exclude improper modules 
plotMods = !(modColors %in% c("grey", "gold")); 

# Create numeric labels for each module 
labs = match(modColors[plotMods], standardColors(50)); 

# Start the plot: open a suitably sized graphical window and set sectioning and margins. Alternatively, # plot into a pdf file. 
pdf(file=paste(pathOut, "PreservationZStatistics.pdf", sep = "/"), w=12.5, h=9)
par(mfrow = c(4,5)) 
par(mar = c(3,3,2,1)) 
par(mgp = c(1.6, 0.4, 0)); 
for (s in 1:ncol(statsZ)) { 
	min = min(statsZ[plotMods, s], na.rm = TRUE); 
	max = max(statsZ[plotMods, s], na.rm = TRUE); 
	if (min > -max/12) {
		min = -max/12 
	}
	plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21, main = colnames(statsZ)[s], cex = 2.2, ylab = colnames(statsZ)[s], xlab = "Module size", log = "x", ylim = c(min -0.1 * (max-min), max + 0.1 * (max-min)), xlim = c(30, 800), cex.lab = 1.2, cex.axis = 1.2)
	labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 1, offs = 0.06); 
	#text(moduleSizes[-1], statsZ[-c(1:2), s], labels = letter[-c(1:2)], col = "black"); 
	#modColors[-2]); 
	abline(h=0) 
	abline(h=2, col = "blue", lty = 2) 
	abline(h=10, col = "darkgreen", lty = 2) 
} 
dev.off()




