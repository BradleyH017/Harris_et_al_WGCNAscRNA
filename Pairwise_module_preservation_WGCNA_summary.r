# Bradley Feb 2022
# Module preservation of cluster 11 orange (auxiliary) module in all other clusters by differential WGCNA - SUMMARY
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


# Load the results of the other clusters
clusters <- seq(1:11)
clusters <- clusters -1
clusters <- paste("cluster_", clusters, sep = "")

# Set up list
mp_res <- vector("list", length = length(clusters))
names(mp_res) <- clusters
# Extract preservation stats
statsObs <- vector("list", length = length(clusters))
names(statsObs) <- clusters
# Extract cor.KME 
 
# 1. ~~~~~~~ Load all module Preservation
for(c in seq_along(clusters)){
	load(paste("BH_analysis/Seurat/August_2021/Epithelial/WGCNA_within_clusters/MP_cluster11_in_others/", clusters[c], "/", "mp_cluster11_in_", clusters[c], ".Rds", sep = ""))
	mp_res[[c]] <- mp
	rm(mp)
	ref = 1 
	test = 2
	statsObs[[c]] <- cbind(mp_res[[c]]$quality$observed[[ref]][[test]][, -1], mp_res[[c]]$preservation$observed[[ref]][[test]][, -1])
	colnames(statsObs[[c]]) <- paste(clusters[c], colnames(statsObs[[c]]), sep = ".")
	
}


# 2. ~~~~~~~~~ Compare the number of 'orange' genes included in the analysis
# Load the data for cluster 11. To see how this has been generated, look at WGCNA_within_clusters.r script
# Note: For this, we had utilised the 5000 most variable genes from the pseudo-bulked expression 
c=11
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021/Epithelial/WGCNA_within_clusters", sep = "/")
pathOut <- paste(pathOut, "/cluster_", c, sep = "")
load(paste(pathOut, "objects/UpTo_geneTraitSignificance.Rds", sep = "/"))
# Extract expression data and transpose
datExpr11 <- exprdata
datExpr11 <- as.data.frame(t(datExpr11))
# Extract the module colours
colors11 <- moduleColors
names(colors11) <- colnames(datExpr11)
trans_only <- T
if(trans_only == T){
	# Put trans genes in a module called orange - EXCLUDE PLCG2, PSTPIP2, CAMP as these weren't significant themselves
	trans <- c("LRMP", "SH2D6", "HTR3E", "ALOX5", "OGDHL", "BMX", "MATK", "SH2D7", "PTGS1", "IL17RB", "AZGP1", "GNG13")
	indx <- grep("black", colors11)
	colors11[names(colors11) %in% trans] <- "orange"
}
orange11 <- colors11[colors11 == "orange"]

# Load in expression from each cluster
orange_list <- list()
for(c in seq_along(clusters)){
	pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021/Epithelial/WGCNA_within_clusters", sep = "/")
	pathOut <- paste(pathOut, clusters[c], sep = "/")
	load(paste(pathOut, "objects/UpTo_geneTraitSignificance.Rds", sep = "/"))
	orange_list[[c]] <- intersect(rownames(exprdata), names(orange11))
	rm(pathOut, exprdata)
}
names(orange_list) <- clusters

# 3. ~~~~~~~~~~~~~ Compare kME of interesting genes across clusters
# Extract kME
kme_list <- vector("list", length = length(clusters))
names(kme_list) <- clusters
test <- as.data.frame(mp_res[[1]]$preservation$observed[[1]][[2]])
indx_cor <- which(colnames(test) == "cor.kME")
test2 <- as.data.frame(mp_res[[1]]$preservation$log.pBonf[[1]][[2]])
indx_b <-  which(colnames(test2) == "log.p.Bonf.cor.kME")
test3 <- as.data.frame(mp_res[[1]]$preservation$log.p[[1]][[2]])
indx_p <- which(colnames(test3) == "log.p.cor.kME")
modules <- rownames(test)
for(c in seq_along(clusters)){
	kme_list[[c]] <- data.frame(cor.kME = mp_res[[c]]$preservation$observed[[1]][[2]][,indx_cor],
									log.p.cor.kME = mp_res[[c]]$preservation$log.p[[1]][[2]][,indx_p],
									log.p.Bonf.cor.kME = mp_res[[c]]$preservation$log.pBonf[[1]][[2]][,indx_b])
	rownames(kme_list[[c]]) <- modules
	kme_list[[c]] <- kme_list[[c]][grep("orange", rownames(kme_list[[c]])),]
}
kme_all <- do.call(rbind, kme_list)
kme_all$adjPval <- kme_all$log.p.Bonf.cor.kME*-1
kme_all$rawPval <- 1-(10^(kme_all$log.p.cor.kME))
kme_all$cluster <- rownames(kme_all)
kme_all$overlap <- sapply(orange_list, length)


## Check colour palette
library(scales)
#show_col(hue_pal()(12))

# Plot this - colour by the scRNAseq cluster
temp <- levels(factor(kme_all$cluster))
temp <- gsub("cluster_", "", temp)
temp <- as.numeric(temp)
order <- order(temp)
ordered <- paste("cluster_", temp[order], sep = "")
kme_all$cluster <- factor(kme_all$cluster, levels = ordered)
ggplot(kme_all, aes(x=rawPval, y= cor.kME, color=cluster)) + 
	geom_point(aes(size=overlap)) + 
	theme_bw() + 
	xlim(c(0,1)) + 
	ylim(c(0,1)) + 
	xlab("1-cor.kME pvalue") + 
	ylab("cor.kME") + 
	geom_vline(xintercept = 0.95, linetype="dotted", color = "red", size=1) + 
	scale_color_manual(values=hue_pal()(12)) + 
	theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size = 20),  axis.title.y = element_text(size = 24), axis.title.x = element_text(size = 24), plot.title=element_text(face="bold", size = 20), legend.title=element_text(size=15), legend.text=element_text(size=12)) + 
	ggtitle("Preservation of cluster 11 black trans-eQTL \ntarget correlations across clusters")



# Extract kIM
kim_list <- vector("list", length = length(clusters))
names(kim_list) <- clusters
test <- as.data.frame(mp_res[[1]]$preservation$observed[[1]][[2]])
indx_cor <- which(colnames(test) == "cor.kIM")
test2 <- as.data.frame(mp_res[[1]]$preservation$log.pBonf[[1]][[2]])
indx_b <-  which(colnames(test2) == "log.p.Bonf.cor.kIM")
test3 <- as.data.frame(mp_res[[1]]$preservation$log.p[[1]][[2]])
indx_p <- which(colnames(test3) == "log.p.cor.kIM")
modules <- rownames(test)
for(c in seq_along(clusters)){
	kim_list[[c]] <- data.frame(cor.kIM = mp_res[[c]]$preservation$observed[[1]][[2]][,indx_cor],
									log.p.cor.kIM = mp_res[[c]]$preservation$log.p[[1]][[2]][,indx_p],
									log.p.Bonf.cor.kIM = mp_res[[c]]$preservation$log.pBonf[[1]][[2]][,indx_b])
	rownames(kim_list[[c]]) <- modules
	kim_list[[c]] <- kim_list[[c]][grep("orange", rownames(kme_list[[c]])),]
}
kim_all <- do.call(rbind, kim_list)
kim_all$adjPval <- 10^(kim_all$log.p.Bonf.cor.kIM)
kim_all$rawPval <- 10^(kim_all$log.p.cor.kIM)
kim_all$cluster <- rownames(kim_all)
kim_all$overlap <- sapply(orange_list, length)


# Plot this
library(scales)
temp <- levels(factor(kim_all$cluster))
temp <- gsub("cluster_", "", temp)
temp <- as.numeric(temp)
order <- order(temp)
ordered <- paste("cluster_", temp[order], sep = "")
kim_all$cluster <- factor(kim_all$cluster, levels = ordered)
ggplot(kim_all, aes(x=rawPval, y= cor.kIM, color=cluster)) + 
	geom_point(aes(size=overlap)) + 
	theme_bw() + 
	xlim(c(0,1)) + 
	ylim(c(0,1)) + 
	xlab("cor.kIM pvalue") + 
	ylab("cor.kIM") + 
	geom_vline(xintercept = 0.05, linetype="dotted", color = "red", size=1) + 
	scale_color_manual(values=hue_pal()(12)) + 
	theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size = 20),  axis.title.y = element_text(size = 24), axis.title.x = element_text(size = 24), plot.title=element_text(face="bold", size = 20), legend.title=element_text(size=15), legend.text=element_text(size=12)) + 
	ggtitle("Preservation of cluster 11 black trans-eQTL \ntarget correlations across clusters")



#### Calculate the MM of genes within the orange module for each cluster (including cluster 11)
# This is to compare the absolute module membership with the relative preservation of this for the paper
# Extract the genes we want - Use the expression from cluster 11 as a reference
datExpr_orange <- datExpr11[,colnames(datExpr11) %in% trans]
MEList = moduleEigengenes(datExpr_orange, colors = rep("orange", ncol(datExpr_orange)))
MEs = MEList$eigengenes
geneModuleMembership = as.data.frame(cor(datExpr_orange, MEs, use = "p"));
nSamples <- nrow(datExpr_orange)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = "MM.orange";
names(MMPvalue) = "p.MM.orange"
# Calculate the mean, significant connectivity of genes in the module
MM_tog <- cbind(geneModuleMembership, MMPvalue)
ave_sig <- mean(MM_tog[MM_tog$p.MM.orange < 0.05,]$MM.orange)
ave_all <- mean(MM_tog$MM.orange)


# Now do this for these genes in every other cluster
# Load the expression for every cluster
cluster_index <- seq(1:12)
cluster_index <- cluster_index -1 

mmres <- as.data.frame(matrix(nrow=length(cluster_index), ncol = 3))
colnames(mmres) <- c("average.MM.orange.all", "average.MM.orange.sig", "size")
rownames(mmres) <- c(clusters, "cluster_11")
for(c in seq_along(cluster_index)){
	cnum <- as.numeric(cluster_index)[c]
	pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021/Epithelial/WGCNA_within_clusters", sep = "/")
	pathOut <- paste(pathOut, "/cluster_", cnum, sep = "")
	load(paste(pathOut, "objects/UpTo_geneTraitSignificance.Rds", sep = "/"))
	# Extract expression data and transpose
	datExprtemp <- exprdata
	datExprtemp <- as.data.frame(t(datExprtemp))
	if(length(intersect(colnames(datExprtemp), trans)) > 2){
		datExpr_orange <- datExpr11[,colnames(datExprtemp) %in% trans]
		MEList = moduleEigengenes(datExpr_orange, colors = rep("orange", ncol(datExpr_orange)))
		MEs = MEList$eigengenes
		geneModuleMembership = as.data.frame(cor(datExpr_orange, MEs, use = "p"));
		nSamples <- nrow(datExpr_orange)
		MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
		names(geneModuleMembership) = "MM.orange";
		names(MMPvalue) = "p.MM.orange"
		# Calculate the mean, significant connectivity of genes in the module
		MM_tog <- cbind(geneModuleMembership, MMPvalue)
		ave_sig <- mean(MM_tog[MM_tog$p.MM.orange < 0.05,]$MM.orange)
		ave_all <- mean(MM_tog$MM.orange)
		mmres[which(rownames(mmres) == paste("cluster_", cnum, sep = "")),1] <- ave_all
		mmres[which(rownames(mmres) == paste("cluster_", cnum, sep = "")),2] <- ave_sig
	}
}
mmres$size <- c(sapply(orange_list, length), length(trans))
mmres$cluster <- rownames(mmres)


# Plot the absolute scores of this with respect to size, color like the scRNAseq
temp <- levels(factor(mmres$cluster))
temp <- gsub("cluster_", "", temp)
temp <- as.numeric(temp)
order <- order(temp)
ordered <- paste("cluster_", temp[order], sep = "")
mmres$cluster <- factor(mmres$cluster, levels = ordered)
ggplot(mmres, aes(x=size, y=average.MM.orange.all, color=cluster)) + 
	geom_point(size=5) + 
	theme_bw() + 
	xlim(0,12) +
	ylim(0,1) + 
	xlab("Overlapping trans-eQTL targets") + 
	ylab("Average MM") + 
	scale_color_manual(values=hue_pal()(12)) + 
	theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size = 20),  axis.title.y = element_text(size = 24), axis.title.x = element_text(size = 24), plot.title=element_text(face="bold", size = 20), legend.title=element_text(size=15), legend.text=element_text(size=12)) + 
	ggtitle("Average MM of trans-eQTL targets in \ncluster 11 black across clusters")


# Merge the absolute preservation and correlation with cluster 11 plots
mmres$cor.kme <- c(kme_all$cor.kME, 1)
ggplot(mmres, aes(x=cor.kme, y= average.MM.orange.all, color=cluster)) + 
	geom_point(aes(size=size)) + 
	theme_bw() + 
	xlim(c(0,1)) + 
	ylim(c(0,1)) + 
	xlab("cor.kME with cluster 11") + 
	ylab("average.ME") + 
	scale_color_manual(values=hue_pal()(12)) + 
	theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size = 20),  axis.title.y = element_text(size = 24), axis.title.x = element_text(size = 24), plot.title=element_text(face="bold", size = 20), legend.title=element_text(size=15), legend.text=element_text(size=12)) + 
	ggtitle("Preservation of cluster 11 black trans-eQTL \ntarget correlations across clusters")











#########################


#### Extract the cor.cor
cor_list <- vector("list", length = length(clusters))
names(cor_list) <- clusters
test <- as.data.frame(mp_res[[1]]$preservation$observed[[1]][[2]])
indx_cor <- which(colnames(test) == "cor.cor")
test2 <- as.data.frame(mp_res[[1]]$preservation$log.pBonf[[1]][[2]])
indx_p <-  which(colnames(test2) == "log.p.Bonf.cor.cor")
modules <- rownames(test)
for(c in seq_along(clusters)){
	cor_list[[c]] <- data.frame(cor.cor.all = mp_res[[c]]$preservation$observed[[1]][[2]][,indx_cor],
									log.p.Bonf.cor.cor = mp_res[[c]]$preservation$log.pBonf[[1]][[2]][,indx_p])
	rownames(cor_list[[c]]) <- modules
	cor_list[[c]] <- cor_list[[c]][grep("orange", rownames(cor_list[[c]])),]
}
cor_all <- do.call(rbind, cor_list)
cor_all$adjPval <- cor_all$log.p.Bonf.cor.cor*-1
plot(cor_all$adjPval, cor_all$cor.cor.all)




