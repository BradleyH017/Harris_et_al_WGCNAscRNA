# Bradley Jan 2022
# Analysis of the variability of gene expression within clusters at the single cell level
# Want to do this for each gene across clusters within datsets and for each gene within interesting clusters, across datasets

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library('Seurat')
library(ggsci)
library(rstatix)
library(optparse)
library(patchwork)
library(gage)
library(fgsea)
library(data.table)
library(limma)
library(qusage)
library(farver)
library(labeling)


# ~~~~ 1. Smillie data
# Choose set and define path - Editting to use te unmerged clusters for epithelial (as not clusters were merged), but merged for Immune
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, sep = "/")
genes.to.test <- c("C11orf53", "COLCA1", "COLCA2", "LRMP", "SH2D6", "PSTPIP2", "HTR3C", "ALOX5", "OGDHL", "MATK", "PIK3CG", "PLCG2", "PTGS1", "IL17RB", "GNG13", "CAMP", "ANKHD1", "GIN1", "SPAG6", "SH2D7", "BMX", "HTR3E", "AZGP1", "TRPM5")

# Load the variability results - This is from interesting cluster analysis script
load(paste(pathOut, "objects/Processed_seur.integrated.clusteredRds", sep = "/"))
genes.to.test <- intersect(genes.to.test, rownames(seur.integrated@assays$RNA))
var_results_smil <- vector("list", length = length(genes.to.test))
names(var_results_smil) <- genes.to.test
for(gene in seq_along(genes.to.test)){
	tryCatch({
		var_results_smil[[gene]] <- read.csv(paste(pathOut, "/tables/int_gene_variability_by_cluster/", genes.to.test[gene], "_variability_analysis_seur_single_cells.csv", sep = ""), row.names=1)
		} , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# Put results together 
var_results_smil_all <- do.call(rbind, var_results_smil)
var_results_smil_all$gene <- unlist(strsplit(rownames(var_results_smil_all), "\\."))[c(T,F,F)]
rownames(var_results_smil_all) <- NULL

# A.  Make into a matrix of results - NOTE: THis is using the raw variance - not standardized
var_mat_smil <- matrix(nrow=length(levels(factor(var_results_smil_all$cluster))), ncol = length(genes.to.test))
rownames(var_mat_smil) <- paste("cluster_", levels(factor(var_results_smil_all$cluster)))
colnames(var_mat_smil) <- genes.to.test
for(r in 1:nrow(var_mat_smil)){
	for(c in 1:ncol(var_mat_smil)){
		var_mat_smil[r,c] <- var_results_smil_all[var_results_smil_all$cluster == levels(factor(var_results_smil_all$cluster))[r], ]$variance[c]
	}
}

# A.1 plot this
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(var_mat_smil, cluster_rows=F, cluster_columns=F, 
#	col=col_fun, 
	name = "variance", column_title = "rs3087967 trans-eQTLs targets across Smillie et al., clusters", column_names_rot = 45, column_title_gp = gpar(fontsize = 15, fontface = "bold"))
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))

# A.2 For visualisation, try z-scoring this within genes
var_mat_smil_z <- apply(var_mat_smil, 2, FUN=scale, center = T, scale = T)
rownames(var_mat_smil_z) <- rownames(var_mat_smil)
col_fun = colorRamp2(c(-(max(var_mat_smil_z)), 0, max(-(max(var_mat_smil_z)))), c("blue", "white", "red"))
ht <- Heatmap(var_mat_smil_z, cluster_rows=F, cluster_columns=F, 
#	col=col_fun, 
	name = "variance", column_title = "rs3087967 trans-eQTLs targets across Smillie et al., clusters", column_names_rot = 45, column_title_gp = gpar(fontsize = 15, fontface = "bold"))
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))


# B  Make into a matrix of results - NOTE: THis is using the standardized variance
var_mat_smil_stand <- matrix(nrow=length(levels(factor(var_results_smil_all$cluster))), ncol = length(genes.to.test))
rownames(var_mat_smil_stand) <- paste("cluster_", levels(factor(var_results_smil_all$cluster)))
colnames(var_mat_smil_stand) <- genes.to.test
for(r in 1:nrow(var_mat_smil_stand)){
	for(c in 1:ncol(var_mat_smil_stand)){
		var_mat_smil_stand[r,c] <- var_results_smil_all[var_results_smil_all$cluster == levels(factor(var_results_smil_all$cluster))[r], ]$variance.standardized[c]
	}
}

# B.1 plot this
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(var_mat_smil_stand, cluster_rows=F, cluster_columns=F, 
#	col=col_fun, 
	name = "variance.standardised", column_title = "rs3087967 trans-eQTLs targets across Smillie et al., clusters", column_names_rot = 45, column_title_gp = gpar(fontsize = 15, fontface = "bold"))
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))

# B.2 For visualisation, try z-scoring this within genes
var_mat_smil_stand_z <- apply(var_mat_smil_stand, 2, FUN=scale, center = T, scale = T)
rownames(var_mat_smil_stand_z) <- rownames(var_mat_smil_stand)
col_fun = colorRamp2(c(-(max(var_mat_smil_stand_z)), 0, max(-(max(var_mat_smil_stand_z)))), c("blue", "white", "red"))
ht <- Heatmap(var_mat_smil_stand_z, cluster_rows=F, cluster_columns=F, 
#	col=col_fun, 
	name = "variance.standardised", column_title = "rs3087967 trans-eQTLs targets across Smillie et al., clusters", column_names_rot = 45, column_title_gp = gpar(fontsize = 15, fontface = "bold"))
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))



# ~~~~ 2. Elmentaite data
var_results_elm <- vector("list", length = length(genes.to.test))
names(var_results_elm) <- genes.to.test
# Some genes missing: 
chuck <- c("HTR3C", "GNG13", "CAMP", "HTR3E")
genes.to.test <- genes.to.test[-grep(paste(chuck, collapse="|"), genes.to.test)]
for(gene in seq_along(genes.to.test)){
	tryCatch({
		var_results_elm[[gene]] <- read.csv(paste("Elmentaite/results/tables/",  genes.to.test[gene], "_variability_analysis_seur_single_cells.csv", sep = ""), row.names=1)
		} , error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# A.  Put results together - NOTE: This is not standardised
var_results_elm_all <- do.call(rbind, var_results_elm)
var_mat_elm <- matrix(nrow=length(levels(factor(var_results_elm_all$cluster))), ncol = length(genes.to.test))
rownames(var_mat_elm) <- paste("cluster_", levels(factor(var_results_elm_all$cluster)))
colnames(var_mat_elm) <- genes.to.test
for(r in 1:nrow(var_mat_elm)){
	for(c in 1:ncol(var_mat_elm)){
		var_mat_elm[r,c] <- var_results_elm_all[var_results_elm_all$cluster == levels(factor(var_results_elm_all$cluster))[r], ]$variance[c]
	}
}

# A. 1 Plot this
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(var_mat_elm, cluster_rows=F, cluster_columns=F, 
#	col=col_fun, 
	name = "variance", column_title = "rs3087967 trans-eQTLs targets across Elmentaite., clusters", column_names_rot = 45, column_title_gp = gpar(fontsize = 15, fontface = "bold"))
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))

# A.2 For visualisation, try z-scoring this within genes
var_mat_elm_z <- apply(var_mat_elm, 2, FUN=scale, center = T, scale = T)
rownames(var_mat_elm_z) <- rownames(var_mat_elm)
col_fun = colorRamp2(c(-(max(var_mat_elm_z)), 0, max(-(max(var_mat_elm_z)))), c("blue", "white", "red"))
ht <- Heatmap(var_mat_elm_z, cluster_rows=F, cluster_columns=F, 
#	col=col_fun, 
	name = "variance", column_title = "rs3087967 trans-eQTLs targets across Elmentaite et al., clusters", column_names_rot = 45, column_title_gp = gpar(fontsize = 15, fontface = "bold"))
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))

# B - USing standardized variance
var_mat_elm_stand <- matrix(nrow=length(levels(factor(var_results_elm_all$cluster))), ncol = length(genes.to.test))
rownames(var_mat_elm_stand) <- paste("cluster_", levels(factor(var_results_elm_all$cluster)))
colnames(var_mat_elm_stand) <- genes.to.test
for(r in 1:nrow(var_mat_elm)){
	for(c in 1:ncol(var_mat_elm)){
		var_mat_elm_stand[r,c] <- var_results_elm_all[var_results_elm_all$cluster == levels(factor(var_results_elm_all$cluster))[r], ]$variance.standardized[c]
	}
}

# B.1 plot this
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(var_mat_elm_stand, cluster_rows=F, cluster_columns=F, 
#	col=col_fun, 
	name = "variance.standardized", column_title = "rs3087967 trans-eQTLs targets across Elmentaite., clusters", column_names_rot = 45, column_title_gp = gpar(fontsize = 15, fontface = "bold"))
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))

# A.2 For visualisation, try z-scoring this within genes
var_mat_elm_stand_z <- apply(var_mat_elm_stand, 2, FUN=scale, center = T, scale = T)
rownames(var_mat_elm_stand_z) <- rownames(var_mat_elm_stand)
test <- var_mat_elm_stand_z
test[is.na(test)] <- 0
col_fun = colorRamp2(c(-(max(test)), 0, max(-(max(test)))), c("blue", "white", "red"))
ht <- Heatmap(var_mat_elm_stand_z, cluster_rows=F, cluster_columns=F, 
#	col=col_fun, 
	name = "variance.standardized", column_title = "rs3087967 trans-eQTLs targets across Elmentaite et al., clusters", column_names_rot = 45, column_title_gp = gpar(fontsize = 15, fontface = "bold"))
draw(ht, padding = unit(c(2, 15, 2, 2), "mm"))



# ~~~~~~ Elmentaite vs Smillie data
# looking at non-standardized variance first
# From Smillie want cluster 11
# From Elmentaite want cluster 18

# A. Using raw variation
int_vars_list <- vector("list", length=2)
names(int_vars_list) <- c("Smillie", "Elmentaite")
int_vars_list[[1]] <- var_mat_smil[grep("11", rownames(var_mat_smil)),]
int_vars_list[[2]] <- var_mat_elm[grep("18", rownames(var_mat_elm)),]

# Function to combine
rbind.named.fill <- function(x) {
    nam <- sapply(x, names)
    unam <- unique(unlist(nam))
    len <- sapply(x, length)
    out <- vector("list", length(len))
    for (i in seq_along(len)) {
        out[[i]] <- unname(x[[i]])[match(unam, nam[[i]])]
    }
    setNames(as.data.frame(do.call(rbind, out), stringsAsFactors=FALSE), unam)
}

# Combine
vars_comp <- rbind.named.fill(int_vars_list)
rownames(vars_comp) <- names(int_vars_list)
vars_comp <- as.data.frame(t(vars_comp))

# Calculate fold change
vars_comp$FC <- vars_comp$Smillie / vars_comp$Elmentaite

# Save this
write.csv(vars_comp, "BH_analysis/Seurat/August_2021/Epithelial/tables/Smillie_vs_Elmentaite_raw_variantion_from_seur_single_cells.csv")

# B. Using variation.standardized
int_vars_list_stand <- vector("list", length=2)
names(int_vars_list_stand) <- c("Smillie", "Elmentaite")
int_vars_list_stand[[1]] <- var_mat_smil_stand[grep("11", rownames(var_mat_smil_stand)),]
int_vars_list_stand[[2]] <- var_mat_elm_stand[grep("18", rownames(var_mat_elm_stand)),]

# Combine
vars_comp_stand <- rbind.named.fill(int_vars_list_stand)
rownames(vars_comp_stand) <- names(int_vars_list_stand)
vars_comp_stand <- as.data.frame(t(vars_comp_stand))

# Calculate fold change
vars_comp_stand$FC <- vars_comp_stand$Smillie / vars_comp_stand$Elmentaite

# Save this
write.csv(vars_comp_stand, "BH_analysis/Seurat/August_2021/Epithelial/tables/Smillie_vs_Elmentaite_variantion_standardized_from_seur_single_cells.csv")


