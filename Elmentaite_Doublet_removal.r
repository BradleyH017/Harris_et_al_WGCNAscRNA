# Bradley October 2021
# scRNASeq preprocessing of Elmentaite (Epithelial) cells, aligned using 10X 
# Doublet_detection.r

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
setwd("Elmentaite/results")
library(DoubletFinder)


# Pre-process the Seurat Object using SCTransform, but make sure this is per-sample - don't do all samples together.
# So will split into each sample
# Also need to ensure data is cleared of low-quality clusters

load("Preprocessing/raw_seur_GeneCellQC.RData")
seur <- seur_final
rm(seur_final)


# Split Object
seur_list <- SplitObject(seur, split.by = "orig.ident")

# Apply SCTransofmr processing using usual parameters
for( x in seq(1:length(seur_list))){
        seur_list[[x]] <- SCTransform(seur_list[[x]]);
        seur_list[[x]] <- RunPCA(seur_list[[x]]);
        seur_list[[x]] <- RunUMAP(seur_list[[x]], reduction = "pca", dims = 1:30)
}


# Save this so don't have to re-compute
pathOut = "Doublet_detection"
save(seur_list, file = "seur_list_SCTransformed_per_sample.Rds")


# As have no ground-truth, will make use of pK identification (no groung-truth). Using 20 PCS based on elbow plot ------------------------------------



########### Perform and summarize parameter sweep 

# test paramters over a range to optomize
sweep.res.list <- vector("list", length = length(seur_list))
for( x in  seq(1:length(seur_list))){
        sweep.res.list[[x]] <- paramSweep_v3(seu = seur_list[[x]], PCs = 1:20, sct=TRUE)
}

save(sweep.res.list, file = "Doublet_detection/seur_sweep.res.list.Rds")


# Summarize
sweep.stats <- vector("list", length = length(seur_list))
for( x in seq(1:length(seur_list))){
        sweep.stats[[x]] <- summarizeSweep(sweep.res.list[[x]], GT = F)
}

bcmvn.list <- vector("list", length = length(seur_list))
for(x in seq(1:length(seur_list))){
        bcmvn.list[[x]] <- find.pK(sweep.stats[[x]])
}

###### Plotting this: 
# Generate sample variable to facet by
for(x in seq(1:length(bcmvn.list))){
        bcmvn.list[[x]]$Sample <- rep(levels(factor(seur_list[[x]]@meta.data$orig.ident)), nrow(bcmvn.list[[x]]))
}

bcmvn.plot.data <- do.call(rbind, bcmvn.list)



# Generate plot
y <- ggplot(data =  bcmvn.plot.data, aes(x = pK, y = BCmetric, group = 1)) +
        geom_line(color = "#41b6c4") +
        geom_point() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45)) +
        facet_wrap(~ Sample, scales = "free", ncol = 4)
        
# Caluclate pK at max BCmetric for each
bcmvn.plot.data$pK <- as.numeric(bcmvn.plot.data$pK)
max_vals <- vector("list", length = length(bcmvn.list))
for( x in seq(1:length(bcmvn.list))){
        max_vals[[x]] <- bcmvn.list[[x]][bcmvn.list[[x]]$BCmetric == max(bcmvn.list[[x]]$BCmetric),]
}
max_labels <- do.call(rbind, max_vals)
max_labels$pK <- paste("pK=", max_labels$pK, sep = "")

# Add text to facet plot
y <- y + geom_text(data = max_labels,
        mapping = aes(x = -Inf, y = -Inf, label = pK),
        hjust   = -0.5,
        vjust   = -4)
        
pdf(file = "Doublet_detection/pK_vs_BCmetric_perSample.pdf")
y
dev.off()

save.image("Doublet_detection/seur_pK_pN_parameterSwept_DoubletDetection_per_sample.Rds")

########### Model homotypic doublet proportions, set pANN thresholds

# Keep the max labels per samples
pK_list <- vector("list", length = nrow(max_labels))
for(x in seq(1:nrow(max_labels))){
        pK_list[[x]] <- max_labels$pK[x];
        pK_list[[x]] <- gsub("pK=", "", pK_list[[x]])
}

# Do this per annotations (cluster) - want to do this for my clusters, so need to load in my final seurat object from DimRed after kArray
# For Elmentaite et al kArray, k.param = 20
k.param = "20"

load(paste("kArray/seur.integrated.clusterd.k", k.param, ".Rds", sep=""))
seur.final.list <- SplitObject(seur.integrated, split.by = "orig.ident")
rm(seur.integrated)
annotations <- vector("list", length = length(pK_list))
for(x in seq(1:length(annotations))){
        annotations[[x]] <- seur.final.list[[x]]@meta.data$seurat_clusters
}

homotypic.prop.list <- vector("list", length = length(annotations))
for( x in seq(1:length(homotypic.prop.list))){
        homotypic.prop.list[[x]] <- modelHomotypic(annotations[[x]])
}


# Setting pANN thresholds
# For 10X systems, claims doublet rate of 0.9% per 1000 cells
# https://www.kumc.edu/Documents/gsf/10X%20Single%20Cell%20Expression_PlanningGuide.pdf
# Increases with the number of cells loaded (8000) per sample, so estimated doublet rate is 8*0.9
# As all the cells were loaded at the same density (8000 cells), this puts doublet detection estimation at a rate of 0.072% of all cells
DR = 0.009*8
nExp.list <- vector("list", length = length(annotations))
for(x in seq(1:length(nExp.list))){
        nExp.list[[x]] <- round(DR*(ncol(seur_list[[x]])))
}

# Adjust this for homotypic proportion
nExp.list.adj <- vector("list", length = length(annotations))
for(x in seq(1:length(nExp.list.adj))){
        nExp.list.adj[[x]] <- round(nExp.list[[x]]*(1-homotypic.prop.list[[x]]))
}


#################### Run Doublet Finder with varying classification stringencies
# pN = The proportion of artifical droplet to include - this is used at 0.25 by authors so will keep at this
# pK  = The maxima of BCmvn distributions per sample
# PCs = The number of statistically significant PCs to include: This is inferred from the elbow plot diags (ElbowPlot(x)).
# nExp = The number of expected doublets per sample, will use the adjusted
# reuse.pANN = whether or not to reuse the pANN to generate list of doublet classifications from existing vector to save time - Will set at F.
# sct = Has this data been sct normalised (for me, yes) 

pK_list <- as.numeric(pK_list)
for(x in seq(1:length(seur_list))){
        seur_list[[x]] <- doubletFinder_v3(seur_list[[x]], PCs = 1:20, pN = 0.25, pK = pK_list[[x]], nExp = nExp.list.adj[[x]], reuse.pANN = F, sct = T)
}

save(seur_list, file = "Doublet_detection/seur_post_doublet_detection_per_sample_final.Rds", sep = "/")

############ For analysis itself, want to remove doublets from the original object, then repeat the clustering, dimensionality reduction, enrichment vs authors etc.
load("Preprocessing/raw_seur_GeneCellQC.RData")
seur <- seur_final
rm(seur_final)

# plot the UMAPs of these samples, labelling by their doublets
doublet_plots <- vector("list", length = length(seur_list))
for(x in seq(1:length(doublet_plots))){
       colnames(seur_list[[x]]@meta.data)[grep("DF.classifications", colnames(seur_list[[13]]@meta.data))] <- "DF.classifications";
       doublet_plots[[x]] <- DimPlot(seur_list[[x]], group.by="DF.classifications", reduction = "umap") + NoLegend();
       doublet_plots[[x]] <- doublet_plots[[x]] + labs(title = levels(factor(seur_list[[x]]@meta.data$orig.ident)))
}

library(patchwork)
pdf(file = paste(pathOut, "Doublet_detection/doublets_per_sample.pdf", sep = "/"))
wrap_plots(doublet_plots[1:8], ncol = 2)
dev.off()

## plot the relationship between logUMI and doublets
library(rstatix)
metas <- vector("list", length = length(seur_list))
for(x in seq(1:length(metas))){
       metas[[x]] <- seur_list[[x]]@meta.data;
       colnames(metas[[x]])[grep("pANN", colnames(metas[[x]]))] <- "pANN"
}


meta.all <- do.call(rbind, metas)
meta.all$lognUMI <- log10(as.numeric(meta.all$nCount_RNA))
meta.all$DF.classifications <- as.factor(meta.all$DF.classifications)
stat.test.list <- vector("list", length = length(seur_list))
for(x in seq(1:length(stat.test.list))){
       stat.test.list[[x]] <- meta.all[meta.all$orig.ident == levels(factor(meta.all$orig.ident))[x],] %>% wilcox_test(lognUMI ~ DF.classifications, paired = F);
       stat.test.list[[x]]$y.position <-  max(meta.all[meta.all$orig.ident == levels(factor(meta.all$orig.ident))[x],]$lognUMI)*1.1
}
stat.test <- do.call(rbind, stat.test.list)

bxp1 <- ggboxplot(meta.all, x = "DF.classifications", y = "lognUMI", fill = "DF.classifications", outlier.shape = NA, palette = c("#F8766D", "#00BFC4")) + theme(legend.position="none") + theme(strip.background = element_blank()) + theme_classic() + facet_wrap( ~ orig.ident, ncol = 4)


pdf(file="Doublet_detection/doublets_vs_lognUMI.pdf")
bxp1
dev.off()

rownames(stat.test) <- levels(factor(meta.all$orig.ident))

# Now make the vector to add to the seur object from the previous attempt.
doublet_meta <- meta.all$DF.classifications
names(doublet_meta) <- rownames(meta.all)


## Now plot them on the previous attempt (after clustering etc)
load(paste("kArray/seur.integrated.clusterd.k", k.param, ".Rds", sep=""))
seur.integrated <- AddMetaData(object = seur.integrated, metadata = doublet_meta, col.name = "DF.classifications")

pdf(file = paste(pathOut, "Doublet_detection/Doublets_on_UMAP_pre_removal.pdf", sep = "/"))
DimPlot(seur.integrated, group.by = "DF.classifications")
dev.off()

# Remove the doublets from the original seur object read into this script (this is the once produced after the preprocessing)
seur <- AddMetaData(object = seur, metadata = doublet_meta, col.name = "DF.classifications")
seur <- subset(seur, subset = DF.classifications == "Singlet")

print(paste("The final dimensions of the matrix are:", dim(seur)))

# Save, this can now be used as the fully processed set for dimensionality reduction, clustering (again), statistical assessment of clusters and further analysis
save(seur, file = "Doublet_detection/seur_raw_post_doublet_removal.Rds")

