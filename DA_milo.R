# Bradley October 2021
# Differential abundance testing of Smillie clusters (mine and theirs), divided by the turquoise hub cluster grouping defined by WGCNS in pseudo-bulked sc-RNASeq from all samples
# DA_milo.r
# Following the tutorial https://rdrr.io/github/MarioniLab/miloR/f/vignettes/milo_demo.Rmd and https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html#1_Load_data

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library(miloR)
set <- "Epithelial"
pathOut = paste("BH_analysis/Seurat/August_2021", set, sep = "/")

# read in the final object
load(paste(pathOut, "objects/Processed_seur.integrated.clusteredRds", sep = "/"))

# Append the site, batch and turquoise hub cluster grouping to the seurat object metadata before conversion into a single cell experiment object
meta <- data.frame(batch = rep("", length(levels(factor(seur.integrated@meta.data$Sample)))),
                Sex = c("F","F", "M","M", "F","F",  "M","M",  "M","M",  "M","M",  "F","F",  "F","F",  "F","F",  "M","M",  "F","F"),
                Site = c("R","R", "R","R", "R","R", "R","R", "R","R", "R","R", "R","R",  "","", "L","L", "R","R", "",""),
                Sample = levels(factor(seur.integrated@meta.data$Sample)))
        meta$Subject <- unlist(strsplit(meta$Sample, "\\."))[c(T,F)]
        for(x in seq(1:nrow(meta))){
               if(as.numeric(gsub("N", "", meta$Subject[x])) < 40) {
                       meta$batch[x] <- "train"
               } else {
                       meta$batch[x] <- "valid"
               }
        }
        rownames(meta) <- NULL
        rownames(meta) <- meta$Sample;
        # Binarise / make variables in the meta continous
#        meta <- meta %>% mutate(batch=recode(batch, "train" = 0, "valid" = 1));
#        meta <- meta %>% mutate(Sex=recode(Sex, "F" = 0, "M" = 1))
        meta$Site.R <- rep(0, nrow(meta))
        for(r in seq(1:nrow(meta))){
                if(meta$Site[r] == "R"){
                        meta$Site.R[r] <- 1;
                }
                if(meta$Site[r] == "T"){
                        meta$Site.R[r] <- 1
                }
                if(meta$Site[r] == ""){
                        meta$Site.R[r] <- NA
                }
                if(meta$Site[r] == "L"){
                        meta$Site.R[r] <- 0
                }
        }
        meta <- meta[,-which(colnames(meta) == "Site")]
        meta <- meta[,-which(colnames(meta) == "Subject")]
cellmeta <- seur.integrated@meta.data
cellmeta$Sex <- ifelse(cellmeta$Sample %in% meta[meta$Sex == "M", ]$Sample, "M", "F")
cellmeta$batch <- ifelse(cellmeta$Sample %in% meta[meta$batch == "train",]$Sample, "train", "valid")
cellmeta$Site <- ifelse(cellmeta$Sample %in% meta[meta$Site == "1",]$Sample, "R", "L")
Low <- c("N20.EpiA", "N20.EpiB", "N8.EpiB", "N13.EpiB", "N13.EpiA")
cellmeta$Blue_hub_gene_group <- ifelse(cellmeta$Sample %in% Low, "Low", "High")


# Convert this into a single cell experiment object 
# Originally using log2TP10K
DefaultAssay(seur.integrated) <- "integrated"
sce <- as.SingleCellExperiment(seur.integrated)

# Convert into a milo object, using the nearest neighbour graph I have already computed (kNN=250, d = nPCs = 50)
seur.milo <- Milo(sce)


# For current best results, use_seurat_graph=F
use_seurat_graph=F

if(use_seurat_graph == T){
	pathOut <- paste(pathOut, "DA/milo/seurat_kNN", sep = "/")
	# using the nearest neighbour graph I have already computed (kNN=250, d = nPCs = 50)
	k=250
	d=50
	milo <- buildFromAdjacency(seur.integrated@graphs$integrated_snn, k = k, d=d, is.binary = F)
	graph(seur.milo) <- graph(milo)
	# Define the representative nearest neighbours graph
	# This is done to reduce computation, and not testing every neighbourhood but those with representative cells
	# As well as d and k, for sampling we need to define a few additional parameters:
	#prop: the proportion of cells to randomly sample to start with. We suggest using prop=0.1 for datasets of less than 30k cells. For bigger datasets using prop=0.05 should be sufficient (and makes computation faster). - Using prop=0.1
	#refined: indicates whether you want to use the sampling refinement algorith, or just pick cells at random. The default and recommended way to go is to use refinement. The only situation in which you might consider using random instead, is if you have batch corrected your data with a graph based correction algorithm, such as BBKNN, but the results of DA testing will be suboptimal. 
	seur.milo <- makeNhoods(seur.milo, prop = 0.1, k=k, d=d, refined = TRUE, reduced_dims = "PCA")
	pdf(file = paste(pathOut, "plots/Nhood_histogram.pdf", sep = "/"))
	plotNhoodSizeHist(seur.milo)
	dev.off()
	# Most neighbourhoods are extremely large (~1000 cells), which greatly exceeds the ~60 used in the tutorials
	}


if(use_seurat_graph == F){
	pathOut <- paste(pathOut, "DA/milo/reconstruct_kNN", sep = "/")
	# See how many PCs to consider
	library(SingleCellExperiment)
	library(scater)
	library(scran)
	k=250
	sce <- runPCA(sce, ncomponents=50)
	percent.var <- attr(reducedDim(sce), "percentVar")
	plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
	chosen.elbow <- PCAtools::findElbowPoint(percent.var)
	chosen.elbow
	# This is shown to equal 4
	# Using the same k as have done for clustering
	d = as.numeric(chosen.elbow)
	set.seed(123)
	seur.milo <- buildGraph(seur.milo, k=k, d = d, reduced.dim = "PCA")
	# Now make the neighbourhoods
	seur.milo <- makeNhoods(seur.milo, prop = 0.1, k=k, d=d, refined = TRUE, reduced_dims = "PCA")
	pdf(file = paste(pathOut, "plots/Nhood_histogram.pdf", sep = "/"))
	plotNhoodSizeHist(seur.milo) 
	dev.off()
	# This brings the average Nhood size to around 400, still large though
}

# Count the cells within each neighbourhood
seur.milo <- countCells(seur.milo, meta.data = as.data.frame(colData(seur.milo)), sample="Sample")
head(nhoodCounts(seur.milo))

# Defining experimental design
# Uses negative Binomial GLM in edgeR
design <- cellmeta[,c("Sample", "batch", "Sex", "Site", "Blue_hub_gene_group")]
# Convert batch info from integer to factor
design$batch <- as.factor(design$batch)
design$Sex <- as.factor(design$Sex)
design$Site <- as.factor(design$Site)
design <- distinct(design)
rownames(design) <- design$Sample
design_model <- as.data.frame(model.matrix(~batch + Sex + Site + Blue_hub_gene_group,data=design))

# Caluclate the nighbourhood distances
#Milo uses an adaptation of the Spatial FDR correction introduced by cydar, where we correct p-values accounting for the amount of overlap between neighbourhoods. Specifically, each hypothesis test P-value is weighted by the reciprocal of the kth nearest neighbour distance. To use this statistic we first need to store the distances between nearest neighbors in the Milo object. This is done by the calcNhoodDistance function (N.B. this step is the most time consuming of the analysis workflow and might take a couple of minutes for large datasets).
seur.milo <- calcNhoodDistance(seur.milo, d=d, reduced.dim = "PCA")
# Now do the DA testing
da_results <- testNhoods(seur.milo, design = ~ Blue_hub_gene_groupLow, design.df = design_model)
head(da_results)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 
nrow(da_results[da_results$SpatialFDR < 0.05,])

# Inspecting the DA results
# Histogram of nominal pvals - Looks okay
pdf(file=paste(pathOut, "plots/NomP_histogram.pdf", sep = "/"))
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
dev.off()

# Volcano
pdf(file=paste(pathOut, "plots/Volcano_DE_nhoods.pdf", sep = "/"))
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
dev.off()

# Plot on the Umap
seur.milo <- buildNhoodGraph(seur.milo)

## Plot single-cell UMAP
seur.milo@metadata <- cellmeta
group.add <- cellmeta$Blue_hub_gene_group
names(group.add) <- rownames(cellmeta)
seur.integrated <- AddMetaData(seur.integrated, group.add, col = "Blue_hub_gene_group")
umap_pl <- DimPlot(seur.integrated, group.by = "Blue_hub_gene_group", label=F)
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(seur.milo, da_results, layout="UMAP",alpha=0.1) + ggtitle("Differential abundance of neighbourhoods") + theme(plot.title=element_text(face="bold"))

pdf(file = paste(pathOut, "plots/umap_and_nhoods.pdf", sep = "/"))
umap_pl + nh_graph_pl
dev.off()

# Annotate the neighbourhoods in the graph
da_results <- annotateNhoods(seur.milo, da_results, coldata_col = "seurat_clusters")
head(da_results)

write.csv(paste(pathOut, "tables/da_results.csv", sep = "/"))

# Exclude neighbourhoods that are a mix of fractions - Shows that lots are a mix
pdf(file = paste(pathOut, "plots/Nhood_fractions_per_seurat_cluster.pdf", sep = "/"))
ggplot(da_results, aes(seurat_clusters_fraction)) + geom_histogram(bins=50) + geom_vline(xintercept = 0.8, linetype="dotted", 
                color = "red", size=1.5)
dev.off()

da_results$celltype <- ifelse(da_results$seurat_clusters_fraction < 0.8, "Mixed", da_results$seurat_clusters)

# Now look at DE, generalising to cluster annotations
old <- da_results
da_results$seurat_clusters <- factor(da_results$seurat_clusters, levels=as.character(11:0))

pdf(file = paste(pathOut, "plots/DAbeeswarm_per_cluster.pdf", sep = "/"))
plotDAbeeswarm(da_results, group.by = "seurat_clusters", alpha=0.01) + ggtitle("Neighbourhood differential abundance by cluster") + xlab("Cluster") + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size = 20),  axis.title.y = element_text(size = 24), axis.title.x = element_text(size = 24), plot.title=element_text(face="bold", size = 28)) 
dev.off()

# Test without mixed neighbourhoods
nomix <- da_results[da_results$celltype != "Mixed",]
pdf(file = paste(pathOut, "plots/DAbeeswarm_per_cluster_no_mixed.pdf", sep = "/"))
plotDAbeeswarm(nomix, group.by = "seurat_clusters", alpha=0.01) + ggtitle("Neighbourhood differential abundance by cluster") + xlab("Cluster") + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size = 20),  axis.title.y = element_text(size = 24), axis.title.x = element_text(size = 24), plot.title=element_text(face="bold", size = 28)) 
dev.off()


save.image(paste(pathOut, "After_beeswarm.Rds", sep = "/"))


# Count proportion of DA neighbourhoods for each cluster, exluding the mixed
clusters <- levels(factor(nomix$seurat_clusters))
prop_list <- vector("list", length = length(clusters))
names(prop_list) <- clusters
for(cluster in seq_along(clusters)){
	prop_list[[cluster]] <- nomix[nomix$seurat_clusters == clusters[cluster],]
	prop_list[[cluster]] <- prop_list[[cluster]][prop_list[[cluster]]$celltype != "Mixed" & prop_list[[cluster]]$SpatialFDR < 0.01,]
}
sig <- sapply(prop_list, nrow)
tot <- table(nomix$seurat_clusters)
frac <- rep(0, length(sig))
names(frac) <- names(sig)
for(c in seq_along(frac)){
	frac[c] <- (sig[which(names(sig) == names(frac)[c])] / tot[which(names(tot) == names(frac)[c])] * 100)
}

######### Plot the abundance of clusters across the blue module grouping 
freq <- data.frame(sample = factor(seur.integrated@meta.data$Sample), seurat_clusters = factor(seur.integrated@meta.data$seurat_clusters))
library(plyr)
tally <- ddply(freq, .(freq$sample, freq$seurat_cluster), nrow, .drop=F)
colnames(tally)<- c("orig.ident", "seurat_cluster", "freq")

#~~~ Reformat
samps <- levels(factor(tally$orig.ident))
indx <- seq(1:length(samps))
per_samp <- vector("list", length = length(samps))
names(per_samp) <- samps
for(x in indx){per_samp[[x]] = tally[tally$orig.ident == samps[x],c(2:3)];
rownames(per_samp[[x]]) <- per_samp[[x]][,1];
per_samp[[x]] <- per_samp[[x]][,-1];
colnames(per_samp[[x]]) <= levels(factor(tally$seurat_cluster))}
freq <- data.frame(matrix(unlist(per_samp), nrow = length(per_samp), byrow = T))
rownames(freq) <- names(per_samp)
colnames(freq) <- levels(factor(tally$seurat_cluster))

# Now convert to proportions
for(r in seq(1:nrow(freq))){
        freq[r,] <- (freq[r,]/rowSums(freq[r,]))*100
}

meta_ab <- cellmeta[,c("Sample", "batch", "Sex", "Site", "Blue_hub_gene_group")]
meta_ab <- distinct(meta_ab)
freq <- freq[match(meta_ab$Sample, rownames(freq)), ]
#colnames(freq) <- paste("Cluster_", colnames(freq), sep = "")
freq$Blue_hub_gene_group <- meta_ab$Blue_hub_gene_group

library(rstatix)
library(ggpubr)

nclust <- ncol(freq) -1

stat.test <- vector("list", length = nclust)
bxp <- vector("list", length = nclust)

for(clust in 1:nclust){
        clust_name <- colnames(freq)[clust]
        temp <- freq[, colnames(freq) %in% c(clust_name, "Blue_hub_gene_group")]
        colnames(temp) <- c("Cluster", "Blue_hub_gene_group")
        stat.test[[clust]] <- temp %>% wilcox_test(Cluster ~ Blue_hub_gene_group, paired = F)
        stat.test[[clust]]$y.position = max(temp$Cluster*1.1)
        stat.test[[clust]]$Cluster <- clust_name
        rm(temp)
}

stat.test <- do.call(rbind, stat.test)
stat.test$padj <- p.adjust(stat.test$p, method = "BH") %>% signif(digits = 2)

# Mine
melt <- melt(freq)
melt$Blue_hub_gene_group <- factor(melt$Blue_hub_gene_group)

ggplot(data = melt, aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill = Blue_hub_gene_group), width = 0.8) + theme_bw() + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size = 20),  axis.title.y = element_text(size = 24), axis.title.x = element_text(size = 24), plot.title = element_text(face = "bold", size = 28), legend.position="bottom") + xlab("Cluster") + ylab("Proportion (%)") + ggtitle("Cluster proportion across grouping") + scale_fill_manual(values=c("#00A087FF", "#DC0000FF"))



