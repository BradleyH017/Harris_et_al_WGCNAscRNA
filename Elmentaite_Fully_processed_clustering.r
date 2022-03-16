# Bradley October 2021
# scRNASeq preprocessing of Elmentaite (Epithelial) cells, aligned using 10X 
# Fully_processed_clustering.r

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library(ggsci)
library(rstatix)
library(optparse)
library(patchwork)
library(Rcpp)
library(rlang)
library(dplyr)
library(gage)
library(fgsea)
library(data.table)
library(limma)
library(qusage)
library(farver)
library(labeling)
setwd("Elmentaite/results")

# Read in the fully processed, post doublet object
load("Doublet_detection/seur_raw_post_doublet_removal.Rds")

# Perform the integration of data across samples
# ReciprocalPCA integration step
# Want to regress out MT
#https://satijalab.org/seurat/articles/integration_rpca.html
split.by = "orig.ident"
vars.to.regress = "MT"


# Doing the integration with the doublet retains - Only need to do this once, then do the kArray to see which k value best atches the authors clusters and is likley tobe the most useful
seur.list <- SplitObject(seur, split.by = split.by)
seur.list <- lapply(X = seur.list, FUN = SCTransform, vars.to.regress = vars.to.regress)
features <- SelectIntegrationFeatures(object.list = seur.list, nfeatures = 3000)
seur.list <- PrepSCTIntegration(object.list = seur.list, anchor.features = features)
seur.list <- lapply(X = seur.list, FUN = RunPCA, features = features)
# Using k.anchor = 5, PCs 1:30 as default
seur.anchors <- FindIntegrationAnchors(object.list = seur.list, normalization.method = "SCT",
anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 5)
seur.integrated <- IntegrateData(anchorset = seur.anchors, normalization.method = "SCT", dims = 1:50)
# Now run PCA
seur.integrated <- RunPCA(seur.integrated, verbose = FALSE)
# Now run UMAP
seur.integrated <- RunUMAP(seur.integrated, reduction = "pca", dims = 1:50)
seur.integrated@meta.data$nCount_RNA <- as.numeric(seur.integrated@meta.data$nCount_RNA)
seur.integrated@meta.data$nFeature_RNA <- as.numeric(seur.integrated@meta.data$nFeature_RNA)



# Before saving, also calculate the log2TP10K matrix, so that this can be used in Authors cluster marker analysis (here) and seurat cluster marker analysis (next in kArray)
## Suggested in https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAST-interoperability.html
##https://github.com/satijalab/seurat/discussions/4032
## That this is also scale normalized, could use cpm or tmm for this. Initially using cpm (as feel like TMM might not perform well on single cell data, cpm is a simple converson)
## Split into 1000 cell chunks for computation
# Load the seur object
counts_list <- vector("list", length = ceiling(ncol(seur.integrated@assays$RNA@counts)/1000))
for(x in seq(1:(length(counts_list)-1))){
       counts_list[[x]] <- as.data.frame(seur.integrated@assays$RNA@counts[,(((1000*x)-999):(1000*x))])
}
counts_list[[length(counts_list)]] <- seur.integrated@assays$RNA@counts[,((floor(ncol(seur.integrated@assays$RNA@counts)/1000)*1000)+1):ncol(seur.integrated@assays$RNA@counts)]
# Now perform the logcpm normalisation of each 1000. These are too large to do this in R, so will save and do this from the command line
for(x in seq(1:length(counts_list))){
        y <- DGEList(counts = counts_list[[x]]);
        cpm <- edgeR::cpm(y, log=F, normalized.lib.sizes=F);
        cpm <- cpm/1000; #So that is per 10k
        logcpm <- log2(cpm +1);
        counts_list[[x]] <- logcpm;
        counts_list[[x]] <- as(counts_list[[x]], "sparseMatrix");
}
log2TPM <- do.call(cbind, counts_list)

seur.integrated[['log2TP10K']] <- CreateAssayObject(counts = log2TPM)


# From kArray, want to use k=20
k.param = "20"
DefaultAssay(seur.integrated) <- "integrated"
seur.integrated <- FindNeighbors(seur.integrated, reduction = "pca", dims = 1:50, k.param = as.numeric(k.param))
# Now find clusters (using usual resolution as this seems to work well)
seur.integrated <- FindClusters(seur.integrated, resolution = 0.6)


# Now save 
save(seur.integrated, file = "objects/Processed_seur.integrated.clusteredRds")


library('MAST')
packageVersion('MAST')
DefaultAssay(seur.integrated) <- "log2TP10K"
Idents(seur.integrated) = "seurat_clusters"
# Using the same thresholds as have done before
markers <- FindAllMarkers(seur.integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
# Save
write.csv(markers, "tables/markers_unmerged_clusters/fully_processed_markers.csv")

#### Calculate the markers against the authors
gene_list_dir <- "/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021/Epithelial/tables/markers/Authors"
sets <- list.files(gene_list_dir)
sets <- paste(gene_list_dir, sets, sep = "/")
smil_sets <- lapply(sets, read.gmt)
# reformat
for(x in seq(1:length(smil_sets))){
names(smil_sets)[x] <- names(smil_sets[[x]]);
smil_sets[[x]] <- smil_sets[[x]][[1]]
}
# Subset for the significant only
markers <- markers[markers$p_val_adj < 0.05,]
seur_clusters <- levels(factor(markers$cluster))
seur_markers <- vector("list", length = length(seur_clusters))
seur_geneList <- vector("list", length = length(seur_clusters))
names(seur_geneList) <- seur_clusters
# Enrichment vs authors
for(b in seq(1:length(seur_markers))){
        seur_markers[[b]] <- markers[markers$cluster == seur_clusters[b],];
        seur_geneList[[b]] <- seur_markers[[b]]$avg_log2FC;
        names(seur_geneList[[b]]) <- seur_markers[[b]]$gene
}
# Now do the enrichment
fgsea_Res <- vector("list", length = length(seur_geneList))
names(fgsea_Res) <- names(seur_geneList)
for(q in seq(1:length(fgsea_Res))){
        print(paste("testing cluster ", names(fgsea_Res)[q], sep = ""));
        fgsea_Res[[q]] <-  fgseaMultilevel(smil_sets, seur_geneList[[q]], minSize = 0, maxSize = 500, scoreType = "pos", eps = 0)
        fgsea_Res[[q]]$cluster <- rep(names(fgsea_Res)[q], nrow(fgsea_Res[[q]]));
        fgsea_Res[[q]] <- fgsea_Res[[q]][order(fgsea_Res[[q]]$padj),];
        fgsea_Res[[q]]$BH_p.adj <- p.adjust(fgsea_Res[[q]]$pval, method = "BH")
}
fgsea_Res <- do.call(rbind, fgsea_Res)
# Save this
fgsea_Res <- data.frame(fgsea_Res)
fgsea_Res <- fgsea_Res[,-which(colnames(fgsea_Res) == "leadingEdge")]
write.csv(fgsea_Res, "tables/markers_unmerged_clusters/fgsea_against_authors.csv")

# Finally, make a table of the comparisons for DE testing between clusters 
clusters <- levels(factor(seur.integrated@meta.data$seurat_clusters))
comparisons <- expand.grid(clusters, clusters)
comparisons <- comparisons[comparisons$Var1 != comparisons$Var2,]
comparisons$comb <- ""
for(r in seq(1:nrow(comparisons))){
       if(as.numeric(comparisons$Var1)[r] < as.numeric(comparisons$Var2)[r]){
               comparisons[r, 3] = paste(comparisons$Var1[r], comparisons$Var2[r], sep = ".")
       } else {
               comparisons[r, 3] = paste(comparisons$Var2[r], comparisons$Var1[r],sep = ".")
       }
}

comp <- as.data.frame(comparisons)
comp <- comp[!duplicated(comp$comb),-3]

# To submi this as an array job, want to do this once per test. So save the list as seperate files for each condition and pass it as an argument
write.csv(comp[,1], "tables/comparisons_for_MAST_test_using_seurat_clusters_cond1.txt", row.names = F, quote=F)
write.csv(comp[,2], "tables/comparisons_for_MAST_test_using_seurat_clusters_cond2.txt", row.names = F, quote=F)





