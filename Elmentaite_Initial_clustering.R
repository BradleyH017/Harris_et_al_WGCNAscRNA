# Bradley October 2021
# scRNASeq preprocessing of Elmentaite (Epithelial) cells, aligned using 10X 
# Initial_clustering

# Set up and load libraries
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
setwd("Elmentaite/results")

# Read in the processed object from the relevant set
load("Preprocessing/raw_seur_GeneCellQC.RData")
seur <- seur_final
rm(seur_final)


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
save(seur.integrated, file = "Initial_clustering/seur_integrated_UMAPd.Rds")


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
save(seur.integrated, file = "Initial_clustering/seur_integrated_log2TP10K.Rds")
 






