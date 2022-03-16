# Bradley October 2021
# scRNASeq preprocessing of Elmentaite (Epithelial) cells, aligned using 10X 
# kArray.r

# Set up and load libraries
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
setwd("Elmentaite/results")
library('Seurat')
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

# Load object
load("Initial_clustering/seur_integrated_UMAPd.Rds")


# Read in the value of k to be used here
 option_list = list(make_option(c("-k", "--kparam"), action = "store", default = NA, type ="character",
                help="k.parameter"))
opt = parse_args(OptionParser(option_list=option_list))
k.param = opt$k;
print(paste("The option for k is", k.param, sep = " "))


# Define the options for this script. What needs to be done:
do_clustering = T
generate_markers = T
enrichment_vs_authors = T
pathOut = "kArray"


# Do clustering if wanted
if(do_clustering == T){
        # load data
        load("Initial_clustering/seur_integrated_log2TP10K.Rds")
        # Do clustering
        DefaultAssay(seur.integrated) <- "integrated"
        seur.integrated <- FindNeighbors(seur.integrated, reduction = "pca", dims = 1:50, k.param = as.numeric(k.param))
        # Now find clusters (using usual resolution as this seems to work well)
        seur.integrated <- FindClusters(seur.integrated, resolution = 0.6)
        # Now save this for exploring later
       save(seur.integrated, file = paste(pathOut, "/seur.integrated.clusterd.k", k.param, ".Rds", sep = ""))
}


# Generate the markers if wanted
if(generate_markers == T){
        # Generate markers
        library('MAST')
        packageVersion('MAST')
        DefaultAssay(seur.integrated) <- "log2TP10K"
        Idents(seur.integrated) = "seurat_clusters"
        # Using the same thresholds as have done before
        markers <- FindAllMarkers(seur.integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
        # Save
        write.csv(markers, paste(pathOut, "/k", k.param, "markers.csv", sep = ""))
} else {
        markers <- read.csv(paste(pathOut, "/k", k.param, ".csv", sep = ""))
}

# Do the enrichment against the putative authors clusters
if(enrichment_vs_authors == T){
        # Load in the smillie sets
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
        write.csv(fgsea_Res, paste(pathOut, "/k", k.param, ".fgsea_vs_authors.csv", sep = ""))
}









