# Bradley August 2021
# scRNA-Seq - tidied pipeline, initial clustering, dimensionality reduction, pre-clustering and marker calculation in arrays
# Initial_clustering_dimred.r

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library('Seurat')

# Choose set and define path
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, sep = "/")


# Read in the processed object from the relevant set
load(paste(pathOut, "Preprocessing/raw_HealthyImmune_seur_GeneCellQC_noN51.RData", sep = "/"))
seur <- seur_final
rm(seur_final)


# Perform the integration of data across samples
# ReciprocalPCA integration step
# Want to regress out MT
#https://satijalab.org/seurat/articles/integration_rpca.html

# Doing the integration with the doublet retains - Only need to do this once, then do the kArray to see which k value best atches the authors clusters and is likley tobe the most useful
split.by = "Sample"
vars.to.regress = "MT"

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
seur.integrated@meta.data$nUMI <- as.numeric(seur.integrated@meta.data$nUMI)
seur.integrated@meta.data$nGene <- as.numeric(seur.integrated@meta.data$nGene)


# Before saving, also calculate the log2TP10K matrix, so that this can be used in Authors cluster marker analysis (here) and seurat cluster marker analysis (next in kArray)
## Suggested in https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAST-interoperability.html
##https://github.com/satijalab/seurat/discussions/4032
## That this is also scale normalized, could use cpm or tmm for this. Using cpm
## Split into 34 matrices for computation
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
save(seur.integrated, file = paste(pathOut, "objects/seur_integrated_UMAPd.Rds", sep = "/"))


# Now calculate the markers of the authors clusters, then save this as csv file for each cluster and reformat so that is readable during the gsea
library('MAST')
	packageVersion('MAST')
        DefaultAssay(seur.integrated) <- "log2TP10K"
        Idents(seur.integrated) = "Cluster"
        # Using the same thresholds as have done before
        markers <- FindAllMarkers(seur.integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
        write.csv(markers, paste(pathOut, "kArray/Authors_clusters.allMarkers.csv", sep = "/"))
        # Need to then format as a gmt file
        smil_markers <- markers
        rm(markers)
        smil_markers <- smil_markers[smil_markers$p_val_adj < 0.05,] # Subset for only the significant
        clusters <- levels(factor(smil_markers$cluster))
        cluster.names <- clusters
        #cluster.names <- gsub("\\+", "", clusters)
        cluster.names <- gsub("\\ ", ".", cluster.names)
        smil_list <- vector("list", length = length(clusters))
        smil_vec <- vector("list", length = length(clusters))
        names(smil_vec) <- cluster.names
        for(x in seq(1:length(clusters))){
                smil_list[[x]] <- smil_markers[smil_markers$cluster == clusters[x],];
                smil_vec[[x]] <- smil_list[[x]]$gene;
                smil_vec[[x]] <- c(cluster.names[x], c("NA", smil_vec[[x]]));
                write.table(smil_vec[[x]], paste(pathOut, "/tables/markers/", cluster.names[x], ".list.gmt", sep = ""), sep = "\t", row.names = F, quote = F);


        ## Then need to go and manually edit these from the command line: 
        ## To remove the x
        #for f in *.list.gmt; do
        #tail -n +2 "$f" > "${f}".tmp && mv "${f}".tmp "$f"
        #echo "Processing $f"
        #done

        ## To transpose
        #for f in *.list.gmt; do
        # awk '
        #{ 
        #    for (i=1; i<=NF; i++)  {
        #        a[NR,i] = $i
        #    }
        #}
        #NF>p { p = NF }
        #END {    
        #    for(j=1; j<=p; j++) {
        #        str=a[1,j]
        #        for(i=2; i<=NR; i++){
        #            str=str" "a[i,j];
        #        }
        #        print str
        #    }
        #}' $f > 2$f
        #done

        ## To ensure is tab delim
        #for f in 2*; do
        #sed 's/ /\t/g' $f > 3$f
        #done


        ## Now remove the intermediate files and rename
        #rm 2*

        #for filename in 32*; do 
        #    [ -f "$filename" ] || continue
        #    mv "$filename" "${filename//32/}"
        #       
        #done
        #
        #}



# Then can run the kArray script to calculate the markers of each cluster at a range of k
