# Bradley August 2021
# scRNA-Seq - tidied pipeline, DE between clusters for analysis of whether merging is neccessary
# DE_between_clusters.r

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library('Seurat')
library(optparse)

# Choose set and define path
set <- "Epithelial"
pathOut = paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, sep = "/")


#### Now read in the object 
load(paste(pathOut, "objects/Processed_seur.integrated.clusteredRds", sep = "/"))



# Now set up the options for DE
option_list = list(
       make_option(c("-c", "--cond1"), action = "store", default = NA, type ="character",
               help="What is the first cluster to compare"),
        make_option(c("-d", "--cond2"), action = "store", default = NA, type ="character",
                help="What is the second cluster to compare")
       )
opt = parse_args(OptionParser(option_list=option_list))

print("Option for cluster 1 is")
opt$cond1
print("Option for cluster 2 is")
opt$cond2


print("Doing MAST")
DefaultAssay(seur.integrated) <- "log2TP10K"
test <- FindMarkers(seur.integrated, ident.1 = opt$cond1, ident.2 = opt$cond2, test.use = "MAST")
test$cond1 <- rep(opt$cond1, nrow(test))
test$cond2 <- rep(opt$cond2, nrow(test))
head(test)
write.table(test, paste(pathOut, "/tables/pairwiseDEonSeuratClusters/MAST/MAST_", opt$cond1, "vs", opt$cond2, ".txt",sep = ""), sep = "\t") 
print("Done MAST")

print("Doing ROC")
DefaultAssay(seur.integrated) <- "log2TP10K"
test <- FindMarkers(seur.integrated, ident.1 = opt$cond1, ident.2 = opt$cond2, test.use = "roc", min.pct = 0.1, logfc.threshold = 0.25)
test$cond1 <- rep(opt$cond1, nrow(test))
test$cond2 <- rep(opt$cond2, nrow(test))
head(test)
write.table(test, paste(pathOut, "/tables/pairwiseDEonSeuratClusters/roc/roc_", opt$cond1, "vs", opt$cond2, ".txt", sep = ""), sep = "\t")
print("Done ROC")



##### End of script



