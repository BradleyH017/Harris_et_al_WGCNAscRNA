# Bradley October 2021
# scRNASeq preprocessing of Elmentaite (Epithelial) cells, aligned using 10X 
# Pairwise DE between clusters to detect potential for merging

# Set up
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
library(optparse)
setwd("Elmentaite/results")

# Read in the object
load("objects/Processed_seur.integrated.clusteredRds")


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
write.table(test, paste("tables/pairwiseDEonSeuratClusters/MAST/MAST_", opt$cond1, "vs", opt$cond2, ".txt",sep = ""), sep = "\t")
print("Done MAST")

print("Doing ROC")
DefaultAssay(seur.integrated) <- "log2TP10K"
test <- FindMarkers(seur.integrated, ident.1 = opt$cond1, ident.2 = opt$cond2, test.use = "roc", min.pct = 0.1, logfc.threshold = 0.25)
test$cond1 <- rep(opt$cond1, nrow(test))
test$cond2 <- rep(opt$cond2, nrow(test))
head(test)
write.table(test, paste("tables/pairwiseDEonSeuratClusters/roc/roc_", opt$cond1, "vs", opt$cond2, ".txt", sep = ""), sep = "\t")
print("Done ROC")

# Remember to remove the 'x' from the first line of each of these








