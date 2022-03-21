# Bradley October 2021
# scRNASeq preprocessing of Elmentaite (Epithelial) cells, aligned using 10X 
# Preprocessing.Elmentaite.r

# Set up and load libraries
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')
setwd("Elmentaite/data")

# Read in the h5 file 
# Note: before loading, need to make sure in the shell I have ran: 
# module load igmm/apps/hdf5/1.8.13
# And set the following variable
# export C_INCLUDE_PATH=/exports/igmm/software/pkg/el7/apps/hdf5/1.8.13/include:$C_INCLUDE_PATH
seur <- Read10X_h5("cell_ranger_aggr/Elmentaite_aggr/outs/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
# The filtered matrix is comprised of a total of 33354 genes and 19497 cells


# Add the sample names to the cells
samp_to_num <- read.csv("cell_ranger_aggr/Elmentaite_aggr.csv")
samp_to_num <- data.frame(num = rownames(samp_to_num), Sample = samp_to_num$library_id)
for(samp in 1:nrow(samp_to_num)){
	colnames(seur) <- gsub(paste("-", samp_to_num$num[samp], sep = ""), paste("-", samp_to_num$Sample[samp], sep = ""), colnames(seur))
}

# Need to read in the meta
meta <- read.csv("Elmentaite_meta_usable.csv")

# Now create a seurat object
seur = CreateSeuratObject(counts=seur,  min.cells=0, min.features=0, names.delim='\\.')


# Now reformat the meta
meta1 <- seur@meta.data
meta1$orig.ident <- unlist(strsplit(rownames(meta1), "\\-"))[c(F,T)]
colnames(meta)[which(colnames(meta) == "Source.Name")] <- "orig.ident"
cells <- rownames(meta1)
meta2 <- merge(meta1, meta, by = "orig.ident")
if(all(meta2$nCount_RNA == meta1$nCount_RNA)){
	rownames(meta2) <- cells
}
seur@meta.data <- meta2

###### Now ready for the preprocessing

# Generate the pathOut to store all of the plots
setwd("../results")
pathOut = "Preprocessing"
# Checking the nCount_RNAs per cell
libSizeDf <- seur@meta.data[,c("nFeature_RNA", "nCount_RNA")]
libSizeDf$nFeature_RNA <- as.numeric(libSizeDf$nFeature_RNA)
libSizeDf$nCount_RNA <- as.numeric(libSizeDf$nCount_RNA)


pdf(file = paste(pathOut, "nCount_RNA_histos.pdf", sep = "/"))
p1  = ggplot(libSizeDf, aes(x= nCount_RNA)) + geom_histogram(bins = 50)
p2 = ggplot(libSizeDf, aes(x=log10(nCount_RNA))) + geom_histogram(bins = 50)
p1 + p2
dev.off()

#Using a mixture model to fit 2 normal distibutions to the log UMIs. Identifying a mean for non-cell droplets and one for cell droplets
# NOTE: I don't think this is super neccessary - the previous plot clearly shows that all cells in this cohort have a good number of reads
set.seed(100)
library('mixtools')
# Make lib sizes a log scale
log10_lib_size <- log10(as.data.frame(libSizeDf)$nCount_RNA)
# Fit mixture model
mix <- normalmixEM(log10_lib_size, maxrestarts=50, epsilon = 1e-03)

# plots
pdf(file = paste(pathOut, "mixed_model_nCount_RNA_per_cell.pdf", sep = "/")) 
plot(mix, which=2, xlab2="log(mol per cell)")
p1 <- dnorm(log10_lib_size, mean=mix$mu[1], sd=mix$sigma[1])
p2 <- dnorm(log10_lib_size, mean=mix$mu[2], sd=mix$sigma[2])
# find intersection:
if (mix$mu[1] < mix$mu[2]) {
    split <- min(log10_lib_size[p2 > p1])
} else {
    split <- min(log10_lib_size[p1 > p2])
}
# show split on plot:
#log10_lib_size <- log10(libSizeDf$nCount_RNAs)
abline(v=split, lwd=2)
dev.off()

# Barcode rank plot
libSizeDf <- as.data.frame(libSizeDf)
barcode_rank <- rank(-libSizeDf$nCount_RNA)
#plot(barcode_rank, libSizeDf$nCount_RNA, ylab="library size")

# log
#plot(barcode_rank, log10_lib_size)

o <- order(barcode_rank)
log10_lib_size <- log10_lib_size[o]
barcode_rank <- barcode_rank[o]

rawdiff <- diff(log10_lib_size)/diff(barcode_rank)
inflection <- which(rawdiff == min(rawdiff[100:length(rawdiff)], na.rm=TRUE))

pdf(file = paste(pathOut, "log_barcode_rank_plot_inflection.pdf", sep = "/"))
plot(barcode_rank, log10_lib_size)
abline(v=inflection, col="red", lwd=2)
dev.off()


# CellRanger v1 and v2
# Given an expected number of cells, cellranger used to assume a ~10-fold range of library sizes for real cells and estimate this range (cellranger (v1 and v2). 
# The threshold was defined as the 99th quantile of the library size, divided by 10.

n_cells <- 1000
# CellRanger
totals <- sort(libSizeDf$nCount_RNA, decreasing = TRUE)
# 99th percentile of top n_cells divided by 10
thresh = totals[round(0.01*n_cells)]/10
#plot(log10(totals))
#abline(h=log10(thresh), col="red", lwd=2)
table(libSizeDf$nCount_RNA >= thresh)



# Geerate the inflection/knee graphs on log-nCount_RNA and log-Rank graph
# These represent realistic and overly sensitive estimates of sharp changes in nCount_RNA/cell, as would be expected for whether a cell is present or not
library(DropletUtils)
my.counts <- seur@assays$RNA@counts
br.out <- barcodeRanks(my.counts)

pdf(file = paste(pathOut, "Cell_ranger_knee_inflection_log_barcode_rank.pdf", sep = "/"))
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total UMI count")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
    legend=c("knee", "inflection"))
dev.off()
inflection


# Can assess emptyDroplet using simulation - EDIT: this cannot be done because have not got an unfiltered matrix. 
# Use FDR of 0.1%, so at most 1 in 1000 droplets are called empty
#set.seed(100)
#e.out <- emptyDrops(m=my.counts)
#e.out


# Assuming there are no remaining empty droplets
# Exploring the nCount_RNAs for each gene across all of the cells
# randomly subset 1000 cells.
set.seed(100)
tmpCounts <- my.counts[,sample(1000)]
plot(rowSums(my.counts),
     rowMeans(my.counts == 0),
     log = "x",
     xlab="total number of UMIs",
     ylab="proportion of cells expressing the gene"
)


# Count the genes that are not expressed (detected)
not.expressed <- rowSums(my.counts) == 0
table(not.expressed)
not_expr <- rownames(my.counts)[rowSums(my.counts) == 0]


# Compute the relative expression of each gene per cell
#rel_expression <- t( t(my.counts) / Matrix::colSums(my.counts)) * 100
#rownames(rel_expression) <- rowData(seur)$Symbol
#most_expressed <- sort(Matrix::rowSums( rel_expression ),T)[20:1] / ncol(sce)
rel_expression <- (my.counts/Matrix::colSums(my.counts)) * 100
most_expressed <- sort(rowSums(rel_expression), T)[20:1] / ncol(rel_expression)


pdf(file = paste(pathOut, "Top_expressed_genes.pdf", sep = "/"))
boxplot( as.matrix(t(rel_expression[names(most_expressed),])),cex=.1, las=1, xlab="% total count per cell",col=scales::hue_pal()(20)[20:1],horizontal=TRUE)
dev.off()

# Now doing actual normalisation
# Above analysis allows inference for bad quality droplets. But what about bad quality cells
# High % of mitochondrial genes, high levels of ambient RNA contamination etc is indicative of poor quality cells.

# mitochondrial genes, doing this like is done in the seurat workflow
features <- c("^MT-", "^RP", "^HLA-", "^RNA", "^IG[HJKL]")
feature.names <- c("MT", "RP", "HLA", "RNA", "IG_HJKL")
feat_perc <- vector("list", length = length(features))
names(feat_perc) <- paste("percent", features, sep = ".")
for(x in seq(1:length(features))){
        feat_perc[[x]] <- PercentageFeatureSet(seur, pattern = features[x]);
        colnames(feat_perc[[x]]) <- feature.names[x]
        feat_perc[[x]][,1] <- as.numeric(feat_perc[[x]][,1])
}


qc_h <- function(feat_perc, features){
	ggplot(feat_perc, aes_string(x=features)) + geom_histogram()
}
qc_plots <- vector("list", length = length(features))

for(x in seq(1:length(feat_perc))){
	p <- qc_h(feat_perc = feat_perc[[x]], features = feature.names[x])
	qc_plots[[x]] <- p
	}

pdf(paste(pathOut, "QC_gene_histos.pdf", sep = "/"))
par(mfrow=c(3,2))
qc_plots[[1]] + qc_plots[[2]] + qc_plots[[3]] + qc_plots[[4]] + qc_plots[[5]]
dev.off()



# Calculating the effects of different thresholds for MT
MT_thresh <- data.frame(cutoff = seq(from = 0, to = 100, by = 5))
MT_thresh$cells <- rep(0, nrow(MT_thresh))
MT <- feat_perc[[1]]
for(r in seq(1:nrow(MT_thresh))){
        MT_thresh[r,]$cells <- length(MT[MT[,1] < MT_thresh[r,1],])
        }
write.csv(MT_thresh, paste(pathOut, "MT_cutoffs.csv", sep = "/"))


# Looking to see if MT or RP % is specific to particular samples
# First add these to the meta
seur[["MT"]] <- PercentageFeatureSet(object = seur, pattern = "^MT-")
seur[["RP"]] <- PercentageFeatureSet(object = seur, pattern = "^RP")

pdf(file=paste(pathOut, "MT_perc_by_sample.pdf", sep = "/"))
VlnPlot(seur, features = "MT", group.by = "orig.ident") + theme(legend.position = 'none')
dev.off()

pdf(file=paste(pathOut, "RP_perc_by_sample.pdf", sep = "/"))
VlnPlot(seur, features = "RP", group.by = "orig.ident") + theme(legend.position = 'none')
dev.off()


# Identifying low quality cells based on outlier detection - So this is not hard thresholding, but relative to the experiment
# Generate a function to calculate the MAD for these metrics
metrics <- c("nFeature_RNA", "nCount_RNA", "MT", "RP")
qc.lib <- seur@meta.data[, metrics]
qc.lib$nFeature_RNA <- as.numeric(qc.lib$nFeature_RNA)



# log the nCount_RNA and nFeature_RNA
qc.lib$nCount_RNA <- as.numeric(qc.lib$nCount_RNA)
qc.lib$nFeature_RNA <- log10(qc.lib$nFeature_RNA)
qc.lib$nCount_RNA <- log10(qc.lib$nCount_RNA)

# Calculate the number of mad a cells metric value is from the median
mad_trans <- function(x){
        x <- (x-median(x))/mad(x, constant = 1.4826);
        return(x)
        }


qc.lib2 <- sapply(qc.lib, FUN=mad_trans)
qc.lib2 <- as.data.frame(qc.lib2)
rownames(qc.lib2) <- rownames(qc.lib)

hist_plots <- vector("list", length=ncol(qc.lib2))
for(x in seq(1:ncol(qc.lib2))){
        hist_plots[[x]] <- hist(qc.lib2[,x], plot=F);
        names(hist_plots[[x]]) <- colnames(qc.lib2)[x];
        return(hist_plots[[x]])
}


res <- lapply(qc.lib2, hist, plot = FALSE)

pdf(file = paste(pathOut, "Relative_gene_expression_variability_MAD.pdf", sep = "/"))
par(mfrow=c(2,2))
plot(res[[1]], main = colnames(qc.lib2)[1], xlab = "MAD from median")
plot(res[[2]], main = colnames(qc.lib2)[2], xlab = "MAD from median")
plot(res[[3]], main = colnames(qc.lib2)[3], xlab = "MAD from median")
plot(res[[4]], main = colnames(qc.lib2)[4], xlab = "MAD from median")
dev.off()

#~~~~ Subsetting the cells based on QC metrics 

# If we set the threshold at 2.5 MAD from the median, how many cells does that leave us? NOTE: IF the median + 2.5* MAD is <20%, use this instead - dont want to be too strict
mad_cutoff <- 2.5
# Not removing on the relative basis of MT% as this takes us above 100%
# Even if use a cut off comparable to the other dataset, this leaves us with 1900 cells only
# For now, won't subset on the basis of MT% or RP%

# Old way:
#keep <- qc.lib2[abs(qc.lib2$nCount_RNA) < mad_cutoff & abs(qc.lib2$nFeature_RNA) < mad_cutoff  & abs(qc.lib2$RP) < mad_cutoff, ]

# Not on the basis of MT% or RP% 
keep <- qc.lib2[abs(qc.lib2$nCount_RNA) < mad_cutoff & abs(qc.lib2$nFeature_RNA) < mad_cutoff,]

#### From other script, adjusted to set the relative top threshold based on Smillie analysis
#if(median(qc.lib$MT)+(2.5*mad(qc.lib$MT)) < 20){
#	keep <- keep[rownames(keep) %in% rownames(seur@meta.data[seur@meta.data$MT < 20,]),]
#	} else {
#		if(median(qc.lib$MT)+(2.5*mad(qc.lib$MT)) > 56){
#			keep <- keep[rownames(keep) %in% rownames(seur@meta.data[seur@meta.data$MT < 56,]),]
#			} else {
#				keep <- qc.lib2[abs(qc.lib2$MT) < mad_cutoff,]
#			}
#		}
nrow(keep)


# Diagnostic plots to see the cells we have discarded and retained
keep$retain <- c(rep("TRUE", length = nrow(keep)))
discard_cells <- setdiff(rownames(qc.lib2), rownames(keep))
discard <- qc.lib2[rownames(qc.lib2) %in% discard_cells,]
discard$retain <- c(rep("FALSE", length = nrow(discard)))
together <- rbind(keep, discard)

together <- together[order(rownames(together)),]
seur@meta.data <- seur@meta.data[order(rownames(seur@meta.data)),]

all(rownames(seur@meta.data) == rownames(together))
seur@meta.data$retain <- together$retain


# Plot
#pdf(file = paste(pathOut, "MT_perc_cutoffs_relative_absolute.pdf", sep = "/"))
#hist(as.numeric(seur@meta.data$MT), xlab = "MT%", main = "MT% distribution. Red = median, Blue = hard") + abline(v=c(as.numeric(median(seur@meta.data$MT)), as.numeric(20)), col = c("red", "blue"))
#dev.off()

# Are the distributions similar for each sample?
qc.lib$Sample <-  unlist(strsplit(rownames(qc.lib), "\\-"))[c(F,T)]
qc.lib$Sample.short <-  unlist(strsplit(qc.lib$Sample, "\\_"))[c(F,F,T)]


pdf(file = paste(pathOut, "Sample_nCount_RNA_distribution.pdf", sep = "/"))
qc.lib %>%
        ggplot(aes(x=nCount_RNA, color = Sample.short, fill = Sample.short)) +
        xlab("log10 nCount_RNA") +
        geom_density(alpha = 0.2) +
        theme_classic() + facet_grid(rows = "Sample.short" ) + theme(legend.position = "none")
dev.off()

pdf(file = paste(pathOut, "Sample_MT_distribution.pdf", sep = "/"))
qc.lib %>%
        ggplot(aes(x=MT, color = Sample.short, fill = Sample.short)) +
        xlab("MT %") +
        geom_density(alpha = 0.2) +
        theme_classic() + facet_grid(rows = "Sample.short") + theme(legend.position = "none")
dev.off()


sp <- VlnPlot(seur, features = "MT", group.by = "orig.ident", split.by = "retain")
sp <- ggplot(data.frame(seur@meta.data), aes(x=nFeature_RNA, y=MT, col=retain)) +
  geom_point()
pdf(file = paste("Sample_MT_nFeature_RNA_retain.pdf", sep = "/"))
sp + facet_wrap(~orig.ident)
dev.off()

# Analysing the novelty
# Cells with small library size tend to have higher overall ‘novelty’ i.e. they have not reached saturation for any given gene.
# Outlier cell may have a library with low complexity
# Expected novelty is around 0.8
seur@meta.data$nCount_RNA <- as.numeric(seur@meta.data$nCount_RNA)
seur@meta.data$nFeature_RNA <- as.numeric(seur@meta.data$nFeature_RNA)

pdf(file = paste(pathOut, "Overall_novelty.pdf", sep = "/"))
p <- seur@meta.data %>%
        data.frame() %>%
        ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MT)) +
        geom_point() +
        stat_smooth(method=lm) +
        scale_x_log10() +
        scale_y_log10()
p
dev.off()

pdf(file = paste(pathOut, "Sample_novelty.pdf", sep = "/"))
p <- seur@meta.data %>%
        data.frame() %>%
        ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MT)) +
        geom_point() +
        stat_smooth(method=lm) +
        scale_x_log10() +
        scale_y_log10() +
        facet_wrap(~orig.ident) +
        geom_vline(xintercept = 800)
p
dev.off()

# Add the number of UMIs per gene for each cell to the metadata, then plot this
seur@meta.data$log10GenesPerUMI <- log10(seur@meta.data$nFeature_RNA) / log10(seur@meta.data$nCount_RNA)
mean_log10GenesPerUMI <- seur@meta.data %>% pull(log10GenesPerUMI) %>% mean() %>% signif(6)

#p <- seur@meta.data %>%
#        ggplot(aes(x=log10GenesPerUMI, color = Subject, fill = Sample)) +
#        xlab("log10GenesPerUMI") +
#        geom_density(alpha = 0.2) +
#        theme_classic() +
#        geom_vline(xintercept = mean_log10GenesPerUMI) +
#        facet_grid(Subject ~ replicate) + theme(legend.position = "none")



# Analysing cell sparsity
# Remove genes that are not expressed at all
not.expressed <- rowSums(seur@assays$RNA@counts) == 0
bad_genes <- names(not.expressed[not.expressed != F])
genes.use <- setdiff(rownames(seur@assays$RNA@counts), bad_genes)


seur1 <- subset(seur, features = genes.use)


seur1[["MT"]] <- PercentageFeatureSet(object = seur1, pattern = "^MT-")
seur1[["RP"]] <- PercentageFeatureSet(object = seur1, pattern = "^RP")


cell_sparsity <- apply(seur1@assays$RNA@counts == 0, 2, sum)/nrow(seur1@assays$RNA@counts)
gene_sparsity <- apply(seur1@assays$RNA@counts == 0, 1, sum)/ncol(seur1@assays$RNA@counts)



p1 <- hist(cell_sparsity, breaks=50, col="grey80", xlab="Cell sparsity", main="")
p2 <- hist(gene_sparsity, breaks=50, col="grey80", xlab="Gene sparsity", main="")

pdf(file = paste(pathOut, "Cell_and_gene_sparsity", sep = "/"))
par(mfrow=c(1,2))
plot(p1)
plot(p2)
dev.off()

# Filtering to produce the analysis ready seurat object - have changed this from epithelial, so that if 2.5*mad < 20%, use 20% as the cut off
sparse.cells <- names(cell_sparsity[cell_sparsity > 0.99])

# NOT Subsetting based on MT%
#if(median(qc.lib$MT)+(2.5*mad(qc.lib$MT)) < 20){
#	mito.cells <- rownames(seur1@meta.data[seur1@meta.data$MT > 20,])
#} else {
#	mito.cells <- rownames(seur1@meta.data[seur1@meta.data$MT > median(qc.lib$MT)+(2.5*mad(qc.lib$MT)),])
#}


length(sparse.cells)
bad_cells <- c(sparse.cells)


good_cells <- setdiff(colnames(seur1@assays$RNA@counts), bad_cells)

#Finally, remove genes with < 20 genes expressing 
min.cells <- 1 - (20/length(cell_sparsity))
sparse.genes <- gene_sparsity > min.cells
good.genes <- setdiff(rownames(seur1[['RNA']]), names(sparse.genes[sparse.genes == T]))


seur_final <- subset(seur1, cells = good_cells)
seur_final <- subset(seur_final, features = good.genes)

# Before doublet removal, have a total of 
save(seur_final, file = paste(pathOut, "raw_seur_GeneCellQC.RData", sep = "/"))









