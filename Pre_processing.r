# Bradley August 2021
# scRNA-Seq - tidied pipeline, pre-processing
# Pre-processing.r

# Load libraries
setwd("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp")
source('analysis.r')

set <- "Epithelial"

pathOut <- paste("/exports/igmm/eddie/CCGG-tumour-WGS/BradTemp/BH_analysis/Seurat/August_2021", set, "Preprocessing", sep = "/")


# Checking the preprocessing of this data


# First want to make a new seurat object, with no removal of genes etc, using only the healthy epithelial cells
# Load in the correct counts
if(set == "Immune"){
	counts <- readMM('gene_sorted-Imm.matrix.mtx')
	rownames(counts) = readLines('Imm.genes.tsv')
	colnames(counts) = readLines('Imm.barcodes2.tsv')
}

if(set == "Epithelial"){
	counts <- readMM('gene_sorted-Epi.matrix.mtx')
	rownames(counts) = readLines('Epi.genes.tsv')
	colnames(counts) = readLines('Epi.barcodes2.tsv')
}



seur = CreateSeuratObject(counts=counts,  min.cells=0, min.features=0, names.delim='\\.')
meta = read.table('all.meta2.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)
meta1 <- meta[-1,]
meta1 <- meta1[rownames(meta1) %in% rownames(seur@meta.data),]
seur@meta.data <- cbind(seur@meta.data, meta1)

# Subset the seurat object
seur <- subset(seur, subset = Health == "Healthy")
if(set == "Immune"){
	# In this case, want to use the LP cells
	seur <- subset(seur, subset = Location == "LP")
}
if(set == "Epithelial"){
	seur <- subset(seur, subset = Location == "Epi")
}
seur <- subset(x=seur, idents = "N51", invert = T)
	


# Checking the nUMIs per cell
libSizeDf <- seur@meta.data[,c("nGene", "nUMI")]
libSizeDf$nGene <- as.numeric(libSizeDf$nGene)
libSizeDf$nUMI <- as.numeric(libSizeDf$nUMI)
#ggplot(libSizeDf, aes(x=nUMI)) + geom_histogram(bins = 50)
#ggplot(libSizeDf, aes(x=log10(nUMI))) + geom_histogram(bins = 50)



#Using a mixture model to fit 2 normal distibutions to the log UMIs. Identifying a mean for non-cell droplets and one for cell droplets
# NOTE: I don't think this is super neccessary - the previous plot clearly shows that all cells in this cohort have a good number of reads
set.seed(100)
library('mixtools')
# Make lib sizes a log scale
log10_lib_size <- log10(as.data.frame(libSizeDf)$nUMI)
# Fit mixture model
mix <- normalmixEM(log10_lib_size, maxrestarts=50, epsilon = 1e-03)

# plots
pdf(file = paste(pathOut, "mixed_model_nUMI_per_cell.pdf", sep = "/")) 
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
#log10_lib_size <- log10(libSizeDf$nUmis)
abline(v=split, lwd=2)
dev.off()

# Barcode rank plot
libSizeDf <- as.data.frame(libSizeDf)
barcode_rank <- rank(-libSizeDf$nUMI)
#plot(barcode_rank, libSizeDf$nUMI, ylab="library size")

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

n_cells <- 10000
# CellRanger
totals <- sort(libSizeDf$nUMI, decreasing = TRUE)
# 99th percentile of top n_cells divided by 10
thresh = totals[round(0.01*n_cells)]/10
#plot(log10(totals))
#abline(h=log10(thresh), col="red", lwd=2)
table(libSizeDf$nUMI >= thresh)



# Geerate the inflection/knee graphs on log-nUMI and log-Rank graph
# These represent realistic and overly sensitive estimates of sharp changes in nUMI/cell, as would be expected for whether a cell is present or not
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
# Exploring the nUMIs for each gene across all of the cells
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
feat_perc <- vector("list", length = length(features))
names(feat_perc) <- paste("percent", features, sep = ".")
for(x in seq(1:length(features))){
        feat_perc[[x]] <- PercentageFeatureSet(seur, pattern = features[x]);
        colnames(feat_perc[[x]]) <- features[x]
}


qc_h <- function(feat_perc, features){
        hist(feat_perc, breaks = 20, col = "grey80", xlab = paste("Proportion of reads in", features, "genes",  sep = " ") , main=features)
}

qc_plots <- vector("list", length = length(features))
for(x in seq(1:length(feat_perc))){
        qc_plots[[x]] <- qc_h(feat_perc = feat_perc[[x]], features = features[x]);
}


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

VlnPlot(seur, features = "MT", group.by = "Sample", log=T) + theme(legend.position = 'none')
VlnPlot(seur, features = "MT", group.by = "Sample") + theme(legend.position = 'none')

VlnPlot(seur, features = "RP", group.by = "Sample", log=T) + theme(legend.position = 'none')
VlnPlot(seur, features = "RP", group.by = "Sample") + theme(legend.position = 'none')



# Identifying low quality cells based on outlier detection - So this is not hard thresholding, but relative to the experiment
# Generate a function to calculate the MAD for these metrics
metrics <- c("nUMI", "nGene", "MT", "RP")
qc.lib <- seur@meta.data[, metrics]
qc.lib$nGene <- as.numeric(qc.lib$nGene)



# log the nUMI and nGene
qc.lib$nUMI <- as.numeric(qc.lib$nUMI)
qc.lib$nGene <- log10(qc.lib$nGene)
qc.lib$nUMI <- log10(qc.lib$nUMI)

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


# If we set the threshold at 2.5 MAD from the median, how many cells does that leave us? NOTE: IF the median + 2.5* MAD is <20%, use this instead - dont want to be too strict
mad_cutoff <- 2.5
# Not removing on the relative basis of MT% as this looks okay. Instead use the absolute cut off of 20%
keep <- qc.lib2[abs(qc.lib2$nUMI) < mad_cutoff & abs(qc.lib2$nGene) < mad_cutoff  & abs(qc.lib2$RP) < mad_cutoff, ]
if(median(qc.lib$MT)+(2.5*mad(qc.lib$MT)) < 20){
	keep <- keep[rownames(keep) %in% rownames(seur@meta.data[seur@meta.data$MT < 20,]),]
} else {
	keep <- qc.lib2[abs(qc.lib2$MT) < mad_cutoff,]
}
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
pdf(file = paste(pathOut, "MT_perc_cutoffs_relative_absolute.pdf", sep = "/"))
hist(as.numeric(seur@meta.data$MT), xlab = "MT%", main = "MT% distribution. Red = relative, Blue = hard") + abline(v=c(as.numeric(median(seur@meta.data$MT + (2.5*mad(seur@meta.data$MT)))), as.numeric(20)), col = c("red", "blue"))
dev.off()

# Are the distributions similar for each sample?
samps <-  unlist(strsplit(rownames(qc.lib), "\\."))
for(x in seq(1:nrow(qc.lib))){
        qc.lib$Patient[x] <- samps[(3*x)-2];
        qc.lib$Sample[x] <- paste(samps[(3*x)-2], samps[(3*x)-1], sep = ".");
        qc.lib$replicate[x] <-  samps[(3*x)-1]
}

pdf(file = paste(pathOut, "Sample_nUMI_distribution.pdf", sep = "/"))
qc.lib %>%
        ggplot(aes(x=nUMI, color = Patient, fill = Sample)) +
        xlab("log10 nUMI") +
        geom_density(alpha = 0.2) +
        theme_classic() + facet_grid(Patient ~ replicate) + theme(legend.position = "none")
dev.off()

pdf(file = paste(pathOut, "Sample_MT_distribution.pdf", sep = "/"))
qc.lib %>%
        ggplot(aes(x=MT, color = Patient, fill = Sample)) +
        xlab("MT %") +
        geom_density(alpha = 0.2) +
        theme_classic() + facet_grid(Patient ~ replicate) + theme(legend.position = "none")
dev.off()


sp <- VlnPlot(seur, features = "MT", group.by = "Sample", split.by = "retain")
sp <- ggplot(data.frame(seur@meta.data), aes(x=nGene, y=MT, col=retain)) +
  geom_point()
pdf(file = paste("Sample_MT_nGene_retain.pdf", sep = "/"))
sp + facet_wrap(~Sample)
dev.off()

# Analysing the novelty
# Cells with small library size tend to have higher overall ‘novelty’ i.e. they have not reached saturation for any given gene.
# Outlier cell may have a library with low complexity
# Expected novelty is around 0.8
seur@meta.data$nUMI <- as.numeric(seur@meta.data$nUMI)
seur@meta.data$nGene <- as.numeric(seur@meta.data$nGene)
seur@meta.data$replicate <- unlist(strsplit(seur@meta.data$Sample, "\\."))[c(F,T)]

pdf(file = paste(pathOut, "Overall_novelty.pdf", sep = "/"))
p <- seur@meta.data %>%
        data.frame() %>%
        ggplot(aes(x=nUMI, y=nGene, color=MT)) +
        geom_point() +
        stat_smooth(method=lm) +
        scale_x_log10() +
        scale_y_log10()
dev.off()

pdf(file = paste(pathOut, "Sample_novelty.pdf", sep = "/"))
p <- seur@meta.data %>%
        data.frame() %>%
        ggplot(aes(x=nUMI, y=nGene, color=MT)) +
        geom_point() +
        stat_smooth(method=lm) +
        scale_x_log10() +
        scale_y_log10() +
        facet_wrap(Subject ~ replicate) +
        geom_vline(xintercept = 800)
dev.off()

# Add the number of UMIs per gene for each cell to the metadata, then plot this
seur@meta.data$log10GenesPerUMI <- log10(seur@meta.data$nGene) / log10(seur@meta.data$nUMI)
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

if(median(qc.lib$MT)+(2.5*mad(qc.lib$MT)) < 20){
	mito.cells <- rownames(seur1@meta.data[seur1@meta.data$MT > 20,])
} else {
	mito.cells <- rownames(seur1@meta.data[seur1@meta.data$MT > median(qc.lib$MT)+(2.5*mad(qc.lib$MT)),])
}


length(sparse.cells)
length(mito.cells)
bad_cells <- c(mito.cells, sparse.cells)


good_cells <- setdiff(colnames(seur1@assays$RNA@counts), bad_cells)

#Finally, remove genes with < 20 genes expressing 
min.cells <- 1 - (20/length(cell_sparsity))
sparse.genes <- gene_sparsity > min.cells
good.genes <- setdiff(rownames(seur1[['RNA']]), names(sparse.genes[sparse.genes == T]))


seur_final <- subset(seur1, cells = good_cells)
seur_final <- subset(seur_final, features = good.genes)


save(seur_final, file = paste(pathOut, "raw_HealthyImmune_seur_GeneCellQC_noN51.RData", sep = "/"))





