# Harris_et_al_WGCNAscRNA
Repository for code used in the analysis of scRNAseq data using WGCNA


## Pipeline order
To replicate the analysis performed in this paper, download the datasets into an accessible directory, along with the analysis.r script (which will need ot be adjusted for your own R environment). Scripts can then be ran as follows

### Preprocessing and clustering - Smillie
  1. Pre_processing.r
  2. kArray_clustering_gsea_vs_authors.r
  3. kArray_sum.r
  4. Doublet.detection.r
  5. Fully_processed_clustering.r
  6. Cluster_merging.r
  7. Interesting_cluster_analysis.temp.R

### WGCNA - Smillie
  1. WGCNA_across_clusters_Peters.r
  2. Pairwise_module_preservation_WGCNA.r
  3. Pariwise_module_preservation_WGCNA_summary.r
  4. Gene-gene_correlations_intracultuster_from_WGCNA.r

### Abundance testing - Smillie
  1. lm_of_cluster11.r
  2. DA_milo.R

### Preprocessing and clustering - Elmentaite
  1. 

### Genotyping and variability analysis
  1. Freebayes_variant_calling.bash
  2. samtools_variant_calling.bash
  3. Comparing_Elmentaite_grouping_wrt_Smillie_grouping.r

Any questions about the code or analysis please contact bradley.harris@ed.ac.uk
