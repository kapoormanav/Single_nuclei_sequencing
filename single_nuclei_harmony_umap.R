####Combine data from 40 brains (~280,000 nuclei) using Harmony

combined_coga2 = merge(
  x = test_sample,
  y = reference.list2,
  add.cell.ids = sample_name,
  merge.data = FALSE,
  project = "combined_coga"
)


combined_coga2[["percent.mt"]] <- PercentageFeatureSet(combined_coga2, pattern = "^MT-")
	combined_coga2 <- subset(combined_coga2, subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 25)
	combined_coga2 <- NormalizeData(combined_coga2, normalization.method = "LogNormalize", scale.factor = 10000)
	combined_coga2 <- FindVariableFeatures(combined_coga2, selection.method = "vst", nfeatures = 2000)
    combined_coga2 <- ScaleData(object = combined_coga2, features = rownames(combined_coga2))
    
    RunPCA(pc.genes = combined_coga2@var.genes, npcs = 50, verbose = FALSE)
    
    combined_coga2 <- RunPCA(combined_coga2, npcs = 50, verbose = FALSE)
    
    coga_f <- combined_coga2 %>%
            RunHarmony("stim", plot_convergence = F)
            
            combined_coga2 <- FindNeighbors(combined_coga2, reduction = "pca", dims = 1:50, nn.eps = 0.5)
            
            combined_coga2 <- FindClusters(combined_coga2, resolution = 0.9, n.start = 10)
            
            library(ggplot2)
combined_coga2 <- RunTSNE(combined_coga2, dims = 1:50, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)
combined_coga2 <- RunUMAP(combined_coga2, dims = 1:50, min.dist = 0.9)

p1 <- DimPlot(combined_coga2, reduction = "tsne", pt.size = 0.25, label =T) + ggtitle(label = "FIt-SNE")
p2 <- DimPlot(combined_coga2, reduction = "umap", pt.size = 0.25) + ggtitle(label = "UMAP")
p1 <- AugmentPlot(plot = p1)
p2 <- AugmentPlot(plot = p2)


p3 = CombinePlots(plots = list(p1, p2), legend = "none")

pdf("~/integrated_tsne.plot_microglial.pdf", width = 14, height = 10)
print (p4)
dev.off()

coga_f <- RunHarmony(combined_coga2, 'dataset' , plot_convergence = F)

p3 <- FeaturePlot(combined_coga2, features = c("SERPINI1", "VSNL1", "GAP43", "SNAP25"), reduction = "tsne", pt.size = 0.1, combine = FALSE, order =T, min.cutoff = 'q10')
p3 <- lapply(X = p3, FUN = function(x) AugmentPlot(x + DarkTheme() + NoLegend()))
p4 = CombinePlots(plots = p3)

combined_coga2.markers <- FindAllMarkers(combined_coga2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

p3 <- FeaturePlot(combined_coga3, features = c("CD74", "CX3CR1", "MS4A7", "SLC2A5"), reduction = "umap", pt.size = 0.1, combine = FALSE, order =T, min.cutoff = 'q10')
p3 <- lapply(X = p3, FUN = function(x) AugmentPlot(x + DarkTheme() + NoLegend()))
p4 = CombinePlots(plots = p3)

pdf("~/integrated_umap.plot_microglial.pdf", width = 14, height = 10)
print (p4)
dev.off()

p3 <- FeaturePlot(combined_coga3, features = c("ALDH1L1", "F3", "SLC1A4", "CLU"), reduction = "umap", pt.size = 0.1, combine = FALSE, order =T, min.cutoff = 'q50', cols = c("lightgrey", "blue", "darkblue"))
#p3 <- lapply(X = p3, FUN = function(x) AugmentPlot(x + DarkTheme() + NoLegend()))
p4 = CombinePlots(plots = p3)

pdf("~/integrated_umap.plot_astrocytes.pdf", width = 14, height = 10)
print (p4)
dev.off()

p3 <- FeaturePlot(combined_coga3, features = c("MAPT", "SPI1", "CRHR1", "FKBP5"), reduction = "umap", pt.size = 0.1, combine = FALSE, order =T, min.cutoff = 'q50', cols = c("lightgrey", "blue", "darkblue"))
#p3 <- lapply(X = p3, FUN = function(x) AugmentPlot(x + DarkTheme() + NoLegend()))
p4 = CombinePlots(plots = p3)

pdf("~/integrated_umap.plot_alcohol1.pdf", width = 14, height = 10)
print (p4)
dev.off()

p3 <- FeaturePlot(combined_coga3, features = c("ADH1B", "ALDH2", "KCNJ6", "SULT1A1"), reduction = "umap", pt.size = 0.1, combine = FALSE, order =T, min.cutoff = 'q50', cols = c("lightgrey", "blue", "darkblue"))
#p3 <- lapply(X = p3, FUN = function(x) AugmentPlot(x + DarkTheme() + NoLegend()))
p4 = CombinePlots(plots = p3)

pdf("~/integrated_umap.plot_alcohol2.pdf", width = 14, height = 10)
print (p4)
dev.off()

pc1 <- FetchData(object = combined_coga3, vars = c('orig.ident' , 'seurat_clusters', 'nFeature_RNA', 'nCount_RNA', 'percent.mt'), slot="data")

##########
combined_coga3 <- AddMetaData(
  object = combined_coga3,
  metadata = alc_pheno,
  col.name = "status")
  
  table(combined_coga3@meta.data$status)
  
  combined_coga3 <- AddMetaData(
  object = combined_coga3,
  metadata = sex_pheno,
  col.name = "sex")
  
  combined_coga3 <- AddMetaData(
  object = combined_coga3,
  metadata = age_pheno,
  col.name = "age")
  
  combined_coga3 <- AddMetaData(
  object = combined_coga3,
  metadata = batch_pheno,
  col.name = "batch")
  
  table(combined_coga3@meta.data$batch)

p1 <- DimPlot(combined_coga3, reduction = "umap", group.by = "status")

astro_df_ex <- FindMarkers(combined_coga3, subset.ident = "2", group.by = 'status', test.use = "negbinom", min.pct = 0.1, min.cells.group = 1, logfc.threshold=0,  pseudocount.use= 0.001, max.cells.per.ident= 6000, random.seed = 12, latent.vars = "")



#############################################################################
#Annotate the data with phenotypes
#############################################################################
> library(Seurat)
> library(harmony)
pheno = read.csv ("~/pheno_for_s_nuc_with_cell_id.txt", header=T)
load ("harmony_clusters_seurat.RData")
> rownames (pheno) = pheno$cell_ident
> alc_pheno = pheno[c(10)]
> sex_pheno = pheno[c(8)]
> age_pheno = pheno[c(7)]
> pmi_pheno = pheno[c(13)]
> batch_pheno = pheno[c(18)]

combined_coga3 <- AddMetaData(
  object = combined_coga3,
  metadata = alc_pheno,
  col.name = "status")
  
combined_coga3 <- AddMetaData(
  object = combined_coga3,
  metadata = sex_pheno,
  col.name = "sex")
  
combined_coga3 <- AddMetaData(
  object = combined_coga3,
  metadata = age_pheno,
  col.name = "age")
  
combined_coga3 <- AddMetaData(
  object = combined_coga3,
  metadata = batch_pheno,
  col.name = "batch")
  
combined_coga3 <- AddMetaData(
  object = combined_coga3,
  metadata = pmi_pheno,
  col.name = "pmi")
  
  table(combined_coga3@meta.data$status)
  
  astro_df_ex <- FindMarkers(combined_coga3, ident.1 = "Alcoholic", subset.ident = "2", group.by = 'status', test.use = "negbinom", features = c("SPI1"), logfc.threshold=0,  pseudocount.use= 0.001, max.cells.per.ident= 6000, min.pct = 0.05, random.seed = 12)


p3 <- FeaturePlot(combined_coga3, features = c("LRRC37A", "LRRC37A2"), reduction = "umap", pt.size = 0.1, combine = FALSE, order =T, min.cutoff = 'q50', cols = c("lightgrey", "blue", "darkblue"))