## scATAC-seq integration
## Integrating multiple chomatin accessibility data
# Reference: https://stuartlab.org/signac/articles/integrate_atac

library(Signac)
library(Seurat)
library(ggplot2)

# load the pre-processed multiome data
pbmc.multi <- readRDS("pbmc_multiome/pbmc_multiomic.rds")

# load the pre-processed atac data
pbmc.atac <- readRDS("pbmc_vignette/pbmc.rds")

# quantify multiome peaks in the scATAC-seq dataset
counts <- FeatureMatrix(
  fragments = Fragments(pbmc.atac),
  features = granges(pbmc.multi),
  cells = colnames(pbmc.atac)
)

# add new assay with multiome peaks
pbmc.atac[['ATAC']] <- CreateChromatinAssay(
  counts = counts,
  fragments = Fragments(pbmc.atac)
)

# compute LSI
DefaultAssay(pbmc.atac) <- "ATAC"
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 10)
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- RunSVD(pbmc.atac)

# first add dataset-identifying metadata
pbmc.atac$dataset <- "ATAC"
pbmc.multi$dataset <- "Multiome"

# merge
pbmc.combined <- merge(pbmc.atac, pbmc.multi)

# process the combined dataset
pbmc.combined <- FindTopFeatures(pbmc.combined, min.cutoff = 10)
pbmc.combined <- RunTFIDF(pbmc.combined)
pbmc.combined <- RunSVD(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(pbmc.combined, group.by = "dataset")

## Integration
# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(pbmc.multi, pbmc.atac),
  anchor.features = rownames(pbmc.multi),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = pbmc.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(integrated, group.by = "dataset")

(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))

## Reference mapping
## For integrating large, high quality datasets

# compute UMAP and store the UMAP model
pbmc.multi <- RunUMAP(pbmc.multi, reduction = "lsi", dims = 2:30, return.model = TRUE)

# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = pbmc.multi,
  query = pbmc.atac,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

# map query onto the reference dataset
pbmc.atac <- MapQuery(
  anchorset = transfer.anchors,
  reference = pbmc.multi,
  query = pbmc.atac,
  refdata = pbmc.multi$predicted.id,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

p1 <- DimPlot(pbmc.multi, reduction = "umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference")
p2 <- DimPlot(pbmc.atac, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Query")

p1 | p2

## RNA Imputation
# predict gene expression values
rna <- TransferData(
  anchorset = transfer.anchors,
  refdata = LayerData(pbmc.multi, assay = "SCT", layer = "data"),
  weight.reduction = pbmc.atac[["lsi"]],
  dims = 2:30
)

# add predicted values as a new assay
pbmc.atac[["predicted"]] <- rna

DefaultAssay(pbmc.atac) <- "predicted"

FeaturePlot(
  object = pbmc.atac,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  reduction = "ref.umap",
  ncol = 3
)

