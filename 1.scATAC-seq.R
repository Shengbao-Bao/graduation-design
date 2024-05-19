##1.set env
setwd("sci/repeat/") 
set.seed(2024)
rm(list=ls())
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)  
library(ggplot2)
library(patchwork)


##2.pre-processing workflow
#peak/cell
counts <- Read10X_h5(filename = "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
head(counts)
#meatadata
metadata <- read.csv(
  file = "atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
) 
names(metadata)
head(metadata)[,1:5]

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = 'atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
pbmc
pbmc[['peaks']]
granges(pbmc)
#extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg19"
# add the gene information to the object
Annotation(pbmc) <- annotations


##3.computing QC metrics
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pdf('1.pdf')
DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 3, 'High', 'Low')
pdf('2.pdf')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
dev.off()

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
pdf('3.pdf')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
dev.off()

pdf('4.pdf')
VlnPlot(
  object = pbmc,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
pbmc

##4.Normalization and linear dimensional reduction
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pdf('5.pdf')
DepthCor(pbmc)
dev.off()

##5.Non-linear dimension reduction and clustering
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

pdf('6.pdf')
DimPlot(object = pbmc, label = TRUE) + NoLegend()
dev.off()

##6.Create a gene activity matrix
gene.activities <- GeneActivity(pbmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'

pdf('7.pdf')
FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
dev.off()

#7.Integrating with scRNA-seq data
# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("scRNA/pbmc_10k_v3.rds")
pbmc_rna <- UpdateSeuratObject(pbmc_rna)
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)
pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
pdf('8.pdf')
plot1 + plot2
dev.off()

# replace each label with its most likely prediction
for(i in levels(pbmc)) {
  cells_to_reid <- WhichCells(pbmc, idents = i)
  newid <- names(which.max(table(pbmc$predicted.id[cells_to_reid])))
  Idents(pbmc, cells = cells_to_reid) <- newid
}
View(pbmc@meta.data)
sort(table(pbmc@active.ident))
pbmc = StashIdent(pbmc, "celltype")
CD14monocytes=subset(pbmc, celltype== 'CD14+ Monocytes')
CD4memory=subset(pbmc, celltype== 'CD4 Memory')
CD4native=subset(pbmc,celltype=='CD4 Naive')
CD8native=subset(pbmc,celltype=='CD8 Naive')
save(CD14monocytes,file='CD14monocytes.Rdata')
save(CD4memory,file='CD4memory.Rdata')
save(CD4native,file='CD4native.Rdata')
save(CD8native,file='CD8native.Rdata')

