## Load required packages
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(MASS)
library(viridis)
library(dplyr)
library(sctransform)
library(glmGamPoi)
library(stringr)
library(RColorBrewer)
library(patchwork)
library(future)
set.seed(100)

## Creat output directory 
out_dir <- "~/Cd36_OSNs/output/scATAC/seurat/integrated_analysis/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

# Load scATAC object
integrated_scATAC <- readRDS(str_c(out_dir,"integrated_OE_scATAC.rds"))

# call peaks for mature ORNs using MACS2
mature_ORNs <- rownames(integrated_scATAC@meta.data)[which(integrated_scATAC$cell_type=="Mature ORNs")]
mature_ORNs_scATAC <- subset(integrated_scATAC,cells=mature_ORNs)
mature_ORNs_peaks <- CallPeaks(mature_ORNs_scATAC, macs2.path = "/public/home/shipy3/miniconda3/conda_software/bin/macs2",effective.genome.size=2.3e+09,group.by="cell_type",idents="Mature ORNs")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
mature_ORNs_peaks <- keepStandardChromosomes(mature_ORNs_peaks, pruning.mode = "coarse")
mm10_blacklist_bed <- read.table("/public/home/shipy3/DB/mm10/annotation/mm10-blacklist.v2.bed",sep="\t")
mm10_blacklist.gr <- GRanges(mm10_blacklist_bed$V1,
  IRanges(start = mm10_blacklist_bed$V2+1, end = mm10_blacklist_bed$V3))
mature_ORNs_peaks <- subsetByOverlaps(x = mature_ORNs_peaks, ranges = mm10_blacklist.gr, invert = TRUE)

# Load data and create Seurat objects
samples <- c("19d_0818","8w_0720")
dirs <- c("~/Cd36_OSNs/output/scATAC/cellranger/19d_0818/outs/","~/Cd36_OSNs/output/scATAC/cellranger/8w_0720/outs/")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

kept_mature_ORNs.ls <- lapply(1:length(samples),function(i){
  cells <- grep(samples[i],mature_ORNs,value=TRUE)
  cells <- gsub(str_c(samples[i],"_"),"",cells)
  cells
  })

# Create the objects
mature_ORNs_obj.ls <- lapply(1:length(samples),function(i){
  metadata <- read.csv(
    file = str_c(dirs[i],"singlecell.csv"),
    header = TRUE,
    row.names = 1)[-1,] # remove the first row

  # create fragment object for high-quality cells
  frag <- CreateFragmentObject(
    path = str_c(dirs[i],"fragments.tsv.gz"),
    cells = kept_mature_ORNs.ls[[i]]
  )
  
  # create a matrix of peaks (common peaks) x cell
  counts <- FeatureMatrix(
    fragments = frag,
    features = mature_ORNs_peaks,
    cells = kept_mature_ORNs.ls[[i]]
  )

  # Create a 'ChromatinAssay' object from re-quantified matrices
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    annotation = annotation,
    fragments = frag)

  # Create Seurat object
  obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = 'ATAC',
    project = samples[i],
    meta.data = metadata)
  obj
  })

# add information to identify dataset of origin
mature_ORNs_obj.ls[[1]]$dataset <- "19d"
mature_ORNs_obj.ls[[2]]$dataset <- "8w"

# merge all datasets, adding a cell ID to make sure cell names are unique
combined_mature_ORNs_scATAC <- merge(
  x = mature_ORNs_obj.ls[[1]],
  y = list(mature_ORNs_obj.ls[[2]]),
  add.cell.ids = samples
)

# process merged object
combined_mature_ORNs_scATAC <- RunTFIDF(combined_mature_ORNs_scATAC)
combined_mature_ORNs_scATAC <- FindTopFeatures(combined_mature_ORNs_scATAC, min.cutoff = 20)
combined_mature_ORNs_scATAC <- RunSVD(object=combined_mature_ORNs_scATAC)
pdf(str_c(out_dir,"combined_mature_ORNs_scATAC_total_conuts_vs_PC_corr.pdf"),width=9,height=5)
DepthCor(combined_mature_ORNs_scATAC)
dev.off()
combined_mature_ORNs_scATAC <- RunUMAP(combined_mature_ORNs_scATAC, dims = 2:30,reduction = 'lsi')

# split objects
add_prefix_mature_ORNs_obj.ls <- SplitObject(combined_mature_ORNs_scATAC,split.by = "dataset")

# preprocess each object
for (i in 1:length(add_prefix_mature_ORNs_obj.ls)){
  add_prefix_mature_ORNs_obj.ls[[i]] <- RunTFIDF(add_prefix_mature_ORNs_obj.ls[[i]])
  add_prefix_mature_ORNs_obj.ls[[i]] <- FindTopFeatures(add_prefix_mature_ORNs_obj.ls[[i]], min.cutoff = 10)
  add_prefix_mature_ORNs_obj.ls[[i]] <- RunSVD(object=add_prefix_mature_ORNs_obj.ls[[i]])
}


# To find integration anchors across datasets, we need to project them into a shared low-dimensional space. To do this, we’ll use reciprocal LSI projection (projecting each dataset into the others LSI space) by setting reduction="rlsi"
integration.anchors <- FindIntegrationAnchors(
  object.list = list(add_prefix_mature_ORNs_obj.ls[[1]], add_prefix_mature_ORNs_obj.ls[[2]]),
  anchor.features = rownames(add_prefix_mature_ORNs_obj.ls[[1]]),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings,perform dataset integration using a pre-computed Anchorset of specified low dimensional representations
integrated_mature_ORNs_scATAC <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined_mature_ORNs_scATAC[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated_mature_ORNs_scATAC <- RunUMAP(integrated_mature_ORNs_scATAC, reduction = "integrated_lsi", dims = 2:30)

# Create a gene activity matrix
gene.activities <- GeneActivity(integrated_mature_ORNs_scATAC)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
integrated_mature_ORNs_scATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated_mature_ORNs_scATAC <- NormalizeData(
  object = integrated_mature_ORNs_scATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated_mature_ORNs_scATAC$nCount_RNA)
)

# visualize the activities of canonical marker genes
markers <- c("Syt1", #神经元
              "Omp",# Mature ORNs
              "Nqo1","Acsm4",
              "Ncam2","Nfix","Nfib",
              "Cd36","Cd55",
              "Gucy1b2","Trpc2",
              "Calb2","Pcp4l1",
              "Avp","Rprm","Fxyd2",
              "Ighm","Igkc","Cd74")
DefaultAssay(integrated_mature_ORNs_scATAC) <- 'RNA'
pdf(str_c(out_dir,"integrated_mature_ORNs_scATAC_markers_GeneScore_UMAP.pdf"))
FeaturePlot(
  object = integrated_mature_ORNs_scATAC,
  features = markers,
  cols=c("lightgrey","red"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  combine=FALSE
)
dev.off()

DefaultAssay(filtered_mature_ORNs_scATAC) <- 'RNA'
pdf(str_c(out_dir,"integrated_mature_OSNs_scATAC_A_P_markers_GeneScore_UMAP.pdf"))
FeaturePlot(
  object = filtered_mature_ORNs_scATAC,
  features = c("Nrp1","Plxna1","Sema3a"),
  cols=c("lightgrey","red"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  combine=FALSE
)
dev.off()

# tune the resolution
DefaultAssay(integrated_mature_ORNs_scATAC) <- "ATAC"
integrated_mature_ORNs_scATAC <- FindNeighbors(object = integrated_mature_ORNs_scATAC, reduction = 'integrated_lsi', dims = 2:30)
for ( i in seq(0.2,1,0.2)){
  integrated_mature_ORNs_scATAC <- FindClusters(integrated_mature_ORNs_scATAC, resolution = i)
  pdf(str_c(out_dir,"integrated_mature_ORNs_scATAC_PC30_resolution",i,"_UMAP.pdf"))
  p <- DimPlot(integrated_mature_ORNs_scATAC, reduction = "umap",label=TRUE,group.by="seurat_clusters")
  print(p)
  dev.off()
}

# select 0.3 resolution 
integrated_mature_ORNs_scATAC <- FindClusters(integrated_mature_ORNs_scATAC, resolution = 0.3)

# quality metrices for each cluster
DefaultAssay(integrated_mature_ORNs_scATAC) <- "ATAC"
integrated_mature_ORNs_scATAC <- NucleosomeSignal(integrated_mature_ORNs_scATAC)
integrated_mature_ORNs_scATAC <- TSSEnrichment(integrated_mature_ORNs_scATAC)
integrated_mature_ORNs_scATAC$pct_reads_in_peaks <- integrated_mature_ORNs_scATAC$peak_region_fragments / integrated_mature_ORNs_scATAC$passed_filters * 100

pdf(str_c(out_dir,"/md01/shipy3/Projects/mouse_ORs/output/scATAC/two_sample_OE_add/integrated_mature_ORNs_scATAC_cluster_qc_VlnPlot.pdf"))
VlnPlot(integrated_mature_ORNs_scATAC,features=c("passed_filters","TSS.enrichment","nucleosome_signal","pct_reads_in_peaks"),combine=FALSE)
dev.off()

# Find differentially expressed genes between clusters
DefaultAssay(integrated_mature_ORNs_scATAC) <- "RNA"
clusters <- levels(Idents(integrated_mature_ORNs_scATAC))
cluster_markers.ls <- lapply(1:length(clusters),function(i){
  df <- FindMarkers(integrated_mature_ORNs_scATAC, ident.1 =clusters[i],assay="RNA",only.pos=TRUE)
  if ( nrow(df)==0 ){
    return(df)
  } else {
    df <- df %>%
      filter(p_val<0.05) %>%
      arrange(desc(avg_log2FC))
    return(df)
  }
  })

# annotate cell types
new.cluster.ids <- c("Ventral ORNs","Dorsal ORNs",rep("Ventral ORNs",2),"Cd36+ ORNs","Gucy1b2+ cells")
names(new.cluster.ids) <- levels(Idents(integrated_mature_ORNs_scATAC))
integrated_mature_ORNs_scATAC <- RenameIdents(integrated_mature_ORNs_scATAC, new.cluster.ids)
cell_type_df <- data.frame(row.names=names(Idents(integrated_mature_ORNs_scATAC)),cell_type=as.vector(Idents(integrated_mature_ORNs_scATAC)))
integrated_mature_ORNs_scATAC <- AddMetaData(integrated_mature_ORNs_scATAC,cell_type_df)

saveRDS(integrated_mature_ORNs_scATAC,str_c(out_dir,"integrated_mature_ORNs_scATAC.rds"))
