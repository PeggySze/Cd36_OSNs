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
library(ggpubr)
set.seed(100)

## Creat output directory 
out_dir <- "~/Cd36_OSNs/output/scATAC/Signac/integrated_analysis/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

## Load data and create Seurat objects
samples <- c("19d_0818","8w_0720")
dirs <- c("~/Cd36_OSNs/output/scATAC/cellranger/19d_0818/outs/","~/Cd36_OSNs/output/scATAC/cellranger/8w_0720/outs/")

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

obj.ls <- lapply(1:length(dirs),function(i){
  counts <- Read10X_h5(str_c(dirs[i],"filtered_peak_bc_matrix.h5"))
  metadata <- read.csv(
    file = str_c(dirs[i],"singlecell.csv"),
    header = TRUE,
    row.names = 1)[-1,] # remove the first row
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    annotation = annotation,
    fragments = str_c(dirs[i],"fragments.tsv.gz"))
  obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = 'peaks',
    project = samples[i],
    meta.data = metadata)
  obj
  })

## Quality control
for (i in 1:length(samples)){
  # Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to fragments < 147 bp (nucleosome-free)
  obj.ls[[i]] <- NucleosomeSignal(obj.ls[[i]]) 
  # compute TSS enrichment score per cell
  obj.ls[[i]] <- TSSEnrichment(obj.ls[[i]],fast = FALSE)
  # add blacklist ratio and fraction of reads in peaks
  obj.ls[[i]]$pct_reads_in_peaks <- obj.ls[[i]]$peak_region_fragments / obj.ls[[i]]$passed_filters * 100
  obj.ls[[i]]$blacklist_ratio <- obj.ls[[i]]$blacklist_region_fragments / obj.ls[[i]]$peak_region_fragments
  pdf(str_c(out_dir,samples[i],"_OE_scATAC_qc.pdf"))
  p1 <- VlnPlot(obj.ls[[i]],features ="passed_filters",pt.size = 0)+ NoLegend()
  p2 <- ggplot(data=obj.ls[[i]][["passed_filters"]],aes(x=passed_filters))+geom_density()
  p3 <- VlnPlot(obj.ls[[i]],features ="TSS.enrichment",pt.size = 0)+ NoLegend()
  p4 <- ggplot(data=obj.ls[[i]][["TSS.enrichment"]],aes(x=TSS.enrichment))+geom_density()
  p5 <- VlnPlot(obj.ls[[i]],features ="nucleosome_signal",pt.size = 0)+ NoLegend()
  p6 <- ggplot(data=obj.ls[[i]][["nucleosome_signal"]],aes(x=nucleosome_signal))+geom_density()
  p7 <- VlnPlot(obj.ls[[i]],features ="pct_reads_in_peaks",pt.size = 0)+ NoLegend()
  p8 <- ggplot(data=obj.ls[[i]][["pct_reads_in_peaks"]],aes(x=pct_reads_in_peaks))+geom_density()
  p9 <- TSSPlot(obj.ls[[i]]) + NoLegend()
  p10 <- FragmentHistogram(object=obj.ls[[i]],region = "chr1-1-100000000")
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  print(p7)
  print(p8)
  print(p9)
  print(p10)
  dev.off()
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
for (i in 1:length(samples)){
  df <- obj.ls[[i]]@meta.data
  df$log10_fragments <- log10(df$passed_filters)
  df$ATAC_density <- get_density(df$log10_fragments,df$TSS.enrichment,n=200)
  pdf(str_c(out_dir,samples[i],"_OE_scATAC_qc_02.pdf"),width=8)
  p <- ggplot(df,aes(x=log10_fragments,y=TSS.enrichment,color=ATAC_density))+
    geom_point(size=0.5) +
    theme_bw() +
    labs(x="Log10(Unique fragments)",y="TSS Enrichment")+
    scale_color_viridis() +
    xlim(2.7,5.3) +
    ylim(0,30)+
    #geom_vline(xintercept=c(3.5,4.4),linetype="longdash")+
    #geom_hline(yintercept=5,linetype="longdash")+
    theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),axis.title = element_text(face = "bold",size = rel(1.2)),axis.title.y = element_text(angle=90,vjust =2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(color="black",size=rel(1.2)), 
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        legend.title = element_text(face="italic"))
  print(p)
  dev.off()
}

## Figure S4D
library(gridExtra)
library(grid)
# https://community.rstudio.com/t/common-axis-title-in-grid-arrange/96353/2
d19_df <- obj.ls[[1]]@meta.data
d19_df$log10_fragments <- log10(d19_df$passed_filters)
d19_df$ATAC_density <- get_density(d19_df$log10_fragments,d19_df$TSS.enrichment,n=200)
w8_df <- obj.ls[[2]]@meta.data
w8_df$log10_fragments <- log10(w8_df$passed_filters)
w8_df$ATAC_density <- get_density(w8_df$log10_fragments,w8_df$TSS.enrichment,n=200)
pdf(str_c(out_dir,"combined_OE_scATAC_qc.pdf"),width=13)
p1 <- ggplot(d19_df,aes(x=log10_fragments,y=TSS.enrichment,color=ATAC_density))+
  geom_point(size=0.5) +
  theme_bw() +
  labs(x="",y="",title="19d")+
  scale_color_viridis() +
  xlim(2.7,5.3) +
  ylim(0,30)+
  geom_vline(xintercept=c(log10(5000),log10(40000)),linetype="longdash")+
  geom_hline(yintercept=4,linetype="longdash")+
  theme(plot.title = element_text(face = "bold",
                                         size = rel(2.5), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_blank(),
               axis.text = element_text(color="black",size=rel(2)), 
               axis.line = element_line(colour="black",size=1),
               axis.ticks = element_line(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(0, 0.5, 0, 0.5), "cm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold"),
               legend.position = " none"
          )
p2 <- ggplot(w8_df,aes(x=log10_fragments,y=TSS.enrichment,color=ATAC_density))+
  geom_point(size=0.5) +
  theme_bw() +
  labs(x="",y="",title="8w",color="Density")+
  scale_color_continuous(type = "viridis", breaks=c(min(w8_df$ATAC_density),max(w8_df$ATAC_density)),labels = c("Low", "High"),guide = guide_colorbar(frame.colour = "black"))+
  xlim(2.7,5.3) +
  ylim(0,30)+
  geom_vline(xintercept=c(log10(5000),log10(40000)),linetype="longdash")+
  geom_hline(yintercept=4,linetype="longdash")+
  theme(plot.title = element_text(face = "bold",
                                         size = rel(2.5), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_blank(),
               axis.text = element_text(color="black",size=rel(1.8)), 
               axis.line = element_line(colour="black",size=1),
               axis.ticks = element_line(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = "black",size=1.5, linetype="solid"),
               legend.title = element_text(size=15),
               legend.text=element_text(size=12),
               plot.margin=unit(c(0, 0.5, 0, 0.5), "cm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold"))
left <- textGrob("TSS Enrichment",rot = 90, gp = gpar(fontsize = 20))
bottom <- textGrob("Log10 (# fragments)", gp = gpar(fontsize = 20))
result <- p1 + p2
gt <- patchwork::patchworkGrob(result)
gridExtra::grid.arrange(gt, left = left, bottom = bottom)
dev.off()

## Figure S4C
df <- rbind(obj.ls[[1]]@meta.data, obj.ls[[2]]@meta.data)
df$orig.ident <- factor(df$orig.ident,levels=c("19d_0818","8w_0720"),labels=c("19d","8w"))
df$log10_fragments <- log10(df$passed_filters)
newpalette <- c("#F39B7FFF","#E64B35FF")
pdf(str_c(out_dir,"combined_OE_scATAC_qc_violin.pdf"),height=10)
p1 <- ggplot(data=df,aes(x=orig.ident,y=log10_fragments,fill=orig.ident))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white",outlier.size=0.8)+
  theme_classic()+
  scale_fill_manual(values=newpalette) +
  labs(y="Log10\n(# fragments)") + 
  geom_hline(yintercept=c(log10(5000),log10(40000)),linetype="longdash")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=rel(2),vjust =2),
        axis.text.y = element_text(color="black",size=rel(1.8)), 
        axis.text.x = element_blank(),
        axis.line = element_line(colour="black",size = 1),
        axis.ticks = element_line(),
        plot.margin=unit(c(0, 0.5, 0.5, 0.5), "cm"),
        legend.position = " none")
p2 <- ggplot(data=df,aes(x=orig.ident,y=pct_reads_in_peaks,fill=orig.ident))+
  geom_violin(trim=FALSE,)+
  geom_boxplot(width=0.1, fill="white",outlier.size=0.8)+
  theme_classic()+
  scale_fill_manual(values=newpalette) +
  labs(y="FRiP(%)") + 
  geom_hline(yintercept=40,linetype="longdash")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=rel(2),vjust =2),
        axis.text.y = element_text(color="black",size=rel(1.8)), 
        axis.text.x = element_blank(),
        axis.line = element_line(colour="black",size = 1),
        axis.ticks = element_line(),
        plot.margin=unit(c(0, 0.5, 0.5, 0.5), "cm"),
        legend.position = " none")
p3 <- ggplot(data=df,aes(x=orig.ident,y=nucleosome_signal,fill=orig.ident))+
  geom_violin(trim=FALSE,)+
  geom_boxplot(width=0.1, fill="white",outlier.size=0.8)+
  theme_classic()+
  scale_fill_manual(values=newpalette) +
  labs(y="Nucleosome\nsignal") + 
  geom_hline(yintercept=2,linetype="longdash")+
  coord_cartesian(ylim=c(0,5))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=rel(2),vjust =2),
        axis.text = element_text(color="black",size=rel(1.8)), 
        axis.line = element_line(colour="black",size = 1),
        axis.ticks = element_line(),
        plot.margin=unit(c(0, 0.5, 0.5, 0.5), "cm"),
        legend.position = " none")
p1 / p2 / p3 
dev.off()

# criteria of filtering 
# filtered out cells with less than 5000 or more than 40000 sequencing fragments
# discarded cells with a TSS enrichment less than 4
# subset = 
#  passed_filters >= 5000 &
#  passed_filters <= 40000 &
#  TSS.enrichment >= 4 &
#  nucleosome_signal < 2 &
#  pct_reads_in_peaks >= 40 
for (i in 1:length(samples)){
  obj.ls[[i]] <- subset(obj.ls[[i]],
     subset= passed_filters >= 5000 &
      passed_filters <= 40000 &
      TSS.enrichment >= 4 &
      nucleosome_signal < 2 &
      pct_reads_in_peaks >= 40 
      )
}
kept_cells.ls <- lapply(1:length(samples),function(i){
  colnames(obj.ls[[i]])
  })

## integrate two scATAC-seq datasets
# Creating a common peak set
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
peaks_gr.ls <- lapply(1:length(samples),function(i){
  df <- read.table(file = str_c(dirs[i],"peaks.bed"),col.names = c("chr", "start", "end"))
  gr <- makeGRangesFromDataFrame(df)
  gr
  })
combined.peaks <- reduce(x = c(peaks_gr.ls[[1]], peaks_gr.ls[[2]]))

# Create the objects
obj.ls <- lapply(1:length(samples),function(i){
  metadata <- read.csv(
    file = str_c(dirs[i],"singlecell.csv"),
    header = TRUE,
    row.names = 1)[-1,] # remove the first row

  # create fragment object for high-quality cells
  frag <- CreateFragmentObject(
    path = str_c(dirs[i],"fragments.tsv.gz"),
    cells = kept_cells.ls[[i]]
  )
  
  # create a matrix of peaks (common peaks) x cell
  counts <- FeatureMatrix(
    fragments = frag,
    features = combined.peaks,
    cells = kept_cells.ls[[i]]
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
obj.ls[[1]]$dataset <- "19d"
obj.ls[[2]]$dataset <- "8w"

# merge all datasets, adding a cell ID to make sure cell names are unique
combined_scATAC <- merge(
  x = obj.ls[[1]],
  y = list(obj.ls[[2]]),
  add.cell.ids = samples
)

# compute TSS enrichment score per cell
combined_scATAC <- TSSEnrichment(combined_scATAC, fast = FALSE)


## Figure S4E,S4F
combined_scATAC$dataset <- factor(combined_scATAC$dataset,levels=c("19d","8w"))
newpalette <- c("#F39B7FFF","#E64B35FF")
positionEnrichment <- GetAssayData(object = combined_scATAC, assay = "ATAC", slot = "positionEnrichment")
enrichment.matrix <- positionEnrichment[["TSS"]]
enrichment.matrix <- as.data.frame(as.matrix(enrichment.matrix[1:(nrow(x =enrichment.matrix)- 2), ]))
TSS_profile_df <- c()
for (i in 1:length(samples)){
  dataset_enrichment.matrix <- enrichment.matrix[grepl(samples[i],rownames(enrichment.matrix)),]
  TSS_profile_df <- data.frame(pos=names(colMeans(dataset_enrichment.matrix)),mean=colMeans(dataset_enrichment.matrix),dataset=samples[i]) %>% rbind(TSS_profile_df,.)
}
TSS_profile_df$dataset <- factor(TSS_profile_df$dataset,levels=samples,labels=c("19d","8w"))
TSS_profile_df$pos <- as.numeric(TSS_profile_df$pos)
pdf(str_c(out_dir,"combined_OE_scATAC_TSS_profile&fragment_size_distribution.pdf"),width=13)
p1 <- ggplot(data=TSS_profile_df,aes(x=pos,y=mean,color=dataset))+
  geom_line(stat = "identity", size = 0.2) +
  labs(x="Distance from TSS (bp)",y="Mean TSS enrichment score") +
  scale_color_manual(values=newpalette) +
  theme_classic()+
  theme(axis.title.x = element_text(size=rel(2)),
        axis.title.y = element_text(size=rel(2),vjust =2),
        axis.text = element_text(color="black",size=rel(1.8)), 
        axis.line = element_line(colour="black",size = 1),
        axis.ticks = element_line(),
        plot.margin=unit(c(0, 0.5, 0.5, 0.5), "cm"),
        legend.position = " none")
p2 <- FragmentHistogram(object=combined_scATAC,region = "chr1-1-100000000",group.by="dataset")
p2$data$group <- factor(p2$data$group,levels=c("19d","8w"))
p3 <- ggplot(data=p2$data,aes(x=length,color=group))+
  geom_density(size=1) +
  labs(x="Fragment length (bp)",y="Density") +
  scale_color_manual(values=newpalette) +
  coord_cartesian(xlim=c(0,1000))+
  theme_classic()+
  theme(axis.title.x = element_text(size=rel(2)),
        axis.title.y = element_text(size=rel(2),vjust =2),
        axis.text = element_text(color="black",size=rel(1.8)), 
        axis.line = element_line(colour="black",size = 1),
        axis.ticks = element_line(),
        plot.margin=unit(c(0, 0.5, 0.5, 0.5), "cm"),
        legend.position = " none")
p1 + p3 
dev.off()


# process merged object
combined_scATAC <- RunTFIDF(combined_scATAC)
combined_scATAC <- FindTopFeatures(combined_scATAC, min.cutoff = 50)
combined_scATAC <- RunSVD(object=combined_scATAC)

# assess the correlation between each LSI component and sequencing depth using the DepthCor() function, because The first LSI component often captures sequencing depth (technical variation) rather than biological variation.  If this is the case, the component should be removed from downstream analysis
pdf(str_c(out_dir,"combined_scATAC_OE_total_conuts_vs_PC_corr.pdf"),width=9,height=5)
DepthCor(combined_scATAC)
dev.off()
# there is a very strong correlation between the first LSI component and the total number of counts for the cell, so we will perform downstream steps without this component

combined_scATAC <- RunUMAP(combined_scATAC, dims = 2:30,reduction = 'lsi')

# split objects
add_prefix_obj.ls <- SplitObject(combined_scATAC,split.by = "dataset")

# preprocess each object
for (i in 1:length(add_prefix_obj.ls)){
  add_prefix_obj.ls[[i]] <- RunTFIDF(add_prefix_obj.ls[[i]])
  add_prefix_obj.ls[[i]] <- FindTopFeatures(add_prefix_obj.ls[[i]], min.cutoff = 20)
  add_prefix_obj.ls[[i]] <- RunSVD(object=add_prefix_obj.ls[[i]])
}

# To find integration anchors across datasets, we need to project them into a shared low-dimensional space. To do this, we’ll use reciprocal LSI projection (projecting each dataset into the others LSI space) by setting reduction="rlsi"
integration.anchors <- FindIntegrationAnchors(
  object.list = list(add_prefix_obj.ls[[1]], add_prefix_obj.ls[[2]]),
  anchor.features = rownames(add_prefix_obj.ls[[1]]),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings,perform dataset integration using a pre-computed Anchorset of specified low dimensional representations
integrated_scATAC <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined_scATAC[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated_scATAC <- RunUMAP(integrated_scATAC, reduction = "integrated_lsi", dims = 2:30)

# combined vs integrated
p1 <- DimPlot(combined_scATAC, group.by = "dataset")
p2 <- DimPlot(integrated_scATAC, group.by = "dataset")
pdf(str_c(out_dir,"combined_vs_integrated_OE_scATAC_UMAP.pdf"),width=15)
(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
dev.off()

# Create a gene activity matrix
# To create a gene activity matrix, we extract gene coordinates and extend them to include the 2 kb upstream region (as promoter accessibility is often correlated with gene expression). We then count the number of fragments for each cell that map to each of these regions, using the using the FeatureMatrix() function
gene.activities <- GeneActivity(integrated_scATAC)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
integrated_scATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated_scATAC <- NormalizeData(
  object = integrated_scATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated_scATAC$nCount_RNA)
)

# visualize the activities of canonical marker genes
markers <- c("Syt1", #神经元
              "Omp",# Mature ORNs
              "Nqo1","Ncam2","Cd36",
              "Gap43","Gng8",#Immature ORNs
              "Sox11","Neurod1","Neurog1",#INPs,immediate neuronal precursors 
              "Ascl1","Kit" ,#GBCs
              "Cebpd","Krt5","Trp63","Krt14",#HBCs
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#支持细胞  Sustentacular cells
              "Atp1a2","Fabp7", #Ensheathing glia
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Ascl3","Cftr",#Microvillar cells
              "Krt18","Trpm5", #Brush cells
              "Col1a1","Bglap",#Osteogenic cells 成骨细胞
              "Eng","Sox17",#Pericytes
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#中性粒细胞 Neutrophils
              "Hmgb2","Top2a",#Late activated neural stem cells
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Mcpt8","Ccl4","Ptprc","Itga2","Fcer1a", #Basophils
              "C1qa","Ms4a7"#Macrophages
)
DefaultAssay(integrated_scATAC) <- 'RNA'
pdf(str_c(out_dir,"integrated_OE_scATAC_markers_GeneScore_UMAP.pdf"))
FeaturePlot(
  object = integrated_scATAC,
  features = markers,
  cols=c("lightgrey","red"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  combine=FALSE
)
dev.off()

# tune the resolution
DefaultAssay(integrated_scATAC) <- "ATAC"
integrated_scATAC <- FindNeighbors(object = integrated_scATAC, reduction = 'integrated_lsi', dims = 2:30)
for ( i in seq(2.2,3,0.2)){
  integrated_scATAC <- FindClusters(integrated_scATAC, resolution = i)
  pdf(str_c(out_dir,"integrated_OE_scATAC_PC30_resolution",i,"_UMAP.pdf"))
  p <- DimPlot(integrated_scATAC, reduction = "umap",label=TRUE,group.by="seurat_clusters")
  print(p)
  dev.off()
}

# select 0.8 resolution
integrated_scATAC <- FindClusters(integrated_scATAC, resolution = 0.8)

# To help interpret the scATAC-seq data, we can classify cells based on an scRNA-seq experiment from the same biological system
# Load the pre-processed scRNA-seq data
integrated_scRNA <- readRDS("~/Cd36_OSNs/output/scRNA/seurat/integrated_analysis/integrated_OE_scRNA.rds")
Idents(integrated_scRNA) <- integrated_scRNA$cell_type
DefaultAssay(integrated_scRNA) <- "integrated"

Idents(integrated_scATAC) <- integrated_scATAC$ATAC_snn_res.0.8
DefaultAssay(integrated_scATAC) <- "ATAC"
transfer.anchors <- FindTransferAnchors(
  reference = integrated_scRNA,
  query = integrated_scATAC,
  features = VariableFeatures(object = integrated_scRNA), 
  reference.assay = "integrated",
  query.assay = "RNA", 
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = Idents(integrated_scRNA),
  weight.reduction = integrated_scATAC[['integrated_lsi']],
  dims = 2:30
)
predicted.labels$cluster <- as.vector(Idents(integrated_scATAC))
predicted.labels$cluster <- factor(predicted.labels$cluster,levels=levels(Idents(integrated_scATAC)))
table(predicted.labels$cluster,predicted.labels$predicted.id)

pdf(str_c(out_dir,"integrated_OE_prediction_scores_distribution.pdf"))
hist(predicted.labels$prediction.score.max)
dev.off()

# Find differentially accessible peaks between clusters
DefaultAssay(integrated_scATAC) <- "ATAC"
clusters <- levels(Idents(integrated_scATAC))
da_peaks.ls <- lapply(1:length(clusters),function(i){
  df <- FindMarkers(
    object = integrated_scATAC,
    ident.1 = clusters[i],
    min.pct = 0.05,
    test.use = 'LR',
    only.pos = TRUE,
    latent.vars = 'nCount_ATAC')
  if ( nrow(df)==0 ){
    return(df)
  } else {
    df <- df %>%
      filter(p_val<0.05) %>%
      arrange(desc(avg_log2FC))
    ann_df <- ClosestFeature(object=integrated_scATAC,
   regions=rownames(df))
    return(ann_df)
  }
  })

# Find differentially expressed genes between clusters
DefaultAssay(integrated_scATAC) <- "RNA"
clusters <- levels(Idents(integrated_scATAC))
cluster_markers.ls <- lapply(1:length(clusters),function(i){
  df <- FindMarkers(integrated_scATAC, ident.1 =clusters[i],assay="RNA",only.pos=TRUE)
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
# C0,1,2,3,4,5,7,9,21 : Mature ORNs  C21: Gucy1b2+
# C6 : Immature ORNs,INPs,Mature ORNs
# C8 : Sustentacular cells 
# C11,15: Neutrophils
# C12,13 : B cells
# C14 : INPs
# C16 : Brush cells
# C17 : Periglomerular cells
# C18 : Basophils
# C19 : Monocytes
# C20 : Erythrocytes

new.cluster.ids <- c(rep("HBC_lineage_cells",8),"Sustentacular cells",rep("HBC_lineage_cells",2),"Neutrophils","B cells","B cells","HBC_lineage_cells","Neutrophils","Brush cells","Periglomerular cells","Basophils","Monocytes","Erythrocytes","HBC_lineage_cells")
names(new.cluster.ids) <- levels(integrated_scATAC)
integrated_scATAC <- RenameIdents(integrated_scATAC, new.cluster.ids)
cell_type_df <- data.frame(row.names=names(Idents(integrated_scATAC)),rough_annotation=as.vector(Idents(integrated_scATAC)))
integrated_scATAC <- AddMetaData(integrated_scATAC,cell_type_df)

## subset HBC lineage cells for further annotation
HBC_lineage_cells <- rownames(integrated_scATAC@meta.data)[which(integrated_scATAC$rough_annotation =="HBC_lineage_cells")]

# retain HBC lineage cells
HBC_lineage_obj.ls <- lapply(1:length(obj.ls),function(i){
  cells <- grep(samples[i],HBC_lineage_cells,value=TRUE)
  cells <- gsub(str_c(samples[i],"_"),"",cells)
  HBC_lineage_obj <- subset(obj.ls[[i]],cells=cells)
  HBC_lineage_obj
})

# merge datasets, adding a cell ID to make sure cell names are unique
combined_HBC_lineage_scATAC <- merge(
  x = HBC_lineage_obj.ls[[1]],
  y = list(HBC_lineage_obj.ls[[2]]),
  add.cell.ids = samples
)

# process merged object
combined_HBC_lineage_scATAC <- RunTFIDF(combined_HBC_lineage_scATAC)
combined_HBC_lineage_scATAC <- FindTopFeatures(combined_HBC_lineage_scATAC, min.cutoff = 20)
combined_HBC_lineage_scATAC <- RunSVD(object=combined_HBC_lineage_scATAC)
pdf(str_c(out_dir,"combined_HBC_lineage_scATAC_total_conuts_vs_PC_corr.pdf"),width=9,height=5)
DepthCor(combined_HBC_lineage_scATAC)
dev.off()
combined_HBC_lineage_scATAC <- RunUMAP(combined_HBC_lineage_scATAC, dims = 2:30,reduction = 'lsi')

# split objects
add_prefix_HBC_lineage_obj.ls <- SplitObject(combined_HBC_lineage_scATAC,split.by = "dataset")

# preprocess each object
for (i in 1:length(add_prefix_HBC_lineage_obj.ls)){
  add_prefix_HBC_lineage_obj.ls[[i]] <- RunTFIDF(add_prefix_HBC_lineage_obj.ls[[i]])
  add_prefix_HBC_lineage_obj.ls[[i]] <- FindTopFeatures(add_prefix_HBC_lineage_obj.ls[[i]], min.cutoff = 10)
  add_prefix_HBC_lineage_obj.ls[[i]] <- RunSVD(object=add_prefix_HBC_lineage_obj.ls[[i]])
}

# integrate scATAC of HBC lineage cells
HBC_lineage_integration.anchors <- FindIntegrationAnchors(
  object.list = list(add_prefix_HBC_lineage_obj.ls[[1]], add_prefix_HBC_lineage_obj.ls[[2]]),
  anchor.features = rownames(add_prefix_HBC_lineage_obj.ls[[1]]),
  reduction = "rlsi",
  dims = 2:30
)

integrated_HBC_lineage_scATAC <- IntegrateEmbeddings(
  anchorset = HBC_lineage_integration.anchors,
  reductions = combined_HBC_lineage_scATAC[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated_HBC_lineage_scATAC <- RunUMAP(integrated_HBC_lineage_scATAC, reduction = "integrated_lsi", dims = 2:30)

# Create a gene activity matrix
gene.activities <- GeneActivity(integrated_HBC_lineage_scATAC)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
integrated_HBC_lineage_scATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated_HBC_lineage_scATAC <- NormalizeData(
  object = integrated_HBC_lineage_scATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated_HBC_lineage_scATAC$nCount_RNA)
)

# visualize the activities of canonical marker genes
markers <- c("Syt1", #神经元
              "Omp",# Mature ORNs
              "Nqo1","Ncam2","Cd36",
              "Gap43","Gng8",#Immature ORNs
              "Sox11","Neurod1","Neurog1",#INPs,immediate neuronal precursors 
              "Ascl1","Kit" ,#GBCs
              "Cebpd","Krt5","Trp63","Krt14",#HBCs
              "Gucy1b2",
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#支持细胞  Sustentacular cells
              "Atp1a2","Fabp7", #Ensheathing glia
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Ascl3","Cftr",#Microvillar cells
              "Krt18","Trpm5", #Brush cells
              "Col1a1","Bglap",#Osteogenic cells 成骨细胞
              "Eng","Sox17",#Pericytes
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#中性粒细胞 Neutrophils
              "Hmgb2","Top2a",#Late activated neural stem cells
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Mcpt8","Ccl4", #Basophils
              "C1qa","Ms4a7","Ptprc","Itga2","Fcer1a"#Macrophages
)
DefaultAssay(integrated_HBC_lineage_scATAC) <- 'RNA'
pdf(str_c(out_dir,"integrated_HBC_lineage_scATAC_markers_GeneScore_UMAP.pdf"))
FeaturePlot(
  object = integrated_HBC_lineage_scATAC,
  features = markers,
  cols=c("lightgrey","red"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  combine=FALSE
)
dev.off()

# highlight cells in UMAP
highlight_clusters <- c(10,14,6,21)
pdf(str_c(out_dir,"integrated_HBC_lineage_scATAC_highlight_cells_UMAP.pdf"))
for (i in 1:length(highlight_clusters)){
  cells <- rownames(integrated_scATAC@meta.data)[which(integrated_scATAC$ATAC_snn_res.0.8 == highlight_clusters[i])]
  p <- DimPlot(integrated_HBC_lineage_scATAC, reduction = "umap",label=FALSE,cells.highlight=cells) + labs(title=str_c("C",highlight_clusters[i]))
  print(p)
}
dev.off()


# tune the resolution
DefaultAssay(integrated_HBC_lineage_scATAC) <- "ATAC"
integrated_HBC_lineage_scATAC <- FindNeighbors(object = integrated_HBC_lineage_scATAC, reduction = 'integrated_lsi', dims = 2:30)
for ( i in seq(1.2,2,0.2)){
  integrated_HBC_lineage_scATAC <- FindClusters(integrated_HBC_lineage_scATAC, resolution = i)
  pdf(str_c(out_dir,"integrated_HBC_lineage_scATAC_PC30_resolution",i,"_UMAP.pdf"))
  p <- DimPlot(integrated_HBC_lineage_scATAC, reduction = "umap",label=TRUE)
  print(p)
  dev.off()
}

# Select 1.4 resolution 
integrated_HBC_lineage_scATAC <- FindClusters(integrated_HBC_lineage_scATAC, resolution = 1.4)

# classify cells based on an scRNA-seq experiment from the same biological system
transfer.anchors <- FindTransferAnchors(
  reference = integrated_scRNA,
  query = integrated_HBC_lineage_scATAC,
  features = VariableFeatures(object = integrated_scRNA), 
  reference.assay = "integrated",
  query.assay = "RNA", 
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = Idents(integrated_scRNA),
  weight.reduction = integrated_HBC_lineage_scATAC[['integrated_lsi']],
  dims = 2:30
)
predicted.labels$cluster <- as.vector(Idents(integrated_HBC_lineage_scATAC))
predicted.labels$cluster <- factor(predicted.labels$cluster,levels=levels(Idents(integrated_HBC_lineage_scATAC)))
table(predicted.labels$cluster,predicted.labels$predicted.id)

colnames(table(predicted.labels$cluster,predicted.labels$predicted.id))[apply(table(predicted.labels$cluster,predicted.labels$predicted.id),1,which.max)]

# annotate cell types
Idents(integrated_HBC_lineage_scATAC) <- integrated_HBC_lineage_scATAC$ATAC_snn_res.1.4
new.cluster.ids <- c(rep("Mature ORNs",14),"Immature ORNs","INPs","GBCs","Neutrophils","Ensheathing glia","Bowman's gland","HBCs","Mature ORNs","Immature ORNs")
names(new.cluster.ids) <- levels(integrated_HBC_lineage_scATAC)
integrated_HBC_lineage_scATAC <- RenameIdents(integrated_HBC_lineage_scATAC, new.cluster.ids)
cell_type_df <- data.frame(row.names=names(Idents(integrated_HBC_lineage_scATAC)),cell_type=as.vector(Idents(integrated_HBC_lineage_scATAC)))
integrated_HBC_lineage_scATAC <- AddMetaData(integrated_HBC_lineage_scATAC,cell_type_df)

integrated_scATAC$cell_type <- integrated_scATAC$rough_annotation
reannotated_cell_types <- unique(integrated_HBC_lineage_scATAC$cell_type)
for (i in 1:length(reannotated_cell_types)){
  cells <- rownames(integrated_HBC_lineage_scATAC@meta.data)[which(integrated_HBC_lineage_scATAC$cell_type==reannotated_cell_types[i])]
  integrated_scATAC$cell_type[which(rownames(integrated_scATAC@meta.data) %in% cells)] <- reannotated_cell_types[i]
}


## Figure S4A - UMAP with cell type labels
integrated_scATAC$cell_type <- gsub("ORNs","OSNs",integrated_scATAC$cell_type)
Idents(integrated_scATAC) <- integrated_scATAC$cell_type
Idents(integrated_scATAC) <- factor(Idents(integrated_scATAC),levels=c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs","Sustentacular cells","Ensheathing glia","Bowman's gland","Periglomerular cells","Brush cells","B cells","Erythrocytes","Basophils","Neutrophils","Monocytes"))
library(colorspace)
scRNA_newpalette <- c(brewer.pal(9,"YlGn")[c(9,7,5,4,2)],brewer.pal(12,"Set3")[8],brewer.pal(9,"RdPu")[c(2,4,5,7)],brewer.pal(9,"Blues")[c(2,6,8)],brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7])
scATAC_newpalette <- darken(scRNA_newpalette, 0.2)
cell_number_df <- integrated_scATAC@meta.data[,"cell_type",drop=FALSE] %>%
  group_by(cell_type) %>%
  summarise(counts=n())
cell_number_df$log10_counts <- log10(cell_number_df$counts)
cell_number_df$cell_type <- factor(cell_number_df$cell_type,levels=rev(c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs","Sustentacular cells","Ensheathing glia","Bowman's gland","Periglomerular cells","Brush cells","B cells","Erythrocytes","Basophils","Neutrophils","Monocytes")))
cell_type_df <- cell_number_df[,"cell_type",drop=FALSE]
cell_type_df$x <- 1
pdf("/md01/shipy3/Projects/mouse_ORs/output/scATAC/two_sample_OE_add/integrated_OE_scATAC_UMAP_with_cell_types.pdf",width=9.5)
p1 <- DimPlot(integrated_scATAC, reduction = "umap",label=TRUE,cols=scATAC_newpalette,repel=TRUE,label.size=5)+
  labs(title="scATAC-seq")+
  #geom_segment(aes(x = -14.641, y = -16.31245, xend = -12.2, yend = -16.31245),arrow = arrow(length = unit(0.3, "cm")))+
  #geom_segment(aes(x = -14.641, y = -16.31245, xend = -14.641, yend = -13.31245),arrow = arrow(length = unit(0.3, "cm")))+
  #annotate(geom = "text", x = -13.5, y = -17, label = "UMAP_1", color = "black") +
  #annotate(geom = "text", x = -15.2, y = -15, label = "UMAP_2", color = "black",angle = 90) +
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_text(size=rel(1)),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0.3), "cm"))
p2 <- ggplot(data=cell_number_df,aes(x=log10_counts,y=cell_type,fill=cell_type))+
  geom_bar(stat = "identity",width=0.8)+
  labs(x="Log10 (# cells)")+
  #scale_x_continuous(limits=c(0,11100),breaks = c(1000,11000))+
  #scale_y_discrete(position = "right")+
  scale_fill_manual(values=rev(scATAC_newpalette))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0), "cm"),axis.text.x=element_text(size=rel(1),color="black"),axis.title.x=element_text(size=rel(1)))
p3 <- ggplot(data=cell_type_df,aes(x=x,y=cell_type))+
  geom_tile(aes(fill=cell_type),color="black")+
  scale_y_discrete(position = "right",expand=c(0,0))+
  scale_fill_manual(values=rev(scATAC_newpalette))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0.2), "cm"),axis.text.y=element_text(size=rel(1.5),color="black"),axis.title.y=element_blank())
p1 + p2 + p3 + plot_layout(widths=c(8,1.5,0.2))
dev.off()

## marker peaks in HBC lineage
DefaultAssay(integrated_scATAC) <- "ATAC"

markers <- c("Krt5","Trp63","Krt14",#HBCs
  "Ascl1","Kit" ,#GBCs
  "Sox11","Neurod1","Neurog1",#INPs
  "Gap43","Gng8", #Immature ORNs
  "Omp" # Mature ORNs
  )
pdf(str_c(out_dir,"scATAC_HBC_lineage_cell_type_markers_CoveragePlot.pdf"))
for (i in 1:length(markers)){
  p <- CoveragePlot(
    object = integrated_scATAC,
    region = markers[i],
    extend.upstream = 2000,
    extend.downstream = 2000,
    peaks=FALSE
  )
  print(p)
}
dev.off()


Krt5_cov_plot <- CoveragePlot(object = integrated_scATAC,region="chr15-101711800-101715000",peaks=FALSE,annotation = FALSE)
Krt5_cov_df <- Krt5_cov_plot$data %>%
  filter(group %in% c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs"))
Krt5_cov_plot <- ggplot(data = Krt5_cov_df,mapping = aes(x = position, y = coverage, fill = group)) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "left", ncol = 1) +
    ylim(c(0, signif(max(Krt5_cov_df$coverage, na.rm = TRUE), digits = 2))) +
    theme_bw() + 
    labs(y="Normalized accessibility",title="Krt5") +
    #theme_browser(legend = FALSE) +
    #theme_classic()+
    theme(plot.title=element_text(face = "bold",size = rel(1),hjust = 0.5),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",strip.text.y.left = element_text(angle = 0,color="black"), strip.background = element_blank(),axis.text.y = element_blank(),panel.spacing.y = unit(x = 0, units = "line"),plot.margin=unit(c(0,0,0,0), "cm"),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks=element_blank()) +
    scale_fill_manual(values = scATAC_newpalette[1:6])
Krt5_gene_plot <- AnnotationPlot(object = integrated_scATAC,region = "chr15-101711800-101715000")
Krt5_gene_plot <- Krt5_gene_plot + theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid"),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),plot.margin=unit(c(0.2,0.03,0,0), "cm"))

# Krt5 : chr15-101711800-101715000
# Kit : chr5:75574200-75575500
# Neurod1 : chr2:79456000-79457300
# Gap43 : chr16:42274500-42275800
# Omp : chr7-98144500-98146500

regions <- c(
  "chr5-75574200-75575500",
  "chr2-79456000-79457300",
  "chr16-42277000-42279100",
  "chr7-98144500-98146500")
names <- c("Kit","Neurod1","Gap43","Omp")
cov_plot.ls <- lapply(1:length(regions),function(i){
  cov_plot <- CoveragePlot(object = integrated_scATAC,region=regions[i],peaks=FALSE,annotation = FALSE)
  cov_df <- cov_plot$data %>%
    filter(group %in% c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs"))
  cov_plot <- ggplot(data = cov_df,mapping = aes(x = position, y = coverage, fill = group)) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "left", ncol = 1) +
    ylim(c(0, signif(max(cov_df$coverage, na.rm = TRUE), digits = 2))) +
    theme_bw() + 
    labs(title=names[i]) +
    theme(plot.title=element_text(face = "bold",size = rel(1),hjust = 0.5),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",strip.text.y.left = element_blank(), strip.background = element_blank(),axis.text.y = element_blank(),panel.spacing.y = unit(x = 0, units = "line"),plot.margin=unit(c(0, 0,0,0), "cm"),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks=element_blank()) +
    scale_fill_manual(values = scATAC_newpalette[1:6])
  cov_plot
  })

gene_plot.ls <- lapply(1:length(regions),function(i){
  gene_plot <- AnnotationPlot(object = integrated_scATAC,region = regions[i])
  gene_plot <- gene_plot + scale_y_discrete(position = "right")+ theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid"),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y=element_blank(),plot.margin=unit(c(0.2,0.03,0,0), "cm"))
  gene_plot
  })

pdf(str_c(out_dir,"integrated_OE_scATAC_cell_type_markers_combined_CoveragePlot.pdf"),width=4.5,height=3.5)
((Krt5_cov_plot + Krt5_gene_plot + plot_layout(heights=c(3,0.6))) | (cov_plot.ls[[1]] + gene_plot.ls[[1]] + plot_layout(heights=c(3,0.6))) | (cov_plot.ls[[2]] + gene_plot.ls[[2]] + plot_layout(heights=c(3,0.6))) | (cov_plot.ls[[3]] + gene_plot.ls[[3]] + plot_layout(heights=c(3,0.6))) | (cov_plot.ls[[4]] + gene_plot.ls[[4]] + plot_layout(heights=c(3,0.6))))  +  plot_layout(widths=c(1.1,1,1,1,1))
dev.off()

## Figure S4G - marker peaks of all cell types
markers <- c("Syt1", #神经元
              "Omp",# Mature ORNs
              "Nqo1","Ncam2",
              "Gap43","Gng8",#Immature ORNs
              "Sox11","Neurod1","Neurog1",#INPs,immediate neuronal precursors 
              "Ascl1","Kit" ,#GBCs
              "Krt5","Trp63","Krt14",#HBCs
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#支持细胞  Sustentacular cells
              "Atp1a2","Fabp7", #Ensheathing glia
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #球周细胞Periglomerular cells
              "Ascl3","Cftr",#Microvillar cells
              "Krt18","Trpm5", #Brush cells
              "Col1a1","Bglap",#Osteogenic cells
              "Eng","Sox17",#Pericytes
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#中性粒细胞 Neutrophils
              "Hmgb2","Top2a",#Late activated neural stem cells
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
              "Mcpt8","Ccl4","Itga2", #Basophils
              "C1qa","Ms4a7"#Macrophages
)
pdf(str_c(out_dir,"integrated_OE_scATAC_cell_type_markers_CoveragePlot.pdf"))
for (i in 1:length(markers)){
  p <- CoveragePlot(
    object = integrated_scATAC,
    region = markers[i],
    extend.upstream = 1000,
    extend.downstream = 1000,
    peaks=FALSE
  )
  print(p)
}
dev.off()


scRNA_newpalette <- c(brewer.pal(9,"YlGn")[c(9,7,5,4,2)],brewer.pal(12,"Set3")[8],brewer.pal(9,"RdPu")[c(2,4,5,7)],brewer.pal(9,"Blues")[c(2,6,8)],brewer.pal(9,"BuPu")[5],brewer.pal(9,"BuPu")[7])
scATAC_newpalette <- darken(scRNA_newpalette, 0.2)

Krt5_cov_plot <- CoveragePlot(object = integrated_scATAC,region="chr15-101711800-101715000",peaks=FALSE,annotation = FALSE)
Krt5_cov_df <- Krt5_cov_plot$data
Krt5_cov_plot <- ggplot(data = Krt5_cov_df,mapping = aes(x = position, y = coverage, fill = group)) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "left", ncol = 1) +
    ylim(c(0, signif(max(Krt5_cov_df$coverage, na.rm = TRUE), digits = 2))) +
    theme_bw() + 
    labs(y="Normalized accessibility",title="Krt5") +
    #theme_browser(legend = FALSE) +
    #theme_classic()+
    theme(plot.title=element_text(face = "bold",size = rel(1.5),hjust = 0.5),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",strip.text.y.left = element_text(angle = 0,size = rel(1.5)), strip.background = element_rect(colour="black",fill="#FFFFFF",linetype="solid"),axis.text.y = element_blank(),panel.spacing.y = unit(x = 0, units = "line"),plot.margin=unit(c(0,0,0,0), "cm"),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks=element_blank()) +
    scale_fill_manual(values = scATAC_newpalette)
Krt5_gene_plot <- AnnotationPlot(object = integrated_scATAC,region = "chr15-101712400-101714500")
Krt5_gene_plot <- Krt5_gene_plot + theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid"),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),plot.margin=unit(c(0.2,0.03,0,0), "cm"))


regions <- c(
  "chr5-75574200-75575500",# Kit
  "chr2-79456000-79457300",# Neurod1
  "chr16-42277000-42279100",# Gap43
  "chr7-98144500-98146500",# Omp
  "chr2-58052364-58053864",# Ermn
  "chr1-172297064-172298564",# Atp1a2
  "chr15-79164000-79166000",# Sox10
  "chr8-110167000-110168500",# Calb2
  "chr7-143093500-143094250",# Trpm5
  "chr7-24897000-24897700",# Cd79a
  "chr7-103827500-103828300",# Hbb−bs
  #"chr14-56085000-56085500",# Mcpt8 
  "chr13-114844000-114848000", # Itga2
  "chr3-90695400-90696200",# S100a9
  "chr3-90603300-90604000"# S100a4
  )
names <- c("Kit","Neurod1","Gap43","Omp","Ermn","Atp1a2","Sox10","Calb2","Trpm5","Cd79a","Hbb-bs","Itga2","S100a9","S100a4")
cov_plot.ls <- lapply(1:length(regions),function(i){
  cov_plot <- CoveragePlot(object = integrated_scATAC,region=regions[i],peaks=FALSE,annotation = FALSE)
  cov_df <- cov_plot$data 
  cov_plot <- ggplot(data = cov_df,mapping = aes(x = position, y = coverage, fill = group)) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    facet_wrap(facets = ~group, strip.position = "left", ncol = 1) +
    ylim(c(0, signif(max(cov_df$coverage, na.rm = TRUE), digits = 2))) +
    theme_bw() + 
    labs(title=names[i]) +
    theme(plot.title=element_text(face = "bold",size = rel(1.5),hjust = 0.5),panel.background=element_rect(fill='transparent', color='black',linetype="solid"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",strip.text.y.left = element_blank(), strip.background = element_blank(),axis.text.y = element_blank(),panel.spacing.y = unit(x = 0, units = "line"),plot.margin=unit(c(0, 0,0,0), "cm"),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks=element_blank()) +
    scale_fill_manual(values = scATAC_newpalette)
  cov_plot
  })

gene_plot.ls <- lapply(1:length(regions),function(i){
  gene_plot <- AnnotationPlot(object = integrated_scATAC,region = regions[i])
  gene_plot <- gene_plot + theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid"),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y=element_blank(),plot.margin=unit(c(0.2,0.03,0,0), "cm"))
  gene_plot
  })

pdf(str_c(out_dir,"integrated_OE_scATAC_all_cell_type_markers_combined_CoveragePlot.pdf"),width=18,height=6.5)
((Krt5_cov_plot + Krt5_gene_plot + plot_layout(heights=c(7,0.6))) | 
(cov_plot.ls[[1]] + gene_plot.ls[[1]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[2]] + gene_plot.ls[[2]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[3]] + gene_plot.ls[[3]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[4]] + gene_plot.ls[[4]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[5]] + gene_plot.ls[[5]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[6]] + gene_plot.ls[[6]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[7]] + gene_plot.ls[[7]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[8]] + gene_plot.ls[[8]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[9]] + gene_plot.ls[[9]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[10]] + gene_plot.ls[[10]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[11]] + gene_plot.ls[[11]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[12]] + gene_plot.ls[[12]] + plot_layout(heights=c(7,0.6))) |
(cov_plot.ls[[13]] + gene_plot.ls[[13]] + plot_layout(heights=c(7,0.6))) | 
(cov_plot.ls[[14]] + gene_plot.ls[[14]] + plot_layout(heights=c(7,0.6)))) +  plot_layout(widths=c(1.1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
dev.off()

## Figure S4B - UMAP split by sample
newpalette <- c("#F39B7FFF","#E64B35FF")
pdf(str_c(out_dir,"integrated_OE_scATAC_split_by_sample_UMAP.pdf"),width=7,height=12)
DimPlot(integrated_scATAC, reduction = "umap",label=FALSE,cols=newpalette,group.by="dataset",split.by="dataset",ncol=1)+labs(title="scATAC-seq")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"))
dev.off()

saveRDS(integrated_scATAC,str_c(out_dir,"integrated_OE_scATAC.rds"))
saveRDS(integrated_HBC_lineage_scATAC,str_c(out_dir,"integrated_HBC_lineage_scATAC.rds"))
