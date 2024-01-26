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
out_dir <- "~/Cd36_OSNs/output/scATAC/Signac/HBC_lineage_scATAC/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}


## Load data
integrated_scATAC <- readRDS("~/Cd36_OSNs/output/scATAC/Signac/integrated_analysis/integrated_OE_scATAC.rds")

# call peaks 
HBC_lineage_cells <- rownames(integrated_scATAC@meta.data)[which(integrated_scATAC$cell_type %in% c("HBCs","GBCs","INPs","Immature ORNs","Mature ORNs"))]
HBC_lineage_scATAC <- subset(integrated_scATAC,cells=HBC_lineage_cells)
HBC_lineage_peaks <- CallPeaks(HBC_lineage_scATAC, macs2.path = "/md01/shipy3/ori/miniconda3/conda_software/bin/macs2",effective.genome.size=2.3e+09,group.by="cell_type",idents=c("HBCs","GBCs","INPs","Immature ORNs","Mature ORNs"))

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
HBC_lineage_peaks <- keepStandardChromosomes(HBC_lineage_peaks, pruning.mode = "coarse")
mm10_blacklist_bed <- read.table("/md01/shipy3/ori/DB/mm10/annotation/mm10-blacklist.v2.bed",sep="\t")
mm10_blacklist.gr <- GRanges(mm10_blacklist_bed$V1,
  IRanges(start = mm10_blacklist_bed$V2+1, end = mm10_blacklist_bed$V3))
HBC_lineage_peaks <- subsetByOverlaps(x = HBC_lineage_peaks, ranges = mm10_blacklist.gr, invert = TRUE)

# Load data and create Seurat objects
samples <- c("19d_0818","8w_0720")
dirs <- c("~/Cd36_OSNs/output/scATAC/cellranger/19d_0818/outs/","~/Cd36_OSNs/output/scATAC/cellranger/8w_0720/outs/")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"


kept_cells.ls <- lapply(1:length(samples),function(i){
  cells <- grep(samples[i],HBC_lineage_cells,value=TRUE)
  cells <- gsub(str_c(samples[i],"_"),"",cells)
  cells
  })

# Create the objects
HBC_lineage_obj.ls <- lapply(1:length(samples),function(i){
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
    features = HBC_lineage_peaks,
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
HBC_lineage_obj.ls[[1]]$dataset <- "19d"
HBC_lineage_obj.ls[[2]]$dataset <- "8w"

# merge all datasets, adding a cell ID to make sure cell names are unique
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


# To find integration anchors across datasets, we need to project them into a shared low-dimensional space. To do this, we’ll use reciprocal LSI projection (projecting each dataset into the others LSI space) by setting reduction="rlsi"
integration.anchors <- FindIntegrationAnchors(
  object.list = list(add_prefix_HBC_lineage_obj.ls[[1]], add_prefix_HBC_lineage_obj.ls[[2]]),
  anchor.features = rownames(add_prefix_HBC_lineage_obj.ls[[1]]),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings,perform dataset integration using a pre-computed Anchorset of specified low dimensional representations
integrated_HBC_lineage_scATAC <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined_HBC_lineage_scATAC[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated_HBC_lineage_scATAC <- RunUMAP(integrated_HBC_lineage_scATAC, reduction = "integrated_lsi", dims = 2:30)

pdf(str_c(out_dir,"integrated_HBC_lineage_UMAP.pdf"))
DimPlot(integrated_HBC_lineage_scATAC)
dev.off()

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


# add annotation
integrated_mature_ORNs_scATAC <- readRDS("~/Cd36_OSNs/output/scATAC/Signac/integrated_analysis/integrated_mature_ORNs_scATAC.rds")
integrated_mature_ORNs_scATAC$Cd36_or_not <- integrated_mature_ORNs_scATAC$cell_type
integrated_mature_ORNs_scATAC$Cd36_or_not[which(integrated_mature_ORNs_scATAC$Cd36_or_not %in% c("Dorsal ORNs","Ventral ORNs"))] <- "Non-Cd36+ ORNs"
cell_type_df <- integrated_scATAC@meta.data[HBC_lineage_cells,"cell_type",drop=FALSE]
Cd36_or_not <- c("Cd36+ ORNs","Non-Cd36+ ORNs","Gucy1b2+ cells")
for (i in 1:length(Cd36_or_not)){
  cells <- rownames(integrated_mature_ORNs_scATAC@meta.data)[which(integrated_mature_ORNs_scATAC$Cd36_or_not==Cd36_or_not[i])]
  cell_type_df[cells,"cell_type"] <- Cd36_or_not[i]
}
integrated_HBC_lineage_scATAC <- AddMetaData(integrated_HBC_lineage_scATAC,cell_type_df)
Idents(integrated_HBC_lineage_scATAC) <- integrated_HBC_lineage_scATAC$cell_type
Idents(integrated_HBC_lineage_scATAC) <- factor(Idents(integrated_HBC_lineage_scATAC),levels=c("HBCs","GBCs","INPs","Immature ORNs","Non-Cd36+ ORNs","Cd36+ ORNs","Gucy1b2+ cells"))


# calculate chromVar scores
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# Tshz1 motif1
Tshz1_meme_01 <- read.table("~/Cd36_OSNs/input/motifs/Tshz1_RTGACTCA.meme",sep=" ",skip=5,header=FALSE)
colnames(Tshz1_meme_01) <- c("A","C","G","T")
Tshz1_meme_01 <- as.matrix(t(round(Tshz1_meme_01*2672)))
Tshz1_pfm_01 <- PFMatrix(ID="ENSMUSG00000046982_Tshz1_01", name="Tshz1_01", matrixClass="Unknown",
                       strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                       tags=list(ensembl="ENSMUSG00000046982",species="Homo sapiens"),
                       profileMatrix=Tshz1_meme_01)
pfm[["ENSMUSG00000046982_Tshz1_01"]] <- Tshz1_pfm_01

# Tshz1 motif2
Tshz1_meme_02 <- read.table("~/Cd36_OSNs/input/motifs/Tshz1_ATGACMTCATY.meme",sep=" ",skip=5,header=FALSE)
colnames(Tshz1_meme_02) <- c("A","C","G","T")
Tshz1_meme_02 <- as.matrix(t(round(Tshz1_meme_02*1881)))
Tshz1_pfm_02 <- PFMatrix(ID="ENSMUSG00000046982_Tshz1_02", name="Tshz1_02", matrixClass="Unknown",
                       strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                       tags=list(ensembl="ENSMUSG00000046982",species="Homo sapiens"),
                       profileMatrix=Tshz1_meme_02)
pfm[["ENSMUSG00000046982_Tshz1_02"]] <- Tshz1_pfm_02

# Tshz2 motif1
Tshz2_meme_01 <- read.table("~/Cd36_OSNs/input/motifs/Tshz2_TCTGBTWTHTRTCA_01.meme",sep=" ",skip=5,header=FALSE)
colnames(Tshz2_meme_01) <- c("A","C","G","T")
Tshz2_meme_01 <- as.matrix(t(round(Tshz2_meme_01*4655)))
Tshz2_pfm_01 <- PFMatrix(ID="ENSMUSG00000047907_Tshz2_01", name="Tshz2_01", matrixClass="Unknown",
                       strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                       tags=list(ensembl="ENSMUSG00000047907",species="Homo sapiens"),
                       profileMatrix=Tshz2_meme_01)
pfm[["ENSMUSG00000047907_Tshz2_01"]] <- Tshz2_pfm_01


# Tshz2 motif2
Tshz2_meme_02 <- read.table("~/Cd36_OSNs/input/motifs/Tshz2_TRATTTRATTW_02.meme",sep=" ",skip=5,header=FALSE)
colnames(Tshz2_meme_02) <- c("A","C","G","T")
Tshz2_meme_02 <- as.matrix(t(round(Tshz2_meme_02*3234)))
Tshz2_pfm_02 <- PFMatrix(ID="ENSMUSG00000047907_Tshz2_02", name="Tshz2_02", matrixClass="Unknown",
                       strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                       tags=list(ensembl="ENSMUSG00000047907",species="Homo sapiens"),
                       profileMatrix=Tshz2_meme_02)
pfm[["ENSMUSG00000047907_Tshz2_02"]] <- Tshz2_pfm_02



# Tshz2 motif3
Tshz2_meme_03 <- read.table("~/Cd36_OSNs/input/motifs/Tshz2_DGVCAKCTGG_03.meme",sep=" ",skip=5,header=FALSE)
colnames(Tshz2_meme_03) <- c("A","C","G","T")
Tshz2_meme_03 <- as.matrix(t(round(Tshz2_meme_03*5102)))
Tshz2_pfm_03 <- PFMatrix(ID="ENSMUSG00000047907_Tshz2_03", name="Tshz2_03", matrixClass="Unknown",
                       strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                       tags=list(ensembl="ENSMUSG00000047907",species="Homo sapiens"),
                       profileMatrix=Tshz2_meme_03)
pfm[["ENSMUSG00000047907_Tshz2_03"]] <- Tshz2_pfm_03


# Trps1 motif
Trps1_meme <- read.table("~/Cd36_OSNs/input/motifs/Trps1_MA1970.1.meme",sep="",skip=11,header=FALSE)
colnames(Trps1_meme) <- c("A","C","G","T")
Trps1_meme <- as.matrix(t(round(Trps1_meme*18283)))
Trps1_pfm_df <- PFMatrix(ID="ENSMUSG00000038679_Trps1", name="Trps1", matrixClass="Unknown",
                       strand="*", bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                       tags=list(ensembl="ENSMUSG00000038679",species="Homo sapiens"),
                       profileMatrix=Trps1_meme)
pfm[["ENSMUSG00000038679_Trps1"]] <- Trps1_pfm_df


pfm_df <- c()
for (i in 1:length(pfm)){
  pfm_df <- data.frame(ID=pfm[[i]]@ID,name=pfm[[i]]@name,species=pfm[[i]]@tags$species) %>% rbind(pfm_df,.)
}

integrated_HBC_lineage_scATAC <- AddMotifs(
  object = integrated_HBC_lineage_scATAC,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

DefaultAssay(integrated_HBC_lineage_scATAC) <- "ATAC"
integrated_HBC_lineage_scATAC <- RunChromVAR(
  object = integrated_HBC_lineage_scATAC,
  genome = BSgenome.Mmusculus.UCSC.mm10
)


## chromVar z-score of Lhx2, Mef2a, Tshz1 for each cell type
Idents(integrated_HBC_lineage_scATAC) <- integrated_HBC_lineage_scATAC$cell_type

chromvar_df <- as.data.frame(integrated_HBC_lineage_scATAC@assays$chromvar@data) %>%
    tibble::rownames_to_column("motif")  %>%
    dplyr::filter(motif %in% motifs)
chromvar_df <- reshape2::melt(chromvar_df,id="motif",variable.name="cell",value.name ="chromvar_zscore")
meta_df <- FetchData(integrated_HBC_lineage_scATAC,vars="cell_type")
chromvar_df <- merge(chromvar_df,meta_df,by.x="cell",by.y="row.names")

average_chromvar_df <- chromvar_df %>%
  dplyr::group_by(motif,cell_type) %>%
  dplyr::summarise(mean_chromvar=mean(chromvar_zscore)) %>%
  filter(cell_type %in% c("GBCs","INPs","Immature ORNs","Cd36+ ORNs","Non-Cd36+ ORNs") )
average_chromvar_df$cell_type <- factor(average_chromvar_df$cell_type,levels=c("GBCs","INPs","Immature ORNs","Cd36+ ORNs","Non-Cd36+ ORNs"),labels=c("GBCs","INPs","Immature OSNs","Cd36+ OSNs","Cd36- OSNs"))
average_chromvar_df$motif <- factor(average_chromvar_df$motif,levels=motifs,labels=motif_names)


exclude_HBCs_scATAC <- subset(integrated_HBC_lineage_scATAC,subset=cell_type %in% c("GBCs","INPs","Immature ORNs","Cd36+ ORNs","Non-Cd36+ ORNs"))
p <- DotPlot(exclude_HBCs_scATAC,assay="chromvar",features=motifs)
pct_df <- p$data[,c('pct.exp','id')]
pct_df$motif <- rep(motif_names,5)
pct_df$motif <- factor(pct_df$motif,levels=motif_names)
pct_df$id <- factor(pct_df$id,levels=c("GBCs","INPs","Immature ORNs","Cd36+ ORNs","Non-Cd36+ ORNs"),labels=c("GBCs","INPs","Immature OSNs","Cd36+ OSNs","Cd36- OSNs"))

merged_df <- merge(average_chromvar_df,pct_df,by.x=c("motif","cell_type"),by.y=c("motif","id"))
merged_df$motif <- factor(merged_df$motif,levels=rev(motif_names))

## Figure 4H
pdf(str_c(out_dir,"Lhx2_Mef2a_Tshz1_chromvar_zcore_dotplot.pdf"),width=8,height=2.8)
ggplot(data=merged_df,aes(x=cell_type,y=motif))+
  geom_point(aes(size=pct.exp,fill=mean_chromvar),pch=21,color="black")+
  scale_size_continuous(range=c(1,10),breaks = c(0,25,50,75),labels=c("0","25","50","75"),limits=c(0,100))+
  scale_fill_gradientn(colours =as.vector(ArchRPalettes[["solarExtra"]])[c(2,5:8)] ,limits=c(-1,2),oob = scales::squish,,guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.hjust=0.5,title.position = "top"))+
  labs(x="",y="",size="Percent positive",color="Average chromvar z-score")+
  theme_bw()+
  theme(panel.background=element_rect(color='black',linetype="solid",size=1.6),axis.text.x = element_text(color='black',size=11),axis.text.y=element_text(color='black',size=11),plot.title = element_text(hjust = 0.5,size=13),legend.position="bottom")
dev.off()

## DARs between cell types
filtered_HBC_lineage_scATAC <- subset(integrated_HBC_lineage_scATAC,subset=rough_cell_type %in% c("HBCs","GBCs","INPs","Immature ORNs","Mature ORNs"))
DefaultAssay(filtered_HBC_lineage_scATAC) <- "ATAC"
cell_types <- c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs")
cell_type_DARs.ls <- lapply(1:length(cell_types),function(i){
  df <- FindMarkers(
    object = filtered_HBC_lineage_scATAC,
    ident.1 = cell_types[i],
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'nCount_ATAC'
    )
  df <- df %>%
    dplyr::filter(p_val_adj < 0.05,avg_log2FC>0) %>%
    dplyr::arrange(desc(avg_log2FC)) 
  df
  })
cell_type_DARs <- unlist(sapply(1:length(cell_types),function(i){
  rownames(cell_type_DARs.ls[[i]])
  }))
cell_type_DARs_ann.ls <- lapply(1:length(cell_types),function(i){
  df <- ClosestFeature(object=filtered_HBC_lineage_scATAC,
  regions=rownames(cell_type_DARs.ls[[i]]))
  df
  })
cell_type_DARs <- cell_type_DARs[-which(duplicated(cell_type_DARs))]
avg_df <- AverageExpression(filtered_HBC_lineage_scATAC,assays="ATAC",features=cell_type_DARs)$ATAC
scaled_df <- apply(avg_df,1,scale)
rownames(scaled_df) <- colnames(avg_df)

library(ComplexHeatmap)
pdf(str_c(out_dir,"filtered_HBC_lineage_scATAC_each_cell_type_marker_peaks_heatmap.pdf"),width=12,height=5)
Heatmap(scaled_df,
  name="ATAC z-score",
  col=as.vector(ArchRPalettes[["blueYellow"]]),
  cluster_columns=FALSE,
  column_title_gp = gpar(fontsize = 15),
  column_title="Marker peaks for each cell type",
  show_column_names=FALSE,
  show_row_names=TRUE,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 11),
  cluster_rows=FALSE,
  cluster_column_slices=FALSE,
  width = unit(14, "cm"),
  height=unit(5.5, "cm"),
  border_gp = gpar(col = "black")
)
dev.off()


## motif enrichment of DARs for each cell type
cell_type_DARs_enriched_motifs.ls <- lapply(1:length(cell_types),function(i){
  df <- FindMotifs(
    object = filtered_HBC_lineage_scATAC,
    features = rownames(cell_type_DARs.ls[[i]])
  ) 
  #df <- df %>%
  #  filter(pvalue<0.05,fold.enrichment>=1.5)
  df
  })

saveRDS(integrated_HBC_lineage_scATAC,str_c(out_dir,"integrated_HBC_lineage_scATAC.rds"))

