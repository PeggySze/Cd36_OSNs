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
integrated_scATAC <- readRDS(str_c(out_dir,"integrated_mature_ORNs_scATAC.rds"))

# remove Gucy1b2+ cells from scATAC-seq dataset
filtered_mature_ORNs_scATAC <- subset(integrated_mature_ORNs_scATAC,subset=cell_type %in% c("Dorsal ORNs","Ventral ORNs","Cd36+ ORNs"))

# Load scRNA object
integrated_mature_ORNs.scRNA <- readRDS("~/Cd36_OSNs/output/scRNA/seurat/integrated_analysis/integrate_Tsukahara_Brann/remove_homecage_4_integrated_mature_ORNs_Seurat_obj.rds")
filtered_Ifi27_cells_mature_ORNs <- subset(integrated_mature_ORNs.scRNA,subset=cell_type!="Ifi27+ cells")
filtered_Ifi27_cells_mature_ORNs$Cd36_or_not <- ifelse(filtered_Ifi27_cells_mature_ORNs$cell_type=="Cd36+ ORNs","Cd36+ ORNs","Cd36- ORNs")
Idents(filtered_Ifi27_cells_mature_ORNs) <- filtered_Ifi27_cells_mature_ORNs$Cd36_or_not

## Figure 3A - mature ORNs scATAC-seq UMAP showing Cd36 genescore 
DefaultAssay(filtered_mature_ORNs_scATAC) <- "RNA"
Cd36_GeneScore_df <- FetchData(filtered_mature_ORNs_scATAC,vars = c("Cd36","UMAP_1","UMAP_2"),slot = "data")
pdf(str_c(out_dir,"mature_ORNs_scATAC_UMAP_highlight_Cd36+_ORNs.pdf"),height=9,width=8)
ggplot(data=Cd36_GeneScore_df,aes(x=UMAP_1,y=UMAP_2,color=Cd36))+
  geom_point(size=0.4)+
  scale_colour_gradientn(colours = palettes,limits=c(0,quantile(Cd36_GeneScore_df$Cd36,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(Cd36_GeneScore_df$Cd36,seq(0,1,0.01))[96]),labels = c("0", "1.1"),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.position = "top")) + 
  labs(title="Mature ORNs scATAC-seq",color="Cd36 GeneScore") +
  geom_segment(aes(x = -7, y = -3.5, xend = -5, yend = -3.5),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  geom_segment(aes(x = -7, y = -3.5, xend = -7, yend = -2.5),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  annotate(geom = "text", x = -6, y = -3.7, label = "UMAP 1", color = "black",size=4.5) +
  annotate(geom = "text", x = -7.5, y = -3, label = "UMAP 2", color = "black",angle = 90,size=4.5) +
  annotate(geom = "text", x = 6, y = 2.2, label = "Cd36+ ORNs", color = "black",size=5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(2.5)),axis.title=element_blank(),axis.ticks=element_blank(),legend.key.size = unit(0.8, 'cm'),legend.text=element_text(size=rel(1.2)),legend.position="bottom",legend.title=element_text(size=rel(1.4)),axis.text=element_blank())
dev.off()

## Find Cd36+ vs Cd36- OSNs differentially accessible peaks 
DefaultAssay(filtered_mature_ORNs_scATAC) <- "ATAC"
filtered_mature_ORNs_scATAC$Cd36_or_not <- ifelse(filtered_mature_ORNs_scATAC$cell_type=="Cd36+ ORNs","Cd36+ ORNs","Cd36- ORNs")
filtered_mature_ORNs_scATAC$Cd36_or_not <- factor(filtered_mature_ORNs_scATAC$Cd36_or_not,levels=c("Cd36+ ORNs","Cd36- ORNs"))
Idents(filtered_mature_ORNs_scATAC) <- filtered_mature_ORNs_scATAC$Cd36_or_not

Cd36_vs_other_da_peaks <- FindMarkers(
  object = filtered_mature_ORNs_scATAC,
  ident.1 = "Cd36+ ORNs",
  ident.2 = "Cd36- ORNs",
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

Cd36_up_peaks <- Cd36_vs_other_da_peaks %>%
  dplyr::filter(p_val_adj < 0.05,avg_log2FC>=log2(1.5)) %>%
  dplyr::arrange(desc(avg_log2FC)) 
Cd36_up_peaks_ann <-  ClosestFeature(object=filtered_mature_ORNs_scATAC,
  regions=rownames(Cd36_up_peaks))
strig2GRanges <- function(string) {
  chr <- sapply(1:length(string),function(i){
    unlist(strsplit(string[i],"-"))[1]
    })
  start <- sapply(1:length(string),function(i){
    as.numeric(unlist(strsplit(string[i],"-"))[2])
    })
  end <- sapply(1:length(string),function(i){
    as.numeric(unlist(strsplit(string[i],"-"))[3])
    })
  gr <- GRanges(chr,IRanges(start = start, end = end), strand = "*")
  gr
}
Cd36_up_peaks.gr <- strig2GRanges(rownames(Cd36_up_peaks))
Cd36_up_peaks_ann_df <- merge(Cd36_up_peaks,Cd36_up_peaks_ann[,c("gene_name","query_region","distance")],by.x="row.names",by.y="query_region",all.x=TRUE,sort=FALSE)
write.table(Cd36_up_peaks_ann_df,str_c(out_dir,"Cd36_ORNs_up_peaks_with_annotation.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)


Cd36_down_peaks <- Cd36_vs_other_da_peaks %>%
  dplyr::filter(p_val_adj < 0.05,avg_log2FC<=-log2(1.5)) %>%
  dplyr::arrange(avg_log2FC) 
Cd36_down_peaks_ann <- ClosestFeature(object=filtered_mature_ORNs_scATAC,
  regions=rownames(Cd36_down_peaks))
Cd36_down_peaks_ann_df <-  merge(Cd36_down_peaks,Cd36_down_peaks_ann[,c("gene_name","query_region","distance")],by.x="row.names",by.y="query_region",all.x=TRUE,sort=FALSE)
write.table(Cd36_down_peaks_ann_df,str_c(out_dir,"Cd36_ORNs_down_peaks_with_annotation.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)


## Figure 3B - Heatmap showing DARs
DefaultAssay(filtered_mature_ORNs_scATAC) <- "ATAC"
set.seed(100)
cell_types <- c("Cd36+ ORNs","Non-Cd36+ ORNs")
sampled_df <- c()
for (i in cell_types){
  cells <- rownames(filtered_mature_ORNs_scATAC@meta.data)[which(filtered_mature_ORNs_scATAC$Cd36_or_not==i)]
  if (length(cells)>=500){
    sampled_cells <- sample(cells,500,replace=FALSE)
  } else{
    sampled_cells <- sample(cells,500,replace=TRUE)
  }
  sampled_df <- data.frame(idx=1:500,cell_type=i,sampled_cell=sampled_cells) %>% rbind(sampled_df ,.)
}
scATAC_df <- FetchData(object = filtered_mature_ORNs_scATAC, vars = unique(c("Cd36_or_not",rownames(Cd36_up_peaks),rownames(Cd36_down_peaks))),slot = "data",cells=sampled_df$sampled_cell)
scaled_scATAC_df <- t(apply(scATAC_df[,-1],2,scale))
colnames(scaled_scATAC_df) <- rownames(scATAC_df)
rescaled_scATAC_df <- scaled_scATAC_df
rescaled_scATAC_df[which(rescaled_scATAC_df>=0.5)] <- 0.5
rescaled_scATAC_df[which(rescaled_scATAC_df<=-0.5)] <- -0.5

cell_type_cols <- c("#00A087FF",brewer.pal(8,"Set2")[8])
names(cell_type_cols) <- c("Cd36+ ORNs","Other ORNs")
library(ComplexHeatmap)
pdf(str_c(out_dir,"Cd36_vs_other_ORNs_DARs_heatmap.pdf"),width=3,height=6)
ha2 <- HeatmapAnnotation(
  df=data.frame(row.names=rownames(scATAC_df),cell_type=scATAC_df$Cd36_or_not),
  col=list(cell_type=cell_type_cols),
  annotation_legend_param = list(cell_type=list(title="Cell type")),
  gp = gpar(col = "black"),
  show_legend=FALSE,
  show_annotation_name=FALSE
  )
h2 <- Heatmap(rescaled_scATAC_df,
  name="ATAC z-score",
  col=as.vector(ArchRPalettes[["blueYellow"]]),
  cluster_columns=FALSE,
  column_split=factor(as.vector(scATAC_df$Cd36_or_not),levels=c("Cd36+ ORNs","Non-Cd36+ ORNs")),
  column_title_gp = gpar(fontsize = 11),
  show_column_names=FALSE,
  show_row_names=FALSE,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 11),
  column_gap = unit(0, "mm"),
  cluster_rows=FALSE,
  cluster_column_slices=FALSE,
  top_annotation=ha2,
  heatmap_legend_param = list(direction = "horizontal",legend_width = unit(4, "cm"), title_position = "topcenter",border="black")
)
draw(h2, heatmap_legend_side = "bottom")
dev.off()

## motif enrichment of DARs
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motifs
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

filtered_mature_ORNs_scATAC <- AddMotifs(
  object = filtered_mature_ORNs_scATAC,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

Cd36_up_enriched.motifs <- FindMotifs(
  object = filtered_mature_ORNs_scATAC,
  features = Cd36_up_peaks_ann_df$Row.names
) 
Cd36_up_enriched.motifs$rank <- 1:nrow(Cd36_up_enriched.motifs)

Cd36_down_enriched.motifs <- FindMotifs(
  object = filtered_mature_ORNs_scATAC,
  features = Cd36_down_peaks_ann_df$Row.names
)
Cd36_down_enriched.motifs$rank <- 1:nrow(Cd36_down_enriched.motifs)

## Figure S5D - ranking plot of Cd36+ OSNs down peaks enriched motifs
library(ggrepel)
Cd36_down_enriched_motifs_df <- Cd36_down_enriched.motifs
Cd36_down_enriched_motifs_df$motif.name <- gsub("\\(var\\.\\d\\)","",Cd36_down_enriched_motifs_df$motif.name)
Cd36_down_enriched_motifs_df$motif.name <- gsub("\\_02","",Cd36_down_enriched_motifs_df$motif.name)
Cd36_down_enriched_motifs_df <- Cd36_down_enriched_motifs_df[-grep("\\:\\:",Cd36_down_enriched_motifs_df$motif.name),]
convert2_mouse_symbol <- function(symbol){
  first_char <- toupper(substr(symbol,1,1))
  other_chars <- tolower(substr(symbol,2,nchar(symbol)))
  coverted_symbol <- str_c(first_char,other_chars)
  return(coverted_symbol)
}
Cd36_down_enriched_motifs_df$mouse_symbol <- sapply(1:nrow(Cd36_down_enriched_motifs_df),function(i){
  mouse_symbol <- convert2_mouse_symbol(Cd36_down_enriched_motifs_df$motif.name[i])
  mouse_symbol
  })
Cd36_down_enriched_motifs_df$label <- ifelse(Cd36_down_enriched_motifs_df$mouse_symbol %in% c(Cd36_down_enriched_motifs_df$mouse_symbol[1:10],"Lhx2"),Cd36_down_enriched_motifs_df$mouse_symbol,"")
Cd36_down_enriched_motifs_df$rank <- 1:nrow(Cd36_down_enriched_motifs_df)
Cd36_down_enriched_motifs_df$mlog10pvalue <- -log10(Cd36_down_enriched_motifs_df$pvalue)
pdf(str_c(out_dir,"Cd36_OSNs_down_peaks_enriched_motifs_ranking_plot.pdf"),height=7.5)
ggplot(data=Cd36_down_enriched_motifs_df[1:50,],aes(x=rank,y=mlog10pvalue))+
  geom_point(aes(color=fold.enrichment,size=percent.observed))+
  labs(x="Motif rank",y="-log10(p_value)",color="Enrichment",size="Percentage of peaks with corresponding motif",title="Enriched motifs of down-regulated peaks in Cd36+ OSNs")+
  geom_text_repel(aes(label=label),size=5,point.padding = 0.2,
    nudge_x = .15,
    nudge_y = .5,
    segment.curvature = -1e-20,
    arrow = arrow(length = unit(0.015, "npc")))+
  scale_colour_gradient(low = "yellow", high = "red")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position="bottom")
dev.off()
write.table(Cd36_down_enriched_motifs_df,str_c(out_dir,"Cd36_down_peaks_enriched_motifs.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

# Figure S5C - ranking plot of Cd36+ OSNs up peaks enriched motifs
Cd36_up_enriched_motifs_df <- Cd36_up_enriched.motifs
Cd36_up_enriched_motifs_df$motif.name <- gsub("\\(var\\.\\d\\)","",Cd36_up_enriched_motifs_df$motif.name)
Cd36_up_enriched_motifs_df$motif.name <- gsub("\\_02","",Cd36_up_enriched_motifs_df$motif.name)
Cd36_up_enriched_motifs_df <- Cd36_up_enriched_motifs_df[-grep("\\:\\:",Cd36_up_enriched_motifs_df$motif.name),]
convert2_mouse_symbol <- function(symbol){
  first_char <- toupper(substr(symbol,1,1))
  other_chars <- tolower(substr(symbol,2,nchar(symbol)))
  coverted_symbol <- str_c(first_char,other_chars)
  return(coverted_symbol)
}
Cd36_up_enriched_motifs_df$mouse_symbol <- sapply(1:nrow(Cd36_up_enriched_motifs_df),function(i){
  mouse_symbol <- convert2_mouse_symbol(Cd36_up_enriched_motifs_df$motif.name[i])
  mouse_symbol
  })
Cd36_up_enriched_motifs_df$label <- ifelse(Cd36_up_enriched_motifs_df$mouse_symbol %in% c(Cd36_up_enriched_motifs_df$mouse_symbol[1:10],"Tshz1","Atf4"),Cd36_up_enriched_motifs_df$mouse_symbol,"")
Cd36_up_enriched_motifs_df$rank <- 1:nrow(Cd36_up_enriched_motifs_df)
Cd36_up_enriched_motifs_df$mlog10pvalue <- -log10(Cd36_up_enriched_motifs_df$pvalue)
pdf(str_c(out_dir,"Cd36_OSNs_up_peaks_enriched_motifs_ranking_plot.pdf"),height=7.5)
ggplot(data=Cd36_up_enriched_motifs_df[1:50,],aes(x=rank,y=mlog10pvalue))+
  geom_point(aes(color=fold.enrichment,size=percent.observed))+
  labs(x="Motif rank",y="-log10(p_value)",color="Enrichment",size="Percentage of peaks with corresponding motif",title="Enriched motifs of up-regulated peaks in Cd36+ OSNs")+
  geom_text_repel(aes(label=label),size=5,point.padding = 0.2,
    nudge_x = .15,
    nudge_y = .5,
    segment.curvature = -1e-20,
    arrow = arrow(length = unit(0.015, "npc")),max.overlaps=500)+
  scale_colour_gradient(low = "yellow", high = "red")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1),legend.position="bottom")
dev.off()
write.table(Cd36_up_enriched_motifs_df,str_c(out_dir,"Cd36_up_peaks_enriched_motifs.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)


## Figure 3D - expression levels vs motif enrichment levels for top enriched TFs (split into two figures)
shown_TF_motifs <- c(Cd36_up_enriched_motifs_df$motif[c(1:8,10)],"ENSMUSG00000046982_Tshz1_02")
shown_TF_names <- c(Cd36_up_enriched_motifs_df$mouse_symbol[c(1:8,10)],"Tshz1")

Cd36_down_enrichment_df <- Cd36_down_enriched_motifs_df %>%
  filter(motif %in% shown_TF_motifs) %>%
  select(motif,fold.enrichment,pvalue,p.adjust,mouse_symbol) %>%
  mutate(cell_type="Cd36- ORNs")
Cd36_up_enrichment_df <- Cd36_up_enriched_motifs_df %>%
  filter(motif %in% shown_TF_motifs) %>%
  select(motif,fold.enrichment,pvalue,p.adjust,mouse_symbol) %>%
  mutate(cell_type="Cd36+ ORNs")
enrichment_df <- rbind(Cd36_down_enrichment_df,Cd36_up_enrichment_df)

DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"
exp_df <- as.data.frame(AverageExpression(filtered_Ifi27_cells_mature_ORNs,features=shown_TF_names,assays="RNA")$RNA)
exp_df$mouse_symbol <- rownames(exp_df)
exp_df <- reshape2::melt(exp_df,id="mouse_symbol",variable.name="cell_type",value.name ="expression")

combined_df <- merge(enrichment_df,exp_df,by=c("mouse_symbol","cell_type"))
combined_df$mouse_symbol <- factor(combined_df$mouse_symbol,levels=rev(c("Mef2a","Mef2b","Mef2c","Mef2d",shown_TF_names[5:10])))
combined_df$cell_type <- factor(combined_df$cell_type,levels=rev(c("Cd36- ORNs","Cd36+ ORNs")),labels=rev(c("Cd36- OSNs","Cd36+ OSNs")))
combined_df$mlog10pvalue <- -log10(combined_df$pvalue)

newpalette <- brewer.pal(9,"Greens")[1:8]
p1 <- ggplot(data=combined_df,aes(y=mouse_symbol,x=cell_type)) + 
  geom_tile(aes(fill=expression),colour = "black")+
  scale_fill_gradientn(colours=newpalette,guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.position="top"),limits=c(0,10),oob = scales::squish)+
  geom_point(aes(size=mlog10pvalue))+
  scale_size_continuous(range=c(1,14),breaks = seq(0,25,10),limits=c(0,25))+
  labs(size="-log10(p_value)",fill="gene expression")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks=element_blank(),axis.text.x=element_text(angle =45,hjust=1,color="black",size=rel(1.2)),axis.title.x=element_blank(),plot.margin=unit(c(0.2, 0, 0, 0.2), "cm"),axis.text.y=element_text(color="black",size=rel(1.4)),axis.title.y=element_blank(),legend.position="none")


Cd36_vs_other_exp_df <- FindMarkers(filtered_Ifi27_cells_mature_ORNs,ident.1="Cd36+ ORNs",features=shown_TF_names,logfc.threshold=-Inf,min.pct=-Inf)
Cd36_vs_other_exp_df$TF <- rownames(Cd36_vs_other_exp_df)
Cd36_vs_other_exp_df$x <- 1
Cd36_vs_other_exp_df$TF <- factor(Cd36_vs_other_exp_df$TF,levels=rev(c("Mef2a","Mef2b","Mef2c","Mef2d",shown_TF_names[5:10])))
p2 <- ggplot(data=Cd36_vs_other_exp_df,aes(x=x,y=TF)) + 
  geom_tile(aes(fill=avg_log2FC),colour = "black")+
  scale_fill_gradientn(colours =c(brewer.pal(9,"Blues")[2],"#FFFFFF",brewer.pal(9,"Reds")[c(4,6)]) ,limits=c(-0.2,0.6),oob = scales::squish,guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.position="top"))+
  labs(fill="Cd36+ vs Cd36- OSNs expression Log2FoldChange")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),plot.margin=unit(c(0.2, 0, 0, 0.2), "cm"),axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position="none")

pdf(str_c(out_dir,"Cd36_vs_other_OSNs_enriched_motifs_enrichment&expression_partI.pdf"),height=8,width=3)
p1 + p2 + plot_layout(widths=c(2,0.4))
dev.off()

# part II
shown_TF_motifs <- Cd36_down_enriched_motifs_df$motif[c(1:10)]
shown_TF_names <- Cd36_down_enriched_motifs_df$mouse_symbol[c(1:10)]

Cd36_down_enrichment_df <- Cd36_down_enriched_motifs_df %>%
  filter(motif %in% shown_TF_motifs) %>%
  select(motif,fold.enrichment,pvalue,p.adjust,mouse_symbol) %>%
  mutate(cell_type="Cd36- ORNs")
Cd36_up_enrichment_df <- Cd36_up_enriched_motifs_df %>%
  filter(motif %in% shown_TF_motifs) %>%
  select(motif,fold.enrichment,pvalue,p.adjust,mouse_symbol) %>%
  mutate(cell_type="Cd36+ ORNs")
enrichment_df <- rbind(Cd36_down_enrichment_df,Cd36_up_enrichment_df)

DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"
exp_df <- as.data.frame(AverageExpression(filtered_Ifi27_cells_mature_ORNs,features=shown_TF_names,assays="RNA")$RNA)
exp_df$mouse_symbol <- rownames(exp_df)
exp_df <- reshape2::melt(exp_df,id="mouse_symbol",variable.name="cell_type",value.name ="expression")

combined_df <- merge(enrichment_df,exp_df,by=c("mouse_symbol","cell_type"))
combined_df$mouse_symbol <- factor(combined_df$mouse_symbol,levels=rev(shown_TF_names))
combined_df$cell_type <- factor(combined_df$cell_type,levels=rev(c("Cd36- ORNs","Cd36+ ORNs")),labels=rev(c("Cd36- OSNs","Cd36+ OSNs")))
combined_df$mlog10pvalue <- -log10(combined_df$pvalue)

newpalette <- brewer.pal(9,"Greens")[1:8]
p1 <- ggplot(data=combined_df,aes(y=mouse_symbol,x=cell_type)) + 
  geom_tile(aes(fill=expression),colour = "black")+
  scale_fill_gradientn(colours=newpalette,guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.position="top"),,limits=c(0,10),oob = scales::squish)+
  geom_point(aes(size=mlog10pvalue))+
  scale_size_continuous(range=c(1,14),breaks = seq(0,25,10),limits=c(0,25))+
  labs(size="-log10(p_value)",fill="gene expression")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks=element_blank(),axis.text.x=element_text(angle =45,hjust=1,color="black",size=rel(1.2)),axis.title.x=element_blank(),plot.margin=unit(c(0.2, 0, 0, 0.2), "cm"),axis.text.y=element_text(color="black",size=rel(1.4)),axis.title.y=element_blank(),legend.position="right")

Cd36_vs_other_exp_df <- FindMarkers(filtered_Ifi27_cells_mature_ORNs,ident.1="Cd36+ ORNs",features=shown_TF_names,logfc.threshold=-Inf,min.pct=-Inf)
Cd36_vs_other_exp_df$TF <- rownames(Cd36_vs_other_exp_df)
Cd36_vs_other_exp_df$x <- 1
Cd36_vs_other_exp_df$TF <- factor(Cd36_vs_other_exp_df$TF,levels=rev(shown_TF_names))
p2 <- ggplot(data=Cd36_vs_other_exp_df,aes(x=x,y=TF)) + 
  geom_tile(aes(fill=avg_log2FC),colour = "black")+
  scale_fill_gradientn(colours =c(brewer.pal(9,"Blues")[2],"#FFFFFF",brewer.pal(9,"Reds")[c(4,6)]) ,limits=c(-0.2,0.6),oob = scales::squish,guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.position="top"))+
  labs(fill="Cd36+ vs Cd36- OSNs expression Log2FoldChange")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),plot.margin=unit(c(0.2, 0, 0, 0.2), "cm"),axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position="right")

pdf(str_c(out_dir,"Cd36_vs_other_OSNs_enriched_motifs_enrichment&expression_partII.pdf",height=8,width=8.4)
p1 + p2 + plot_layout(widths=c(8.2,1.5))
dev.off()

## chromVAR
# Compute a per-cell motif activity score by running chromVAR
DefaultAssay(filtered_mature_ORNs_scATAC) <- "ATAC"
filtered_mature_ORNs_scATAC <- RunChromVAR(
  object = filtered_mature_ORNs_scATAC,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# differential gene expression vs differential TF activity
DefaultAssay(filtered_mature_ORNs_scATAC) <- 'chromvar'
diff_chromvar_df <- FindMarkers(
  object = filtered_mature_ORNs_scATAC,
  ident.1 = 'Cd36+ ORNs',
  ident.2 = 'Cd36- ORNs',
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  logfc.threshold=-Inf,
  min.pct=-Inf
)
diff_chromvar_df <- merge(diff_chromvar_df,pfm_df,by.x="row.names",by.y="ID",all.x=TRUE)
diff_chromvar_df <- diff_chromvar_df[-grep("Tshz1-01",diff_chromvar_df$Row.names),]
diff_chromvar_df$name[1:5] <- c("Trps1","Tshz1","Tshz2","Tshz2","Tshz2")
diff_chromvar_df$name <- gsub("\\(var\\.\\d\\)","",diff_chromvar_df$name)
diff_chromvar_df <- diff_chromvar_df[-grep("\\:\\:",diff_chromvar_df$name),]
convert2_mouse_symbol <- function(symbol){
  first_char <- toupper(substr(symbol,1,1))
  other_chars <- tolower(substr(symbol,2,nchar(symbol)))
  coverted_symbol <- str_c(first_char,other_chars)
  return(coverted_symbol)
}
diff_chromvar_df$mouse_symbol <- sapply(1:nrow(diff_chromvar_df),function(i){
  mouse_symbol <- convert2_mouse_symbol(diff_chromvar_df$name[i])
  mouse_symbol
  })

DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"
genes <- unique(diff_chromvar_df$mouse_symbol[which(diff_chromvar_df$mouse_symbol %in% rownames(filtered_Ifi27_cells_mature_ORNs))])
diff_exp_df <- FindMarkers(filtered_Ifi27_cells_mature_ORNs,ident.1="Cd36+ ORNs",features=genes,logfc.threshold=-Inf,min.pct=-Inf)

merged_df <- merge(diff_exp_df,diff_chromvar_df,by.x="row.names",by.y="mouse_symbol",all.x=TRUE)
merged_df <- merged_df %>% 
  arrange(abs(avg_diff))
merged_df$label <- ifelse(merged_df$Row.names %in% c("Mef2a","Mef2b","Mef2c","Mef2d","Cebpg","Atf4","Tshz1","Pknox2","Pbx3","Pou6f2","Lhx2","Trps1","Tshz2"),merged_df$Row.names,"")
merged_df$label[which(merged_df$label=="Cebpg")[1]] <- ""
merged_df$label[which(merged_df$label=="Tshz2")[1:2]] <- ""
merged_df$chromvar_mlog10pvalue <- -log10(merged_df$p_val.y)
merged_df_01 <- merged_df %>%
  filter(label=="")
merged_df_02 <- merged_df %>%
  filter(label!="")

df_01 <- merged_df %>%
  filter(avg_log2FC<=-log2(1.5),p_val.x<0.05,avg_diff<0,p_val.y<0.05)
# Tshz2
df_02 <- merged_df %>%
  filter(avg_log2FC<0,p_val.x<0.05,avg_diff>0,p_val.y<0.05)
df_03 <- merged_df %>%
  filter(avg_log2FC>0,avg_diff<0,p_val.y<0.05)

# Figure 3C
pdf(str_c(out_dir,"Cd36_vs_other_OSNs_TF_expression_vs_activity.pdf"),width=10)
ggplot()+
  geom_rect(data=data.frame(xmin=0,xmax=Inf,ymin=0,ymax=Inf),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Reds")[3] , alpha = 0.2)+
  geom_rect(data=data.frame(xmin=-Inf,xmax=0,ymin=-Inf,ymax=0),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Reds")[3] , alpha = 0.2)+
  geom_rect(data=data.frame(xmin=0,xmax=Inf,ymin=-Inf,ymax=0),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Blues")[3] , alpha = 0.2)+
  geom_rect(data=data.frame(xmin=-Inf,xmax=0,ymin=0,ymax=Inf),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill = brewer.pal(9,"Blues")[3] , alpha = 0.2)+
  geom_point(data=merged_df_01,aes(x=avg_diff,y=avg_log2FC,size=pmax(pmin(chromvar_mlog10pvalue, 60), 0)),color="lightgrey")+
  geom_point(data=merged_df_02,aes(x=avg_diff,y=avg_log2FC,size=pmax(pmin(chromvar_mlog10pvalue, 60), 0),fill=avg_log2FC),pch=21,color="black")+
  scale_size_continuous(range=c(1,8),breaks = seq(0,60,20),limits=c(0,60))+
  scale_fill_gradientn(colours =c(brewer.pal(9,"Blues")[3],"#FFFFFF",brewer.pal(9,"Reds")[c(4,5,6,8)]) ,limits=c(-0.2,0.8),oob = scales::squish,guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.position="top",title.hjust = 0.5))+
  labs(y="Differential expression log2FoldChange",x="Differential TF chromVAR activity",title="Cd36+ vs Cd36- OSNs",size="-log10(Differential TF chromVAR activity p value)",fill="log2(expression FoldChange)")+
  geom_text_repel(data=merged_df_02,aes(x=avg_diff,y=avg_log2FC,label=label),size=5,point.padding=unit(1.6, "lines"),arrow = arrow(length=unit(0.01, "npc")),max.overlaps=500,force=3,segment.color = "#cccccc")+
  ylim(c(-3.1,3.1))+
  xlim(c(-1.5,1.5))+
  geom_hline(yintercept=0,linetype="longdash")+
  geom_vline(xintercept=0,linetype="longdash")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1))
dev.off()


## Figure 3E,G,I - gene expression of Mef2a,Lhx2 and Tshz1 (UMAP + Violin plot)
markers <- c("Mef2a","Lhx2","Tshz1")
gene_exp_palettes <- c(brewer.pal(8,"Pastel2")[8],brewer.pal(9,"Greens")[3:8])
DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"
gene_exp_UMAP_plot.ls <- lapply(1:length(markers),function(i){
  gene_exp_df <- FetchData(filtered_Ifi27_cells_mature_ORNs,vars = c(markers[i],"UMAP_1","UMAP_2"),slot = "data")
  colnames(gene_exp_df)[1] <- "gene"
  p <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.1)+
    scale_colour_gradientn(colours = gene_exp_palettes,limits=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),labels = c("0", round(quantile(gene_exp_df$gene,seq(0,1,0.01))[96],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1,title.position="top", title.hjust = 0.5)) + 
    scale_x_continuous(limits=c(-12.5,7.5),expand = c(0,0))+
    scale_y_continuous(limits=c(-7.5,8),expand = c(0,0))+
    labs(title=markers[i],color="Gene expression",x="UMAP 1",y="UMAP 2") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=3),plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_text(size=rel(1.5),vjust=0.5),axis.title.y=element_text(size=rel(1.5),vjust=0.5),plot.margin=unit(c(0.3, 0.2, 0.7, 0.5), "cm"),legend.title = element_text(size=15),legend.text=element_text(size=10),legend.position=c(0.2,0.15),legend.direction="horizontal",legend.key.size = unit(0.8, 'cm'))
    p
  })

violin_palettes <- c(brewer.pal(9,"Greens")[7],brewer.pal(8,"Set2")[8])
#violin_palettes <- colorspace::lighten(violin_palettes,amount = 0.3)
fmt_dcimals <- function(decimals=0){
    function(x) format(x,nsmall = decimals,scientific = FALSE)
  }
library(ggpubr)
Cd36_vs_other_exp_df <- FindMarkers(filtered_Ifi27_cells_mature_ORNs,features=markers,ident.1="Cd36+ ORNs",logfc.threshold=-Inf,min.pct=-Inf)
gene_exp_violin_plot.ls <- lapply(1:length(markers),function(i){
  gene_exp_df <- FetchData(filtered_Ifi27_cells_mature_ORNs,vars = c(markers[i],"UMAP_1","UMAP_2","Cd36_or_not"),slot = "data")
  colnames(gene_exp_df)[1] <- "gene" 
  gene_exp_df$Cd36_or_not <- factor(gene_exp_df$Cd36_or_not,levels=c("Cd36+ ORNs","Cd36- ORNs"),labels=c("Cd36+ OSNs","Cd36- OSNs"))
  Cd36_vs_other_p_val_adj <- signif(Cd36_vs_other_exp_df$p_val_adj[which(rownames(Cd36_vs_other_exp_df)==markers[i])],3)
  stat_df <- tibble::as_tibble(data.frame(group1="Cd36+ OSNs",group2="Cd36- OSNs",p.adj=Cd36_vs_other_p_val_adj,y.position=round(max(gene_exp_df$gene),1)+0.7))
  stat_df$p.signif <- ifelse(stat_df$p.adj<0.05,stat_df$p.adj,"ns")
  p <- ggplot() + 
    geom_violin(data=gene_exp_df,aes(x=Cd36_or_not,y=gene,fill=Cd36_or_not),scale="width",trim=FALSE)+
    geom_boxplot(data=gene_exp_df,aes(x=Cd36_or_not,y=gene,fill=Cd36_or_not),width=0.05, fill="white",outlier.size=0.8)+
    scale_fill_manual(values=violin_palettes)+
    labs(y="Gene expression",title=markers[i]) + 
    scale_y_continuous(labels = fmt_dcimals(0),limits=c(round(min(gene_exp_df$gene),1),round(max(gene_exp_df$gene),1)+1))+
    stat_pvalue_manual(stat_df, label = "p = {p.adj}",label.size=5) + 
    theme_classic()+
    theme(
          plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black",size=rel(1.8)),
          axis.text.y = element_text(color="black",size=rel(1.5)), 
          axis.text.x = element_text(size=rel(1.5),color="black",vjust=1),
          axis.line = element_line(colour="black",size=1.5),
          axis.ticks = element_line(),
          plot.margin=unit(c(0, 0.2, 0.7, 0.2), "cm"),
          legend.position = " none")
  p
  })

pdf(str_c(out_dir,"Mef2a_Tshz1_Lhx2_gene_expression_UMAP_VlnPlot.pdf"),width=12,height=21)
(gene_exp_UMAP_plot.ls[[1]] + gene_exp_violin_plot.ls[[1]] + plot_layout(widths=c(5.5,5))) / (gene_exp_UMAP_plot.ls[[2]] + gene_exp_violin_plot.ls[[2]]+plot_layout(widths=c(5.5,5))) / (gene_exp_UMAP_plot.ls[[3]] + gene_exp_violin_plot.ls[[3]] + plot_layout(widths=c(5.5,5))) + plot_annotation(title="Gene expression of TFs",theme = theme(plot.title = element_text(hjust = 0.5,size = rel(2.8))))
dev.off()


## Figure 3F,H,J - ChromVar Z-score for Mef2a, Lhx2 and Tshz1 (UMAP + Violin plot)
motifs <- c("MA0052.4","MA0700.2","ENSMUSG00000046982-Tshz1-02",)
motif_names <- c("Mef2a","Lhx2","Tshz1")

chromvar_UMAP_plot.ls <- lapply(1:length(motifs),function(i){
  chromvar_df <- as.data.frame(filtered_mature_ORNs_scATAC@assays$chromvar@data) %>%
    tibble::rownames_to_column("motif")  %>%
    dplyr::filter(motif %in% motifs[i])
  chromvar_df <- reshape2::melt(chromvar_df,id="motif",variable.name="cell",value.name ="chromvar_zscore")
  meta_df <- FetchData(filtered_mature_ORNs_scATAC,vars=c("UMAP_1","UMAP_2"))
  chromvar_df <- merge(chromvar_df,meta_df,by.x="cell",by.y="row.names")
  p <- ggplot(data=chromvar_df,aes(x=UMAP_1,y=UMAP_2,color=chromvar_zscore))+
    geom_point(size=0.4)+
    scale_colour_gradientn(colours = as.vector(ArchRPalettes[["solarExtra"]]),limits=c(-1.5,1.5),oob = scales::squish,breaks=c(-1.5,0,1.5),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1,title.position="top", title.hjust = 0.5)) + 
    labs(title=motif_names[i],x="UMAP 1",y="UMAP 2",color="ChromVAR z-score")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_text(hjust = 0.3,size=rel(3),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_text(size=rel(1.5),vjust=0.5),axis.title.y=element_text(size=rel(1.5),vjust=0.5),legend.direction="horizontal",legend.position=c(0.8,0.1),legend.key.size = unit(0.65, 'cm'),legend.title=element_text(size=rel(1.1))) 
  p
  })

DefaultAssay(filtered_mature_ORNs_scATAC) <- "chromvar"
newpalette <- c(brewer.pal(9,"Greens")[5],brewer.pal(9,"Pastel1")[9])
chromvar_violin_plot.ls <- lapply(1:length(motifs),function(i){
  chromvar_df <- FetchData(filtered_mature_ORNs_scATAC,vars=c(motifs[i],"Cd36_or_not"))
  colnames(chromvar_df)[1] <- "motif"
  chromvar_df$Cd36_or_not <- factor(chromvar_df$Cd36_or_not,levels=c("Cd36+ ORNs","Cd36- ORNs"),labels=c("Cd36+ OSNs","Cd36- OSNs"))
  Cd36_vs_other_p_val_adj <- signif(Cd36_vs_other_differential.activity$p_val_adj[which(Cd36_vs_other_differential.activity$Row.names==motifs[i])],3)
  stat_df <- tibble::as_tibble(data.frame(group1="Cd36+ OSNs",group2="Cd36- OSNs",p.adj=Cd36_vs_other_p_val_adj,y.position=round(max(chromvar_df$motif),1)+1.5))
  p <- ggplot() + 
    geom_violin(data=chromvar_df,aes(x=Cd36_or_not,y=motif,fill=Cd36_or_not),scale="width",trim=FALSE)+
    geom_boxplot(data=chromvar_df,aes(x=Cd36_or_not,y=motif,fill=Cd36_or_not),width=0.05, fill="white",outlier.size=0.8)+
    scale_fill_manual(values=newpalette)+
    scale_y_continuous(labels = fmt_dcimals(1),limits=c(round(min(chromvar_df$motif),1),round(max(chromvar_df$motif),1)+2)) + 
    labs(y="ChromVAR z-score",title=motif_names[i]) +
    stat_pvalue_manual(stat_df, label = "p = {p.adj}",label.size=4.5) +
    theme_classic()+
    theme(
          plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color="black",size=rel(1.8)),
          axis.text.y = element_text(color="black",size=rel(1.5)), 
          axis.text.x = element_text(size=rel(1.5),color="black",vjust=1),
          axis.line = element_line(colour="black",size=1.5),
          axis.ticks = element_line(),
          plot.margin=unit(c(0, 0.2, 0.7, 0.2), "cm"),
          legend.position = " none") 
  p 
  }) 


pdf(str_c(out_dir,"Mef2a_Tshz1_Lhx2_chromVAR_zscore_UMAP_VlnPlot.pdf"),width=12,,height=21)
(chromvar_UMAP_plot.ls[[1]] + chromvar_violin_plot.ls[[1]]+plot_layout(widths=c(5.5,5))) / (chromvar_UMAP_plot.ls[[2]] + chromvar_violin_plot.ls[[2]] + plot_layout(widths=c(5.5,5))) / (chromvar_UMAP_plot.ls[[3]] + chromvar_violin_plot.ls[[3]]+plot_layout(widths=c(5.5,5))) + plot_annotation(title="Accessibility of cis-elements",theme = theme(plot.title = element_text(hjust = 0.5,size = rel(2.8))))
dev.off()
