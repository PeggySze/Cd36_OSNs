## Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)

## Creat output directory 
out_dir <- "~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

## Load the dataset
integrated_scRNA <- readRDS("~/Cd36_OSNs/output/scRNA/seurat/integrated_analysis/integrated_OE_scRNA.rds")


## subset 19d nosort scRNA-seq data
nonsort_scRNA <- subset(integrated_scRNA,subset=orig.ident=="19d_rep2"& cell_type %in% c("HBCs","GBCs","INPs","Immature ORNs","Mature ORNs" ))
nonsort_scRNA <- DietSeurat(
  object=nonsort_scRNA,
  counts=TRUE,
  data=TRUE,
  scale.data=FALSE,
  assays="RNA"
  )

# perform standard preprocessing (log-normalization),and and identify variable features
nonsort_scRNA <- NormalizeData(nonsort_scRNA,verbose = FALSE)
nonsort_scRNA <- FindVariableFeatures(nonsort_scRNA,selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Scaling the data
nonsort_scRNA <- ScaleData(nonsort_scRNA, features =rownames(nonsort_scRNA))

# Perform linear dimensional reduction
nonsort_scRNA <- RunPCA(nonsort_scRNA, npcs = 50, verbose = FALSE)
pdf("/data/R03/shipy3/Projects/mouse_ORs/output/scRNA/HBC_lineage_trajectory/nonsort_HBC_lineage_scRNA_pc.pdf")
# ‘Elbow plot’
ElbowPlot(nonsort_scRNA,ndims=50)
dev.off()

# tune the number of PCs
markers <- c(
  "Syt1", # neurons
  "Omp","Stoml3","Cnga2", "Adcy3",# Mature ORNs
  "Nqo1","Acsm4",
  "Ncam2","Nfix","Nfib","Bcl11b",
  "Cd36","Tshz1",
  "Gap43","Gng8",#Immature ORNs
  "Sox11","Neurod1","Neurog1",#INPs,immediate neuronal precursors Neurog1
  "Ascl1","Kit" ,#GBCs
  "Cebpd","Krt5","Trp63","Krt14",#HBCs
  "Sox2","Ermn","Cyp2g1","Cyp1a2",#支持细胞  Sustentacular cells
  "Gucy2d","Gucy1b2")
for ( nPCs in seq(25,40,5)){
  nonsort_scRNA <- FindNeighbors(nonsort_scRNA, dims = 1:nPCs)
  nonsort_scRNA <- FindClusters(nonsort_scRNA, resolution = 0.2)
  nonsort_scRNA <- RunUMAP(nonsort_scRNA, dims = 1:nPCs)
  pdf(str_c(out_dir,"nonsort_HBC_lineage_scRNA_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p <- DimPlot(nonsort_scRNA, reduction = "umap",label=TRUE)
  print(p)
  dev.off()
  pdf(str_c(out_dir,"nonsort_HBC_lineage_scRNA_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(nonsort_scRNA,features=markers,combine=FALSE,cols=c("lightgrey", "red"),order=TRUE)
  print(p)
  dev.off()
}

# select 35 PCs
nonsort_scRNA <- FindNeighbors(nonsort_scRNA, dims = 1:35)
nonsort_scRNA <- RunUMAP(nonsort_scRNA, dims = 1:35)

# tune the resolution
for (i in seq(0.4,1.2,0.2)){
  nonsort_scRNA <- FindClusters(nonsort_scRNA, resolution = i)
  pdf(str_c(out_dir,"nonsort_HBC_lineage_scRNA_PC35_resolution",i,"_UMAP.pdf"))
  p1 <- DimPlot(nonsort_scRNA, reduction = "umap",label=TRUE,group.by="seurat_clusters")
  print(p1)
  dev.off()
}

# select 0.6 resolution
nonsort_scRNA <- FindClusters(nonsort_scRNA, resolution = 0.6)

# Constructs a phylogenetic tree based on a distance matrix constructed in PCA space
# Reordering identity classes according to position on the tree
nonsort_scRNA <- BuildClusterTree(nonsort_scRNA,dims =1:35,reorder=TRUE,reorder.numeric=TRUE)
pdf(str_c(out_dir,"nonsort_HBC_lineage_scRNA_PC35_resolution0.6_ClusterTree.pdf"))
PlotClusterTree(nonsort_scRNA,use.edge.length=FALSE)
dev.off()

pdf(str_c(out_dir,"nonsort_HBC_lineage_scRNA_PC35_resolution0.6_reorder_UMAP.pdf"))
DimPlot(nonsort_scRNA, reduction = "umap",label=TRUE)
dev.off()

pdf(str_c(out_dir,"nonsort_HBC_lineage_scRNA_metric_for_each_cluster.pdf"),width=12)
VlnPlot(nonsort_scRNA,features =c("nCount_RNA","nFeature_RNA","percent.mt"),combine=FALSE)
dev.off()

# markers for each cluster
clusters <- levels(nonsort_scRNA)
cluster_markers.ls <- lapply(1:length(clusters),function(i){
  df <- FindMarkers(nonsort_scRNA,ident.1=clusters[i],only.pos = TRUE)
  df <- df %>%
    filter(p_val_adj<0.05) %>%
    arrange(desc(avg_log2FC))
  df
  })

# annotate cell types
# 1 : HBCs
# 2 : 中性粒细胞 Neutrophils
# 3 : GBCs
# 4,9 : erythrocyte
# 14,15 : Ighm+ cells
new.cluster.ids <- c("HBCs","Neutrophils","GBCs","Erythrocytes","INPs","Immature ORNs","Immature ORNs","Non-Cd36+ ORNs","Erythrocytes","Immature ORNs","Cd36+ ORNs","Non-Cd36+ ORNs","Non-Cd36+ ORNs","Ighm+ cells","Ighm+ cells")
names(new.cluster.ids) <- levels(nonsort_scRNA)
nonsort_scRNA <- RenameIdents(nonsort_scRNA, new.cluster.ids)
cell_type_df <- data.frame(row.names=names(Idents(nonsort_scRNA)),cell_type=as.vector(Idents(nonsort_scRNA)))
nonsort_scRNA <- AddMetaData(nonsort_scRNA,cell_type_df)

pdf(str_c(out_dir,"nonsort_HBC_lineage_scRNA_UMAP_with_cell_type.pdf"))
DimPlot(nonsort_scRNA, reduction = "umap",label=TRUE)+NoLegend()
dev.off()


# retain HBC lineage
filtered_nonsort_scRNA <- subset(nonsort_scRNA,subset=cell_type %in% c("HBCs","GBCs","INPs","Immature ORNs","Cd36+ ORNs","Non-Cd36+ ORNs"))
filtered_nonsort_scRNA <-  DietSeurat(
  object=filtered_nonsort_scRNA,
  counts=TRUE,
  data=TRUE,
  scale.data=FALSE,
  assays="RNA"
  )
filtered_nonsort_scRNA <- NormalizeData(filtered_nonsort_scRNA,verbose = FALSE)
filtered_nonsort_scRNA <- FindVariableFeatures(filtered_nonsort_scRNA,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
filtered_nonsort_scRNA <- ScaleData(filtered_nonsort_scRNA, features =rownames(filtered_nonsort_scRNA))

filtered_nonsort_scRNA <- RunPCA(filtered_nonsort_scRNA, npcs = 50, verbose = FALSE)
pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_pc.pdf"))
# ‘Elbow plot’
ElbowPlot(filtered_nonsort_scRNA,ndims=50)
dev.off()


# tune the number of PCs
for ( nPCs in seq(25,40,5)){
  filtered_nonsort_scRNA <- FindNeighbors(filtered_nonsort_scRNA, dims = 1:nPCs)
  filtered_nonsort_scRNA <- FindClusters(filtered_nonsort_scRNA, resolution = 0.2)
  filtered_nonsort_scRNA <- RunUMAP(filtered_nonsort_scRNA, dims = 1:nPCs)
  pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p <- DimPlot(filtered_nonsort_scRNA, reduction = "umap",label=TRUE)
  print(p)
  dev.off()
  pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(filtered_nonsort_scRNA,features=markers,combine=FALSE,cols=c("lightgrey", "red"),order=TRUE)
  print(p)
  dev.off()
}


# select 40 PCs
filtered_nonsort_scRNA <- FindNeighbors(filtered_nonsort_scRNA, dims = 1:40)
filtered_nonsort_scRNA <- RunUMAP(filtered_nonsort_scRNA, dims = 1:40)
pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_PC40_markers_UMAP.pdf"))
FeaturePlot(filtered_nonsort_scRNA,features=markers,combine=FALSE,cols=c("lightgrey", "red"),order=TRUE)
dev.off()

# tune the resolution
for (i in seq(1.4,2,0.2)){
  filtered_nonsort_scRNA <- FindClusters(filtered_nonsort_scRNA, resolution = i)
  pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_PC40_resolution",i,"_UMAP.pdf"))
  p1 <- DimPlot(filtered_nonsort_scRNA, reduction = "umap",label=TRUE,group.by="seurat_clusters")
  print(p1)
  dev.off()
}

# select 0.6 resolution
filtered_nonsort_scRNA <- FindClusters(filtered_nonsort_scRNA, resolution = 0.6)

# Constructs a phylogenetic tree based on a distance matrix constructed in PCA space
# Reordering identity classes according to position on the tree
filtered_nonsort_scRNA <- BuildClusterTree(filtered_nonsort_scRNA,dims =1:40,reorder=TRUE,reorder.numeric=TRUE)
pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_PC40_resolution0.6_ClusterTree.pdf"))
PlotClusterTree(filtered_nonsort_scRNA,use.edge.length=FALSE)
dev.off()

pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_PC40_resolution0.6_reorder_UMAP.pdf"))
DimPlot(filtered_nonsort_scRNA, reduction = "umap",label=TRUE)
dev.off()

# annotate cell types
new.cluster.ids <- c("HBCs","GBCs","Non-Cd36+ ORNs","Nearly mature ORNs","Cd36+ ORNs","Non-Cd36+ ORNs","Non-Cd36+ ORNs","INPs","Immature ORNs","Immature ORNs")
names(new.cluster.ids) <- levels(filtered_nonsort_scRNA)
filtered_nonsort_scRNA <- RenameIdents(filtered_nonsort_scRNA, new.cluster.ids)
cell_type_df <- data.frame(row.names=names(Idents(filtered_nonsort_scRNA)),detail_cell_type=as.vector(Idents(filtered_nonsort_scRNA)))
filtered_nonsort_scRNA <- AddMetaData(filtered_nonsort_scRNA,cell_type_df)


pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_UMAP_with_cell_type.pdf"))
DimPlot(filtered_nonsort_scRNA, reduction = "umap",label=TRUE)+NoLegend()
dev.off()


## Figure S7A - UMAP with clusters
library(ArchR)
newpalette <- ArchRPalettes$stallion[c(1,19,3,7,5,6,16,9,8,13)]
Idents(filtered_nonsort_scRNA) <- filtered_nonsort_scRNA$tree.ident
Idents(filtered_nonsort_scRNA) <- factor(Idents(filtered_nonsort_scRNA),levels=1:10)
pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_UMAP_with_clusters.pdf"),width=7.5,height=7)
DimPlot(filtered_nonsort_scRNA, reduction = "umap",label=TRUE,label.size=6,cols=as.vector(newpalette),repel=FALSE,pt.size=1.2)+ NoAxes()+
  labs(title="Neuronal lineage scRNA-seq")+
  geom_segment(aes(x = -12, y = -7, xend = -9.5, yend = -7),arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = -12, y = -7, xend = -12, yend = -4),arrow = arrow(length = unit(0.3, "cm")))+
  annotate(geom = "text", x = -10.5, y = -7.5, label = "UMAP_1", color = "black") +
  annotate(geom = "text", x = -12.5, y = -5.3, label = "UMAP_2", color = "black",angle = 90) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),legend.position="none")
dev.off()


## Figure S7C - Gene expression of Omp, Adcy3
markers <- c("Omp","Adcy3")
gene_exp_palettes <- as.vector(ArchRPalettes[["solarExtra"]])
gene_exp_UMAP_plot.ls <- lapply(1:length(markers),function(i){
  gene_exp_df <- FetchData(filtered_nonsort_scRNA,vars = c(markers[i],"UMAP_1","UMAP_2"),slot = "data")
  colnames(gene_exp_df)[1] <- "gene"
  p <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.1)+
    scale_colour_gradientn(colours = gene_exp_palettes,limits=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),labels = c("0", round(quantile(gene_exp_df$gene,seq(0,1,0.01))[96],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.position = "top",title.hjust=0.5)) + 
    labs(title=markers[i],y="UMAP 2",x="UMAP 1",color="Gene expression") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_text(hjust = 0.5,size=rel(2.5),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_text(size=rel(1.2),vjust=0.5),axis.title.y=element_text(size=rel(1.2),vjust=0.5),legend.position = c(0.25,0.15),plot.margin=unit(c(0, 0.2, 0.7, 0.2), "cm"),legend.direction="horizontal",legend.key.size = unit(0.5, 'cm'))
  p
  })

pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_mature_ORNs_markers_expression_UMAP_01.pdf"),height=9,width=4)
(gene_exp_UMAP_plot.ls[[1]] / gene_exp_UMAP_plot.ls[[2]])
dev.off()


## Figure S7D - UMAP with cell types
newpalette <- ArchRPalettes$stallion[c(1,19,9,8,7,2,13)]
Idents(filtered_nonsort_scRNA) <- filtered_nonsort_scRNA$detail_cell_type
Idents(filtered_nonsort_scRNA) <- factor(Idents(filtered_nonsort_scRNA),levels=c("HBCs","GBCs","INPs","Immature ORNs","Nearly mature ORNs","Cd36+ ORNs","Non-Cd36+ ORNs"),labels=c("HBCs","GBCs","INPs","Immature OSNs","Adcy3- OSNs","Cd36+ OSNs","Cd36- OSNs"))
pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_UMAP_with_cell_types.pdf"),width=6,height=6.5)
DimPlot(filtered_nonsort_scRNA, reduction = "umap",label=TRUE,label.size=5,cols=as.vector(newpalette),repel=TRUE,pt.size=1)+
  labs(title="ORNs lineage scRNA-seq")+ NoAxes()+ NoLegend()+
  geom_segment(aes(x = -12, y = -7, xend = -9.5, yend = -7),arrow = arrow(length = unit(0.3, "cm")))+
  geom_segment(aes(x = -12, y = -7, xend = -12, yend = -4.5),arrow = arrow(length = unit(0.3, "cm")))+
  annotate(geom = "text", x = -10.5, y = -7.5, label = "UMAP_1", color = "black") +
  annotate(geom = "text", x = -12.6, y = -5.9, label = "UMAP_2", color = "black",angle = 90) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"))
dev.off()

## Figure 4B - Gene expression of Tshz1 (exclude HBCs)
remove_HBCs_nonsort_scRNA <- subset(filtered_nonsort_scRNA,subset=detail_cell_type!="HBCs"& UMAP_1>(-10))
gene_exp_palettes <- as.vector(ArchRPalettes[["solarExtra"]])
gene_exp_df <- FetchData(remove_HBCs_nonsort_scRNA,vars = c("Tshz1","UMAP_1","UMAP_2"),slot = "data")
colnames(gene_exp_df)[1] <- "gene"
plot <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.3)+
    scale_colour_gradientn(colours = gene_exp_palettes,limits=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),labels = c("0", round(quantile(gene_exp_df$gene,seq(0,1,0.01))[96],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.8,title.position = "top",title.hjust=0.5,barheight=0.8)) + 
    labs(title="Tshz1",color="Gene expression") +
    geom_segment(aes(x = -9, y = -7, xend = -6.5, yend = -7),arrow = arrow(length = unit(0.3, "cm")),color="black")+
    geom_segment(aes(x = -9, y = -7, xend = -9, yend = -4.5),arrow = arrow(length = unit(0.3, "cm")),color="black")+
    annotate(geom = "text", x = -7.5, y = -7.5, label = "UMAP_1", color = "black") +
    annotate(geom = "text", x = -9.6, y = -5.9, label = "UMAP_2", color = "black",angle = 90) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),legend.key.size = unit(1, 'cm'),legend.text=element_text(size=rel(1.2)),legend.position="bottom",legend.title=element_text(size=rel(1.4)))
pdf(str_c(out_dir,"remove_HBCs_nonsort_HBC_lineage_scRNA_mature_ORNs_Tshz1_expression_UMAP.pdf"),height=6,width=4)
plot 
dev.off()

## Figure 4D-E, Gene expression of Lhx2 and Mef2a 
gene_exp_palettes <- as.vector(ArchRPalettes[["solarExtra"]])
gene_exp_df <- FetchData(remove_HBCs_nonsort_scRNA,vars = c("Lhx2","UMAP_1","UMAP_2"),slot = "data")
colnames(gene_exp_df)[1] <- "gene"
plot <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.3)+
    scale_colour_gradientn(colours = gene_exp_palettes,limits=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),labels = c("0", round(quantile(gene_exp_df$gene,seq(0,1,0.01))[96],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.8,title.position = "top",title.hjust=0.5,barheight=0.8)) + 
    labs(title="Lhx2",color="Gene expression") +
    geom_segment(aes(x = -9, y = -7, xend = -6.5, yend = -7),arrow = arrow(length = unit(0.3, "cm")),color="black")+
    geom_segment(aes(x = -9, y = -7, xend = -9, yend = -4.5),arrow = arrow(length = unit(0.3, "cm")),color="black")+
    annotate(geom = "text", x = -7.5, y = -7.5, label = "UMAP_1", color = "black") +
    annotate(geom = "text", x = -9.6, y = -5.9, label = "UMAP_2", color = "black",angle = 90) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),legend.key.size = unit(1, 'cm'),legend.text=element_text(size=rel(1.2)),legend.position="bottom",legend.title=element_text(size=rel(1.4)))
pdf("/data/R03/shipy3/Projects/mouse_ORs/output/scRNA/HBC_lineage_trajectory/remove_HBCs_nonsort_HBC_lineage_scRNA_mature_ORNs_Lhx2_expression_UMAP.pdf",height=6,width=4)
plot 
dev.off()

gene_exp_palettes <- as.vector(ArchRPalettes[["solarExtra"]])
gene_exp_df <- FetchData(remove_HBCs_nonsort_scRNA,vars = c("Mef2a","UMAP_1","UMAP_2"),slot = "data")
colnames(gene_exp_df)[1] <- "gene"
plot <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.3)+
    scale_colour_gradientn(colours = gene_exp_palettes,limits=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),labels = c("0", round(quantile(gene_exp_df$gene,seq(0,1,0.01))[96],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.8,title.position = "top",title.hjust=0.5,barheight=0.8)) + 
    labs(title="Mef2a",color="Gene expression") +
    geom_segment(aes(x = -9, y = -7, xend = -6.5, yend = -7),arrow = arrow(length = unit(0.3, "cm")),color="black")+
    geom_segment(aes(x = -9, y = -7, xend = -9, yend = -4.5),arrow = arrow(length = unit(0.3, "cm")),color="black")+
    annotate(geom = "text", x = -7.5, y = -7.5, label = "UMAP_1", color = "black") +
    annotate(geom = "text", x = -9.6, y = -5.9, label = "UMAP_2", color = "black",angle = 90) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),legend.key.size = unit(1, 'cm'),legend.text=element_text(size=rel(1.2)),legend.position="bottom",legend.title=element_text(size=rel(1.4)))
pdf("/data/R03/shipy3/Projects/mouse_ORs/output/scRNA/HBC_lineage_trajectory/remove_HBCs_nonsort_HBC_lineage_scRNA_mature_ORNs_Mef2a_expression_UMAP.pdf",height=6,width=4)
plot 
dev.off()

# save object
saveRDS(filtered_nonsort_scRNA,str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA.rds"))

