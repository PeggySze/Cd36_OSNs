## Load required packages
library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tibble) 
library(sctransform)
set.seed(100)

## Specify output directory 
out_dir <- "~/Cd36_OSNs/output/scRNA/seurat/integrated_analysis/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

# Load scRNA object
integrated_scRNA <- readRDS(str_c(out_dir,"integrated_OE_scRNA.rds"))

# subset mature ORNs
Idents(integrated_scRNA) <- integrated_scRNA$cell_type
mature_ORNs.scRNA <- subset(integrated_scRNA,idents="Mature ORNs")
DefaultAssay(mature_ORNs.scRNA) <- "RNA"
mature_ORNs.scRNA <- DietSeurat(
  object=mature_ORNs.scRNA,
  counts=TRUE,
  data=TRUE,
  scale.data=FALSE,
  assays="RNA"
  )

# identiﬁcation of the OR expressed in each OSN
mm10_genes <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_symbols.txt")
OR_genes <- as.data.frame(readxl::read_excel("/md01/shipy3/Projects/mouse_ORs/input/PMC7055050_mouse_OR_genes.xlsx"))[,c("Gene symbol","Ensembl gene ID")]
OR_genes <- OR_genes[-which(duplicated(OR_genes)),]
merged_OR_genes <- merge(mm10_genes,OR_genes,by.x="V2",by.y="Ensembl gene ID") # 1140 OR genes
Olfr_counts <- FetchData(object=mature_ORNs.scRNA,vars=merged_OR_genes$V1,slot= "counts")

calculate_ORNs_percentage <- function(cutoff){
  if(cutoff==0){
    each_cell_expressed_Olfr_N <- rowSums(Olfr_counts>cutoff)
  } else {
    each_cell_expressed_Olfr_N <- rowSums(Olfr_counts>=cutoff)
  }
  zero_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==0)*100
  one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==1)*100
  more_than_one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N>1)*100
  df <- data.frame(cutoff=cutoff,Olfr_N=c("0","1",">1"),percentage=c(zero_Olfr_percentage,one_Olfr_percentage,more_than_one_Olfr_percentage))
  return(df)
}

cutoffs <- c(0:8,16,32,64,128,256)
df <- c()
for (i in 1:length(cutoffs)){
  df <- calculate_ORNs_percentage(cutoffs[i]) %>% rbind(df, .)
}
df$cutoff <- factor(df$cutoff,levels=cutoffs)
df$Olfr_N <- factor(df$Olfr_N,levels=c("0","1",">1"))
newpalette <- c(brewer.pal(8,"Set2")[8],brewer.pal(9,"Blues")[c(5,3)])

# figure S2C
pdf(str_c(out_dir,"mature_ORNs_UMIs_cutoff_tunning.pdf"),width=8)
ggplot(data=df,aes(x=cutoff,y=percentage,color=Olfr_N,group=Olfr_N))+
  geom_point(size=2) +
  geom_line() +
  scale_color_manual(values=newpalette)+
  labs(x="Threshold(# of UMIs)",y="Percent of OSNs",color="# of ORs")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title=element_text(size=rel(2)),axis.text.y = element_text(color="black",size=rel(1.8)), axis.text.x = element_text(color="black",size=rel(1.8)),axis.line = element_line(colour="black",size = 1))
dev.off()

## an OR was considered expressed in any given OSN if at least 4 UMIs were detected, any OSN expressing either zero or multiple OR genes was not considered further for any downstream analyses

kept_cells <- rownames(Olfr_counts)[which(rowSums(Olfr_counts>=4)==1)]
one_Olfr_mature_ORNs.scRNA <- subset(mature_ORNs.scRNA,cells=kept_cells)

# split the integrated object into a list, with each dataset as an element
mature_ORNs_scRNA.ls <- SplitObject(one_Olfr_mature_ORNs.scRNA,split.by = "orig.ident")

# perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(mature_ORNs_scRNA.ls)){
  mature_ORNs_scRNA.ls[[i]] <- NormalizeData(mature_ORNs_scRNA.ls[[i]],verbose = FALSE)
  mature_ORNs_scRNA.ls[[i]] <- FindVariableFeatures(mature_ORNs_scRNA.ls[[i]],selection.method = "vst", nfeatures = 3000)
}

# exclude OR genes,mitochondrial genes,Malat1,Xist frow variable features
mm10_genes <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_symbols.txt")
OR_genes <- as.data.frame(readxl::read_excel("/md01/shipy3/Projects/mouse_ORs/input/PMC7055050_mouse_OR_genes.xlsx"))[,c("Gene symbol","Ensembl gene ID")]
OR_genes <- OR_genes[-which(duplicated(OR_genes)),]
merged_OR_genes <- merge(mm10_genes,OR_genes,by.x="V2",by.y="Ensembl gene ID") # 1140 OR genes

MT_genes <- grep("^mt-",rownames(mature_ORNs.scRNA),value=TRUE)

for (i in 1:length(mature_ORNs_scRNA.ls)){
  VariableFeatures(mature_ORNs_scRNA.ls[[i]]) <- VariableFeatures(mature_ORNs_scRNA.ls[[i]])[which(!VariableFeatures(mature_ORNs_scRNA.ls[[i]]) %in% c(merged_OR_genes$V1,MT_genes,"Malat1","Xist"))][1:1500]
}


# select integration features and exclude OR genes,mitochondrial genes,Malat1,Xist 
features <- SelectIntegrationFeatures(mature_ORNs_scRNA.ls,nfeatures=3000)
features <- features[which(!features %in% c(merged_OR_genes$V1,MT_genes,"Malat1","Xist"))][1:1500]


# identify anchors using the FindIntegrationAnchors function
mature_ORNs_scRNA.anchors <- FindIntegrationAnchors(object.list = mature_ORNs_scRNA.ls,anchor.features = features,dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
mature_ORNs.scRNA <- IntegrateData(anchorset = mature_ORNs_scRNA.anchors, dims = 1:30,features.to.integrate = rownames(mature_ORNs_scRNA.ls[[1]]))

# switch to integrated assay
DefaultAssay(mature_ORNs.scRNA) <- "integrated"

# scale and center features in the dataset
mature_ORNs.scRNA <- ScaleData(mature_ORNs.scRNA, features =rownames(mature_ORNs.scRNA))

# Perform linear dimensional reduction
mature_ORNs.scRNA <- RunPCA(mature_ORNs.scRNA, npcs = 50,features = VariableFeatures(object = mature_ORNs.scRNA))
# Determine the ‘dimensionality’ of the dataset
pdf(str_c(out_dir,"integrated_one_Olfr_mature_ORNs_scRNA_pc.pdf"))
ElbowPlot(mature_ORNs.scRNA,ndims=50)
dev.off()

# tune the number of PCs
for (nPCs in seq(25,40,5)){
  DefaultAssay(mature_ORNs.scRNA) <- "integrated"
  mature_ORNs.scRNA<- FindNeighbors(mature_ORNs.scRNA, dims = 1:nPCs)
  mature_ORNs.scRNA <- FindClusters(mature_ORNs.scRNA, resolution = 0.2)
  mature_ORNs.scRNA <- RunUMAP(mature_ORNs.scRNA, dims = 1:nPCs)
  pdf(str_c(out_dir,"integrated_one_Olfr_mature_ORNs_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(mature_ORNs.scRNA, reduction = "umap",label=TRUE)
  p2 <- DimPlot(mature_ORNs.scRNA, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  DefaultAssay(mature_ORNs.scRNA) <- "RNA"
  pdf(str_c(out_dir,"integrated_one_Olfr_mature_ORNs_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(mature_ORNs.scRNA,features=c("Nqo1","Acsm4","Ncam2","Nfix","Nfib","Ncam2","Bcl11b","Cd36","Cd55","Cd74","Obp1a","Obp1b","Obp2a","Obp2b","Omp","Cd9","Cartpt","Calb2","S100a5","Pcp4l1","Avp","Rprm","Fxyd2","Malat1","Gucy1b2"),combine=FALSE)
  print(p)
  dev.off()
}

# select 40 PCs
DefaultAssay(mature_ORNs.scRNA) <- "integrated"
mature_ORNs.scRNA<- FindNeighbors(mature_ORNs.scRNA, dims = 1:40)
mature_ORNs.scRNA <- FindClusters(mature_ORNs.scRNA, resolution = 0.2)
mature_ORNs.scRNA <- RunUMAP(mature_ORNs.scRNA, dims = 1:40)

mature_ORNs.scRNA$log10_nCount_RNA <- log10(mature_ORNs.scRNA$nCount_RNA)
pdf(str_c(out_dir,"integrated_one_Olfr_mature_ORNs_PC40_QC_UMAP.pdf"))
FeaturePlot(mature_ORNs.scRNA,features=c("log10_nCount_RNA","nFeature_RNA","percent.mt"),combine=FALSE,cols=c("lightgrey","red"))
dev.off()

pdf(str_c(out_dir,"integrated_one_Olfr_mature_ORNs_PC40_cluster_QC_VlnPlot.pdf"),height=5,width=9)
VlnPlot(mature_ORNs.scRNA,features=c("nCount_RNA","nFeature_RNA","percent.mt"),combine=FALSE)
dev.off()

# markers for each cluster
DefaultAssay(mature_ORNs.scRNA) <- "RNA"
clusters <- levels(Idents(mature_ORNs.scRNA))
cluster_markers.ls <- lapply(1:length(clusters),function(i){
  df <- FindMarkers(mature_ORNs.scRNA,ident.1=clusters[i],assay="RNA",only.pos=TRUE)
  df <- df %>% 
        dplyr::filter(p_val<0.05) %>%
        dplyr::arrange(desc(avg_log2FC))
  df
  })

saveRDS(mature_ORNs.scRNA,str_c(out_dir,"integrated_one_Olfr_mature_ORNs_scRNA.rds"))

# remove cluster 2 and 5 due to a relatively low number of expressed genes
filtered_mature_ORNs_scRNA <- subset(mature_ORNs.scRNA,idents=c("0","1","3","4"))

# subset filtered mature ORNs
DefaultAssay(filtered_mature_ORNs_scRNA) <- "RNA"
filtered_mature_ORNs_scRNA <- DietSeurat(
  object=filtered_mature_ORNs_scRNA,
  counts=TRUE,
  data=TRUE,
  scale.data=FALSE,
  assays="RNA"
  )

# split the integrated object into a list, with each dataset as an element
filtered_mature_ORNs_scRNA.ls <- SplitObject(filtered_mature_ORNs_scRNA,split.by = "orig.ident")

# perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(filtered_mature_ORNs_scRNA.ls)){
  filtered_mature_ORNs_scRNA.ls[[i]] <- NormalizeData(filtered_mature_ORNs_scRNA.ls[[i]],verbose = FALSE)
  filtered_mature_ORNs_scRNA.ls[[i]] <- FindVariableFeatures(filtered_mature_ORNs_scRNA.ls[[i]],selection.method = "vst", nfeatures = 3000)
}

# exclude OR genes,mitochondrial genes,Malat1,Xist frow variable features
mm10_genes <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_symbols.txt")
OR_genes <- as.data.frame(readxl::read_excel("/md01/shipy3/Projects/mouse_ORs/input/PMC7055050_mouse_OR_genes.xlsx"))[,c("Gene symbol","Ensembl gene ID")]
OR_genes <- OR_genes[-which(duplicated(OR_genes)),]
merged_OR_genes <- merge(mm10_genes,OR_genes,by.x="V2",by.y="Ensembl gene ID") # 1140 OR genes

MT_genes <- grep("^mt-",rownames(filtered_mature_ORNs_scRNA.ls[[1]]),value=TRUE)

for (i in 1:length(filtered_mature_ORNs_scRNA.ls)){
  VariableFeatures(filtered_mature_ORNs_scRNA.ls[[i]]) <- VariableFeatures(filtered_mature_ORNs_scRNA.ls[[i]])[which(!VariableFeatures(filtered_mature_ORNs_scRNA.ls[[i]]) %in% c(merged_OR_genes$V1,MT_genes,"Malat1","Xist"))][1:1500]
}


# select integration features and exclude OR genes,mitochondrial genes,Malat1,Xist 
features <- SelectIntegrationFeatures(filtered_mature_ORNs_scRNA.ls,nfeatures=3000)
features <- features[which(!features %in% c(merged_OR_genes$V1,MT_genes,"Malat1","Xist"))][1:1500]


# identify anchors using the FindIntegrationAnchors function
filtered_mature_ORNs_scRNA.anchors <- FindIntegrationAnchors(object.list = filtered_mature_ORNs_scRNA.ls,anchor.features = features,dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
filtered_mature_ORNs.scRNA <- IntegrateData(anchorset = filtered_mature_ORNs_scRNA.anchors, dims = 1:30,features.to.integrate = rownames(filtered_mature_ORNs_scRNA.ls[[1]]))

# switch to integrated assay
DefaultAssay(filtered_mature_ORNs.scRNA) <- "integrated"

# scale and center features in the dataset
filtered_mature_ORNs.scRNA <- ScaleData(filtered_mature_ORNs.scRNA, features =rownames(filtered_mature_ORNs.scRNA))

# Perform linear dimensional reduction
filtered_mature_ORNs.scRNA <- RunPCA(filtered_mature_ORNs.scRNA, npcs = 50,features = VariableFeatures(object = filtered_mature_ORNs.scRNA))
# Determine the ‘dimensionality’ of the dataset
pdf(str_c(out_dir,"integrated_filtered_mature_ORNs_scRNA_pc.pdf"))
ElbowPlot(filtered_mature_ORNs.scRNA,ndims=50)
dev.off()

# tune the number of PCs
for (nPCs in seq(20,40,5)){
  DefaultAssay(filtered_mature_ORNs.scRNA) <- "integrated"
  filtered_mature_ORNs.scRNA<- FindNeighbors(filtered_mature_ORNs.scRNA, dims = 1:nPCs)
  filtered_mature_ORNs.scRNA <- FindClusters(filtered_mature_ORNs.scRNA, resolution = 0.2)
  filtered_mature_ORNs.scRNA <- RunUMAP(filtered_mature_ORNs.scRNA, dims = 1:nPCs)
  pdf(str_c(out_dir,"integrated_filtered_mature_ORNs_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(filtered_mature_ORNs.scRNA, reduction = "umap",label=TRUE)
  p2 <- DimPlot(filtered_mature_ORNs.scRNA, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  DefaultAssay(filtered_mature_ORNs.scRNA) <- "RNA"
  pdf(str_c(out_dir,"integrated_filtered_mature_ORNs_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(filtered_mature_ORNs.scRNA,features=c("Nqo1","Acsm4","Ncam2","Nfix","Nfib","Ncam2","Bcl11b","Cd36","Cd55","Cd74","Obp1a","Obp1b","Obp2a","Obp2b","Omp","Cd9","Cartpt","Calb2","S100a5","Pcp4l1","Avp","Rprm","Fxyd2","Malat1","Gucy1b2","Gucy2d","Hbb-bs","Hbb-bt","Ighm","Igkc"),combine=FALSE,order=TRUE)
  print(p)
  dev.off()
}

# select 20 PCs
DefaultAssay(filtered_mature_ORNs.scRNA) <- "integrated"
filtered_mature_ORNs.scRNA<- FindNeighbors(filtered_mature_ORNs.scRNA, dims = 1:20)
filtered_mature_ORNs.scRNA <- RunUMAP(filtered_mature_ORNs.scRNA, dims = 1:20)

# tune the resolution
for (i in seq(0.4,1,0.2)){
  filtered_mature_ORNs.scRNA <- FindClusters(filtered_mature_ORNs.scRNA, resolution = i)
  pdf(str_c(out_dir,"integrated_filtered_mature_ORNs_PC20_resolution",i,"_UMAP.pdf"))
  p1 <- DimPlot(filtered_mature_ORNs.scRNA, reduction = "umap",label=TRUE)
  p2 <- DimPlot(filtered_mature_ORNs.scRNA, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
}

# select 1 resolution
filtered_mature_ORNs.scRNA <- FindClusters(filtered_mature_ORNs.scRNA, resolution = 1)

pdf(str_c(out_dir,"integrated_filtered_mature_ORNs_PC20_cluster_QC_VlnPlot.pdf"),height=5,width=14)
VlnPlot(filtered_mature_ORNs.scRNA,features=c("nCount_RNA","nFeature_RNA","percent.mt"),combine=FALSE)
dev.off()

# markers for each cluster
DefaultAssay(filtered_mature_ORNs.scRNA) <- "RNA"
clusters <- levels(Idents(filtered_mature_ORNs.scRNA))
cluster_markers.ls <- lapply(1:length(clusters),function(i){
  df <- FindMarkers(filtered_mature_ORNs.scRNA,ident.1=clusters[i],assay="RNA",only.pos=TRUE)
  df <- df %>% 
        dplyr::filter(p_val<0.05) %>%
        dplyr::arrange(desc(avg_log2FC))
  df
  })

# annotate cell types
new.cluster.ids <- c("Dorsal ORNs","Dorsal ORNs","Ventral ORNs","Ventral ORNs","Ventral ORNs","Dorsal ORNs","Ventral ORNs","Ventral ORNs","Ventral ORNs","Ventral ORNs","Ventral ORNs","Cd36+ ORNs","Cd36+ ORNs","Ighm+ cells")
names(new.cluster.ids) <- levels(filtered_mature_ORNs.scRNA)
filtered_mature_ORNs.scRNA <- RenameIdents(filtered_mature_ORNs.scRNA, new.cluster.ids)
cell_type_df <- data.frame(row.names=names(Idents(filtered_mature_ORNs.scRNA)),cell_type=as.vector(Idents(filtered_mature_ORNs.scRNA)))
filtered_mature_ORNs.scRNA <- AddMetaData(filtered_mature_ORNs.scRNA,cell_type_df)

saveRDS(filtered_mature_ORNs.scRNA,str_c(out_dir,"integrated_filtered_mature_ORNs_scRNA.rds"))


