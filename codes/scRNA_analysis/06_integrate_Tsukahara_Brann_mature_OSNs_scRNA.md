## Regenerate gene-barcode matrices for mouse olfactory epithelium scRNA-seq data
- Download BAM files of previously published raw scRNA-seq data (GSM5283254, GSM5283255, GSM5283256, GSM5283257, GSM5283258, GSM5283259) derived from the main olfactory epithelium of adult mice housed in a typical homecage environment reported by Tsukahara et al 
- Convert BAM files into FASTQ files by bamtofastq executable: 04_Tsukahara_Brann_bamtofastq.sh
- Align data to the Mus musculus reference genome mm10 and converted into gene-barcode matrices using 10x Genomics Cell Ranger: 05_Tsukahara_Brann_cellranger.sh

## Integrate scRNA-seq data of Tsukahara mature OSNs with our mature OSNs from this study
```R
## Load required packages
library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tibble) 
library(sctransform)

## Specify output directory 
out_dir <- "~/Cd36_OSNs/output/scRNA/seurat/integrated_analysis/integrate_Tsukahara_Brann/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

## Load the Tsukahara_Brann dataset
folders <- dir('~/Cd36_OSNs/input/Tsukahara_Brann_OSN/data/cellranger',pattern='^homecage_*',full.names = TRUE)
samples <- paste("homecage",1:6,sep="_")

# Retain cells that are regarded as mature ORNs in Tsukahara_Brann paper
Tsukahara_Brann_home_cage_metadata <- read.csv("~/Cd36_OSNs/input/Tsukahara_Brann_OSN/data/raw/GSE173947_home_cage_metadata.csv.gz",row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
orig_idents <- c("baseline-1","baseline-2","baseline-3","baseline-4","baseline-9","baseline-10")
objList <- lapply(1:length(folders),function(i){
  obj <- CreateSeuratObject(counts = Read10X_h5(str_c(folders[i],"/outs/raw_feature_bc_matrix.h5")),project = samples[i],assay = "RNA") 
  cells <- grep(orig_idents[i],rownames(Tsukahara_Brann_home_cage_metadata),value=TRUE)
  cells <- gsub(str_c(orig_idents[i],"_"),"",cells)
  cells <- str_c(cells,"-1")
  obj <- subset(obj,cells=cells)
  obj
})

## QC and filter out low-quality cells
for ( i in 1:length(objList)){
  objList[[i]][["percent.mt"]] = PercentageFeatureSet(objList[[i]],pattern="^mt-")
  pdf(str_c(out_dir,samples[i],"mature_ORNs_qc_plot.pdf"),width=9)
  p1=VlnPlot(objList[[i]],features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2=ggplot(data=objList[[i]][["nFeature_RNA"]],aes(x=nFeature_RNA))+geom_density()
  p3=ggplot(data=objList[[i]][["nCount_RNA"]],aes(x=nCount_RNA))+geom_density()
  p4=ggplot(data=objList[[i]][["percent.mt"]],aes(x=percent.mt))+geom_density()
  p5=FeatureScatter(objList[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  p6=FeatureScatter(objList[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  dev.off()
}


# Simply merge Seurat objects
merged_obj <- merge(objList[[1]],
                    y=c(objList[[2]],objList[[3]],objList[[4]],objList[[5]],objList[[6]]),
                    add.cell.ids = samples)

# split the combined object into a list, with each dataset as an element
obj.ls <- SplitObject(merged_obj,split.by = "orig.ident")

## perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(objList)){
  obj.ls[[i]] <- NormalizeData(obj.ls[[i]],verbose = FALSE)
  obj.ls[[i]] <- FindVariableFeatures(obj.ls[[i]],selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}

# exclude OR genes,mitochondrial genes,Malat1,Xist frow variable features
mm10_genes <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_symbols.txt")
OR_genes <- as.data.frame(readxl::read_excel("/md01/shipy3/Projects/mouse_ORs/input/PMC7055050_mouse_OR_genes.xlsx"))[,c("Gene symbol","Ensembl gene ID")]
OR_genes <- OR_genes[-which(duplicated(OR_genes)),]
merged_OR_genes <- merge(mm10_genes,OR_genes,by.x="V2",by.y="Ensembl gene ID") # 1140 OR genes
MT_genes <- grep("^mt-",rownames(objList[[1]]),value=TRUE)
for (i in 1:length(obj.ls)){
  VariableFeatures(obj.ls[[i]]) <- VariableFeatures(obj.ls[[i]])[which(!VariableFeatures(obj.ls[[i]]) %in% c(as.vector(merged_OR_genes$V1),MT_genes,"Malat1","Xist"))][1:1500]
}

# select integration features and exclude OR genes,mitochondrial genes,Malat1,Xist 
features <- SelectIntegrationFeatures(obj.ls,nfeatures=3000)
features <- features[which(!features %in% c(as.vector(merged_OR_genes$V1),MT_genes,"Malat1","Xist"))][1:1500]

# identify anchors using the FindIntegrationAnchors function
anchors <- FindIntegrationAnchors(object.list = obj.ls,anchor.features = features,dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
Tsukahara_Brann_mature_ORNs.scRNA <- IntegrateData(anchorset = anchors, dims = 1:30,features.to.integrate = rownames(obj.ls[[1]]))

# switch to integrated assay
DefaultAssay(Tsukahara_Brann_mature_ORNs.scRNA) <- "integrated"

# scale and center features in the dataset
Tsukahara_Brann_mature_ORNs.scRNA <- ScaleData(Tsukahara_Brann_mature_ORNs.scRNA, features =rownames(Tsukahara_Brann_mature_ORNs.scRNA))

# Perform linear dimensional reduction
Tsukahara_Brann_mature_ORNs.scRNA <- RunPCA(Tsukahara_Brann_mature_ORNs.scRNA, npcs = 50)
# Determine the ‘dimensionality’ of the dataset
pdf(str_c(out_dir,"integrated_Tsukahara_Brann_mature_ORNs_scRNA_pc.pdf"))
ElbowPlot(Tsukahara_Brann_mature_ORNs.scRNA,ndims=50)
dev.off()

# tune the number of PCs
for (nPCs in seq(25,40,5)){
  DefaultAssay(Tsukahara_Brann_mature_ORNs.scRNA) <- "integrated"
  Tsukahara_Brann_mature_ORNs.scRNA <- FindNeighbors(Tsukahara_Brann_mature_ORNs.scRNA, dims = 1:nPCs)
  Tsukahara_Brann_mature_ORNs.scRNA <- FindClusters(Tsukahara_Brann_mature_ORNs.scRNA, resolution = 0.2)
  Tsukahara_Brann_mature_ORNs.scRNA <- RunUMAP(Tsukahara_Brann_mature_ORNs.scRNA, dims = 1:nPCs)
  pdf(str_c(out_dir,"integrated_Tsukahara_Brann_mature_ORNs_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(Tsukahara_Brann_mature_ORNs.scRNA, reduction = "umap",label=TRUE)
  p2 <- DimPlot(Tsukahara_Brann_mature_ORNs.scRNA, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  DefaultAssay(Tsukahara_Brann_mature_ORNs.scRNA) <- "RNA"
  pdf(str_c(out_dir,"integrated_Tsukahara_Brann_mature_ORNs_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(Tsukahara_Brann_mature_ORNs.scRNA,features=c("Nqo1","Acsm4","Ncam2","Nfix","Nfib","Ncam2","Cd36","Cd55","Cd74","Omp","Malat1","Gucy1b2"),combine=FALSE)
  print(p)
  dev.off()
}

# select 30 PCs
DefaultAssay(Tsukahara_Brann_mature_ORNs.scRNA) <- "integrated"
Tsukahara_Brann_mature_ORNs.scRNA <- FindNeighbors(Tsukahara_Brann_mature_ORNs.scRNA, dims = 1:30)
Tsukahara_Brann_mature_ORNs.scRNA <- FindClusters(Tsukahara_Brann_mature_ORNs.scRNA, resolution = 0.2)
Tsukahara_Brann_mature_ORNs.scRNA <- RunUMAP(Tsukahara_Brann_mature_ORNs.scRNA, dims = 1:30)

# markers of each cluster 
DefaultAssay(Tsukahara_Brann_mature_ORNs.scRNA) <- "RNA"
clusters <- levels(Idents(Tsukahara_Brann_mature_ORNs.scRNA))
cluster_markers.ls <- lapply(1:length(clusters),function(i){
  df <- FindMarkers(Tsukahara_Brann_mature_ORNs.scRNA, ident.1 =clusters[i],assay="RNA",only.pos=TRUE)
  df <- df %>% 
    dplyr::filter(p_val<0.05) %>%
    dplyr::arrange(desc(avg_logFC),desc(pct.1))
  df
  })


## Integrate Tsukahara_Brann data with our data
filtered_mature_ORNs.scRNA <- readRDS("~/Cd36_OSNs/output/scRNA/seurat/integrated_analysis/integrated_filtered_mature_ORNs_scRNA.rds")
filtered_Ighm_cells_mature_ORNs.scRNA <- subset(filtered_mature_ORNs.scRNA,subset=cell_type!="Ighm+ cells")
rm(filtered_mature_ORNs.scRNA)

DefaultAssay(filtered_Ighm_cells_mature_ORNs.scRNA) <- "RNA"
filtered_Ighm_cells_mature_ORNs.scRNA <- DietSeurat(
  object=filtered_Ighm_cells_mature_ORNs.scRNA,
  counts=TRUE,
  data=TRUE,
  scale.data=FALSE,
  assays="RNA"
  )
filtered_Ighm_cells_mature_ORNs.scRNA_obj.ls <- SplitObject(filtered_Ighm_cells_mature_ORNs.scRNA,split.by = "orig.ident")
rm(filtered_Ighm_cells_mature_ORNs.scRNA)

# remove homecage_4 due to low quality
all_obj.ls <- list("homecage_1"=obj.ls[[1]],"homecage_2"=obj.ls[[2]],"homecage_3"=obj.ls[[3]],"homecage_5"=obj.ls[[5]],"homecage_6"=obj.ls[[6]])
for (i in 1:length(filtered_Ighm_cells_mature_ORNs.scRNA_obj.ls)){
  sample <- names(filtered_Ighm_cells_mature_ORNs.scRNA_obj.ls)[i]
  all_obj.ls[[sample]] <- filtered_Ighm_cells_mature_ORNs.scRNA_obj.ls[[i]]
}

# perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(all_obj.ls)){
  all_obj.ls[[i]] <- NormalizeData(all_obj.ls[[i]],verbose = FALSE)
  all_obj.ls[[i]] <- FindVariableFeatures(all_obj.ls[[i]],selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}

# exclude OR genes,mitochondrial genes,Malat1,Xist frow variable features
mm10_genes <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_symbols.txt")
OR_genes <- as.data.frame(readxl::read_excel("/md01/shipy3/Projects/mouse_ORs/input/PMC7055050_mouse_OR_genes.xlsx"))[,c("Gene symbol","Ensembl gene ID")]
OR_genes <- OR_genes[-which(duplicated(OR_genes)),]
merged_OR_genes <- merge(mm10_genes,OR_genes,by.x="V2",by.y="Ensembl gene ID") # 1140 OR genes
MT_genes <- grep("^mt-",rownames(all_obj.ls[[1]]),value=TRUE)
for (i in 1:length(all_obj.ls)){
  VariableFeatures(all_obj.ls[[i]]) <- VariableFeatures(all_obj.ls[[i]])[which(!VariableFeatures(all_obj.ls[[i]]) %in% c(as.vector(merged_OR_genes$V1),MT_genes,"Malat1","Xist"))][1:1500]
}

# select integration features and exclude OR genes,mitochondrial genes,Malat1,Xist 
features <- SelectIntegrationFeatures(all_obj.ls,nfeatures=3000)
features <- features[which(!features %in% c(as.vector(merged_OR_genes$V1),MT_genes,"Malat1","Xist"))][1:1500]

# identify anchors using the FindIntegrationAnchors function
anchors <- FindIntegrationAnchors(object.list = all_obj.ls,anchor.features = features,dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
integrated_mature_ORNs.scRNA <- IntegrateData(anchorset = anchors, dims = 1:30)

# switch to integrated assay
DefaultAssay(integrated_mature_ORNs.scRNA) <- "integrated"

# scale and center features in the dataset
integrated_mature_ORNs.scRNA <- ScaleData(integrated_mature_ORNs.scRNA, features =rownames(integrated_mature_ORNs.scRNA))

# Perform linear dimensional reduction
integrated_mature_ORNs.scRNA <- RunPCA(integrated_mature_ORNs.scRNA, npcs = 50)
# Determine the ‘dimensionality’ of the dataset
pdf(str_c(out_dir,"integrated_mature_ORNs_scRNA_pc.pdf"))
ElbowPlot(integrated_mature_ORNs.scRNA,ndims=50)
dev.off()

# tune the number of PCs
for (nPCs in seq(20,30,5)){
  DefaultAssay(integrated_mature_ORNs.scRNA) <- "integrated"
  integrated_mature_ORNs.scRNA <- FindNeighbors(integrated_mature_ORNs.scRNA, dims = 1:nPCs)
  integrated_mature_ORNs.scRNA <- FindClusters(integrated_mature_ORNs.scRNA, resolution = 0.2)
  integrated_mature_ORNs.scRNA <- RunUMAP(integrated_mature_ORNs.scRNA, dims = 1:nPCs)
  pdf(str_c(out_dir,"integrated_mature_ORNs_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(integrated_mature_ORNs.scRNA, reduction = "umap",label=TRUE)
  p2 <- DimPlot(integrated_mature_ORNs.scRNA, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  DefaultAssay(integrated_mature_ORNs.scRNA) <- "RNA"
  pdf(str_c(out_dir,"integrated_mature_ORNs_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(integrated_mature_ORNs.scRNA,features=c("Nqo1","Acsm4","Ncam2","Nfix","Nfib","Ncam2","Cd36","Cd55","Cd74","Omp","Malat1","Gucy1b2"),combine=FALSE)
  print(p)
  dev.off()
}

# select 25 PCs
DefaultAssay(integrated_mature_ORNs.scRNA) <- "integrated"
integrated_mature_ORNs.scRNA <- FindNeighbors(integrated_mature_ORNs.scRNA, dims = 1:25)
integrated_mature_ORNs.scRNA <- FindClusters(integrated_mature_ORNs.scRNA, resolution = 0.2)
integrated_mature_ORNs.scRNA <- RunUMAP(integrated_mature_ORNs.scRNA, dims = 1:25)

# markers of each cluster 
DefaultAssay(integrated_mature_ORNs.scRNA) <- "RNA"
clusters <- levels(Idents(integrated_mature_ORNs.scRNA))
cluster_markers.ls <- lapply(1:length(clusters),function(i){
  df <- FindMarkers(integrated_mature_ORNs.scRNA, ident.1 =clusters[i],assay="RNA",only.pos=TRUE)
  df <- df %>% 
    dplyr::filter(p_val<0.05) %>%
    dplyr::arrange(desc(avg_log2FC),desc(pct.1))
  df
  })

# tune the resolution
DefaultAssay(integrated_mature_ORNs.scRNA) <- "integrated"
for (i in seq(0.4,1.4,0.2)){
  integrated_mature_ORNs.scRNA <- FindClusters(integrated_mature_ORNs.scRNA, resolution = i)
  pdf(str_c(out_dir,"integrated_mature_ORNs_PC25_resolution",i,"_UMAP.pdf"))
  p1 <- DimPlot(integrated_mature_ORNs.scRNA, reduction = "umap",label=TRUE)
  p2 <- DimPlot(integrated_mature_ORNs.scRNA, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
}


# select 0.2 resolution
integrated_mature_ORNs.scRNA <- FindClusters(integrated_mature_ORNs.scRNA, resolution = 0.2)

# Constructs a phylogenetic tree based on a distance matrix constructed in PCA space
# Reordering identity classes according to position on the tree
Idents(integrated_mature_ORNs.scRNA) <- integrated_mature_ORNs.scRNA$seurat_clusters
integrated_mature_ORNs.scRNA <- BuildClusterTree(integrated_mature_ORNs.scRNA,dims =1:25,reorder=TRUE,reorder.numeric=TRUE)
pdf(str_c(out_dir,"integrated_mature_ORNs_PC25_resolution0.2_ClusterTree.pdf"))
PlotClusterTree(integrated_mature_ORNs.scRNA,use.edge.length=FALSE)
dev.off()

pdf(str_c(out_dir,"integrated_mature_ORNs_PC25_resolution0.2_reorder_UMAP.pdf"))
DimPlot(integrated_mature_ORNs.scRNA, reduction = "umap",label=TRUE)
dev.off()


# annotate each cluster
new.cluster.ids <- c("Ventral ORNs","Ifi27+ cells","Ventral ORNs","Dorsal ORNs","Ventral ORNs","Ventral ORNs","Cd36+ ORNs","Ventral ORNs")
names(new.cluster.ids) <- levels(integrated_mature_ORNs.scRNA)
integrated_mature_ORNs.scRNA <- RenameIdents(integrated_mature_ORNs.scRNA, new.cluster.ids)
cell_type_df <- data.frame(row.names=names(Idents(integrated_mature_ORNs.scRNA)),cell_type=as.vector(Idents(integrated_mature_ORNs.scRNA)))
integrated_mature_ORNs.scRNA <- AddMetaData(integrated_mature_ORNs.scRNA,cell_type_df)

saveRDS(integrated_mature_ORNs.scRNA,str_c(out_dir,"remove_homecage_4_integrated_mature_ORNs_Seurat_obj.rds"))
```