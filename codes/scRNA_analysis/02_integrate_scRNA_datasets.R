## Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(patchwork)
set.seed(100)

## Creat output directory 
out_dir <- "~/Cd36_OSNs/output/scRNA/seurat/integrated_analysis/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

## Load the dataset
samples <- c("19d_rep1","19d_rep2","8w_rep1","8w_rep2")
scRNA_dirs <- str_c("~/Cd36_OSNs/output/scRNA/cellranger/",samples,"/outs")
objList <- lapply(1:length(scRNA_dirs),function(i){
  CreateSeuratObject(counts = Read10X_h5(str_c(scRNA_dirs[i],"filtered_feature_bc_matrix.h5")),
  project = samples[i],assay = "RNA") 
})

## Quality control
for ( i in 1:length(objList)){
  objList[[i]][["percent.mt"]] = PercentageFeatureSet(objList[[i]],pattern="^mt-")
  pdf(str_c(out_dir,samples[i],"_qc_plot.pdf"),width=9)
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

# figure S1b
df <- rbind(objList[[1]]@meta.data, objList[[2]]@meta.data,objList[[3]]@meta.data,objList[[4]]@meta.data)
df$orig.ident <- factor(df$orig.ident,levels=samples)
library(colorspace)
newpalette <- c("#F39B7FFF",darken("#F39B7FFF",0.2),"#E64B35FF",darken("#E64B35FF",0.2))
pdf(out_dir,"combined_OE_scRNA_qc_violin.pdf",height=10)
p1 <- ggplot(data=df,aes(x=orig.ident,y=nCount_RNA,fill=orig.ident))+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white",outlier.size=0.8)+
  theme_classic()+
  scale_fill_manual(values=newpalette) +
  coord_cartesian(ylim=c(0,50000)) + 
  labs(y="# UMIs") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=rel(2),vjust =2),
        axis.text.y = element_text(color="black",size=rel(1.8)), 
        axis.text.x = element_blank(),
        axis.line = element_line(colour="black",size = 1),
        axis.ticks = element_line(),
        plot.margin=unit(c(0, 0.5, 0.5, 0.5), "cm"),
        legend.position = " none")
p2 <- ggplot(data=df,aes(x=orig.ident,y=nFeature_RNA,fill=orig.ident))+
  geom_violin(trim=FALSE,)+
  geom_boxplot(width=0.1, fill="white",outlier.size=0.8)+
  theme_classic()+
  scale_fill_manual(values=newpalette) +
  labs(y="# genes") + 
  geom_hline(yintercept=c(400,5000),linetype="longdash")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=rel(2),vjust =2),
        axis.text.y = element_text(color="black",size=rel(1.8)), 
        axis.text.x = element_blank(),
        axis.line = element_line(colour="black",size = 1),
        axis.ticks = element_line(),
        plot.margin=unit(c(0, 0.5, 0.5, 0.5), "cm"),
        legend.position = " none")
p3 <- ggplot(data=df,aes(x=orig.ident,y=percent.mt,fill=orig.ident))+
  geom_violin(trim=FALSE,)+
  geom_boxplot(width=0.1, fill="white",outlier.size=0.8)+
  theme_classic()+
  scale_fill_manual(values=newpalette) +
  labs(y="MT reads\npercentage") + 
  geom_hline(yintercept=10,linetype="longdash")+
  coord_cartesian(ylim=c(0,20))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=rel(2),vjust =2),
        axis.text.y = element_text(color="black",size=rel(1.8)), 
        axis.text.x = element_text(color="black",size=rel(1.8),angle =45,hjust=1),
        axis.line = element_line(colour="black",size = 1),
        axis.ticks = element_line(),
        plot.margin=unit(c(0, 0.5, 0.5, 0.5), "cm"),
        legend.position = " none")
p1 / p2 / p3 
dev.off()

## Retain cells that express more than 400 genes and less than 5000 genes and had mitochondrial content <10%
for (i in 1:length(objList)){
  objList[[i]] <- subset(objList[[i]],subset = nFeature_RNA >=400 & nFeature_RNA <= 5000 &percent.mt < 10 )
}

## Integrate single-cell RNA-seq datasets across different replicates and time points
# Simply merge Seurat objects
merged_obj <- merge(objList[[1]],
                    y=c(objList[[2]],objList[[3]],objList[[4]]),
                    add.cell.ids = samples)

# split the combined object into a list, with each dataset as an element
obj.ls <- SplitObject(merged_obj,split.by = "orig.ident")

## perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(obj.ls)){
  obj.ls[[i]] <- NormalizeData(obj.ls[[i]],verbose = FALSE)
  obj.ls[[i]] <- FindVariableFeatures(obj.ls[[i]],selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}

## identify anchors using the FindIntegrationAnchors function
anchors <- FindIntegrationAnchors(object.list = obj.ls,anchor.features = 3000,dims = 1:30)

## pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
integrated_scRNA <- IntegrateData(anchorset = anchors, dims = 1:30,features.to.integrate = rownames(obj.ls[[1]]))


## switch to integrated assay
DefaultAssay(integrated_scRNA) <- "integrated"

## scale and center features in the dataset
integrated_scRNA <- ScaleData(integrated_scRNA, features =rownames(integrated_scRNA))

## Perform linear dimensional reduction
integrated_scRNA <- RunPCA(integrated_scRNA, npcs = 50, verbose = FALSE)
pdf(str_c(out_dir,"integrated_OE_scRNA_pc.pdf"))
# ‘Elbow plot’
ElbowPlot(integrated_scRNA,ndims=50)
dev.off()

# tune the number of PCs
markers <- c("Syt1", # neurons
              "Omp",# Mature OSNs
              "Nqo1","Ncam2", # dorsal/ventral OSNs
              "Gap43","Gng8",#Immature OSNs
              "Sox11","Neurod1","Neurog1",#INPs,immediate neuronal precursors 
              "Ascl1","Kit" ,#GBCs
              "Cebpd","Krt5","Trp63","Krt14",#HBCs
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#Sustentacular cells
              "Atp1a2","Fabp7", #Ensheathing glia
              "Sox9","Sox10",#Bowman's gland
              "Pebp1","Calb2", #Periglomerular cells
              "Ascl3","Cftr",#Microvillar cells
              "Krt18","Trpm5", #Brush cells
              "Col1a1","Bglap",#Osteogenic cells
              "Eng","Sox17",#Pericytes
              "Cd37","Cd79a",#B cells
              "S100a9","S100a8",#Neutrophils
              "Hmgb2","Top2a",#Late activated neural stem cells
              "Lyz2","S100a4", #Monocytes
              "Hbb-bs","Hbb-bt",#Erythrocytes
              "Mcpt8","Ccl4", #Basophils
              "C1qa","Ms4a7"#Macrophages
)

for ( nPCs in seq(20,40,5)){
  DefaultAssay(integrated_scRNA) <- "integrated"
  integrated_scRNA <- FindNeighbors(integrated_scRNA, dims = 1:nPCs)
  integrated_scRNA <- FindClusters(integrated_scRNA, resolution = 0.2)
  integrated_scRNA <- RunUMAP(integrated_scRNA, dims = 1:nPCs)
  pdf(str_c(out_dir,"integrated_OE_scRNA_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(integrated_scRNA, reduction = "umap",label=TRUE)
  p2 <- DimPlot(integrated_scRNA, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  DefaultAssay(integrated_scRNA) <- "RNA"
  pdf(str_c(out_dir,"integrated_OE_scRNA_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(integrated_scRNA,features=markers,combine=FALSE,cols=c("lightgrey", "red"),order=TRUE)
  print(p)
  dev.off()
}

# select 40 PCs
DefaultAssay(integrated_scRNA) <- "integrated"
integrated_scRNA <- FindNeighbors(integrated_scRNA, dims = 1:40)
integrated_scRNA <- RunUMAP(integrated_scRNA, dims = 1:40)

# tune the resolution
for (i in seq(0.4,1.2,0.2)){
  integrated_scRNA <- FindClusters(integrated_scRNA, resolution = i)
  pdf(str_c(out_dir,"integrated_OE_scRNA_PC40_resolution",i,"_UMAP.pdf"))
  p1 <- DimPlot(integrated_scRNA, reduction = "umap",label=TRUE,group.by="seurat_clusters")
  p2 <- DimPlot(integrated_scRNA, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
}

# select 0.8 resolution
integrated_scRNA <- FindClusters(integrated_scRNA, resolution = 0.8)

## roughly annotate cell identity
# C0,1,2,3,4,5,10 : Mature OSNs
# C9 : Immature OSNs
# C17 : INPs/GBCs
# C23 : HBCs
# C12 : Sustentacular cells
# C22 : Ensheathing glia
# C21 : Bowman's gland
# c19,30 : Periglomerular cells
# C25 : Microvillar cells
# C28 : Brush cells
# C26 : Osteogenic cells
# C29 : Pericytes
# C6,11 : B cells
# C8,20,24 : Neutrophils
# C7,18 : Late activated neural stem cells
# C15 : Monocytes
# C13,14 : Erythrocytes
# C16 : Basophils
# C27 : Macrophages
new.cluster.ids <- c("Mature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","B cells","Late activated neural stem cells","Neutrophils","Immature OSNs","Mature OSNs","B cells","Sustentacular cells","Erythrocytes","Erythrocytes","Monocytes","Basophils","INPs/GBCs","Late activated neural stem cells","Periglomerular cells","Neutrophils","Bowman's gland","Ensheathing glia","HBCs","Neutrophils","Microvillar cells","Osteogenic cells","Macrophages","Brush cells","Pericytes","Periglomerular cells")
names(new.cluster.ids) <- levels(Idents(integrated_scRNA))
integrated_scRNA <- RenameIdents(integrated_scRNA,new.cluster.ids)
rough_cell_type_df <- data.frame(rough_cell_type=as.vector(Idents(integrated_scRNA)),row.names=names(Idents(integrated_scRNA)))
integrated_scRNA <- AddMetaData(integrated_scRNA,rough_cell_type_df)

## In order to further annotate cell types, subset HBC lineage cells for subclustering 
HBC_lineage_cells <- rownames(integrated_scRNA@meta.data)[which(integrated_scRNA$rough_cell_type %in% c("HBCs","INPs/GBCs","Immature OSNs","Mature OSNs"))]

# retain HBC lineage cells
HBC_lineage_objList <- lapply(1:length(objList),function(i){
  cells <- grep(samples[i],HBC_lineage_cells,value=TRUE)
  cells <- gsub(str_c(samples[i],"_"),"",cells)
  HBC_lineage_obj <- subset(objList[[i]],cells=cells)
  HBC_lineage_obj
})

# Simply merge Seurat objects
HBC_lineage_merged_obj <- merge(HBC_lineage_objList[[1]],
                    y=c(HBC_lineage_objList[[2]],HBC_lineage_objList[[3]],HBC_lineage_objList[[4]]),
                    add.cell.ids = samples)

# split the combined object into a list, with each dataset as an element
HBC_lineage_obj.ls <- SplitObject(HBC_lineage_merged_obj,split.by = "orig.ident")

# perform standard preprocessing (log-normalization), and identify variable features individually for each
for (i in 1:length(HBC_lineage_obj.ls)){
  HBC_lineage_obj.ls[[i]] <- NormalizeData(HBC_lineage_obj.ls[[i]],verbose = FALSE)
  HBC_lineage_obj.ls[[i]] <- FindVariableFeatures(HBC_lineage_obj.ls[[i]],selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

# identify anchors using the FindIntegrationAnchors function
HBC_lineage_anchors <- FindIntegrationAnchors(object.list = HBC_lineage_obj.ls,anchor.features = 2000,dims = 1:30)

# pass these anchors to the IntegrateData function, which returns a Seurat object that contains a new Assay and holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.
integrated_HBC_lineage_scRNA <- IntegrateData(anchorset = HBC_lineage_anchors, dims = 1:30,features.to.integrate = rownames(HBC_lineage_obj.ls[[1]]))

# switch to integrated assay
DefaultAssay(integrated_HBC_lineage_scRNA) <- "integrated"

## scale and center features in the dataset
integrated_HBC_lineage_scRNA <- ScaleData(integrated_HBC_lineage_scRNA, features =rownames(integrated_HBC_lineage_scRNA))

## Perform linear dimensional reduction
integrated_HBC_lineage_scRNA <- RunPCA(integrated_HBC_lineage_scRNA, npcs = 50, verbose = FALSE)
pdf(str_c(out_dir,"integrated_HBC_lineage_scRNA_pc.pdf")
# ‘Elbow plot’
ElbowPlot(integrated_HBC_lineage_scRNA,ndims=50)
dev.off()

# tune the number of PCs
markers <- c("Syt1", #neurons
              "Omp",# Mature OSNs
              "Nqo1","Ncam2",
              "Gap43","Gng8",#Immature OSNs
              "Sox11","Neurod1","Neurog1",#INPs,immediate neuronal precursors 
              "Ascl1","Kit" ,#GBCs
              "Cebpd","Krt5","Trp63","Krt14",#HBCs
              "Sox2","Ermn","Cyp2g1","Cyp1a2",#Sustentacular cells
              )

for ( nPCs in seq(25,40,5)){
  DefaultAssay(integrated_HBC_lineage_scRNA) <- "integrated"
  integrated_HBC_lineage_scRNA <- FindNeighbors(integrated_HBC_lineage_scRNA, dims = 1:nPCs)
  integrated_HBC_lineage_scRNA <- FindClusters(integrated_HBC_lineage_scRNA, resolution = 0.2)
  integrated_HBC_lineage_scRNA <- RunUMAP(integrated_HBC_lineage_scRNA, dims = 1:nPCs)
  pdf(str_c(out_dir,"integrated_HBC_lineage_scRNA_PC",nPCs,"_resolution0.2_UMAP.pdf"))
  p1 <- DimPlot(integrated_HBC_lineage_scRNA, reduction = "umap",label=TRUE)
  p2 <- DimPlot(integrated_HBC_lineage_scRNA, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
  DefaultAssay(integrated_HBC_lineage_scRNA) <- "RNA"
  pdf(str_c(out_dir,"integrated_HBC_lineage_scRNA_PC",nPCs,"_markers_UMAP.pdf"))
  p <- FeaturePlot(integrated_HBC_lineage_scRNA,features=markers,combine=FALSE,cols=c("lightgrey", "red"),order=TRUE)
  print(p)
  dev.off()
}

# select 35 PCs
DefaultAssay(integrated_HBC_lineage_scRNA) <- "integrated"
integrated_HBC_lineage_scRNA <- FindNeighbors(integrated_HBC_lineage_scRNA, dims = 1:35)
integrated_HBC_lineage_scRNA <- RunUMAP(integrated_HBC_lineage_scRNA, dims = 1:35)


# tune the resolution
for (i in seq(0.4,1,0.2)){
  integrated_HBC_lineage_scRNA <- FindClusters(integrated_HBC_lineage_scRNA, resolution = i)
  pdf(str_c(out_dir,"integrated_HBC_lineage_scRNA_PC35_resolution",i,"_UMAP.pdf"))
  p1 <- DimPlot(integrated_HBC_lineage_scRNA, reduction = "umap",label=TRUE,group.by="seurat_clusters")
  p2 <- DimPlot(integrated_HBC_lineage_scRNA, reduction = "umap",group.by="orig.ident")
  print(p1)
  print(p2)
  dev.off()
}

# select 0.8 resolution
integrated_HBC_lineage_scRNA <- FindClusters(integrated_HBC_lineage_scRNA, resolution = 0.8)


## further annotate cell types in the whole dataset
GBCs <- rownames(integrated_HBC_lineage_scRNA@meta.data)[which(integrated_HBC_lineage_scRNA$integrated_snn_res.0.8=="16")]
integrated_scRNA$cell_type <- integrated_scRNA$rough_cell_type
integrated_scRNA$cell_type[which(integrated_scRNA$rough_cell_type == "INPs/GBCs" & rownames(integrated_scRNA@meta.data) %in% GBCs)] <- "GBCs"
integrated_scRNA$cell_type[which(integrated_scRNA$cell_type == "INPs/GBCs")]  <- "INPs"

# figure 1b - UMAP for the integrated scRNA-seq data with cell number
Idents(integrated_scRNA) <- integrated_scRNA$cell_type
Idents(integrated_scRNA) <- factor(Idents(integrated_scRNA),levels=c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs","Sustentacular cells","Ensheathing glia","Microvillar cells","Bowman's gland","Periglomerular cells","Pericytes","Brush cells","Osteogenic cells","B cells","Late activated neural stem cells","Erythrocytes","Basophils","Neutrophils","Monocytes","Macrophages"))
newpalette <- c(brewer.pal(9,"YlGn")[7:3],brewer.pal(12,"Set3")[8],brewer.pal(9,"RdPu")[2:8],brewer.pal(9,"Blues")[c(2,4,6,8)],brewer.pal(9,"GnBu")[5],brewer.pal(9,"BuPu")[c(5,7)])
cell_number_df <- integrated_scRNA@meta.data[,"cell_type",drop=FALSE] %>%
  group_by(cell_type) %>%
  summarise(counts=n())
cell_number_df$log10_counts <- log10(cell_number_df$counts)
cell_number_df$cell_type <- factor(cell_number_df$cell_type,levels=rev(c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs","Sustentacular cells","Ensheathing glia","Microvillar cells","Bowman's gland","Periglomerular cells","Pericytes","Brush cells","Osteogenic cells","B cells","Late activated neural stem cells","Erythrocytes","Basophils","Neutrophils","Monocytes","Macrophages")))
cell_type_df <- cell_number_df[,"cell_type",drop=FALSE]
cell_type_df$x <- 1
pdf(str_c(out_dir,"integrated_OE_scRNA_UMAP_with_cell_types.pdf"),width=10.5)
p1 <- DimPlot(integrated_scRNA, reduction = "umap",label=TRUE,cols=newpalette,repel=TRUE,label.size=4)+
  labs(title="scRNA-seq",x="UMAP_1",y="UMAP_2")+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_text(size=rel(1)),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0.3), "cm"))
p2 <- ggplot(data=cell_number_df,aes(x=log10_counts,y=cell_type,fill=cell_type))+
  geom_bar(stat = "identity",width=0.8)+
  labs(x="Log10 (# cells)")+
  #scale_x_continuous(limits=c(1,),breaks = c(0,5000,25000))+
  #scale_y_discrete(position = "right")+
  scale_fill_manual(values=rev(newpalette))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0), "cm"),axis.text.x=element_text(size=rel(1),color="black"),axis.title.x=element_text(size=rel(1)))
p3 <- ggplot(data=cell_type_df,aes(x=x,y=cell_type))+
  geom_tile(aes(fill=cell_type),color="black")+
  scale_y_discrete(position = "right",expand=c(0,0))+
  scale_fill_manual(values=rev(newpalette))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0.2, 0.2), "cm"),axis.text.y=element_text(size=rel(1.5),color="black"),axis.title.y=element_blank())
p1 + p2 + p3 + plot_layout(widths=c(8,1.5,0.2))
dev.off()

# figure S1a - UMAP split by sample
newpalette <- c("#F39B7FFF",darken("#F39B7FFF",0.2),"#E64B35FF",darken("#E64B35FF",0.2))
pdf(str_c(out_dir,"integrated_OE_scRNA_split_by_sample_UMAP.pdf"),height=9,width=10)
DimPlot(integrated_scRNA, reduction = "umap",label=FALSE,cols=newpalette,group.by="orig.ident",split.by="orig.ident",ncol=2) + labs(title="scRNA-seq")+theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"))
dev.off()

# figure S1c - Dot plot showing the expression of marker genes for all cell types
markers <- c(
    "Cebpd","Krt5",#HBCs
    "Ascl1","Kit" ,#GBCs
    "Neurod1","Sox11", #INPs
    "Gng8","Gap43",#Immature OSNs
    "Omp","Syt1", #Mmature OSNs
    "Sox2","Ermn",#支持细胞  Sustentacular cells
    "Atp1a2","Fabp7", #Ensheathing glia
    "Ascl3","Cftr",#Microvillar cells
    "Sox9","Sox10",#Bowman's gland
    "Pebp1","Calb2", #球周细胞Periglomerular cells
    "Eng","Sox17",#Pericytes
    "Krt18","Trpm5", #Brush cells
    "Col1a1","Bglap",#Osteogenic cells
    "Cd37","Cd79a",#B cells
    "Hmgb2","Top2a",#Late activated neural stem cells
    "Hbb-bs","Hbb-bt",#红细胞,Erythrocytes
    "Mcpt8",#"Ccl4", #Basophils
    "S100a9","S100a8",#中性粒细胞 Neutrophils
    "Lyz2","S100a4", #Monocytes
    "C1qa","Ms4a7"#Macrophages
    )
Idents(integrated_scRNA) <- factor(Idents(integrated_scRNA),levels=rev(c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs","Sustentacular cells","Ensheathing glia","Microvillar cells","Bowman's gland","Periglomerular cells","Pericytes","Brush cells","Osteogenic cells","B cells","Late activated neural stem cells","Erythrocytes","Basophils","Neutrophils","Monocytes","Macrophages")))
DefaultAssay(integrated_scRNA) <- "RNA"
p <- DotPlot(integrated_scRNA,features=markers,cols =c("lightblue","red"))+RotatedAxis()
DotPlot_df <- p$data
pdf(str_c(out_dir,"integrated_OE_scRNA_cell_type_markers_dotplot.pdf"),width=15,height=6)
ggplot(data=DotPlot_df,aes(x=features.plot,y=id))+
  geom_point(aes(size=pct.exp,color=avg.exp.scaled),)+
  scale_size_continuous(range=c(1,6),breaks = c(0,25,50,75),labels=c("0","25","50","75"),limits=c(0,100))+
  scale_color_gradientn(colours =c(brewer.pal(9,"Blues")[2],brewer.pal(9,"Reds")[c(2,4,6)]) ,limits=c(-1,2.5),oob = scales::squish)+
  labs(x="",y="",size="Percent Expressed",color="Scaled average expression")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.text.x = element_text(angle =45,hjust=1,color='black',size=11),axis.text.y=element_text(color="black",size=11),axis.title=element_blank(),legend.position="bottom",legend.direction="horizontal")
dev.off()


# figure 1c - Dot plot showing the expression of marker genes for cell types in the neuronal lineage
markers <- rev(c(
    "Cebpd","Krt5",#HBCs
    "Ascl1","Kit" ,#GBCs
    "Neurod1","Sox11", #INPs
    "Gng8","Gap43",#Immature OSNs
    "Omp","Syt1" #Mmature OSNs
    ))
Idents(integrated_scRNA) <- factor(Idents(integrated_scRNA),levels=rev(c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs","Sustentacular cells","Ensheathing glia","Microvillar cells","Bowman's gland","Periglomerular cells","Pericytes","Brush cells","Osteogenic cells","B cells","Late activated neural stem cells","Erythrocytes","Basophils","Neutrophils","Monocytes","Macrophages")))
DefaultAssay(integrated_scRNA) <- "RNA"
p <- DotPlot(integrated_scRNA,features=markers,cols =c("lightblue","red"),idents=c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs"))+RotatedAxis()
DotPlot_df <- p$data 
DotPlot_df$id <- factor(DotPlot_df$id,levels=c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs"))
cell_type_df <- data.frame(cell_type=c("HBCs","GBCs","INPs","Immature OSNs","Mature OSNs"),y=1)
pdf(str_c(out_dir,"integrated_OE_scRNA_HBC_lineage_cell_type_markers_dotplot.pdf"),width=6.6,height=6)
p1 <- ggplot(data=DotPlot_df,aes(y=features.plot,x=id))+
  geom_point(aes(size=pct.exp,fill=avg.exp.scaled),pch=21,color="black")+
  scale_size_continuous(range=c(1,8),breaks = c(0,25,50,75),labels=c("0","25","50","75"),limits=c(0,100))+
  scale_fill_gradientn(colours =c(brewer.pal(9,"Blues")[2],brewer.pal(9,"Reds")[c(2,4,6)]) ,limits=c(-1,2.5),oob = scales::squish)+
  labs(x="",y="",size="Percent Expressed",fill="Scaled average\nexpression")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.text.x = element_text(angle =45,hjust=1,color='black',size=11),axis.text.y=element_text(color="black",size=11),legend.position="right",plot.margin=unit(c(0, 0, 0.2, 0.2), "cm"))
p2 <- ggplot(data=cell_type_df,aes(x=cell_type,y=y))+
  geom_tile(aes(fill=cell_type),color="black",height=0.6)+
  labs(y="Cell type")+
  scale_fill_manual(values=brewer.pal(9,"YlGn")[7:3])+
  scale_x_discrete(expand=c(0,0))+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),axis.ticks=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position="none",plot.margin=unit(c(0.2, 0, 0, 0.2), "cm"),axis.text.y=element_blank(),axis.title.y=element_text(size=rel(1),color="black",angle=0,hjust=1))
(p2 / p1) + plot_layout(heights=c(0.2,6))
dev.off()

## save Seurat object
saveRDS(integrated_scRNA,str_c(out_dir,"integrated_OE_scRNA.rds"))















