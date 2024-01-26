## Generate spliced and unspliced matrices by Velocyto
- script : 01_run_RNA_Velocity.sh


## Convert data from Seurat to Python / anndata
- reference : https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
```R
## Load required packages
library(Seurat)
library(SeuratDisk)

# Load 19d nosort data 
filtered_nonsort_scRNA <- readRDS("~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/filtered_nonsort_HBC_lineage_scRNA.rds")
filtered_nonsort_scRNA$detail_cell_type <- gsub("ORNs","OSNs",filtered_nonsort_scRNA$detail_cell_type)
filtered_nonsort_scRNA$detail_cell_type[which(filtered_nonsort_scRNA$detail_cell_type=="Non-Cd36+ OSNs")] <- "Cd36- OSNs"
filtered_nonsort_scRNA$detail_cell_type[which(filtered_nonsort_scRNA$detail_cell_type=="Nearly mature OSNs")] <- "Adcy3- OSNs"

# rename barcode
new_names <- sapply(1:length(colnames(filtered_nonsort_scRNA)),function(i){
	unlist(strsplit(colnames(filtered_nonsort_scRNA)[i],"_"))[3]
	})
filtered_nonsort_scRNA <- RenameCells(filtered_nonsort_scRNA,new.names =new_names)

# remove HBCs
remove_HBCs_nonsort_scRNA <- subset(filtered_nonsort_scRNA,subset=detail_cell_type!="HBCs"& UMAP_1>(-10))

SaveH5Seurat(remove_HBCs_nonsort_scRNA, filename = "~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/remove_HBCs_nonsort_HBC_lineage_scRNA.h5Seurat",overwrite=TRUE)
Convert("~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/remove_HBCs_nonsort_HBC_lineage_scRNA.h5Seurat", dest = "h5ad",overwrite=TRUE)
```

## Estimating RNA velocity by scvelo
- Jupyter Notebook: 
```python
# import modules
import numpy as np
import pandas as pd
import anndata
import scvelo as scv
import scanpy as sc
import os

scv.settings.verbosity = 3 # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False) # for beautified visualization

## set working directory
os.chdir("~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/scvelo/")

## Load data
adata = sc.read_h5ad("~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/remove_HBCs_nonsort_HBC_lineage_scRNA.h5ad")

adata.obs['detail_cell_type'] = adata.obs['detail_cell_type'].astype("category")
adata.obs['detail_cell_type'].cat.reorder_categories(["GBCs","INPs","Immature OSNs","Adcy3- OSNs","Cd36+ OSNs","Cd36- OSNs"], inplace=True)
adata.obs['tree.ident'] = adata.obs['tree.ident'].astype("category")

## Figure S7B 
marker_genes_dict = {
    'HBCs': ["Cebpd","Krt5"],
    'GBCs': ["Ascl1","Kit"],
    'INPs': ["Neurod1","Neurog1"],
    'Immature OSNs' : ["Gng8","Gap43"],
    'Mature OSNs' : ["Omp","Adcy3"],
    'Cd36+ OSNs' : ["Cd36"],                     
}

sc.pl.dotplot(adata, marker_genes_dict, 'tree.ident', dendrogram=False,categories_order=[1,2,8,9,10,4,6,7,3,5],save="filtered_nonsort_HBC_lineage_scRNA_cluster_markers_dotplot.pdf",figsize=(7,4))

# load loom files for spliced/unspliced matrix
ldata=scv.read("~/Cd36_OSNs/output/scRNA/Velocyto/19d_nosort/possorted_genome_bam_RFKAB.loom",cache=True)

# rename barcodes in order to merge
barcodes = [bc.split(':')[1] for bc in ldata.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1' for bc in barcodes]
ldata.obs.index = barcodes

# make variable names unique
ldata.var_names_make_unique()

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

## Computing RNA velocity using scVelo
# Plot pie chart of spliced/unspliced proprtions
scv.pl.proportions(adata, groupby='detail_cell_type')
# the proportions of unspliced counts are higher than expected

## Preprocess the Data
# Preprocessing requisites consist of gene selection by detection (with a minimum number of counts) and high variability (dispersion), normalizing every cell by its total size and logarithmizing X.
# All of this is summarized in a single function scv.pp.filter_and_normalize
scv.pp.filter_and_normalize(adata)

## Computes moments for velocity estimation
# First-/second-order moments are computed for each cell across its nearest neighbors, where the neighbor graph is obtained from euclidean distances in PCA space
scv.pp.moments(adata)

## Estimate RNA velocity
# Positive velocity indicates that a gene is up-regulated, which occurs for cells that show higher abundance of unspliced mRNA for that gene than expected in steady state. Conversely, negative velocity indicates that a gene is down-regulated.
scv.tl.velocity(adata, mode='stochastic')

# Computes velocity graph based on cosine similarities
scv.tl.velocity_graph(adata)

## Visualization
# Figure 4A - visualize as streamlines
scv.pl.velocity_embedding_stream(adata, basis='umap',color='detail_cell_type',density=3,figsize=(5.5,7),title="RNA velocities",fontsize=18,legend_fontsize=13,save="remove_HBCs_nonsort_HBC_lineage_scRNA_velocity_embedding_stream_with_cell_types.pdf")


## Downstream analysis
## Identify important genes
# test which genes have cluster-specific differential velocity expression, being siginificantly higher/lower compared to the remaining population.
scv.tl.rank_velocity_genes(adata, groupby='detail_cell_type', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head(15)

## Speed: length of the velocity vector
## Coherence: how well a velocity vector correlates to its neighbors
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

df = adata.obs.groupby('detail_cell_type')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)


# visualize the velocity graph to portray all velocity-inferred cell-to-cell connections/transitions
scv.pl.velocity_graph(adata, threshold=.1,color='detail_cell_type')

# based on the velocity graph, a velocity pseudotime can be computed
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')

df = adata.obs
df.to_csv("remove_HBCs_nonsort_HBC_lineage_scRNA_velocity_pseudotime.csv",sep='\t', encoding='utf-8')

## PAGA velocity graph
# PAGA provides a graph-like map of the data topology with weighted edges corresponding to the connectivity between two clusters

# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
scv.tl.paga(adata, groups='detail_cell_type')

# Figure S7E
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5,figsize=(5,7),title="PAGA velocity graph",fontsize=20,legend_loc="right margin",save="remove_HBCs_nonsort_HBC_lineage_scRNA_paga_velocity_graph.png")
```

## Pseudotime analysis
```R
## Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)


## Specify output directory 
out_dir <- "~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}


## Load data
filtered_nonsort_scRNA <- readRDS(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA.rds"))


## split nearly mature ORNs (Adcy3- OSNs) into Tshz1+ cells and Tshz1- cells
cells <- rownames(filtered_nonsort_scRNA@meta.data)[which(filtered_nonsort_scRNA$detail_cell_type %in% c("Nearly mature ORNs"))]
Tshz1_exp_df <- FetchData(object=filtered_nonsort_scRNA,vars="Tshz1",cells=cells,slot="data")

## Figure S7F
pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_nearly_mature_ORNs_Tshz1_expression_distribution.pdf"))
ggplot(data=Tshz1_exp_df,aes(x=Tshz1))+
  geom_density() + 
  labs(x="Tshz1 expression",y="Density",title="Tshz1 expression in nearly mature ORNs") +
  geom_vline(xintercept=quantile(Tshz1_exp_df$Tshz1,seq(0,1,0.05))[19],linetype="dashed",size=0.5,color="black")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title=element_text(size=rel(2),),axis.text.y = element_text(color="black",size=rel(1.8)), axis.text.x = element_text(color="black",size=rel(1.8)),axis.line = element_line(colour="black",size = 1),plot.title = element_text(hjust = 0.3,size=rel(2),face="bold"))
dev.off()


Tshz1_cells <- rownames(Tshz1_exp_df)[which(Tshz1_exp_df$Tshz1>= quantile(Tshz1_exp_df$Tshz1,seq(0,1,0.05))[19])]
non_Tshz1_cells <- rownames(Tshz1_exp_df)[which(Tshz1_exp_df$Tshz1< quantile(Tshz1_exp_df$Tshz1,seq(0,1,0.05))[19])]


## Figure 4C - velocity pseudotime for Cd36+ ORNs and Non-Cd36+ ORNs lineage
pseudotime_df <- read.table("~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/scvelo/remove_HBCs_nonsort_HBC_lineage_scRNA_velocity_pseudotime.csv",sep="\t",header=TRUE)
pseudotime_df$X <- str_c(pseudotime_df$orig.ident,"_",pseudotime_df$X,"-1")
remove_HBCs_nonsort_scRNA <- subset(filtered_nonsort_scRNA,subset=detail_cell_type!="HBCs"& UMAP_1>(-10))
UMAP_df <- FetchData(remove_HBCs_nonsort_scRNA,vars=c("UMAP_1","UMAP_2"))
pseudotime_df <- merge(pseudotime_df,UMAP_df,by.x="X",by.y="row.names")
newpalette <- as.vector(ArchRPalettes[["horizon"]])[1:7]
Cd36_pseudotime_df <- pseudotime_df
Cd36_pseudotime_df$velocity_pseudotime[which(Cd36_pseudotime_df$X %in% non_Tshz1_cells | Cd36_pseudotime_df$detail_cell_type == "Cd36- OSNs")] <- NA
Non_Cd36_pseudotime_df <- pseudotime_df
Non_Cd36_pseudotime_df$velocity_pseudotime[which(Cd36_pseudotime_df$X %in% Tshz1_cells | Cd36_pseudotime_df$detail_cell_type == "Cd36+ OSNs")] <- NA

pdf(str_c(out_dir,"remove_HBCs_nonsort_HBC_lineage_scRNA_lineage_pseudotime_UMAP.pdf"),height=14,width=7.5)
p1 <- ggplot(data=Cd36_pseudotime_df,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=velocity_pseudotime),size=1)+
  scale_colour_gradientn(colours=newpalette,na.value="lightgrey",breaks=seq(0,1,0.2),limits=c(0,1),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.vjust=0.5,title.position = "left"))+
  labs(title="Trajectory of Cd36+ OSNs",color="Pseudotime")+
  geom_segment(aes(x = -9, y = -7, xend = -6.5, yend = -7),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  geom_segment(aes(x = -9, y = -7, xend = -9, yend = -4.5),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  annotate(geom = "text", x = -7.5, y = -7.5, label = "UMAP_1", color = "black") +
  annotate(geom = "text", x = -9.6, y = -5.9, label = "UMAP_2", color = "black",angle = 90)  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(2),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),legend.title = element_text(angle = 90,size=20),legend.text=element_text(size=15),legend.key.size = unit(0.8, 'cm'))
p2 <-  ggplot(data=Non_Cd36_pseudotime_df,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=velocity_pseudotime),size=1)+
  scale_colour_gradientn(colours=newpalette,na.value="lightgrey",breaks=seq(0,1,0.2),limits=c(0,1),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.vjust=0.5,title.position = "left"))+
  labs(title="Trajectory of Cd36- OSNs",color="Pseudotime")+
  geom_segment(aes(x = -9, y = -7, xend = -6.5, yend = -7),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  geom_segment(aes(x = -9, y = -7, xend = -9, yend = -4.5),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  annotate(geom = "text", x = -7.5, y = -7.5, label = "UMAP_1", color = "black") +
  annotate(geom = "text", x = -9.6, y = -5.9, label = "UMAP_2", color = "black",angle = 90) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(2),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),legend.title = element_text(angle = 90,size=20),legend.text=element_text(size=15),legend.key.size = unit(0.8, 'cm'))
p1 / p2 + plot_layout(guides = 'collect') 
dev.off()

## Figure 4F - gene expression along pseudotime
genes <- c("Lhx2","Mef2a","Tshz1","Cd36")
newpalette <- as.vector(ArchRPalettes$stallion[c(19,9,8,7,2,13)])
library(ggnewscale)
pseudotime_df <- read.table("~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/scvelo/remove_HBCs_nonsort_HBC_lineage_scRNA_velocity_pseudotime.csv",sep="\t",header=TRUE)
pseudotime_df$X <- str_c(pseudotime_df$orig.ident,"_",pseudotime_df$X,"-1")
remove_HBCs_nonsort_scRNA <- subset(filtered_nonsort_scRNA,subset=detail_cell_type!="HBCs"& UMAP_1>(-10))

plot_01.ls <- lapply(1:3,function(i){
  exp_df <- FetchData(remove_HBCs_nonsort_scRNA,vars=genes[i],slot = "data")
  colnames(exp_df) <- "gene"
  Cd36_pseudotime_df <- pseudotime_df
  Cd36_pseudotime_df <- Cd36_pseudotime_df[-which(Cd36_pseudotime_df$X %in% non_Tshz1_cells | Cd36_pseudotime_df$detail_cell_type == "Cd36- OSNs"),]
  Cd36_pseudotime_df <- merge(Cd36_pseudotime_df,exp_df,by.y="row.names",by.x="X")
  Cd36_pseudotime_df$lineage <- "Cd36+ ORNs linegae"
  Non_Cd36_pseudotime_df <- pseudotime_df
  Non_Cd36_pseudotime_df <- Non_Cd36_pseudotime_df[-which(Non_Cd36_pseudotime_df$X %in% Tshz1_cells | Non_Cd36_pseudotime_df$detail_cell_type == "Cd36+ OSNs"),]
  Non_Cd36_pseudotime_df <- merge(Non_Cd36_pseudotime_df,exp_df,by.y="row.names",by.x="X")
  Non_Cd36_pseudotime_df$lineage <- "Non-Cd36+ ORNs lineage"
  df <- rbind(Cd36_pseudotime_df,Non_Cd36_pseudotime_df)
  df$detail_cell_type <- factor(df$detail_cell_type,levels=c("GBCs","INPs","Immature OSNs","Nearly mature OSNs","Cd36+ OSNs","Cd36- OSNs"),labels=c("GBCs","INPs","Immature OSNs","Adcy3- OSNs","Cd36+ OSNs","Cd36- OSNs"))
  df$lineage <- factor(df$lineage,levels=c("Cd36+ ORNs linegae","Non-Cd36+ ORNs lineage"),labels=c("Cd36+ OSNs linegae","Cd36- OSNs lineage"))
  df <- df %>%
    dplyr::rename(`Cell type`=detail_cell_type,Lineage=lineage)
  p <- ggplot(data=df,aes(x=velocity_pseudotime,y=gene)) + 
    geom_point(data=df,aes(x=velocity_pseudotime,y=gene,color=`Cell type`))+
    #scale_colour_gradientn(colours=newpalette,breaks=seq(0,1,0.2),limits=c(0,1),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.hjust=0.5,title.position = "top"))+
    scale_color_manual(values=newpalette)+
    new_scale_color() +
    geom_smooth(data=df,aes(x=velocity_pseudotime,y=gene,color=Lineage),method='loess',se = FALSE,linewidth=1.5) + 
    scale_color_manual(values=c("#D43F3AFF","#357EBDFF"))+
    labs(x="Pseudotime",y="Gene expression",title=genes[i]) +
    scale_x_continuous(breaks=seq(0,1,0.2))+
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.6),plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),axis.line = element_line(colour="black",size=1),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_text(size=rel(1.3),vjust=0.5),legend.text=element_text(size=rel(1)),legend.title=element_text(size=rel(1.2)),legend.position="none",axis.text.y=element_text(color="black",size=rel(1)))
  p
})

pdf(str_c(out_dir,"filtered_nonsort_HBC_lineage_scRNA_gene_expression_along_pseudotime_01.pdf"),height=11.5,width=9)
plot_01.ls[[1]] / plot_01.ls[[2]] / plot_01.ls[[3]] / plot_02.ls[[1]] + plot_layout(heights=c(2.5,2.5,2.5,2.6))
dev.off()


## Figure 4G upper panel - sequential expression of Lhx2, Mef2a, Tshz1 and Cd36 along the trajectory of Cd36+ OSNs
mm10_genes <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_symbols.txt")
OR_genes <- as.data.frame(readxl::read_excel("/md01/shipy3/Projects/mouse_ORs/input/PMC7055050_mouse_OR_genes.xlsx"))[,c("Gene symbol","Ensembl gene ID")]
OR_genes <- OR_genes[-which(duplicated(OR_genes)),]
merged_OR_genes <- merge(mm10_genes,OR_genes,by.x="V2",by.y="Ensembl gene ID")

remove_HBCs_nonsort_scRNA <- subset(filtered_nonsort_scRNA,subset=detail_cell_type!="HBCs"& UMAP_1>(-10))
Olfr_counts <- FetchData(object=remove_HBCs_nonsort_scRNA,vars=merged_OR_genes$V1,slot="counts")
Olfr_counts$OR_n <- rowSums(Olfr_counts[,1:1140]>=4)
top_OR_idx <- apply(Olfr_counts[,1:1140],1,which.max)
Olfr_counts$top_OR <- colnames(Olfr_counts[,1:1140])[as.vector(top_OR_idx)]
Olfr_counts$top_OR[which(Olfr_counts$OR_n==0)] <- NA
Cd36_ORNs_biased_ORs <- read.table("/md01/shipy3/Projects/mouse_ORs/output/scRNA/integrate_Tsukahara_Brann/Cd36_ORNs_biased_ORs.txt")$V1
Olfr_counts$top_OR_type <- "else"
Olfr_counts$top_OR_type[which(Olfr_counts$top_OR %in% Cd36_ORNs_biased_ORs)] <- "Cd36+ OSNs-biased ORs"
pseudotime_df <- read.table("~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/scvelo/remove_HBCs_nonsort_HBC_lineage_scRNA_velocity_pseudotime.csv",sep="\t",header=TRUE)
pseudotime_df$X <- str_c(pseudotime_df$orig.ident,"_",pseudotime_df$X,"-1")
Cd36_pseudotime_df <- pseudotime_df
Cd36_pseudotime_df <- Cd36_pseudotime_df[-which(Cd36_pseudotime_df$X %in% non_Tshz1_cells | Cd36_pseudotime_df$detail_cell_type == "Cd36- OSNs"),]
Olfr_counts <- merge(Olfr_counts,Cd36_pseudotime_df[,c("X","detail_cell_type","velocity_pseudotime")],by.x="row.names",by.y="X")
Olfr_counts <- Olfr_counts %>%
  arrange(velocity_pseudotime)
Olfr_counts$Cd36_no_not <- factor(Olfr_counts$top_OR_type,levels=c("Cd36+ OSNs-biased ORs","else"),labels=c("Yes","No"))


exp_df <- FetchData(object=remove_HBCs_nonsort_scRNA,vars=c("Lhx2","Mef2a","Tshz1","Cd36","detail_cell_type"),slot="data")
exp_df$detail_cell_type <- factor(exp_df$detail_cell_type,levels=c("GBCs","INPs","Immature ORNs","Nearly mature ORNs","Cd36+ ORNs","Non-Cd36+ ORNs"),labels=c("GBCs","INPs","Immature OSNs","Adcy3- OSNs","Cd36+ OSNs","Cd36- OSNs"))
exp_df <- exp_df[Olfr_counts$Row.names,]
rescale_exp_df <- t((apply(exp_df[,c("Lhx2","Mef2a","Tshz1","Cd36")],2,scale)))
colnames(rescale_exp_df) <- rownames(exp_df)
rescale_exp_df[which(rescale_exp_df>=1.5)] <- 1.5
rescale_exp_df[which(rescale_exp_df<=-1.5)] <- -1.5

pdf(str_c(out_dir,"Tshz1_Cd36_sequential_expression_along_Cd36_pos_pseudotime.pdf"),height=6,width=8.5)
ha1 <- HeatmapAnnotation(
  cell_type=exp_df$detail_cell_type,
  col=list(cell_type=c("GBCs"="#D8A767","INPs"="#E6C2DC","Immature OSNs"="#C06CAB","Adcy3- OSNs"="#8A9FD1","Cd36+ OSNs"="#272E6A")),
  annotation_legend_param = list(cell_type=list(title="Cell type",border="black")),
  show_legend=TRUE,
  show_annotation_name=TRUE,
  border = TRUE,
  annotation_label="Cell type",
  annotation_name_side="left"
  )
ha2 <- HeatmapAnnotation(
  expressed_OR_type=Olfr_counts$Cd36_no_not,
  col=list(expressed_OR_type=c("No"="white","Yes"="firebrick3")),
  annotation_legend_param = list(expressed_OR_type=list(title="Cd36+ ORNs-biased OR\nexpressed or not",border="black")),
  show_legend=TRUE,
  show_annotation_name=TRUE,
  border = TRUE,
  annotation_label="Cd36+ ORNs-biased OR\nexpressed or not",
  annotation_name_side="left"
  ) 
ht <- Heatmap(rescale_exp_df,
  name="Scaled expression",
  col=as.vector(ArchRPalettes[["solarExtra"]]),
  row_names_side="left",
  height=unit(17, "mm"),
  cluster_columns=FALSE,
  cluster_rows=FALSE,
  show_column_names=FALSE,
  show_row_names=TRUE,
  top_annotation=ha1,
  bottom_annotation=ha2,
  heatmap_legend_param = list(direction = "horizontal",legend_width = unit(4, "cm"), title_position = "topcenter",border="black",legend_height=unit(2, "cm")),
  border_gp = gpar(col = "black",lwd=1.5))
draw(ht, heatmap_legend_side = "bottom")
dev.off()


## Figure 4G bottom panel - sequential expression of Lhx2, Mef2a, Tshz1 and Hdac9 along the trajectory of Cd36- OSNs
mm10_genes <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_symbols.txt")
OR_genes <- as.data.frame(readxl::read_excel("/md01/shipy3/Projects/mouse_ORs/input/PMC7055050_mouse_OR_genes.xlsx"))[,c("Gene symbol","Ensembl gene ID")]
OR_genes <- OR_genes[-which(duplicated(OR_genes)),]
merged_OR_genes <- merge(mm10_genes,OR_genes,by.x="V2",by.y="Ensembl gene ID")

remove_HBCs_nonsort_scRNA <- subset(filtered_nonsort_scRNA,subset=detail_cell_type!="HBCs"& UMAP_1>(-10))
Olfr_counts <- FetchData(object=remove_HBCs_nonsort_scRNA,vars=merged_OR_genes$V1,slot="counts")
Olfr_counts$OR_n <- rowSums(Olfr_counts[,1:1140]>=4)
top_OR_idx <- apply(Olfr_counts[,1:1140],1,which.max)
Olfr_counts$top_OR <- colnames(Olfr_counts[,1:1140])[as.vector(top_OR_idx)]
Olfr_counts$top_OR[which(Olfr_counts$OR_n==0)] <- NA
Cd36_ORNs_biased_ORs <- read.table("/md01/shipy3/Projects/mouse_ORs/output/scRNA/integrate_Tsukahara_Brann/Cd36_ORNs_biased_ORs.txt")$V1
Olfr_counts$top_OR_type <- ""
Olfr_counts$top_OR_type[which((!is.na(Olfr_counts$top_OR)) & !(Olfr_counts$top_OR %in% Cd36_ORNs_biased_ORs))] <- "Non Cd36+ OSNs-biased ORs"

pseudotime_df <- read.table("~/Cd36_OSNs/output/scRNA/seurat/HBC_lineage_trajectory/scvelo/remove_HBCs_nonsort_HBC_lineage_scRNA_velocity_pseudotime.csv",sep="\t",header=TRUE)
pseudotime_df$X <- str_c(pseudotime_df$orig.ident,"_",pseudotime_df$X,"-1")

Non_Cd36_pseudotime_df <- pseudotime_df
Non_Cd36_pseudotime_df <- Non_Cd36_pseudotime_df[-which(Non_Cd36_pseudotime_df$X %in% Tshz1_cells | Non_Cd36_pseudotime_df$detail_cell_type %in% c("HBCs","Cd36+ OSNs")),]
Olfr_counts <- merge(Olfr_counts,Non_Cd36_pseudotime_df[,c("X","detail_cell_type","velocity_pseudotime")],by.x="row.names",by.y="X")
Olfr_counts <- Olfr_counts %>%
  arrange(velocity_pseudotime)
Olfr_counts$Cd36_deplted_no_not <- factor(Olfr_counts$top_OR_type,levels=c("Non Cd36+ OSNs-biased ORs",""),labels=c("Yes","No"))

exp_df <- FetchData(object=remove_HBCs_nonsort_scRNA,vars=c("Lhx2","Mef2a","Tshz1","Hdac9","detail_cell_type"),slot="data")
exp_df$detail_cell_type <- factor(exp_df$detail_cell_type,levels=c("GBCs","INPs","Immature ORNs","Nearly mature ORNs","Cd36+ ORNs","Non-Cd36+ ORNs"),labels=c("GBCs","INPs","Immature OSNs","Adcy3- OSNs","Cd36+ OSNs","Cd36- OSNs"))
exp_df <- exp_df[Olfr_counts$Row.names,]
rescale_exp_df <- t((apply(exp_df[,c("Lhx2","Mef2a","Tshz1","Hdac9")],2,scale)))
colnames(rescale_exp_df) <- rownames(exp_df)
rescale_exp_df[which(rescale_exp_df>=1.5)] <- 1.5
rescale_exp_df[which(rescale_exp_df<=-1.5)] <- -1.5

pdf(str_c("Lhx2_Mef2a_Tshz1_Hdac9_sequential_expression_along_Cd36_neg_pseudotime.pdf"),height=6,width=8.5)
ha1 <- HeatmapAnnotation(
  cell_type=exp_df$detail_cell_type,
  col=list(cell_type=c("GBCs"="#D8A767","INPs"="#E6C2DC","Immature OSNs"="#C06CAB","Adcy3- OSNs"="#8A9FD1","Cd36- OSNs"="#9983BD")),
  annotation_legend_param = list(cell_type=list(title="Cell type",border="black")),
  show_legend=TRUE,
  show_annotation_name=TRUE,
  border = TRUE,
  annotation_label="Cell type",
  annotation_name_side="left"
  )
ha2 <- HeatmapAnnotation(
  expressed_OR_type=Olfr_counts$Cd36_deplted_no_not,
  col=list(expressed_OR_type=c("No"="white","Yes"="firebrick3")),
  annotation_legend_param = list(expressed_OR_type=list(title="Non-Cd36+ ORNs-biased OR\nexpressed or not",border="black")),
  show_legend=TRUE,
  show_annotation_name=TRUE,
  border = TRUE,
  annotation_label="Non-Cd36+ ORNs-biased OR\nexpressed or not",
  annotation_name_side="left"
  ) 
ht <- Heatmap(rescale_exp_df,
  name="Scaled expression",
  col=as.vector(ArchRPalettes[["solarExtra"]]),
  row_names_side="left",
  height=unit(17, "mm"),
  cluster_columns=FALSE,
  cluster_rows=FALSE,
  show_column_names=FALSE,
  show_row_names=TRUE,
  top_annotation=ha1,
  bottom_annotation=ha2,
  heatmap_legend_param = list(direction = "horizontal",legend_width = unit(4, "cm"), title_position = "topcenter",border="black",legend_height=unit(2, "cm")),
  border_gp = gpar(col = "black",lwd=1.5))
draw(ht, heatmap_legend_side = "bottom")
dev.off()


```