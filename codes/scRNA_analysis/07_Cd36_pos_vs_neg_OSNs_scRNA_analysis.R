## Load required packages
library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tibble) 
library(sctransform)
library(patchwork)

## Specify output directory 
out_dir <- "~/Cd36_OSNs/output/scRNA/seurat/integrated_analysis/integrate_Tsukahara_Brann/"
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

## Load data
integrated_mature_ORNs.scRNA <- readRDS(str_c(out_dir,"remove_homecage_4_integrated_mature_ORNs_Seurat_obj.rds"))

# exclude a small cluster expressing Ifi27 
filtered_Ifi27_cells_mature_ORNs <- subset(integrated_mature_ORNs.scRNA,subset=cell_type!="Ifi27+ cells")


## Figure 1A - mature ORNs scRNA-seq UMAP showing Cd36 expression 
palettes <- c(brewer.pal(8,"Pastel2")[8],brewer.pal(9,"Greens")[3:8])
DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"
Cd36_gene_exp_df <- FetchData(filtered_Ifi27_cells_mature_ORNs,vars = c("Cd36","UMAP_1","UMAP_2"),slot = "data")
pdf(str_c(out_dir,"integrated_Tsukahara_Brann_mature_ORNs_scRNA_UMAP_highlight_Cd36+_ORNs.pdf"),height=9)
ggplot(data=Cd36_gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=Cd36))+
  geom_point(size=0.1)+
  scale_colour_gradientn(colours = palettes,limits=c(0,quantile(Cd36_gene_exp_df$Cd36,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(Cd36_gene_exp_df$Cd36,seq(0,1,0.01))[96]),labels = c("0", "1.5"),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.position = "top",title.hjust=0.5)) + 
  labs(title="Mature ORNs scRNA-seq",color="Cd36 expression") +
  geom_segment(aes(x = -9, y = -5, xend = -7, yend = -5),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  geom_segment(aes(x = -9, y = -5, xend = -9, yend = -3.5),arrow = arrow(length = unit(0.3, "cm")),color="black")+
  annotate(geom = "text", x = -7.7, y = -5.4, label = "UMAP 1", color = "black",size=4.5) +
  annotate(geom = "text", x = -9.6, y = -4.2, label = "UMAP 2", color = "black",angle = 90,size=4.5) +
  annotate(geom = "text", x = -9, y = 2, label = "Cd36+ ORNs", color = "black",size=5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(2.5)),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.size = unit(0.8, 'cm'),legend.text=element_text(size=rel(1.2)),legend.position="bottom",legend.title=element_text(size=rel(1.4)))
dev.off()


## Figure 1B - Heatmap showing DEGs between Cd36+ ORNs and non-Cd36+ ORNs
filtered_Ifi27_cells_mature_ORNs$Cd36_or_not <- ifelse(filtered_Ifi27_cells_mature_ORNs$cell_type=="Cd36+ ORNs","Cd36+ ORNs","Cd36- ORNs")
Idents(filtered_Ifi27_cells_mature_ORNs) <- filtered_Ifi27_cells_mature_ORNs$Cd36_or_not
DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"

Cd36_ORNs_markers <- FindMarkers(filtered_Ifi27_cells_mature_ORNs,ident.1="Cd36+ ORNs",assay="RNA")
Cd36_ORNs_markers <- Cd36_ORNs_markers %>%
  dplyr::filter(p_val_adj<0.05,avg_log2FC >= log2(1.5)) %>%
  dplyr::arrange(desc(avg_log2FC))
other_ORNs_markers <- FindMarkers(filtered_Ifi27_cells_mature_ORNs,ident.1="Cd36- ORNs",assay="RNA")
other_ORNs_markers <- other_ORNs_markers %>%
  dplyr::filter(p_val_adj<0.05,avg_log2FC >= log2(1.5)) %>%
  dplyr::arrange(desc(avg_log2FC))

set.seed(100)
cell_types <- c("Cd36+ ORNs","Cd36- ORNs")
scRNA_sampled_df <- c()
for (i in cell_types){
  cells <- rownames(filtered_Ifi27_cells_mature_ORNs@meta.data)[which(filtered_Ifi27_cells_mature_ORNs$Cd36_or_not==i)]
  if (length(cells)>=500){
    sampled_cells <- sample(cells,500,replace=FALSE)
  } else{
    sampled_cells <- sample(cells,500,replace=TRUE)
  }
  scRNA_sampled_df <- data.frame(idx=1:500,cell_type=i,sampled_cell=sampled_cells) %>% rbind(scRNA_sampled_df ,.)
}
DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"
scRNA_df <- FetchData(object = filtered_Ifi27_cells_mature_ORNs, vars = c("Cd36_or_not",rownames(Cd36_ORNs_markers),rownames(other_ORNs_markers)),slot = "data",cells=scRNA_sampled_df$sampled_cell)
scRNA_df$Cd36_or_not <- factor(scRNA_df$Cd36_or_not,levels=c("Cd36+ ORNs","Cd36- ORNs"),labels=c("Cd36+ ORNs","Cd36- ORNs"))
scaled_scRNA_df <- t(apply(scRNA_df[,-1],2,scale))
colnames(scaled_scRNA_df) <- rownames(scRNA_df)
rescaled_scRNA_df <- scaled_scRNA_df
rescaled_scRNA_df[rescaled_scRNA_df>=4] <- 4
rescaled_scRNA_df[rescaled_scRNA_df<=-4] <- -4

library(ComplexHeatmap)
cell_type_cols <- c("#00A087FF",brewer.pal(8,"Set2")[8])
names(cell_type_cols) <- c("Cd36+ ORNs","Cd36- ORNs")
highlight_genes <- c("Cd36","Cyp4v3","Abcd2","Mboat1","Snca","Rab3gap1","Macrod1","Hdac9","Tshz1","Zfp407","Cd9","Robo2","Cd55","Klkb1","Sgms2")
pdf(str_c(out_dir,"Cd36_vs_other_ORNs_DEGs_heatmap.pdf"),width=3.8,height=6)
ha1 <- HeatmapAnnotation(
  df=data.frame(row.names=rownames(scRNA_df),cell_type=scRNA_df$Cd36_or_not),
  col=list(cell_type=cell_type_cols),
  annotation_legend_param = list(cell_type=list(title="Cell type")),
  gp = gpar(col = "black"),
  show_legend=FALSE,
  show_annotation_name=FALSE
  )
ha2 <- rowAnnotation(gene = anno_mark(at = which(rownames(rescaled_scRNA_df) %in% highlight_genes), labels = rownames(rescaled_scRNA_df)[which(rownames(rescaled_scRNA_df) %in% highlight_genes)],labels_gp=gpar(fontsize = 10),side = "left"))
h1 <- Heatmap(rescaled_scRNA_df,
  name="RNA z-score",
  col=as.vector(ArchRPalettes[["solarExtra"]]),
  cluster_columns=FALSE,
  column_split=factor(as.vector(scRNA_df$Cd36_or_not),levels=c("Cd36+ ORNs","Cd36- ORNs")),
  column_title_gp = gpar(fontsize = 11),
  show_column_names=FALSE,
  show_row_names=FALSE,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 11),
  column_gap = unit(0, "mm"),
  cluster_rows=FALSE,
  cluster_column_slices=FALSE,
  top_annotation=ha1,
  left_annotation=ha2,
  heatmap_legend_param = list(direction = "horizontal",legend_width = unit(4, "cm"), title_position = "topcenter",border="black")
)
draw(h1, heatmap_legend_side = "bottom")
dev.off()

## GO enrichment of Cd36+ ORNs markers
library(clusterProfiler)
library(org.Mm.eg.db)

GO_ORA <- function(gene_sets){
  ENTREZID_id <- bitr(gene_sets,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db",drop=TRUE)
  BP_GO_ORA <- enrichGO(gene=ENTREZID_id$ENTREZID,
    OrgDb=org.Mm.eg.db,
    ont= "BP",
    readable= TRUE)
  BP_GO_ORA <- simplify(BP_GO_ORA)
  MF_GO_ORA <- enrichGO(gene=ENTREZID_id$ENTREZID,
    OrgDb=org.Mm.eg.db,
    ont= "MF",
    readable= TRUE)
  MF_GO_ORA <- simplify(MF_GO_ORA)
  CC_GO_ORA <- enrichGO(gene=ENTREZID_id$ENTREZID,
    OrgDb=org.Mm.eg.db,
    ont= "CC",
    readable= TRUE)
  CC_GO_ORA <- simplify(CC_GO_ORA)
  BP_GO_ORA_result <- as.data.frame(BP_GO_ORA)
  MF_GO_ORA_result <- as.data.frame(MF_GO_ORA)
  CC_GO_ORA_result <- as.data.frame(CC_GO_ORA)
  BP_GO_ORA_result <- BP_GO_ORA_result[order(BP_GO_ORA_result$pvalue),]
  MF_GO_ORA_result <- MF_GO_ORA_result[order(MF_GO_ORA_result$pvalue),]
  CC_GO_ORA_result <- CC_GO_ORA_result[order(CC_GO_ORA_result$pvalue),]
  results <- list(BP_GO_ORA=BP_GO_ORA_result,
    MF_GO_ORA=MF_GO_ORA_result,
    CC_GO_ORA=CC_GO_ORA_result)
  return(results)
}

Cd36_ORNs_markers_GO <- GO_ORA(rownames(Cd36_ORNs_markers))
write.table(Cd36_ORNs_markers_GO[[1]],str_c(out_dir,"Cd36_up_DEGs_BP_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(Cd36_ORNs_markers_GO[[2]],str_c(out_dir,"Cd36_up_DEGs_MF_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

other_ORNs_markers_GO <- GO_ORA(rownames(other_ORNs_markers))
write.table(other_ORNs_markers_GO[[1]],str_c(out_dir,"Cd36_down_DEGs_BP_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(other_ORNs_markers_GO[[2]],str_c(out_dir,"Cd36_down_DEGs_MF_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
write.table(other_ORNs_markers_GO[[3]],str_c(out_dir,"Cd36_down_DEGs_CC_GO_ORA.txt"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)


## Figure 1C - gene expression for Cd36+ ORNs specific genes in UMAP
markers <- c("Cyp4v3","Mboat1","Klkb1","Sgms2")
gene_exp_palettes <- c(brewer.pal(9,"Greys")[3],brewer.pal(9,"Greens")[3:9])
DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"
gene_exp_UMAP_plot.ls <- lapply(1:length(markers),function(i){
  gene_exp_df <- FetchData(filtered_Ifi27_cells_mature_ORNs,vars = c(markers[i],"UMAP_1","UMAP_2"),slot = "data")
  colnames(gene_exp_df)[1] <- "gene"
  p <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.1)+
    scale_colour_gradientn(colours = gene_exp_palettes,limits=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[99]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[99]),labels = c("0",round(quantile(gene_exp_df$gene,seq(0,1,0.01))[99],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 0.5,title.hjust=0.5,title.position = "top")) + 
    scale_x_continuous(limits=c(-12.5,7.5),expand = c(0,0))+
    scale_y_continuous(limits=c(-7.5,8),expand = c(0,0))+
    labs(title=markers[i],color="Gene expression") +
    #annotate("rect", xmin = -12.5, xmax = -8, ymin = -7.5, ymax = -6.5,fill = "white",color="black",size=1)+
    #annotate(geom = "text", x = -10.3, y = -7, label = str_c("0-",round(quantile(gene_exp_df$gene,seq(0,1,0.01))[96],1)), color = "black",size=5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),plot.title = element_text(hjust = 0.5,size=rel(2.5),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title=element_blank(),plot.margin=unit(c(0, 0.2, 0.7, 0.5), "cm"),legend.title = element_text(size=15),legend.text=element_text(size=10),legend.position=c(0.2,0.15),legend.direction="horizontal",legend.key.size = unit(0.8, 'cm'))
    p
  })

pdf(str_c(out_dir,"Cd36_ORNs_specific_genes_expression_UMAP_08.pdf"),width=11.5,height=13)
(gene_exp_UMAP_plot.ls[[1]] + gene_exp_UMAP_plot.ls[[2]]) / (gene_exp_UMAP_plot.ls[[3]] + gene_exp_UMAP_plot.ls[[4]])
dev.off()


## Figure 1D - expression of dorsal/ventral markers in UMAP
markers <- c("Nqo1","Ncam2")
gene_exp_palettes1 <- c(brewer.pal(8,"Pastel2")[8],brewer.pal(9,"Greens")[3:8])
gene_exp_palettes2 <- c(brewer.pal(8,"Pastel2")[8],brewer.pal(9,"Greens")[3:8])
DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"
gene_exp_df_plot.ls <- list()
gene_exp_df <- FetchData(filtered_Ifi27_cells_mature_ORNs,vars = c(markers[1],"UMAP_1","UMAP_2"),slot = "data")
colnames(gene_exp_df)[1] <- "gene"
gene_exp_df_plot.ls[[1]] <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.1)+
    scale_colour_gradientn(colours = gene_exp_palettes1,limits=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),labels = c("0", round(quantile(gene_exp_df$gene,seq(0,1,0.01))[96],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1,title.position="top", title.hjust = 0.5)) + 
    labs(title=markers[1],color="Gene expression") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),plot.margin=unit(c(0.2, 0.2, 0.4, 0.2), "cm"),legend.title = element_text(size=15),legend.text=element_text(size=10),legend.position=c(0.2,0.15),legend.direction="horizontal",legend.key.size = unit(0.8, 'cm'))


gene_exp_df <- FetchData(filtered_Ifi27_cells_mature_ORNs,vars = c(markers[2],"UMAP_1","UMAP_2"),slot = "data")
colnames(gene_exp_df)[1] <- "gene"
gene_exp_df_plot.ls[[2]] <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.1)+
    scale_colour_gradientn(colours = gene_exp_palettes2,limits=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$gene,seq(0,1,0.01))[96]),labels = c("0", round(quantile(gene_exp_df$gene,seq(0,1,0.01))[96],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1,title.position="top", title.hjust = 0.5)) + 
    labs(title=markers[2],color="Gene expression") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_text(hjust = 0.5,size=rel(3),face="bold"),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),plot.margin=unit(c(0.2, 0.2, 0.4, 0.2), "cm"),legend.title = element_text(size=15),legend.text=element_text(size=10),legend.position=c(0.2,0.15),legend.direction="horizontal",legend.key.size = unit(0.8, 'cm'))

pdf(str_c(out_dir,"mature_ORNs_dorsal_ventral_markers_UMAP.pdf"),width=14,height=7)
p <- gene_exp_df_plot.ls[[1]] + gene_exp_df_plot.ls[[2]] + plot_layout(guides = "collect")
x.grob <- textGrob("UMAP_1", gp=gpar(fontsize=20))
y.grob <- textGrob("UMAP_2",gp=gpar(fontsize=20),rot=90)
grid.arrange(patchworkGrob(p), left = y.grob, bottom = x.grob)
dev.off()


## Figure 1E - the proportion of dorsal/ventral cells in Cd36+ ORNs & Non-Cd36+ ORNs
markers <- c("Nqo1","Acsm4","Ncam2","Nfix","Nfib")
filtered_Ifi27_cells_mature_ORNs$Cd36_or_not <- ifelse(filtered_Ifi27_cells_mature_ORNs$cell_type=="Cd36+ ORNs","Cd36+ ORNs","Other ORNs")
df <- FetchData(filtered_Ifi27_cells_mature_ORNs,vars = c(markers,"Cd36_or_not"),slot = "data")
df$mean_dorsal_exp <- (df$Nqo1+df$Acsm4)/2
df$mean_ventral_exp <- (df$Ncam2+df$Nfix+df$Nfib)/3
df <- df[-which(df$mean_dorsal_exp==0 & df$mean_ventral_exp==0),]
df$region <- ifelse(df$mean_dorsal_exp>df$mean_ventral_exp,"dorsal","ventral")
count_df <- df %>%
  dplyr::group_by(Cd36_or_not,region) %>%
  dplyr::summarise(counts=n()) %>%
  dplyr::mutate(percentage=counts/sum(counts)*100)
count_df$Cd36_or_not <- factor(count_df$Cd36_or_not,levels=rev(c("Cd36+ ORNs","Other ORNs")),labels=rev(c("Cd36+ OSNs","Cd36- OSNs")))
count_df$region <- factor(count_df$region,levels=c("ventral","dorsal"))
count_df$shown_percentage <- str_c(round(count_df$percentage,1),"%")
newpalette <- c("#D81C38","#1565C0")
percentage_plot <- ggplot(data=count_df,aes(x=percentage,y=Cd36_or_not,fill=region))+
  geom_bar(stat="identity",width=0.7)+
  labs(x="Percentage (%)",y="",fill="Spatial location")+
  geom_text(aes(label=shown_percentage),size = 6,color="white",position=position_stack(vjust=0.5))+
  scale_fill_manual(values=newpalette,breaks=c("dorsal","ventral"))+
  theme_classic()+
  theme(axis.title.x = element_text(color="black",size=rel(2)),axis.text.x=element_text(color="black",size=rel(1.8)),axis.title.y = element_blank(),axis.text.y = element_text(color="black",size=rel(2)),axis.line = element_line(colour="black",size=1.2),legend.text=element_text(size=rel(1.6)),legend.title=element_text(size=rel(1.8)),plot.margin=unit(c(0.5, 0.2, 0.5, 0.2), "cm"),legend.position="bottom")

pdf(str_c(out_dir,"mature_ORNs_dorsal_ventral_markers_UMAP&percentage_plot.pdf"),width=12.5,height=10.5)
(gene_exp_df_plot.ls[[1]] + gene_exp_df_plot.ls[[2]])  / percentage_plot + plot_layout(heights=c(6,3))
dev.off()


## Figure 1F - Cd36 expression in reconstructed map of olfactory bulb
map_pre <- readRDS("~/Cd36_OSNs/input/GSE169012_olfactory_bulb_Slide_seq2/map_pre.rds")

plot_features <- function(df, targets, options = "viridis", directions = -1) {
  require(ggplot2)
  require(viridis)
  require(jpeg)
  require(ggpubr)
  img = readJPEG("/md01/shipy3/Projects/mouse_ORs/input/GSE169012_olfactory_bulb_Slide_seq2/projection3.jpg")
  df = df[order(df[,paste0(targets)], decreasing = F),]
  if (grepl("^[[:digit:]]+", targets)) {
    colnames(df)[which(colnames(df)==paste0(targets))] = paste0("gene_",targets)
    targets = paste0("gene_",targets)
  }
  if (grepl("-", targets)) {
    new_name  = gsub("-","_", targets)
    colnames(df)[which(colnames(df)==paste0(targets))] = new_name
    targets = new_name
  }
  if (is.numeric(df[,paste0(targets)])) {
    ggplot(df, aes_string("x", "y", colour = paste0(targets)))+background_image(img)+geom_point(size=3)+
      xlim(c(0,2450))+ylim(c(-3130,300))+
      labs(title=targets,color="Scaled expression")+
      scale_color_viridis(direction = directions, option = options,guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1),breaks=c(min(df[,paste0(targets)]),max(df[,paste0(targets)])),labels=c("Min","Max"))+
      theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5,size=rel(2.5)),plot.margin=unit(c(0.5, 0, 0, 0), "cm")) + 
      guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))
  } else {
    ggplot(df, aes_string("x", "y", colour = paste0(targets)))+background_image(img)+geom_point()+
      xlim(c(0,2450))+ylim(c(-3130,300))+theme_classic()}
}

pdf(str_c(out_dir,"Cd36_gene_expression_in_olfactory_bulb_reconstructed_map.pdf"),height=8,width=6)
plot_features(map_pre,"Cd36",directions=1) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5,size=rel(2.3)),plot.margin=unit(c(0, 0, 0, 0), "cm"),legend.key.size = unit(0.8,'cm'),legend.key.width = unit(1, 'cm'),legend.title=element_text(size=rel(1.5)),legend.text=element_text(size=rel(1.3)),legend.position="bottom") +NoAxes()
dev.off()

## Figure 1G+H - gene expression of Plxna1 (UMAP + violin plot)
gene_exp_palettes <- c(brewer.pal(8,"Pastel2")[8],brewer.pal(9,"Greens")[3:8])
DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"
gene_exp_df <- FetchData(filtered_Ifi27_cells_mature_ORNs,vars = c("Plxna1","UMAP_1","UMAP_2"),slot = "data")
gene_exp_UMAP_plot <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=Plxna1))+
    geom_point(size=0.1)+
    scale_colour_gradientn(colours = gene_exp_palettes,limits=c(0,quantile(gene_exp_df$Plxna1,seq(0,1,0.01))[96]),oob = scales::squish,breaks=c(0,quantile(gene_exp_df$Plxna1,seq(0,1,0.01))[96]),labels = c("0", round(quantile(gene_exp_df$Plxna1,seq(0,1,0.01))[96],1)),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1,title.position="top", title.hjust = 0.5)) + 
    scale_x_continuous(limits=c(-12.5,7.5),expand = c(0,0))+
    scale_y_continuous(limits=c(-7.5,8),expand = c(0,0))+
    labs(color="Gene expression",x="UMAP 1",y="UMAP 2") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.y=element_text(size=13,color="black"),axis.title.x=element_text(size=13,color="black"),plot.margin=unit(c(0.3, 0.2, 0.7, 0.5), "cm"),legend.title = element_text(size=15),legend.text=element_text(size=10),legend.position=c(0.2,0.15),legend.direction="horizontal",legend.key.size = unit(0.8, 'cm'))


violin_palettes <- c(brewer.pal(9,"Greens")[7],brewer.pal(8,"Set2")[8])
fmt_dcimals <- function(decimals=0){
    function(x) format(x,nsmall = decimals,scientific = FALSE)
  }
library(ggpubr)
Cd36_vs_other_exp_df <- FindMarkers(filtered_Ifi27_cells_mature_ORNs,features="Plxna1",ident.1="Cd36+ ORNs",logfc.threshold=-Inf,min.pct=-Inf)
gene_exp_df <- FetchData(filtered_Ifi27_cells_mature_ORNs,vars = c("Plxna1","UMAP_1","UMAP_2","Cd36_or_not"),slot = "data")
gene_exp_df$Cd36_or_not <- factor(gene_exp_df$Cd36_or_not,levels=c("Cd36+ ORNs","Other ORNs"),labels=c("Cd36+ OSNs","Cd36- OSNs"))
Cd36_vs_other_p_val_adj <- signif(Cd36_vs_other_exp_df$p_val_adj[which(rownames(Cd36_vs_other_exp_df)=="Plxna1")],3)
stat_df <- tibble::as_tibble(data.frame(group1="Cd36+ OSNs",group2="Cd36- OSNs",p.adj=Cd36_vs_other_p_val_adj,y.position=round(max(gene_exp_df$Plxna1),1)+0.7))
gene_exp_violin_plot <- ggplot() + 
    geom_violin(data=gene_exp_df,aes(x=Cd36_or_not,y=Plxna1,fill=Cd36_or_not),scale="width",trim=FALSE)+
    geom_boxplot(data=gene_exp_df,aes(x=Cd36_or_not,y=Plxna1,fill=Cd36_or_not),width=0.05, fill="white",outlier.size=0.8)+
    scale_fill_manual(values=violin_palettes)+
    labs(y="Gene expression") + 
    scale_y_continuous(labels = fmt_dcimals(0),limits=c(round(min(gene_exp_df$Plxna1),1),round(max(gene_exp_df$Plxna1),1)+1))+
    stat_pvalue_manual(stat_df, label = "p = {p.adj}",label.size=5)+ 
    theme_classic()+
    theme(
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=rel(1.5),vjust=0.5),
          axis.text.y = element_text(color="black",size=rel(1.5)), 
          axis.text.x = element_text(color="black",size=rel(1.5)),
          axis.line = element_line(colour="black",size=1.4),
          axis.ticks = element_line(),
          plot.margin=unit(c(0.3, 0.2, 0.7, 0.2), "cm"),
          legend.position = " none")

pdf(str_c(out_dir,"Plxna1_gene_expression_UMAP_VlnPlot.pdf"),width=10,height=6)
gene_exp_UMAP_plot + gene_exp_violin_plot + plot_layout(widths=c(4.5,4)) + plot_annotation(title="Plxna1",theme=theme(plot.title = element_text(hjust = 0.5,size=rel(2.5),face="bold")))
dev.off()






