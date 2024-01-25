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

mm10_genes <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_symbols.txt")
OR_genes <- as.data.frame(readxl::read_excel("/md01/shipy3/Projects/mouse_ORs/input/PMC7055050_mouse_OR_genes.xlsx"))[,c("Gene symbol","Ensembl gene ID")]
OR_genes <- OR_genes[-which(duplicated(OR_genes)),]
merged_OR_genes <- merge(mm10_genes,OR_genes,by.x="V2",by.y="Ensembl gene ID")

# our own data
filtered_mature_ORNs.scRNA <- readRDS("~/Cd36_OSNs/output/scRNA/seurat/integrated_analysis/integrated_filtered_mature_ORNs_scRNA.rds")
filtered_Ighm_cells_mature_ORNs.scRNA <- subset(filtered_mature_ORNs.scRNA,subset=cell_type!="Ighm+ cells")
DefaultAssay(filtered_Ighm_cells_mature_ORNs.scRNA) <- "RNA"
df <- FetchData(object =filtered_Ighm_cells_mature_ORNs.scRNA,vars=merged_OR_genes$V1,slot = "data")
expressed_OR_idx <- apply(df,1,which.max)
df$expressed_OR <- colnames(df)[expressed_OR_idx]
XuLab_meta_df <- data.frame(row.names=rownames(df),expressed_OR=df$expressed_OR,stringsAsFactors=FALSE)
filtered_Ighm_cells_mature_ORNs.scRNA <- AddMetaData(object=filtered_Ighm_cells_mature_ORNs.scRNA,metadata =XuLab_meta_df)

# Tsukahara_Brann data
Tsukahara_Brann_home_cage_metadata <- read.csv("/md01/shipy3/learning/Tsukahara_Brann_OSN/data/raw/GSE173947_home_cage_metadata.csv.gz",row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
orig_idents <- c("baseline-1","baseline-2","baseline-3","baseline-4","baseline-9","baseline-10")
samples <- paste("homecage",1:6,sep="_")
Tsukahara_Brann_home_cage_metadata$sample <- ""
for (i in 1:length(orig_idents)){
  Tsukahara_Brann_home_cage_metadata$sample[which(Tsukahara_Brann_home_cage_metadata$orig_ident==orig_idents[i])] <- samples[i]
}
Tsukahara_Brann_home_cage_metadata$cell <- sapply(1:nrow(Tsukahara_Brann_home_cage_metadata),function(i){
  str_c(Tsukahara_Brann_home_cage_metadata$sample[i],"_",unlist(strsplit(rownames(Tsukahara_Brann_home_cage_metadata)[i],"_"))[2],"-1")
  })

Tsukahara_Brann_meta_df <- data.frame(cell=Tsukahara_Brann_home_cage_metadata$cell,Tsukahara_Brann_top_Olfr=Tsukahara_Brann_home_cage_metadata$top_Olfr) %>%
  filter(cell %in% colnames(filtered_Ifi27_cells_mature_ORNs),Tsukahara_Brann_top_Olfr %in% merged_OR_genes$V1) 
DefaultAssay(filtered_Ifi27_cells_mature_ORNs) <- "RNA"
Tsukahara_Brann_Olfr_counts <- FetchData(object =filtered_Ifi27_cells_mature_ORNs,vars=merged_OR_genes$V1,slot = "data",cells=Tsukahara_Brann_meta_df$cell)
expressed_OR_idx <- apply(Tsukahara_Brann_Olfr_counts,1,which.max)
Tsukahara_Brann_Olfr_counts$expressed_OR <- colnames(Tsukahara_Brann_Olfr_counts)[expressed_OR_idx]
Tsukahara_Brann_meta_df_02 <- Tsukahara_Brann_Olfr_counts[,"expressed_OR",drop=FALSE]
Tsukahara_Brann_meta_df <- merge(Tsukahara_Brann_meta_df,Tsukahara_Brann_meta_df_02,by.x="cell",by.y="row.names")

# remove ORNs that have confilit expressed OR genes
XuLab_meta_df <- XuLab_meta_df[which(rownames(XuLab_meta_df) %in% colnames(filtered_Ifi27_cells_mature_ORNs)),,drop=FALSE]
XuLab_meta_df$Tsukahara_Brann_top_Olfr <- NA
Tsukahara_Brann_meta_df <- Tsukahara_Brann_meta_df[-which(Tsukahara_Brann_meta_df$Tsukahara_Brann_top_Olfr != Tsukahara_Brann_meta_df$expressed_OR),]
rownames(Tsukahara_Brann_meta_df) <- Tsukahara_Brann_meta_df$cell
Tsukahara_Brann_meta_df <- Tsukahara_Brann_meta_df[,c("expressed_OR","Tsukahara_Brann_top_Olfr")]
meta_df <- rbind(XuLab_meta_df,Tsukahara_Brann_meta_df)
meta_df <- meta_df[which(rownames(meta_df) %in% colnames(filtered_Ifi27_cells_mature_ORNs)),]

filtered_confilit_mature_ORNs <- subset(filtered_Ifi27_cells_mature_ORNs,cells=rownames(meta_df))
filtered_confilit_mature_ORNs <- AddMetaData(filtered_confilit_mature_ORNs,metadata=meta_df) 
# 48074 cells 

## Figure S3A - histogram of cell number for each ORN group
counts_df <- meta_df %>%
  dplyr::group_by(expressed_OR) %>%
  dplyr::summarise(counts=n())
pdf(str_c(out_dir,"filtered_confilit_mature_ORNs_each_ORN_group_cell_number_histogram.pdf"))
ggplot(data=counts_df,aes(x=counts))+
  geom_histogram(fill=brewer.pal(9,"Blues")[5],binwidth=10,color="black",alpha = 0.7)+
  labs(x="# of cells per ORN group",y="# of ORN groups",title="945 ORN groups with >= 10 cells")+
  geom_vline(xintercept=10,linetype="dashed",size=0.5,color="black")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5,size=rel(1.6)),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1))
dev.off()

# retain 945 ORN subgroups that at least have 10 cells
all_ORNs_expressed_OR_df <- filtered_confilit_mature_ORNs@meta.data %>%
  dplyr::group_by(expressed_OR) %>%
  dplyr::summarise(all_counts=n()) %>%
  dplyr::arrange(desc(all_counts))
retained_ORs <- all_ORNs_expressed_OR_df$expressed_OR[which(all_ORNs_expressed_OR_df$all_counts>=10)]
more_than_10_cells_mature_ORNs <- subset(filtered_confilit_mature_ORNs,subset=expressed_OR %in% retained_ORs)

all_ORNs_expressed_OR_df <- more_than_10_cells_mature_ORNs@meta.data %>%
  dplyr::group_by(expressed_OR) %>%
  dplyr::summarise(all_counts=n()) %>%
  arrange(desc(all_counts))

more_than_10_cells_mature_ORNs$Cd36_or_not <- ifelse(more_than_10_cells_mature_ORNs$cell_type=="Cd36+ ORNs","Cd36+ ORNs","Non-Cd36+ ORNs")
Cd36_ORNs_expressed_OR_df <- more_than_10_cells_mature_ORNs@meta.data %>%
  dplyr::filter(Cd36_or_not=="Cd36+ ORNs") %>%
  dplyr::group_by(expressed_OR) %>%
  dplyr::summarise(Cd36_ORNs_counts=n()) %>%
  arrange(desc(Cd36_ORNs_counts))

other_ORNs_expressed_OR_df <- more_than_10_cells_mature_ORNs@meta.data %>%
  dplyr::filter(Cd36_or_not=="Non-Cd36+ ORNs") %>%
  dplyr::group_by(expressed_OR) %>%
  dplyr::summarise(other_ORNs_counts=n()) %>%
  arrange(desc(other_ORNs_counts))

df <- merge(all_ORNs_expressed_OR_df,Cd36_ORNs_expressed_OR_df,all.x=TRUE)
df <- merge(df,other_ORNs_expressed_OR_df,all.x=TRUE)
df[is.na(df)] <- 0

table(more_than_10_cells_mature_ORNs$Cd36_or_not)
# Cd36+ ORNs Non-Cd36+ ORNs
#      2573          44829

exp_p <- 2573/(2573+44829)
df$Cd36_biased_pval <- sapply(1:nrow(df),function(i){
  binom.test(x=df$Cd36_ORNs_counts[i],n=df$all_counts[i],p=exp_p,alternative="greater")$p.value
  })
df$Cd36_biased_fdr <- p.adjust(df$Cd36_biased_pval,method="fdr")
df$other_biased_pval <- sapply(1:nrow(df),function(i){
  binom.test(x=df$Cd36_ORNs_counts[i],n=df$all_counts[i],p=exp_p,alternative="less")$p.value
  })
df$other_biased_fdr <- p.adjust(df$other_biased_pval,method="fdr")
df$gene_type <- "Nonsignificant"
df$gene_type[which(df$Cd36_biased_fdr<0.05)] <- "Cd36+_ORNs_biased"
df$gene_type[which(df$other_biased_fdr<0.05)] <- "Other_ORNs_biased"
df$gene_type <- factor(df$gene_type,levels=rev(c("Cd36+_ORNs_biased","Other_ORNs_biased","Nonsignificant")))
df$Cd36_ORNs_percentage <- df$Cd36_ORNs_counts/df$all_counts*100
df$other_ORNs_percentage <- df$other_ORNs_counts/df$all_counts*100
df$adjusted_Cd36_ORNs_ratio <- (df$Cd36_ORNs_counts/2573) / ((df$Cd36_ORNs_counts/2573) + (df$other_ORNs_counts/44829))*100
df <- df %>%
  arrange(gene_type)

df$pval <- sapply(1:nrow(df),function(i){
  binom.test(x=df$Cd36_ORNs_counts[i],n=df$all_counts[i],p=exp_p,alternative="two.sided")$p.value
  })
df$fdr <- p.adjust(df$pval,method="fdr")
df$mlog10pval <- -log10(df$pval)
df$mlog10fdr <- -log10(df$fdr)
df$gene_type <- factor(df$gene_type,levels=c("Cd36+_ORNs_biased","Other_ORNs_biased","Nonsignificant"),labels=c("Cd36+ ORNs-biased","Cd36- ORNs-biased","Unbiased"))


Cd36_ORNs_biased_df <- df[which(df$gene_type=="Cd36+_ORNs_biased"),] %>% arrange(Cd36_biased_fdr)
write.table(Cd36_ORNs_biased_df[,1,drop=FALSE],str_c(out_dir,"Cd36_ORNs_biased_ORs.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE)

other_ORNs_biased_df <- df[which(df$gene_type=="Other_ORNs_biased"),] %>% arrange(other_biased_fdr)
write.table(other_ORNs_biased_df[,1,drop=FALSE],str_c(out_dir,"other_ORNs_biased_ORs.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE)

## Figure 2A - Perentage of Cd36+ OSNs contribution in each ORN subgroup
pdf(str_c(out_dir,"945_ORN_subgroups_adjusted_Cd36_ORNs_contribution_with_pvalue.pdf"),width=7.5,height=3.5)
p1 <- ggplot(data=df,aes(x=adjusted_Cd36_ORNs_ratio,y=1)) +
  geom_jitter(aes(size=pmin(mlog10fdr,3),fill=gene_type),pch=21,height=1) + 
  scale_size_continuous(range=c(0.5,2),breaks = c(0,1.3,2,3),labels=c("0","1.3","2",">=3"),limits=c(0,3))+
  labs(x="Adjusted Cd36+ ORNs contribution (%)",fill="OR gene type",size="-log10(FDR)")+ 
  scale_fill_manual(values = c("firebrick3","dodgerblue3","lightgrey")) +
  geom_vline(xintercept=50,linetype="dotted",size=0.5,color="black")+
  scale_x_continuous(breaks=c(0,25,50,75,100)) + 
  theme_bw() + 
  theme(panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=1.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.title.x = element_text(size=rel(1),vjust=-0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x= element_text(color="black",size=rel(1.2)),
    axis.ticks.y = element_blank(),
    plot.margin=unit(c(0.005,0,0,0),"cm"))
df$gene_type <- factor(df$gene_type,levels=c("Unbiased","Cd36+ ORNs-biased","Cd36- ORNs-biased"))
p2 <- ggplot(data=df,aes(x = adjusted_Cd36_ORNs_ratio, fill = gene_type)) +
  geom_histogram(alpha = 0.7,position = "stack",binwidth=1,color="black") +
  scale_fill_manual(values = c("lightgrey","firebrick3","dodgerblue3")) +
  coord_trans(y = "log1p") +
  scale_y_continuous(breaks=c(0,1,50,100,800))+
  labs(y="# ORs") +
  theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.y=element_text(size=10),axis.text.y = element_text(color="black",size=5), axis.text.x = element_blank(),axis.title.x=element_blank(),axis.ticks.x = element_blank(),plot.margin=unit(c(0,0,0,0),"cm"),legend.position="none")
null <- ggplot()+theme_void()
cowplot::plot_grid(p2,null,p1,null,axis="tblr",align="v"
  ,rel_widths=c(1.6,0.4),rel_heights=c(0.2,0.8))
dev.off()

## Figure2B - visualize Cd36+ ORNs biased / other ORNs-biased / nonbiased ORNs in UMAP
Cd36_ORNs_biased_df <- df %>%
  filter(gene_type=="Cd36+ ORNs-biased")%>% 
  arrange(Cd36_biased_fdr)
other_ORNs_biased_df <- df %>%
  filter(gene_type=="Cd36- ORNs-biased")%>% 
  arrange(other_biased_fdr)

ORs <- c(Cd36_ORNs_biased_df$expressed_OR[c(1,6)],other_ORNs_biased_df$expressed_OR[1:2])
OR_types <- c(rep("Cd36+ ORNs-biased",2),rep("Cd36- ORNs-biased",2))
cols <- c("lightgrey","red")
DefaultAssay(filtered_confilit_mature_ORNs) <- "RNA"
gene_exp_df_plot_01.ls <- lapply(c(1,3),function(i){
  gene_exp_df <- FetchData(filtered_confilit_mature_ORNs,vars = c(ORs[i],"UMAP_1","UMAP_2"),slot = "data")
  colnames(gene_exp_df)[1] <- "gene"
  gene_exp_df$gene[which(gene_exp_df$gene<=quantile(gene_exp_df$gene,seq(0,1,0.01))[100])] <- 0
  gene_exp_df <- gene_exp_df %>% arrange(gene)
  p <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.8)+
    scale_colour_gradientn(colours = cols,breaks=c(0,max(gene_exp_df$gene)),labels = c("0", "Max"),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1)) + 
    labs(title=OR_types[i],y=ORs[i]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_text(hjust = 0.5,size=rel(1.5)),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),legend.position = 'none',axis.title.y=element_text(size=rel(2)))
  p
  })
gene_exp_df_plot_02.ls <- lapply(c(2,4),function(i){
  gene_exp_df <- FetchData(filtered_confilit_mature_ORNs,vars = c(ORs[i],"UMAP_1","UMAP_2"),slot = "data")
  colnames(gene_exp_df)[1] <- "gene"
  gene_exp_df$gene[which(gene_exp_df$gene<=quantile(gene_exp_df$gene,seq(0,1,0.01))[100])] <- 0
  gene_exp_df <- gene_exp_df %>% arrange(gene)
  p <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.8)+
    scale_colour_gradientn(colours = cols,breaks=c(0,max(gene_exp_df$gene)),labels = c("0", "Max"),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1)) + 
    labs(y=ORs[i]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_text(hjust = 0.5,size=rel(2)),axis.text=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),legend.position = 'none',axis.title.y=element_text(size=rel(2)))
  p
  })
legend_plot <- ggplot(data=gene_exp_df,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.8)+
    scale_colour_gradientn(colours = cols,breaks=c(0,max(gene_exp_df$gene)),labels = c("0", "Max"),guide = guide_colorbar(frame.colour = "black",frame.linewidth= 1)) + 
    labs(title=ORs[1],color="Gene expression") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_rect(fill='transparent', color='black',linetype="solid",size=2.5),plot.title = element_text(hjust = 0.5,size=rel(2)),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),legend.key.size = unit(1, 'cm'),legend.text=element_text(size=rel(1.5)),legend.position="bottom",legend.title=element_text(size=rel(1.8)))
library("grid")
library("gridExtra")
library("cowplot")
legend <- get_legend(legend_plot) 

pdf(str_c(out_dir,"biased&unbiased_OR_genes_UMAP.pdf"),width=6.5,height=9)
(gene_exp_df_plot_01.ls[[1]] + gene_exp_df_plot_01.ls[[2]] ) / (gene_exp_df_plot_02.ls[[1]] + gene_exp_df_plot_02.ls[[2]]) + legend
dev.off()

## Figure 2C - visualize the location of Cd36+ ORNs biased / other ORNs-biased / nonbiased OR genes
gene_pos <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_pos.txt")
colnames(gene_pos) <- c("chrom","start","end","symbol","gene_id")
mm10_genes <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_symbols.txt")
OR_genes <- as.data.frame(readxl::read_excel("/md01/shipy3/Projects/mouse_ORs/input/PMC7055050_mouse_OR_genes.xlsx"))[,c("Gene symbol","Ensembl gene ID")]
OR_genes <- OR_genes[-which(duplicated(OR_genes)),]
merged_OR_genes <- merge(mm10_genes,OR_genes,by.x="V2",by.y="Ensembl gene ID")
merged_OR_genes <- merge(merged_OR_genes,gene_pos,by.x="V2",by.y="gene_id",all.x=TRUE)
OR_cluster <- as.data.frame(readxl::read_excel("/md01/shipy3/Projects/mouse_ORs/input/elife_OR_cluster.xlsx",sheet=1))
OR_cluster_chroms <- sapply(1:nrow(OR_cluster),function(i){
  unlist(strsplit(OR_cluster[i,2],":"))[1]
  })
OR_cluster_starts <- sapply(1:nrow(OR_cluster),function(i){
  as.numeric(unlist(strsplit(unlist(strsplit(OR_cluster[i,2],":"))[2],"-"))[1])
  })
OR_cluster_starts[5] <- 78515000
OR_cluster_starts[8] <- 48978000
OR_cluster_ends <- sapply(1:nrow(OR_cluster),function(i){
  as.numeric(unlist(strsplit(unlist(strsplit(OR_cluster[i,2],":"))[2],"-"))[2])
  })
OR_cluster.gr <- GRanges(OR_cluster_chroms,IRanges(start = OR_cluster_starts, end = OR_cluster_ends), strand = "*")
seqlevels(OR_cluster.gr) <- paste("chr",c(1:11,13:17,19,"X"),sep="")
OR_cluster.gr <- sort(OR_cluster.gr, ignore.strand=TRUE)
OR_cluster_df <- as.data.frame(OR_cluster.gr)[,1:3] %>%
  mutate(OR_cluster_region=str_c(seqnames,"-",start,"-",end)) %>% 
  select(OR_cluster_region,start,end,seqnames) %>%
  rename(chr=seqnames)
OR_cluster_pos <- OR_cluster_df$OR_cluster_region

expressed_ORNs_df <- df
expressed_ORNs_df <- merge(expressed_ORNs_df,merged_OR_genes,by.x="expressed_OR",by.y="V1",all.x=TRUE,sort=FALSE)
expressed_ORNs_df$chrom <- factor(expressed_ORNs_df$chrom,levels=paste("chr",c(1:11,13:17,19,"X"),sep=""))
expressed_ORNs_df <- expressed_ORNs_df %>%
  arrange(chrom,start)
expressed_ORNs.gr <- GRanges(expressed_ORNs_df$chrom,IRanges(start = expressed_ORNs_df$start, end = expressed_ORNs_df$end), strand = "*",gene=expressed_ORNs_df$expressed_OR)
ov <- findOverlaps(OR_cluster.gr,expressed_ORNs.gr)
expressed_ORNs.df <- c()
for (i in 1:length(ov)){
  expressed_ORNs.df <- data.frame(chr=OR_cluster_pos[queryHits(ov)[i]],start=start(expressed_ORNs.gr)[subjectHits(ov)[i]],end=end(expressed_ORNs.gr)[subjectHits(ov)[i]],Cd36_ORNs_ratio=expressed_ORNs_df$adjusted_Cd36_ORNs_ratio[subjectHits(ov)[i]],gene_type=expressed_ORNs_df$gene_type[subjectHits(ov)[i]]) %>% rbind(expressed_ORNs.df,.)
}
expressed_ORNs.df <- expressed_ORNs.df %>%
  arrange(gene_type)

correspondance <- OR_cluster_df[c("chr","start","end","OR_cluster_region","start","end")]

library(circlize)
pdf(str_c(out_dir,"945_ORN_subgroups_circos_plot.pdf"))
f1 = function() {
    circos.par(gap.after = 2, start.degree = 90)
    circos.initializeWithIdeogram(species="mm10",plotType = c("ideogram", "labels"), ideogram.height = 0.03)
}

#chr_bg_color = rand_color(21, transparency = 0.8)
chr_bg_color <- c(rep("white",9),adjust_transparency("#eb3528", alpha = 0.2),rep("white",11))
names(chr_bg_color) = paste0("chr",c(1:19,"X","Y"))
choose_color <- function(gene_types){
  colors <- c()
  for (i in 1:length(gene_types)){
    if (gene_types[i]=="Cd36+_ORNs_biased"){
      colors[i] <- "firebrick3"
    } else if (gene_types[i]=="Other_ORNs_biased"){
      colors[i] <- "dodgerblue3"
    } else {
      colors[i] <- "black"
    }
  }
  return(colors)
}

f2 = function() {
    circos.par(cell.padding = c(0, 0, 0, 0), gap.after = c(rep(1, nrow(OR_cluster_df)-1), 10),"track.height" = 0.4)
    circos.genomicInitialize(OR_cluster_df, plotType = NULL)
    circos.genomicTrack(expressed_ORNs.df, ylim = c(-10, 100), 
        panel.fun = function(region, value, ...) {
            for(h in seq(0, 100, by = 20)) {
                circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#AAAAAA")
            }
            circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 3, col = "#888888")
        
            circos.genomicPoints(region, value, 
                col = choose_color(value[[2]]), 
                pch = 20, cex = 0.7)
    }, bg.col = chr_bg_color[OR_cluster_df$chr], track.margin = c(0.02, 0))
    circos.yaxis(side = "left", at = seq(0, 100, by = 20), 
        sector.index = get.all.sector.index()[1], labels.cex = 0.8)
}
circos.nested(f1, f2, correspondance, connection_col = chr_bg_color[correspondance[[1]]])
circos.clear()
dev.off()


## ## test if biased OR genes enriched in a certain chromosome
expressed_ORs_chr_counts_df <- expressed_ORNs_df %>%
  group_by(chrom) %>%
  summarise(all_counts=n())
Cd36_biased_ORs_chr_counts_df <- expressed_ORNs_df %>%
  filter(gene_type=="Cd36+_ORNs_biased") %>%
  group_by(chrom) %>%
  summarise(Cd36_biased_counts=n())
other_biased_ORs_chr_counts_df <- expressed_ORNs_df %>%
  filter(gene_type=="Other_ORNs_biased") %>%
  group_by(chrom) %>%
  summarise(other_biased_counts=n())
unbiased_ORs_chr_counts_df <- expressed_ORNs_df %>%
  filter(gene_type=="Nonsignificant") %>%
  group_by(chrom) %>%
  summarise(unbiased_counts=n())
ORs_chr_counts_df <- merge(expressed_ORs_chr_counts_df,Cd36_biased_ORs_chr_counts_df,all.x=TRUE)
ORs_chr_counts_df <- merge(ORs_chr_counts_df,other_biased_ORs_chr_counts_df,all.x=TRUE)
ORs_chr_counts_df <- merge(ORs_chr_counts_df,unbiased_ORs_chr_counts_df,all.x=TRUE)
ORs_chr_counts_df[is.na(ORs_chr_counts_df)] <- 0
ORs_chr_counts_df$Cd36_biased_p_val <- sapply(1:nrow(ORs_chr_counts_df),function(i){
  exp_p <- ORs_chr_counts_df$all_counts[i]/sum(ORs_chr_counts_df$all_count)
  binom.test(x=ORs_chr_counts_df$Cd36_biased_counts[i],n=sum(ORs_chr_counts_df$Cd36_biased_counts),p=exp_p,alternative="greater")$p.value
    })
ORs_chr_counts_df$other_biased_p_val <- sapply(1:nrow(ORs_chr_counts_df),function(i){
  exp_p <- ORs_chr_counts_df$all_counts[i]/sum(ORs_chr_counts_df$all_count)
  binom.test(x=ORs_chr_counts_df$other_biased_counts[i],n=sum(ORs_chr_counts_df$other_biased_counts),p=exp_p,alternative="greater")$p.value
    })
ORs_chr_counts_df[which(ORs_chr_counts_df$Cd36_biased_p_val<0.05),] # Cd36+ ORNs biased OR genes enriched in chr10
ORs_chr_counts_df$chrom[which(ORs_chr_counts_df$other_biased_p_val<0.05)] # other ORNs biased OR genes enriched in chr14

## Figure 2D - visualize OR number for each chromosome
plot_df <- ORs_chr_counts_df[,c("chrom","all_counts","Cd36_biased_counts")]
plot_df$else_counts <- plot_df$all_counts - plot_df$Cd36_biased_counts
plot_df <- reshape2::melt(plot_df[,c("chrom","Cd36_biased_counts","else_counts")],id="chrom",variable.name="OR_type",value.name="counts")
plot_df$OR_type <- factor(plot_df$OR_type,levels=c("Cd36_biased_counts","else_counts"),labels=c("Cd36+ ORNs-biased","Else"))
ORs_chr_counts_df <- ORs_chr_counts_df %>%
  arrange(desc(all_counts))
plot_df$chrom <- factor(plot_df$chrom,levels=ORs_chr_counts_df$chrom)
pdf(str_c(out_dir,"each_chrom_Cd36_biased_ORs_number_barplot.pdf"),width=10)
ggplot(data=plot_df,aes(x=chrom,y=counts,fill=OR_type))+
  geom_bar(stat="identity",width=0.5,position="stack")+
  scale_fill_manual(values=c("firebrick3","lightgrey"))+
  labs(title="Number of OR genes for each chromosome",y="Number of OR genes",fill="OR type")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(2),face="bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black",size=rel(1.8)),
        axis.text.y = element_text(color="black",size=rel(1.5)), 
        axis.text.x = element_text(size=rel(1.8),color="black",angle = 45, hjust = 1),
        axis.line = element_line(colour="black",size=1.5),
        axis.ticks = element_line(),
        plot.margin=unit(c(0, 0.2, 0.7, 0.2), "cm"),
        legend.position = "bottom") 
dev.off()


## Figure 2E - visualize the position of Cd36-biased ORN groups in olfactory bulb reconstructed spatial map
Cd36_ORNs_biased_ORs <- read.table(str_c(out_dir,"Cd36_ORNs_biased_ORs.txt"),header=FALSE)
Cd36_ORNs_biased_ORs <- Cd36_ORNs_biased_ORs$V1
require(ggplot2)
require(viridis)
require(jpeg)
require(ggpubr)
library(gghighlight)
map_pre <- readRDS("~/Cd36_OSNs/input/GSE169012_olfactory_bulb_Slide_seq2/map_pre.rds") 
img = readJPEG("~/Cd36_OSNs/input/GSE169012_olfactory_bulb_Slide_seq2/projection3.jpg")
visualization_df <- map_pre 
visualization_df$annotation <- ifelse(visualization_df$OR %in% Cd36_ORNs_biased_ORs,"Cd36+ ORNs biased","Else")
visualization_df$annotation <- factor(visualization_df$annotation,levels=c("Cd36+ ORNs biased","Else"))
visualization_df <- visualization_df %>%
  arrange(desc(annotation))
pdf(str_c(out_dir,"highlight_Cd36_biased_ORN_groups_in_olfactory_bulb_2D_reconstructed_spatial_map.pdf"), width = 5.5, height = 7)
ggplot(visualization_df, aes(x = x, y = y,color=annotation))+
  background_image(img)+ 
  geom_point(size=2.5)+
  scale_color_manual(values=c("firebrick3","snow3"))+
  xlim(c(0,2450))+ylim(c(-3130,300))+
  labs(title="Reconstructed OB spatial map\n(Wang IH et al.)",color="OR type")+
  theme_classic()+NoAxes()+
  theme(legend.position="bottom",plot.title=element_text(hjust = 0.5,size=rel(1.5)),legend.text=element_text(size=rel(1.1)),legend.title=element_text(size=rel(1.3)))
dev.off()