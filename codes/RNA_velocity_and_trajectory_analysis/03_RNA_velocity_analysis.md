### Estimating RNA velocity
- Generate spliced and unspliced matrices by Velocyto
	- script : 01_run_RNA_Velocity.sh
- Convert data from Seurat to Python / anndata
	- reference : https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
```R
## Load required packages
library(Seurat)
library(SeuratDisk)

# Load 19d nosort data 
filtered_nonsort_scRNA <- readRDS("/data/R03/shipy3/Projects/mouse_ORs/output/scRNA/HBC_lineage_trajectory/filtered_nonsort_HBC_lineage_scRNA.rds")
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

SaveH5Seurat(remove_HBCs_nonsort_scRNA, filename = "/data/R03/shipy3/Projects/mouse_ORs/output/scRNA/HBC_lineage_trajectory/CellPath/remove_HBCs_nonsort_HBC_lineage_scRNA.h5Seurat",overwrite=TRUE)
Convert("/data/R03/shipy3/Projects/mouse_ORs/output/scRNA/HBC_lineage_trajectory/CellPath/remove_HBCs_nonsort_HBC_lineage_scRNA.h5Seurat", dest = "h5ad",overwrite=TRUE)




```