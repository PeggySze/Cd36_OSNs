## Get fasta of functional OR genes
```R
mm10_genes <- read.table("/public/home/shipy3/DB/10X_reference/cellranger-arc-mm10/gene_symbols.txt")
OR_genes <- as.data.frame(readxl::read_excel("~/Cd36_OSNs/input/PMC7055050_mouse_OR_genes.xlsx"))[,c("Gene symbol","Gene name","Ensembl gene ID")]
OR_genes <- OR_genes[-which(duplicated(OR_genes)),]
merged_OR_genes <- merge(mm10_genes,OR_genes,by.x="V2",by.y="Ensembl gene ID") 
functional_ORs <- merged_OR_genes[-grep("pseudogene",merged_OR_genes$`Gene name`),]
write.table(functional_ORs[,"V2",drop=FALSE],"~/Cd36_OSNs/output/ORs_phylogenetic_tree/functional_ORs_Ensembl_ID.txt",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
```
- get protein sequences for all functional OR genes
```python
import os

# get functional ORs IDs
functional_ORs_file = open("~/Cd36_OSNs/output/ORs_phylogenetic_tree/functional_ORs_Ensembl_ID.txt")
functional_ORs = []
for line in functional_ORs_file:
	functional_ORs.append(line.strip())
functional_ORs_file.close()

# get protein sequences for each functional OR gene
fasta_file = open("/md01/shipy3/ori/DB/mm10/fasta/Mus_musculus.GRCm38.pep.all.fa")
seq = ""
for line in fasta_file:
	if line[0] == ">" and seq == "":
		header = line
	elif line[0] != ">":
		seq = seq + line
	elif line[0] == ">" and seq != "":
		gene_id = header.split()[3].split(":")[1].split(".")[0]
		if gene_id in functional_ORs:
			output_file = open(os.path.join("~/Cd36_OSNs/output/ORs_phylogenetic_tree", gene_id+".fasta"),"a")
			output_file.write(header+seq)
			output_file.close()
		header = line
		seq = ""
gene_id = header.split()[3].split(":")[1].split(".")[0]
if gene_id in functional_ORs:
	output_file = open(os.path.join("~/Cd36_OSNs/output/ORs_phylogenetic_tree", gene_id+".fasta"),"a")
	output_file.write(header+seq)
	output_file.close()
fasta_file.close()


## get protein sequences for all functional OR genes
fasta_file = open("/md01/shipy3/ori/DB/mm10/fasta/Mus_musculus.GRCm38.pep.all.fa")
seq = ""
output_file = open("~/Cd36_OSNs/output/ORs_phylogenetic_tree/functional_ORs.fasta","a")
for line in fasta_file:
	if line[0] == ">" and seq == "":
		header = line
	elif line[0] != ">":
		seq = seq + line
	elif line[0] == ">" and seq != "":
		gene_id = header.split()[3].split(":")[1].split(".")[0]
		if gene_id in functional_ORs:
			output_file.write(header+seq)
		header = line
		seq = ""
gene_id = header.split()[3].split(":")[1].split(".")[0]
if gene_id in functional_ORs:
	output_file.write(header+seq)
output_file.close()
fasta_file.close()
```

-reference : https://www.protocols.io/view/week-5-aligning-with-muscle-and-making-trees-with-e6nvwbo7vmkj
## Alignment
```shell
## download MUSCLE v5 (MUltiple Sequence Comparison by Log- Expectation)
cd /md01/shipy3/software
chmod +x muscle5.1.linux_intel64

## align OR protein sequences by MUSCLE v5 
# script : 
# /md01/shipy3/Projects/mouse_ORs/output/ORs_phylogenetic_tree/MUSCLE_alignment.sh
# /md01/shipy3/Projects/mouse_ORs/output/ORs_phylogenetic_tree/MUSCLE_super5_alignment.sh
```

## Phylogenetic tree construction with RAxML
- phylip file : http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html
```python
## modify header of aligned fasta file (keep protein ID 后6位)
afa_file = open("~/Cd36_OSNs/output/ORs_phylogenetic_tree/aln_functional_ORs.afa")
seq = ""
output_file = open("~/Cd36_OSNs/output/ORs_phylogenetic_tree/modified_aln_functional_ORs.afa","w")
for line in afa_file:
	if line[0] == ">" and seq == "":
		pep_id = line.split()[0].split(">")[1].split(".")[0][12:18]
	elif line[0] != ">":
		seq = seq + line
	elif line[0] == ">" and seq != "":
		output_file.write(">"+pep_id+"\n")
		output_file.write(seq)
		pep_id = line.split()[0].split(">")[1].split(".")[0][12:18]
		seq = ""
output_file.write(">"+pep_id+"\n")
output_file.write(seq)
output_file.close()
afa_file.close()

## modify header of aligned fasta file (keep protein ID)
afa_file = open("~/Cd36_OSNs/output/ORs_phylogenetic_tree/aln_functional_ORs.afa")
seq = ""
output_file = open("~/Cd36_OSNs/output/ORs_phylogenetic_tree/proteinID_aln_functional_ORs.afa","w")
for line in afa_file:
	if line[0] == ">" and seq == "":
		pep_id = line.split()[0].split(">")[1].split(".")[0]
	elif line[0] != ">":
		seq = seq + line
	elif line[0] == ">" and seq != "":
		output_file.write(">"+pep_id+"\n")
		output_file.write(seq)
		pep_id = line.split()[0].split(">")[1].split(".")[0]
		seq = ""
output_file.write(">"+pep_id+"\n")
output_file.write(seq)
output_file.close()
afa_file.close()
```
```shell
## install fast2phy
cd ~/Cd36_OSNs/output/ORs_phylogenetic_tree/
git clone https://github.com/davidmnoriega/fast2phy.git
cd fast2phy
python setup.py install

## convert aligned fasta (.afa) file into Phylip file (.phy)
fast2phy ~/Cd36_OSNs/output/ORs_phylogenetic_tree/modified_aln_functional_ORs.afa -o ~/Cd36_OSNs/output/ORs_phylogenetic_tree/RAxML/functional_ORs.phy 

## construct phylogenetic tree 
cd ~/Cd36_OSNs/output/ORs_phylogenetic_tree/RAxML/
raxmlHPC -f a -# 1000 -m PROTGAMMAAUTO -p 12345 -x 12345 -s ~/Cd36_OSNs/output/ORs_phylogenetic_tree/proteinID_aln_functional_ORs.afa -n functional_ORs.tree -T 20
```


## Visualize phylogenetic tree by iTOL
```R
## load required packages
library(RColorBrewer)

## iTOL color file
iTOL_file <- read.table("~/Cd36_OSNs/output/ORs_phylogenetic_tree/RAxML/iTOL_file.txt",comment.char="&")


Cd36_ORNs_biased_ORs <- read.table("~/Cd36_OSNs/output/scRNA/integrate_Tsukahara_Brann/Cd36_ORNs_biased_ORs.txt")
other_ORNs_biased_ORs <- read.table("~/Cd36_OSNs/output/scRNA/integrate_Tsukahara_Brann/other_ORNs_biased_ORs.txt")

# add color for OR types
iTOL_file$OR <- sapply(1:nrow(iTOL_file),function(i){
	unlist(strsplit(iTOL_file$V1[i],"_"))[1]
	})
iTOL_file$OR_type_color <- ""
iTOL_file$OR_type_color[which(iTOL_file$OR %in% Cd36_ORNs_biased_ORs$V1)] <- "#cd2626"
#iTOL_file$OR_type_color[which(iTOL_file$OR %in% other_ORNs_biased_ORs$V1)] <- "#1874cd"
iTOL_file$OR_type_color[which(iTOL_file$OR_type_color=="")] <- "#C0C0C0"

# add color for chromosome
iTOL_file$chr <- sapply(1:nrow(iTOL_file),function(i){
	unlist(strsplit(iTOL_file$V1[i],"_"))[2]
	})
chrs <- unique(iTOL_file$chr)
colors <- colorRampPalette(brewer.pal(8,"Set2"))(length(chrs))
names(colors) <- chrs
iTOL_file$chr_color <- sapply(1:nrow(iTOL_file),function(i){
	colors[iTOL_file$chr[i]]
	})
write.table(iTOL_file[,c(1,5)],"~/Cd36_OSNs/output/ORs_phylogenetic_tree/RAxML/iTOL_file_01.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(iTOL_file[,c(1,7)],"~/Cd36_OSNs/output/ORs_phylogenetic_tree/RAxML/iTOL_file_02.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
```