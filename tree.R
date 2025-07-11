library(phyloseq)
library(vegan)
library(data.table)
library(readxl)
library(Rcpp)
library(phyloseq)
library(dada2)
library(stringr)
library(readr)
library(readxl)
library(Biostrings)
library(ggtree)
library(ggplot2)
library(DECIPHER)
library(phangorn)
library(phyloseq)
library(ggtree)
setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1")

# 定义文件路径
asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/otu_table_rarefied.csv"
tax_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/Trivals/Desktop/usyd_tax_table-010822.csv"
sample_metadata_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/Trivals/Desktop/Matadata.csv"

# 读取原始ASV表和分类信息的CSV文件,read和fread等函数会自动将第一行识别成列名
asv_data <- fread(asv_file_path, header = TRUE, data.table = FALSE)
tax_data <- fread(tax_file_path, header = TRUE, data.table = FALSE)
sample_data <- fread(sample_metadata_path, header = TRUE, data.table = FALSE)

# 将数据的第一列设置为行名
rownames(asv_data) <- asv_data$SampleID  # 假设第一列没有名字，通常会自动命名为V1
rownames(tax_data) <- tax_data$SampleID
rownames(sample_data) <- sample_data$SampleID  

# 移除第一列
asv_data <- asv_data[,-1]
tax_data <- tax_data[,-1]
sample_data <- sample_data[,-1]

#确保time和storage都是因子
sample_data$Time <- as.factor(sample_data$Time)
sample_data$Storage <- as.factor(sample_data$Storage)
sample_data$Soil <- as.factor(sample_data$Soil)
# 创建phyloseq对象
physeq_obj <- phyloseq::phyloseq(
  otu_table(as.matrix(asv_data), taxa_are_rows = FALSE), 
  tax_table(as.matrix(tax_data)),  
  sample_data(sample_data)
)
taxa_names(physeq_obj) <- names(Seqs)


Seqs <- getSequences(colnames(asv_data))
names(Seqs) <- paste0("ASV",1:length(Seqs))

save(Seqs,file="Seqs.Rdata1")

alignment <- AlignSeqs(DNAStringSet(Seqs), anchor=NA)
alignment

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)

a <- as.matrix(dm)

options(getClass.msg=FALSE)

tree.NJ <- NJ(dm)

fit = pml(tree.NJ, data=phang.align)

fitGTR_0 <- update(fit, k=4, inv=0.2)
fitGTR <-optim.pml(fitGTR_0, model="GTR", optInv=TRUE, optGamma=TRUE, 
                   rearrangement = "stochastic", control = pml.control(trace = 0))

fitGTR$tree

save(fitGTR, file="T1_Bac_fitGTR.RData")

tree <- fitGTR$tree
plot(tree)

Bac_tree <- merge_phyloseq(physeq_obj, phy_tree(tree))
Bac_tree

plot_tree(Bac_tree)

save(Bac_tree,file="Bac_tree2.19.RData")
save.image(file="Image_Bac_Tree.Rdata")


###替换ASV

###替换ASV，之前跑完了

asv_table <- otu_table(Bac_tree)
asv_table <- as.data.frame(as.matrix(asv_table))

tax_table <- tax_table(Bac_tree)
tax_table <- as.data.frame(as.matrix(tax_table))


write.csv(asv_table, file = "asv_table.csv")

write.csv(tax_table, file = "tax_table.csv")





library(picante)
setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/R_practice/TREE_rare")
asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/rara_abundance/rare_number.csv"
# asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/rara_abundance/abundance_number.csv"

# 定义文件路径
# asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/otu_table_rarefied.csv"
tax_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/Trivals/Desktop/usyd_tax_table-010822.csv"
sample_metadata_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/Trivals/Desktop/Matadata.csv"

# 读取原始ASV表和分类信息的CSV文件,read和fread等函数会自动将第一行识别成列名
asv_data <- fread(asv_file_path, header = TRUE, data.table = FALSE)
tax_data <- fread(tax_file_path, header = TRUE, data.table = FALSE)
sample_data <- fread(sample_metadata_path, header = TRUE, data.table = FALSE)

# 将数据的第一列设置为行名
rownames(asv_data) <- asv_data$SampleID  # 假设第一列没有名字，通常会自动命名为V1
rownames(tax_data) <- tax_data$SampleID
rownames(sample_data) <- sample_data$SampleID  

# 移除第一列
asv_data <- asv_data[,-1]
tax_data <- tax_data[,-1]
sample_data <- sample_data[,-1]

dist.mat <- phyloseq::UniFrac(phyloseq_tree, weighted = F)

asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/rara_abundance/T_rare.csv"


# 读取原始ASV表和分类信息的CSV文件,read和fread等函数会自动将第一行识别成列名
asv_data <- fread(asv_file_path, header = TRUE, data.table = FALSE)


# 将数据的第一列设置为行名
rownames(asv_data) <- asv_data$SampleID  # 假设第一列没有名字，通常会自动命名为V1
# 移除第一列
asv_data <- asv_data[,-1]

phylo_tree <- phy_tree(phyloseq_tree)
tree <- prune.sample(asv_data, phylo_tree )
dis <- cophenetic(tree)

set.seed(123)
mntd <- ses.mntd(asv_data , dis, abundance.weighted = TRUE, null.model = 'taxa.labels', runs = 999)
mntd


###########
# 定义文件路径
asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/asv_table.csv"
tax_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/tax_table.csv"
sample_metadata_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/sam_table.csv"

# 读取原始ASV表和分类信息的CSV文件,read和fread等函数会自动将第一行识别成列名
feature_table <- fread(asv_file_path, header = TRUE, data.table = FALSE)
tax_table <- fread(tax_file_path, header = TRUE, data.table = FALSE)
sample_table<- fread(sample_metadata_path, header = TRUE, data.table = FALSE)


# 将数据的第一列设置为行名
rownames(sample_table) <- sample_table$Label # 假设第一列没有名字，通常会自动命名为V1
rownames(tax_table) <- tax_table$Label
rownames(feature_table) <- feature_table$Label
# 移除第一列
sample_table <- sample_table[,-1]
tax_table<- tax_table[,-1]
feature_table<- feature_table[,-1]

# 创建phyloseq对象
physeq_obj <- phyloseq::phyloseq(
  otu_table(as.matrix(feature_table), taxa_are_rows = FALSE), 
  tax_table(as.matrix(tax_table)),  
  sample_data(sample_table)
)


setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1")
tree <- read.tree(file = "tree.nwk")

Bac_tree <- merge_phyloseq(physeq_obj, phy_tree(tree))

dist.mat <- phyloseq::UniFrac(Bac_tree, weighted = T)

write.csv(dist.mat, "uni_abundance.csv")

