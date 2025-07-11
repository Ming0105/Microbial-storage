library(phyloseq)
library(vegan)
library(data.table)
library(readxl)

# 定义文件路径
asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/original/asv_table.csv"
tax_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/original/tax_table.csv"
sample_metadata_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/original/sam_table.csv"

# 读取原始ASV表和分类信息的CSV文件,read和fread等函数会自动将第一行识别成列名
asv_data <- fread(asv_file_path, header = TRUE, data.table = FALSE)
tax_data <- fread(tax_file_path, header = TRUE, data.table = FALSE)
sample_data <- fread(sample_metadata_path, header = TRUE, data.table = FALSE)

# 将数据的第一列设置为行名
rownames(asv_data) <- asv_data$SampleID  # 假设第一列没有名字，通常会自动命名为V1
rownames(tax_data) <- tax_data$SampleID
rownames(sample_data) <- sample_data$SampleID  # 假设第一列名为SampleID

# 移除第一列
asv_data <- asv_data[,-1]
tax_data <- tax_data[,-1]
sample_data <- sample_data[,-1]
# 
# #确保time和storage都是因子
# sample_data$Time <- as.factor(sample_data$Time)
# sample_data$Storage <- as.factor(sample_data$Storage)

# 创建phyloseq对象
physeq_obj <- phyloseq::phyloseq(
  otu_table(as.matrix(asv_data), taxa_are_rows = FALSE), 
  tax_table(as.matrix(tax_data)),  
  sample_data(sample_data)
)


physeq_obj <- physeq_obj %>% 
  subset_taxa(Kingdom=="Bacteria")%>%
  subset_taxa(!(Family %in% "Mitochondria") &
                !(Order %in%  "Chloroplast") )


#  移除总丰度低于10的物种
physeq_obj <- prune_taxa(taxa_sums(physeq_obj) > 15, physeq_obj) 
physeq_obj <- filter_taxa(physeq_obj ,function(x)sum(x>0)>2,TRUE)

summary(taxa_sums(physeq_obj))
physeq_obj 
print(physeq_obj)
sample_sums <- rowSums(otu_table(physeq_obj))
min_read <- min(sample_sums)
min_read


# 进行抽平操作
set.seed(123) # 为了确保结果的可重复性
physeq_rarefied <- rarefy_even_depth(physeq_obj, sample.size = 33000, rngseed = 123, replace = TRUE)

#提取抽平后的OTU表
otu_table_rarefied <- otu_table(physeq_rarefied)

asv_table <- as.data.frame(as.matrix(otu_table_rarefied))


tax_table <- tax_table(physeq_rarefied  )
tax_table <- as.data.frame(as.matrix(tax_table))



setwd("C:/Users/Admin/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy")

write.csv(asv_table, file = "asv_table.csv",row.names = T)
write.csv(tax_table, file = "tax_table.csv",row.names = T)
