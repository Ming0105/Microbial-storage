library(phyloseq)
library(vegan)
library(data.table)
library(readxl)

asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/rara_abundance/abudance_number.csv"

# 定义文件路径
asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/rara_abundance/rare_number.csv"
tax_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/Trivals/Desktop/usyd_tax_table-010822.csv"
sample_metadata_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/Trivals/Desktop/Matadata1.csv"

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
# 
# load("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/R_practice/TREE/Bac_tree.RData")

# 创建phyloseq对象
physeq_obj <- phyloseq::phyloseq(
  otu_table = otu_table(as.matrix(asv_data), taxa_are_rows = FALSE), 
  tax_table = tax_table(as.matrix(tax_data)),  
  sample_data = sample_data
)


########### Tree
library(DECIPHER)
library(phangorn)

# 1. alpha diversity--------------------------------------------------

# 计算alpha多样性指数
alpha_div <- estimate_richness(physeq_obj, measures = c("Observed", "ACE", "Shannon", "Simpson", "Chao1"))
rownames(alpha_div) <- gsub("^X", "", rownames(alpha_div))

rownames(alpha_div) <- gsub("\\.", "_", rownames(alpha_div))


# 定义目标文件路径
file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/rara_abundance/alpha_diversity_abundance.csv"

# 保存为CSV文件
write.csv(alpha_div, file = file_path, row.names = FALSE)


#差异比较————————————————————————重头开始
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(phyloseq)
library(vegan)
library(data.table)
library(readxl)

asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/rara_abundance/rare_number.csv"

# 定义文件路径

asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/otu_table_rarefied.csv"
tax_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/Trivals/Desktop/usyd_tax_table-010822.csv"
sample_metadata_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/Trivals/Desktop/Matadata1.csv"

# 读取原始ASV表和分类信息的CSV文件,read和fread等函数会自动将第一行识别成列名
asv_data <- fread(asv_file_path, header = TRUE, data.table = FALSE)
tax_data <- fread(tax_file_path, header = TRUE, data.table = FALSE)
sample_data <- fread(sample_metadata_path, header = TRUE, data.table = FALSE)

# 假设OTU列名为"OTU"，将asv_data转换为长格式
otu_long <- asv_data %>%
  pivot_longer(cols = -SampleID, names_to = "OTU", values_to = "Abundance")

# 合并数据
merged_data <- otu_long %>%
  left_join(sample_data, by = "SampleID") %>%
  left_join(tax_data, by = c("OTU" = "SampleID"))


#热图（全部）----------------------------------------------------


# 提取所需列并按组求和
summarized_data <- merged_data %>%
  select(Soil, Storage, Time, Phylum, Abundance) %>%
  group_by(Soil, Storage, Time, Phylum) %>%
  summarise(Sum_Abundance = sum(Abundance, na.rm = TRUE))


# # 对求和后的丰度进行对数变换
# summarized_data <- summarized_data %>%
#   mutate(Log2_Sum_Abundance = log2(Sum_Abundance + 1))
# 
# # 移除 Sum_Abundance 列
# final_data <- summarized_data %>%
#   select(Soil, Storage, Time, Phylum, Log2_Sum_Abundance)


# 创建 heatmap_data
heatmap_data <- summarized_data %>%
  unite("Soil_Storage_Time", Soil, Storage,Time, sep = "_") %>%
  spread(key = Phylum, value = Sum_Abundance) %>%
  as.data.frame()  # 将tibble转换为data.frame

# 将 heatmap_data 保存为 CSV 文件
write.csv(heatmap_data, "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/heatmap_data.csv", row.names = FALSE)

heatmap_data1<- read.csv("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/heatmap_data.csv", header = TRUE, row.names = 1)



#热图(0<genus<0.0001))-------------------------------------------

# merged_data <- merged_data %>%
#   mutate(
#     Class = ifelse(is.na(Class), paste("NA", row_number()), Class),
#     Order = ifelse(is.na(Order), paste("NA1", row_number()), Order),
#     Family = ifelse(is.na(Family), paste("NA2", row_number()), Family),
#     Genus = ifelse(is.na(Genus), paste("NA3", row_number()), Genus),
#     R_Genus = paste(Phylum, Class, Order, Family, Genus, sep = "_")
#   )

 merged_data <- merged_data %>%
   mutate(
     Class = ifelse(is.na(Class), paste("NA", row_number()), Class),
     # Order = ifelse(is.na(Order), paste("NA1", row_number()), Order),
     # Family = ifelse(is.na(Family), paste("NA2", row_number()), Family),
     R_Genus = paste(Phylum, Class, sep = "_")
  )


filtered_data <- merged_data %>% 
  filter(!grepl("NA", R_Genus))

# length(unique(filtered_data$R_Genus))
# 
# length(unique(merged_data$Genus))


genus_abundance_sum <- filtered_data %>%
  group_by(SampleID, Soil, Storage, Time, R_Genus) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

# rare genus
final_data <- genus_abundance_sum %>%
  filter(Abundance>0) %>%
  # filter(Abundance<0.0001) %>%
  # filter(Abundance>0.001) %>%
  group_by(Soil, Storage, Time) %>%
 # summarized_data("")
  arrange(SampleID, Abundance) %>% # 按照Sample排序，然后按照abundance升序排列
  # slice_tail(n =50 ) %>%
  select(Soil, Storage, Time, R_Genus, Abundance)



# 将数据转换为适当的格式以用于ComplexHeatmap
heatmaprare_data <- final_data %>%
  unite("Soil_Storage_Time", Soil, Storage, Time, sep = "_") %>%
  spread(key = R_Genus, value = Abundance) %>%
  as.data.frame()

write.csv(heatmaprare_data, "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/heatmaprare_data.csv", row.names = FALSE)

heatmap_data1<- read.csv("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/heatmaprare_data.csv", header = TRUE, row.names = 1)

write.csv(heatmapabundance_data, "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/heatmapabundance_data.csv", row.names = FALSE)

heatmap_data1<- read.csv("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/heatmapabundance_data.csv", header = TRUE, row.names = 1)

#---------做28day与1day相减
# 读取文件
file_a_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/heatmap/1_day_genus.csv"
file_b_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/heatmap/28_day_genus.csv"

data_a <- read.csv(file_a_path, stringsAsFactors = FALSE)
data_b <- read.csv(file_b_path, stringsAsFactors = FALSE)

# 将NA值替换为0
data_a[is.na(data_a)] <- 0
data_b[is.na(data_b)] <- 0


# 将前四列设置为因子类型
data_a[,1] <- as.factor(data_a[,1])
data_b[,1] <- as.factor(data_b[,1])

# 创建一个新的数据框，先仅包含第一列（假设第一列是你想保留的列名）
df_difference <- data.frame(data_a[,1])

# 计算从第二列开始的差值，并将结果添加到df_difference中
df_difference <- cbind(df_difference, data_b[,2:ncol(data_b)] - data_a[,2:ncol(data_a)])

# 确保第一列的列名保持不变，这一步可能不是必须的，因为cbind操作通常会保留列名
colnames(df_difference)[1] <- colnames(data_a)[1]


# 假设 df_difference 是你的数据框

# 从第二列开始的列索引
start_col <- 2

# 对从第二列开始的每一列进行操作
for (i in start_col:ncol(df_difference)) {
  df_difference[, i] <- sapply(df_difference[, i], function(x) {
    if (x > 0) {
      return(1)
    } else if (x == 0) {
      return(0)
    } else {
      return(-1)
    }
  })
}

# 输出修改后的数据框
write.csv(df_difference, "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/heatmap/df_difference.csv", row.names = FALSE)
heatmap_data1<- read.csv("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/heatmap/df_difference.csv", header = TRUE, row.names = 1)



library(ComplexHeatmap)
library(circlize)


# heatmap_data1  <- heatmap_data1[substring(text = rownames(heatmap_data1),first = 1, last = 1)==1,]


# 将数据框转换为矩阵
heatmap_matrix <- as.matrix(heatmap_data1)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# 为注解定义颜色方案
soil_colors <- c("1" = "#d595a7", "2" = "#cddeb7", "3" = "#612a79")
storage_colors <- c("4" = "#ce3d32", "-20" = "#5db1dd",
                    "RTAD" = "#7a65a5", "RTFM" = "#e4af69")
# time_colors <- colorRamp2(c(0,28), 
#                           c("#FFF7BC", "#FEE391"))

#annotation_df$Soil <- factor(annotation_df$Soil ,levels=c("1","2","3"))

# 自定义 Rare Genus 的颜色映射
# rare_genus_colors <- colorRampPalette(c("#4575B4", "#E0F3F8", "#FEE090", "#D73027"))(100)
# rare_genus_colors <- colorRamp2(c(0,0.001,0.5), c("#4575B4", "#FEE090","#D73027"))
rare_genus_colors <- c("#FEE090", "#E0F3F8", "#4575B4")


# # # 设置 PDF 设备，指定文件保存位置
 # pdf(file = "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/PDF/heatmap1.pdf")

# 从行名中提取Soil和Storage的信息，将Time转换为数值
annotation_df <- data.frame(
  Soil = factor(sub("(.*)_(.*)_(.*)", "\\1", colnames(t(heatmap_matrix)))),
  Storage = factor(sub("(.*)_(.*)_(.*)", "\\2", colnames(t(heatmap_matrix)))),
  # Time = as.numeric(sub("(.*)_(.*)_(.*)", "\\3", colnames(t(heatmap_matrix)))),
  row.names = colnames(t(heatmap_matrix))
)

# # 调整Storage类别的顺序
annotation_df$Storage <- factor(annotation_df$Storage, levels = c("-20", "4", "RTAD", "RTFM"))

# identical(rownames(annotation_df),rownames(heatmap_matrix))

# 创建并绘制热图
#data1=t(heatmap_matrix)


ht <- Heatmap(
  t(heatmap_matrix),
  name = "Phylum Abundance",
 # column_km = 3,
 column_split = annotation_df$Soil,
  rect_gp = gpar(col = "white", lwd = 0.5),
  # use_raster = TRUE,
  show_heatmap_legend =TRUE,
  row_dend_reorder = FALSE,
  cluster_rows = FALSE,
  # clustering_distance_rows = "euclidean",
  # clustering_method_rows = "complete",
  # column_names_rot = 90,
  # column_dend_reorder=FALSE,
  show_row_names =  FALSE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  col = rare_genus_colors,  # 使用自定义颜色
  top_annotation = HeatmapAnnotation(
    df = annotation_df,
    # col = list(Soil = soil_colors, Storage = storage_colors, Time = time_colors),
    # col = list(Storage = storage_colors, Time = time_colors),
    col = list(Storage = storage_colors, Soil = soil_colors),
    annotation_name_side = NULL
  )
)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

print(ht)
# 
length(annotation_df$Soil)
# 创建并绘制热图
ht <- Heatmap(
  t(heatmap_matrix), 
  name = "Phylum Abundance",
  column_km = 3,
  rect_gp = gpar(col = "white", lwd = 0.5),
  use_raster = TRUE,
  show_heatmap_legend = FALSE,
  row_dend_reorder = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  column_dend_reorder = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = rare_genus_colors,  # 使用自定义颜色
  top_annotation = HeatmapAnnotation(
    df = annotation_df,
    col = list(Soil = soil_colors, Storage = storage_colors, Time = time_colors),
    annotation_name_side = NULL
  )
)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
# 
# # 关闭 PDF 设备，保存文件
dev.off()
# 