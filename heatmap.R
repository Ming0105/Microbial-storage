#差异比较————————————————————————重头开始

library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(phyloseq)
library(vegan)
library(data.table)
library(readxl)

# #percentage
asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/percentage_all.csv"

# rare
asv_file_path <-"C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/rare.csv"

#abundance
asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/abundance.csv"

# 定义文件路径

asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/otu_table_rarefied.csv"
tax_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/Trivals/Desktop/usyd_tax_table-010822.csv"
sample_metadata_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/Trivals/Desktop/Matadata.csv"

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

#PCoA准备 （全部）----------------------------------------------------


# 提取所需列并按组求和
# summarized_data <- merged_data %>%
  # select(Soil, Storage, Time,  Phylum, Abundance) %>%
  # group_by( Storage,Phylum) %>%
  # summarise(Sum_Abundance = sum(Abundance, na.rm = TRUE))



# 创建 heatmap_data
# heatmap_data <- summarized_data %>%
#   unite("Soil_Storage_Time", Soil, Storage, Time) %>%
#   spread(key = Phylum, value = Sum_Abundance) %>%
#   as.data.frame()  # 将tibble转换为data.frame

# # 将 heatmap_data 保存为 CSV 文件
# write.csv(heatmap_data, "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/pcoa_abundance_phylum.csv", row.names = FALSE)
# 
# heatmap_data1<- read.csv("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/heatmap_data.csv", header = TRUE, row.names = 1)
# 


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
    Order = ifelse(is.na(Order), paste("NA1", row_number()), Order),
    # Family = ifelse(is.na(Family), paste("NA2", row_number()), Family),
    R_Genus = paste(Phylum, Class, Order, sep = "_")
  )


filtered_data <- merged_data %>% 
  filter(!grepl("NA", R_Genus))

#只看Storage
#汇总
 # genus_abundance_sum <- merged_data %>%
 #   group_by(SampleID, Soil, Storage, Time, Phylum) %>%
 #   summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
 #   ungroup()

# #将Storage一样的求平均数
#  genus_abundance_sum <- genus_abundance_sum  %>%
#    group_by(Storage,Phylum) %>%
#    summarise(Abundance = mean(Abundance, na.rm = TRUE)) %>%
#    ungroup()
#
#  final_data <- genus_abundance_sum %>%
#    filter(Abundance > 0) %>%
#    arrange(Storage, Abundance) %>%
#    mutate(Storage = factor(Storage)) %>%
#    group_by(Storage) %>%
#    slice_tail(n = 10) %>%
#    select(Storage, Phylum, Abundance)
 
 #正常热图
genus_abundance_sum <- merged_data %>%
  group_by(Storage, Time, Phylum) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup()



# rare genus
final_data <- genus_abundance_sum %>%
  # filter(Abundance>0) %>%
  # filter(Abundance<0.0001) %>%
  filter(Abundance>=0.001) %>%
  group_by(Storage) %>%
  # summarized_data("")
  arrange(SampleID, Abundance) %>% # 按照Sample排序，然后按照abundance升序排列
  # slice_tail(n =10 ) %>%
  select(Soil, Time, Storage, Phylum, Abundance)



# 将数据转换为适当的格式以用于ComplexHeatmap
heatmap_abundance <- genus_abundance_sum %>%
  unite("Soil_Storage_Time", Storage, Time, sep = "_") %>%
  spread(key = Phylum, value = Abundance) %>%
  as.data.frame()


#Just storage, dont match the name
write.csv(heatmap_abundance, "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/heatmap_abundance.csv", row.names = FALSE)

heatmap_data1<- read.csv("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/heatmap_abundance.csv", header = TRUE, row.names = 1)

library(ComplexHeatmap)
library(circlize)


# heatmap_data1  <- heatmap_data1[substring(text = rownames(heatmap_data1),first = 1, last = 1)==1,]


# 将数据框转换为矩阵
heatmap_matrix <- as.matrix(heatmap_data1)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# 为注解定义颜色方案
# soil_colors <- c("1" = "#4575B4", "2" = "#E0F3F8", "3" = "#FEE090")
storage_colors <- c("-20" = "#FC9272", "4" = "#BCBDDC", "RTAD" = "#C994C7", "RTFM" = "#FDBB84")
time_colors <- colorRamp2(c(0,28),
                          c("#FFF7BC", "#FEE391"))


# # 自定义 Rare Genus 的颜色映射
rare_genus_colors <- colorRampPalette(c("#4575B4", "#E0F3F8", "#FEE090", "#D73027"))(100)
# rare_genus_colors <- colorRamp2(c(0,0.02,0.1,0.2), c("#4575B4","#E0F3F8","#FEE090","#D73027"))


# # # 设置 PDF 设备，指定文件保存位置
 # pdf(file = "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/PDF/heatmap1.pdf")

# 从行名中提取Soil和Storage的信息，将Time转换为数值
annotation_df <- data.frame(
  Storage = factor(sub("(.*)_(.*)", "\\1", colnames(t(heatmap_matrix)))),
  Time = as.numeric(sub("(.*)_(.*)", "\\2", colnames(t(heatmap_matrix)))),
  row.names = colnames(t(heatmap_matrix))
)
# annotation_df$Soil <- factor(annotation_df$Soil ,levels=c("1","2","3"))
# # 调整Storage类别的顺序
annotation_df$Storage <- factor(annotation_df$Storage, levels = c("-20", "4", "RTAD", "RTFM"))

# identical(rownames(annotation_df),rownames(heatmap_matrix))

# 创建并绘制热图



ht <- Heatmap(
  t(heatmap_matrix),
  name = "Class Abundance",
  # column_km = 3,
  column_split = annotation_df$Soil,
  rect_gp = gpar(col = "white", lwd = 0.5),
  #use_raster = TRUE,
  show_heatmap_legend =TRUE,
  row_dend_reorder = FALSE,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  column_names_rot = 90,
  column_dend_reorder=FALSE,
  show_row_names =  FALSE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  col = rare_genus_colors,  # 使用自定义颜色
  top_annotation = HeatmapAnnotation(
    df = annotation_df,
    col = list(Storage = storage_colors, Time = time_colors),
    # col = list(Storage = storage_colors, Time = time_colors),
    annotation_name_side = NULL
  )
)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
# 
#Venn图----------------------------------------------

