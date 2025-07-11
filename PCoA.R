rm(list=ls())
options(stringsAsFactors = F)
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)

# setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/PCOA")
# 定义文件路径
asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/rara_abundance/abundance_number.csv"

# asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/CN_cycle.csv"
sample_metadata_path <- "C:/Users/midu6195//OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/sam_table.csv"

# 读取原始ASV表和分类信息的CSV文件,read和fread等函数会自动将第一行识别成列名
asv_data <- fread(asv_file_path, header = TRUE, data.table = FALSE)
sample_data <- fread(sample_metadata_path, header = TRUE, data.table = FALSE)
# 使用subset()函数筛选Soil列为"A"的行
# sample_data <- subset(sample_data, Soil == "B")

#

# # #fread函数不支持指定row.name
# # # 设置第一列为行名
rownames(asv_data) <- asv_data[[1]]
# # # 删除第一列
asv_data <- asv_data[, -1]
# ## 首先，获取满足条件的列名
# columns_to_keep <- grep("^B", rownames(asv_data), value = TRUE)
# 
# # 然后，根据这些列名筛选数据框
# asv_data <- asv_data[columns_to_keep,]
# 
# # 计算每一行的总和
# row_sums <- colSums(asv_data)
# 
# # 筛选出总和不为0的行
# asv_data <- asv_data[,row_sums > 0 ]

# # # 计算矩阵距离
otu.dist1 <- vegdist(asv_data,method="bray")

# 
# otu.dist1 <-dist.mat

# 进行主坐标分析
otu_pcoa1<- cmdscale(otu.dist1,eig=TRUE)
pc121 <- as.data.frame(otu_pcoa1$points[,1:2])

#处理分组数据，与主坐标结果合并
pc121$samples<-rownames(pc121)
colnames(sample_data)[1]<-'samples'
pc<-round(otu_pcoa1$eig/sum(otu_pcoa1$eig)*100,digits = 2)
pc121<-merge(pc121,sample_data,by='samples')
colnames(pc121)[2:3]<-c('PC1','PC2')
# pc121 <- pc121[pc121$Soil == "C", ] #选择soil

# 处理分类变量
# pc121$Time <- factor(pc121$Time, levels = unique(sample_data$Time))
pc121$Storage <- factor(pc121$Storage, levels = unique(sample_data$Storage))
pc121$Soil <- factor(pc121$Soil, levels = unique(sample_data$Soil))
temp<-as.matrix(otu.dist1)
table(colnames(temp)==sample_data$samples)
# pc121$Time <- as.factor(pc121$Time)
pc121$Storage <- as.factor(pc121$Storage)

# 设置随机数种子
set.seed(123)
ADONIS <- adonis2(otu.dist1 ~ Storage * Time, data=sample_data)
ADONIS


mycol<- colorRampPalette(c("#4575B4",   "#D73027"))(5)


 
myshape <- c(15, 16, 17, 18)



# # 定义土壤类型对应的颜色
# 定义土壤类型对应的颜色
soil_colors <- c("A" = "#FFEBAD", "B" = "#C7E9B4", "C" = "#A7C0DE")




p <- ggplot(data = pc121, aes(PC1, PC2, color=Time, shape = Storage)) +
  geom_point(size=3)+
  scale_fill_manual(values = soil_colors)+
  scale_color_gradientn(colours = mycol) +  # 使用连续颜色映射
  scale_shape_manual(values = myshape)+ # 保持形状映射不变
  geom_hline(yintercept=0,linetype=2) +
  geom_vline(xintercept=0,linetype=2)+
  labs(x=paste0("PCoA1(",round(pc[1],2),"%",")"),y=paste0("PCoA2(",round(pc[2],2),"%)"))+
  stat_ellipse(data=pc121,aes(x = PC1, y =PC2, fill =Soil), linetype = 1,
               geom='polygon',alpha=0.3,
               level = 0.99,show.legend = F,inherit.aes = F) +
  guides(fill = guide_legend(title = 'Time',order = 1,
                             override.aes = list(size=4,fill=mycol,shape=21)),
         shape=guide_legend(title='Storage',override.aes = list(size=3)))+
  # ggtitle(label = paste(  'PERMANOVA:F=', Fvalue,', p=',TEST_adonis,sep = ''))+
  theme_bw()+theme(#panel.grid=element_blank(),
    title=element_text(size=10))
p

P <- p + coord_cartesian(xlim = c(-0.15, 0.15), ylim = c(-0.1, 0.1))
P
# # 调整坐标轴的范围(unifrac)
#A
P1 <- p1 + coord_cartesian(xlim = c(-0.25, -0.15), ylim = c(-0.2, -0.1))
#B
P2 <- p2 + coord_cartesian(xlim = c(-0.15, -0.05), ylim = c(0.15, 0.25))
#C
P3 <- p3 + coord_cartesian(xlim = c(0.25, 0.35), ylim = c(-0.1, 0))
print(P2)

# # 调整坐标轴的范围(bray)
#A
p1 <- p1 + coord_cartesian(xlim = c(-0.2, -0.3), ylim = c(-0.3, -0.2))
#B
p2 <- p2 + coord_cartesian(xlim = c(-0.1, -0.2), ylim = c(0.25, 0.35))
#C
p3 <- p3 + coord_cartesian(xlim = c(0.4, 0.5), ylim = c(0, -0.1))
print(p3)


library(gridExtra)
plot <- grid.arrange(p, P1, P2, P3, ncol = 2)
setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/PCOA")

ggsave("PCoA_function.pdf", p , width = 10, height = 8)
#----------------------------------
time_colors <- c("1" = "#4575B4", "3" = "#5370A9", "7" = "#666594", "14" = "#8C5370", "28" = "#D73027")

p <- ggplot(data = pc121, aes(PC1, PC2, color=as.factor(Time), shape = Storage)) +
  geom_point(size=3)+
  scale_fill_manual(values = soil_colors) +
  scale_color_manual(values = time_colors) +  # 使用自定义颜色
  scale_shape_manual(values = myshape) + # 保持形状映射不变
  geom_hline(yintercept=0, linetype=2) +
  geom_vline(xintercept=0, linetype=2) +
  labs(x=paste0("PCoA1(",round(pc[1],2),"%",")"), y=paste0("PCoA2(",round(pc[2],2),"%)")) +
  stat_ellipse(data=pc121, aes(x = PC1, y =PC2, fill = Soil), linetype = 1,
               geom='polygon', alpha=0.3,
               level = 0.99, show.legend = F, inherit.aes = F) +
  guides(fill = guide_legend(title = 'Time', order = 1,
                             override.aes = list(size=4, fill=time_colors, shape=21)),
         shape = guide_legend(title = 'Storage', override.aes = list(size=3))) +
  theme_bw() + 
  theme(title=element_text(size=10))
p
