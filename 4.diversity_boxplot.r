#test
rm(list=ls())
library(reshape2)
library(tidyverse)
library(ggpubr)
library(rstatix)

df1 <- read.csv("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/asv_table.csv")
mygroup <- read.csv("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/sam_table.csv")

# #取不同天数
# # mygroup <- subset(mygroup,Time==0|Time==1|Time==28)
# mygroup <- subset(mygroup,Time==0|Time==28)


colnames(mygroup)
mygroup$Storage<-factor(mygroup$Storage,levels = unique(mygroup$Storage))
df2<-merge(df1,mygroup,by='Label')
myindex<-'Shannon'

df2$new<-df2[,myindex]

# 首先，进行ANOVA

anova_result <- aov(new ~ Storage * Time, data = df2)


# 然后，使用Tukey HSD进行事后测试
tukey_result <- TukeyHSD(anova_result)
#
# 查看Tukey HSD测试结果
tukey_result

stat.test <- df2 %>%
  t_test(new ~ Storage)

#tukey hd

# Box plots
stat.test <- stat.test %>% 
  add_xy_position(x = "Storage", dodge = 0.8,step.increase = 0.1)

stat.test$y.position
mycol<-c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072","#67A9CF","#FD8D3C")
mycol<-c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072","#67A9CF")
Time_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a","#FD8D3C")  # 假设有四种 Time 类别
# mycol<-c("#8DD3C7", "#FFFFB3", "#BEBADA")
stat.test$p.adj.signif<-cut(stat.test$p,breaks =c(0,0.001,0.01,0.05,1),
                        labels = c('***','**','*','ns')  )

# 首先确保 Time 是因子类型
df2$Storage <- as.factor(df2$Storage)

ggboxplot(df2, x = "Storage", y = "new", fill = "Storage",
          add = c('jitter'),shape=21,
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.25,
          add.params = list(size=3,alpha=0.5,width=0.5),
          outlier.shape = NA)+ 
 stat_pvalue_manual(stat.test,   label = "p.adj.signif", tip.length = 0.02,
                    size = 6,
                     hide.ns = T)+
  scale_fill_manual(values = mycol)+
  labs(y = myindex, x = "Day") +  # 设置 x 轴标签为 "Day"
  # scale_color_manual(values = Storage_colors) +  # 为 Storage 的点设置颜色
  labs(y=myindex,x=NULL)+
  theme_bw()+theme(legend.title = element_blank(),
                     legend.position = 'none',
    axis.title = element_text(face = 'bold',size=12))

unique(df2$Storage)

stat.test2<-as.data.frame(apply(stat.test,2,as.character))
write.table(stat.test2,paste0('step4.',myindex,'.boxplot.pdf.stat.xls'),row.names = F,sep = '\t',quote = F)

ggsave(paste0('step4.',myindex,'.boxplot.pdf'),width = 4,height = 4)

#-------------------------折线图

df <- read.csv("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/sam_table.csv")


df1 <- df[,c(2:4,10)]

# 使用 aggregate 计算均值
df2 <- aggregate(df1[, 4], by = list(df1$Time, df1$Storage), FUN = mean)

# # 计算标准误差（SE）
df2$se <- aggregate(df1[, 4], by = list(df1$Time, df1$Storage), FUN = function(x) sd(x) / sqrt(length(x)))$x


colnames(df2)[1:2] <- c("Time","Storage")

# df4 <- rbind(df4[1,],df4)
# 
# df4[1:4,2] <- c(unique(df4$Storage)[-1])


# ggplot(df4[df4$Time!=14,], aes(x = Time, y = x, group = Storage, color = Storage))
ggplot(df2, aes(x = Time, y = x, group = Storage, color = Storage)) +
  geom_line() +  # 绘制折线图
  geom_point() +  # 绘制散点图
  # geom_errorbar(aes(ymin = x - se, ymax = x + se), width = 0.2) +  # 添加SE的误差线
  theme_minimal() +  # 使用简洁主题
  labs(
       x = "Time",
       y = "Shannon",
       color = "Storage") +
  scale_color_brewer(palette = "Set1")  # 使用色彩美观的调色板


#----------------------折线图，alpha change
ggplot(df2, aes(x = Time, y = x, group = Storage, color = Storage)) +
  geom_boxplot() +  
  theme_minimal() +  # 使用简洁主题
  labs(
    x = "Time",
    y = "Shannon",
    color = "Storage") +
  scale_color_brewer(palette = "Set1")  # 使用色彩美观的调色板





#---------------------------------
library(ggplot2)

# 确保 Time 和 Storage 是因子，并设置顺序
df$Time <- factor(df$Time, levels = c(1, 3, 7, 14, 28))
df$Storage <- factor(df$Storage)
df$Soil <- factor(df$Soil)  # Soil 作为形状变量

# 进行双因素方差分析 (Time 和 Storage)
anova_result <- aov(Shannon ~ Time * Storage, data = df)
anova_summary <- summary(anova_result)

# 从 ANOVA 结果中提取 p 值
p_value_time <- anova_summary[[1]]$`Pr(>F)`[1]  # 提取 Time 的 p 值
p_value_storage <- anova_summary[[1]]$`Pr(>F)`[2]  # 提取 Storage 的 p 值
p_value_interaction <- anova_summary[[1]]$`Pr(>F)`[3]  # 提取 Time:Storage 交互项的 p 值

# 创建箱线图，并添加数据点
p <- ggplot(df, aes(x = Time, y = Shannon, fill = Storage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +  # 让箱线图分组
  geom_jitter(aes(fill = Storage, shape = Soil), size = 2.5, color = "black", alpha = 0.8,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +  # 让点仅出现在自己的箱线图上
  labs(x = "Time", y = "Shannon Diversity Index") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  scale_fill_brewer(palette = "Set3") +  # 颜色匹配 Storage
  scale_shape_manual(values = c(21, 22, 23)) +  # Soil 形状 (可调整为其他形状)
  annotate("text", x = 3, y = max(df$Shannon) * 0.998, 
           label = sprintf("ANOVA: Time p = %.3f\nStorage p = %.3f\nInteraction p = %.3f", 
                           p_value_time, p_value_storage, p_value_interaction), 
           size = 4, hjust = 0)

# 显示图
print(p)


