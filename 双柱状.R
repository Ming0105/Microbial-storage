library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(multcomp)
library(ggpubr)
setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/Rresult/PDF/alpha_diversity")
df <- read.csv("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/alpha_diversity_indices_2.19.csv")

df1 <- select(df, Storage, Rare, Abundance)

# 假设df1是你的数据框
# 首先，将数据转换为长格式
df_long <- pivot_longer(df1, cols = c(Rare, Abundance), names_to = "Variable", values_to = "Value")

# # 分别对Rare和Abundance进行ANOVA分析
# anova_results <- df_long %>%
#   group_by(Variable) %>%
#   do(tidy(aov(Value ~ Storage, data = .)))

# 将字母结果添加到原始数据框中
df_long <- df_long %>%
  left_join(letters_results, by = "Storage")
#颜色
col <- c("Rare" = "#f89f68", "Abundance" = "#4b84b3")

# 计算每个Storage和Variable组合的平均值和标准差
df_stats <- df_long %>%
  group_by(Storage, Variable) %>%
  summarise(Mean = mean(Value), SD = sd(Value), .groups = 'drop')

ggplot(df_long, aes(
  x = factor(Storage, levels = unique(Storage)), # 确保Storage的顺序
  y = ifelse(Variable == "Rare", Value, -Value), # 使用调整后的Value
  fill = Variable))+
  geom_col(position = position_dodge(width = 0.8), width = 0.7)+
   geom_errorbar(data = df_stats, aes(ymin = Mean - SD, ymax = Mean + SD, x = factor(Storage, levels = unique(Storage)), group = Variable),
                 position = position_dodge(width = 0.8), width = 0.25) + # 添加误差线# 调整柱状图的位置和宽度
  coord_flip()+  # 坐标翻转
  scale_fill_manual(values = col)+  # 自定义颜色
  theme_bw()+  # 白色背景主题
  theme(panel.grid = element_blank(),  # 移除网格线
        axis.text.x = element_text(color = "black", size = 13),  # 自定义轴文本
        axis.text.y = element_text(color = "black", size = 14, face = "italic"),
        axis.title = element_text(color = "black", size = 16),
        legend.position = "top",  # 图例位置
        legend.title = element_blank(),  # 移除图例标题
        legend.text = element_text(color = "black", size = 15))+
  labs(x=NULL, y="Value", fill="Type") +
  scale_y_continuous(breaks = seq(-500, 1100, 300), labels = as.character(abs(seq(-500, 1100, 300))), limits = c(-500, 1100)) + #设置y轴合适范围，并将负数转化为正数
  geom_hline(yintercept = 0, size = 0.4)
p_test <- ggplot(data = df_stats, aes(x = factor(Storage, levels = unique(Storage)), group = Variable)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD, group = Variable),
                position = position_dodge(width = 0.8), width = 0.25)

# 显示测试图形
print(p_test)

# 保存图像
ggsave("Rare_abudance.pdf", width = 8, height = 6)
