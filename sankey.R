library(data.table)
library(readxl)
library(dplyr)
library(ggplot2)
library(networkD3)
# 定义文件路径
Retained_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/liuxian/raw.csv"

# 读取数据
df <- fread(Retained_file_path, header = TRUE, data.table = FALSE)

# 将数据从宽格式转换为长格式
df_long <- df %>%
  pivot_longer(cols = c(Retained, Lost), names_to = "Type", values_to = "Value") %>%
  mutate(Type = factor(Type, levels = c("Lost", "Retained")),
         Time = factor(Time, levels = unique(df$Time)))  # 将Time转换为因子，并保证按数据中出现的顺序排列

# 绘制图形
# 绘制图形
p <- ggplot(df_long, aes(x = Time, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") + # 绘制堆叠柱状图
  facet_grid(Storage ~ ., scales = "free_x", space = "free_x") + # 竖直排列Storage
  labs(x = "Time", y = "Total", fill = "Component Type") + # 设置标签
  theme_minimal() + # 使用简洁的主题
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1), # 添加方框以区分Storage
    strip.background = element_rect(fill = "lightblue", colour = "black", size=1), # 美化标签背景
    strip.text.x = element_text(size = 12, face = "bold"), # 调整标签字体大小和加粗
    panel.grid.major = element_blank(), # 去除主要网格线
    panel.grid.minor = element_blank(), # 去除次要网格线
    axis.text.x = element_text(angle = 90, vjust = 0.5) # 如有必要，调整x轴标签的显示角度
  ) +
  scale_x_discrete(drop = FALSE)  # 确保x轴只包含存在于数据集中的时间点

# 输出图形
print(p)


# 输出图形
print(p)









if (!require("ggalluvial")) install.packages("ggalluvial")
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(tidyr)

# 调整示例数据
df_long <- df %>%
  pivot_longer(cols = c(Retained, Lost), names_to = "Type", values_to = "Value") %>%
  mutate(Type = factor(Type, levels = c("Lost", "Retained")))  # 确保Lost在Retained之上

library(ggplot2)
library(ggalluvial)

p <- ggplot(data = df_long,
       aes(axis1 = Time, axis2 = Type, y = Value)) +
  geom_alluvium(aes(fill = Type)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  facet_grid(~Storage) +  # 竖直方向上分面
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    strip.background = element_rect(fill = "grey", colour = "black", size = 1),  # 设置分面标签的背景和边框颜色
    strip.text = element_text(face = "bold", size = 12),  # 设置分面标签的文字样式
    panel.spacing = unit(2, "lines"),  # 增加分面之间的间距
    panel.border = element_rect(colour = "gray", fill = NA, size = 1),  # 设置分面之间边框的颜色和宽度
    axis.text.x = element_blank(),  # 去除x轴刻度标签
    axis.ticks.x = element_blank()  # 去除x轴刻度线
  ) +
  labs(x = "Time", y = "Total", fill = "Component Type")

p
setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/liuxian")
ggsave("Sankey.pdf", p, width = 15, height = 8)

