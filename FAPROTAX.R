library(vegan)
library(randomForest)
library("rfPermute")
library(ggplot2)
library(tidyverse)
library(patchwork)
library(tibble)
library(reshape2) # 或者 library(data.table) 根据你选择使用哪个包
library(ggplot2)
library(microeco)
library(export)
library(ggtern)
library(magrittr)
library(ggdendro)
library(data.table)
library(Tax4Fun)
library(ape)

library("GUniFrac")
install.packages("RJSONIO")
install.packages(system.file("extdata", "biom_0.3.12.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "qiimer_0.9.4.tar.gz", package="microeco"), repos = NULL, type = "source")
install.packages(system.file("extdata", "Tax4Fun_0.3.1.tar.gz", package="microeco"), repos = NULL, type = "source")
setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/3.20") 
# 定义文件路径
asv_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/T_asv.csv"
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

setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1")
tree <- read.tree(file = "tree.nwk")

tax_table %<>% tidy_taxonomy

dataset <- microtable$new(sample_table = sample_table,
                          otu_table = feature_table, 
                          tax_table = tax_table,
                          phylo_tree = tree)

# remember first clone the whole dataset
# see https://chiliubio.github.io/microeco_tutorial/notes.html#clone-function
Soil_A<- clone(dataset)
# select 'CW'
Soil_A$sample_table <- subset(Soil_A$sample_table, Storage %in% c("RTAD", "RTFM"))
# or: group_CW$sample_table <- subset(group_CW$sample_table, grepl("CW", Group))
# use tidy_dataset to trim all the basic files
Soil_A$tidy_dataset()
Soil_A 

dataset <- Soil_A  



# # use R subset function to filter taxa in tax_table
# dataset$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
# # another way with grepl function
# dataset$tax_table %<>% .[grepl("Bacteria|Archaea", .$Kingdom), ]
# dataset
# 
# # This will remove the lines containing the taxa word regardless of taxonomic ranks and ignoring word case in the tax_table.
# # So if you want to filter some taxa not considerd pollutions, please use subset like the previous operation to filter tax_table.
# dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
# dataset

dataset$tidy_dataset()
print(dataset)

# dataset$sample_sums() %>% range
# 
# dataset$rarefy_samples(sample.size = 10000)
# 
# dataset$sample_sums() %>% range
# 
# dataset$save_table(dirpath = "basic_files", sep = ",")

# use default parameters
dataset$cal_abund()
## The result is stored in object$taxa_abund ...
# return dataset$taxa_abund
class(dataset$taxa_abund)
## [1] "list"
# show part of the relative abundance at Phylum level
dataset$taxa_abund$Phylum[1:5, 1:5]

# dataset$save_abund(dirpath = "taxa_abund")

# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = TRUE)

# return dataset$alpha_diversity
class(dataset$alpha_diversity)

# save dataset$alpha_diversity to a directory
dataset$save_alphadiv(dirpath = "alpha_diversity")

# unifrac = FALSE means do not calculate unifrac metric
# require GUniFrac package installed
dataset$cal_betadiv(unifrac = TRUE)
# return dataset$beta_diversity
names(dataset$beta_diversity)
# save dataset$beta_diversity to a directory
dataset$save_betadiv(dirpath = "beta_diversity")



test <- clone(dataset)
ncol(test$tax_table)
test$add_rownames2taxonomy(use_name = "OTU")
ncol(test$tax_table)

# # 获取otu_table中存在的ASV ID
# asv_ids_in_otu_table <- rownames(dataset$otu_table)
# #
# # 筛选tax_table以保留这些存在的ASV
# filtered_tax_table <- dataset$tax_table[asv_ids_in_otu_table, , drop = FALSE]
# 
# #更新dataset对象中的tax_table为筛选后的版本
# dataset$tax_table <- filtered_tax_table

#功能，详细版
# create object of trans_func
t2 <- trans_func$new(dataset)
# mapping the taxonomy to the database
# this can recognize prokaryotes or fungi automatically if the names of taxonomic levels are standard.
# for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data
# default database for prokaryotes is FAPROTAX database
t2$cal_spe_func(prok_database = "FAPROTAX")

# calculate the percentages for communities
# here do not consider the abundance
t2$cal_spe_func_perc(abundance_weighted = TRUE)
t2$res_spe_func_perc[1:5, 1:2]

t2$trans_spe_func_perc()
t2$res_spe_func_perc_trans


write.csv(t2$res_spe_func_perc_trans, file = "function.csv") 

names(t2$res_spe_func_perc_trans)


t2$plot_spe_func_perc(group = "storage")

t1 <- trans_func$new(dataset)
# https://chiliubio.github.io/microeco_tutorial/intro.html#tax4fun for the installation description
# and provide the file path of SILVA123
t1$cal_tax4fun(folderReferenceData = "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/SILVA123")
# return two files: t1$tax4fun_KO: KO file; t1$tax4fun_path: pathway file.
t1$tax4fun_KO$Tax4FunProfile[1:5, 1:2]


# must transpose to taxa row, sample column
pathway_file <- t1$tax4fun_path$Tax4FunProfile %>% t %>% as.data.frame
# filter rownames, only keep ko+number
rownames(pathway_file) %<>% gsub("(^.*);\\s.*", "\\1", .)
# load the pathway hierarchical metadata
data(Tax4Fun2_KEGG)
# further create a microtable object, familiar?
func1 <- microtable$new(otu_table = pathway_file, tax_table = Tax4Fun2_KEGG$ptw_desc, sample_table = t1$sample_table)
print(func1)


func1$tidy_dataset()
# calculate abundance automatically at three levels: Level.1, Level.2, Level.3
func1$cal_abund()


# bar plot at Level.1
func2 <- trans_abund$new(func1, taxrank = "Level.2", groupmean = "Storage")
func2$plot_bar(legend_text_italic = FALSE)

func2 <- trans_diff$new(dataset = func1, method = "lefse", group = "Storage", alpha = 0.05, lefse_subgroup = NULL)
func2$plot_diff_bar(threshold = 3, width = 0.8)

#功能.简略版
t1 <- trans_func$new(dataset)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
# use list to prepare data
tmp <- list()
# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
tmp$func <- as.data.frame(t(t1$res_spe_func_perc), check.names = FALSE)
# assign the list as taxa_abund in your microtable object
dataset$taxa_abund <- tmp
# use trans_diff class to perform differential test
t2 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Storage", taxa_level = "all", p_adjust_method = "none")
t2$plot_diff_abund(add_sig = T) + ggplot2::ylab("Relative abundance (%)")
warnings()
 write.csv(tmp$func, file = "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/tmp_func.csv", row.names = T)
 ggsave("Lefse.pdf", p, width = 10, height = 8)
# # 获取OTU ID列表
# otu_ids <- rownames(dataset$otu_table)
# tax_ids <- rownames(dataset$tax_table)
# 
# # 过滤tax_table以保留仅在otu_table中出现的OTU
# filtered_tax_ids <- tax_ids[tax_ids %in% otu_ids]
# 
# # 使用过滤后的OTU列表更新tax_table
# filtered_tax_table <- dataset$tax_table[filtered_tax_ids, ]
# 
# # 更新dataset对象中的tax_table
# dataset$tax_table <- filtered_tax_table
# 
# # 检查更新后的结果
# print(nrow(dataset$tax_table))

# head(dataset@tax_table)


#abundance
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", high_level = "Phylum", ntaxa = 10)
p <- t1$plot_box(group = "Storage", xtext_angle = 30)

ggsave("abundance_top10.pdf", p, width = 10, height = 8)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 30,  groupmean = "Storage")
t1$data_abund 
t1$plot_heatmap(facet = "Sample", xtext_keep = FALSE, withmargin = FALSE,  low = "#E0F3F8",
                high = "#D73027")

ggsave("rare_genus_airdring.pdf", p , width = 10, height = 8)
# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 50, groupmean = "Storage")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))


t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 20, groupmean = "Storage")
g1 <- t1$plot_bar(coord_flip = TRUE)
g1 <- g1 + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
g1
g1 <- t1$plot_bar(clustering_plot = TRUE)
# In this case, g1 (aplot object) is the combination of different ggplot objects
# to adjust the main plot, please select g1[[1]]
g1[[1]] <- g1[[1]] + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
g1
# save the figure
ggsave("test.png", g1, width = 8, height = 5)



dataset1 <- dataset$merge_samples(use_group = "Storage")

t1 <- trans_venn$new(dataset1, ratio = NULL)
t1$plot_venn()
t1 <- trans_venn$new(dataset1, ratio = "numratio")
t1$plot_venn()
# 设置保存图形的文件路径和名称
pdf_file_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/venn_A.pdf"

# 打开PDF图形设备
pdf(file = pdf_file_path, width = 7, height = 7) # 可以根据需要调整尺寸

# 绘制Venn图
t1$plot_venn()

# 关闭PDF图形设备
dev.off()


t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 10, groupmean = "Storage")
t1$plot_donut(label = FALSE)
t1$plot_donut(label = TRUE)


library(ggradar)
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Storage")
t1$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 0.5)
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "Time")
t1$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 0.5)



t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = NULL, groupmean = "Storage1")
plot <- t1$plot_tern()
plot + theme(legend.position="none")
# UpSet plot.
tmp <- dataset$merge_samples(use_group = "Storage")
tmp
t1 <- trans_venn$new(dataset = tmp)

t1$data_summary
# only show some sets with large intersection numbers
t1$data_summary %<>% .[.[, 1] > 20, ]
g1 <- t1$plot_bar(left_plot = TRUE, bottom_height = 0.5, left_width = 0.15, up_bar_fill = "grey50", left_bar_fill = "grey50", bottom_point_color = "black")
g1
# g1 is aplot class and can be saved with ggplot2::ggsave, aplot::ggsave or cowplot::save_plot function
# as g1 is comprised of several sub-plots, please adjust the details for each sub-plot
g1[[1]]
g1[[2]]

dataset1 <- dataset$merge_samples(use_group = "Storage")
t1 <- trans_venn$new(dataset1)

# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
t2 <- t1$trans_comm(use_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample
class(t2)

# calculate taxa abundance, that is, the frequency
t2$cal_abund()
# transform and plot
t3 <- trans_abund$new(dataset = t2, taxrank = "Genus", ntaxa = 8)

t3$plot_bar(bar_type = "part", legend_text_italic = T, xtext_angle = 30, color_values = RColorBrewer::brewer.pal(8, "Set2"),
            order_x = c("-20", "4", "RTFM", "RTAD")) + ylab("Frequency (%)")

setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/三元")
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = NULL, groupmean = "Storage1",  high_level = "Phylum" )
# 只修改需要隐藏的图例
plot <- t1$plot_tern() + 
  guides(color = guide_none())  # 假设taxonomy是通过color控制的
print(plot)

ggsave("rare_genus_4.18.pdf", plot , width = 10, height = 8)
# PCoA

# create an trans_beta object
# measure parameter must be one of names(dataset$beta_diversity)
t1 <- trans_beta$new(dataset = dataset, group = "Storage", measure = "wei_unifrac")

# # PCoA, PCA, DCA and NMDS are available
# t1$cal_ordination(ordination = "PCoA")
# # t1$res_ordination is the ordination result list
# class(t1$res_ordination)
# # plot the PCoA result with confidence ellipse
# t1$plot_ordination(plot_color = "Storage", plot_shape = "Storage", plot_type = c("point", "ellipse"))
# 
# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "anova")
# plot_group_order parameter can be used to adjust orders in x axis
p <- t1$plot_group_distance(boxplot_add = "mean")

ggsave("distance_Storage_.pdf", p, width = 10, height = 8)

#lefse
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Storage", alpha = 0.05, lefse_subgroup = NULL, taxa_level = "genus")
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 3)
# we show 20 taxa with the highest LDA (log10)
t1$plot_diff_bar(use_number = 1:20, width = 0.8, group_order = c("H", "M", "L"))
t1$plot_diff_abund(use_number = 1:20, width = 0.8, group_order = c("H", "M", "L"))
t1$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, clade_label_level = 5, group_order = c("H", "M", "L"))



#env
t1 <- trans_env$new(dataset = dataset, add_data = sample_table[, 4:21])

t1 <- trans_env$new(dataset = dataset, env_cols = 4:21)
t1$cal_diff(group = "Storage", method = "wilcox")
head(t1$res_diff)

t1$cal_diff(method = "anova", group = "Storage")
# place all the plots into a list
tmp <- list()
for(i in colnames(t1$data_env)){
  tmp[[i]] <- t1$plot_diff(measure = i, add_sig_text_size = 5, xtext_size = 12) + theme(plot.margin = unit(c(0.1, 0, 0, 1), "cm"))
}
plot(gridExtra::arrangeGrob(grobs = tmp, ncol = 3))


# 6.1 trans_diff class
setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/lefse")
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Storage", alpha = 0.05, lefse_subgroup = NULL, p_adjust_method = "none")
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 2.5)
# we show 20 taxa with the highest LDA (log10)
p <- t1$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = c("RTAD", "RTFM"))
t1$plot_diff_abund(use_number = 1:30, group_order = c("RTAD", "RTFM"))

ggsave("rare.pdf", p , width = 10, height = 8)















