#test
###written by 983499199@qq.com; Huimin Zhang
##https://shop113546122.taobao.com/
###video reccord 
rm(list=ls())
library(ggvenn)
library(RColorBrewer)
library(data.table)
#最多4个
library(ggVennDiagram)


asv_file_path <-"C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/2.19/abudance_species.csv"

# 定义文件路径
asv_file_path <-"C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/asv_table.csv"
sample_metadata_path <- "C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy/sam_table.csv"

# 读取原始ASV表和分类信息的CSV文件,read和fread等函数会自动将第一行识别成列名
asv_data <- fread(asv_file_path, header = TRUE, data.table = FALSE)
sample_data <- fread(sample_metadata_path, header = TRUE, data.table = FALSE)


# 转置ASV表
asv_data_transposed <- t(asv_data)
colnames(asv_data_transposed) <- asv_data_transposed[1, ]
asv_data_T <- asv_data_transposed[-1, ]
asv_data_T <- as.data.frame(asv_data_T) #data frame

#筛选表格
#不分soil选，2个venn
sample_data <- subset(sample_data, Time %in% c(1,28) & (Storage == "-20"))
asv_data_T <- asv_data_T[, grep("-20", colnames(asv_data_T))]


# # #筛选表格
# # #分soil选
# sample_data <- subset(sample_data,Time==0&Soil=="A"|Time==28&Soil=="A")
# asv_data_T <- asv_data_T[, grep("B0|B28RTAD|B28-20|B28RTFM|^B284$", colnames(asv_data_T), perl = TRUE)]

# #两两对比
# sample_data <- subset(sample_data,Storage=="4"&Soil=="B"|Storage=="RTFM"&Soil=="B"|Storage=="-20"&Soil=="B")
# asv_data_T <- asv_data_T[, grep("^B.*4|-20|RTFM", colnames(asv_data_T), perl = TRUE)]
# 
# #两两对比
# sample_data <- subset(sample_data,Storage=="RTAD"&Soil=="B"|Storage=="RTFM"&Soil=="B")
# asv_data_T <- asv_data_T[, grep("^B.*RTAD|RTFM", colnames(asv_data_T), perl = TRUE)]

#两两对比
# 选择 sample_data 中名字为 "B1RTAD" 和 "B28 RTAD" 的行


# 选择 asv_data_T 中只带有 "B1RTAD" 和 "B28 RTAD" 的列
# asv_data_T <- asv_data_T[, grep("1-20$|^28-20$", colnames(asv_data_T))]



# #不分
# sample_data <- subset(sample_data,Time==0|Time==28)
# 
# 
# 只选soil
# sample_data <- subset(sample_data, Storage =="RTAD")
# asv_data_T <- asv_data_T[, grep("^RTAD", colnames(asv_data_T))]


df<-asv_data_T[,sample_data$Label]

colnames(df)<-sample_data$Time

#变成numetric, 并且加上列名
df1 <- data.frame(apply(df,MARGIN = 2,FUN = function(x) as.numeric(x)))
row.names(df1) <- row.names(df)
outord1<-sample_data$Time[!duplicated(sample_data$Time)]
df2<-as.data.frame(t(apply(df1, 1, function(x) tapply(x, colnames(df), mean))))
outord1_char <- as.character(outord1)  # 将数字索引转换为字符型
df2_selected <- df2[, outord1_char, drop = F]  # 使用字符型索引选择列
df3<-cbind(c(rownames(df2)),df2)
colnames(df3)[1]<-'tax'

# #变成numetric, 并且加上列名
# df1 <- data.frame(apply(df,MARGIN = 2,FUN = function(x) as.numeric(x)))
# row.names(df1) <- row.names(df)
# outord1<-sample_data$Time[!duplicated(sample_data$Time)]
# df2<-as.data.frame(t(apply(df1, 1, function(x) tapply(x, colnames(df), mean))))
# df2<-df2[,outord1,drop=F]
# df3<-cbind(c(rownames(df2)),df2)
# colnames(df3)[1]<-'tax'


write.table(df3,'step22.merge_for_venn.txt',row.names = F,sep = '\t',quote = F)

###
data<-read.delim("step22.merge_for_venn.txt",header = T)

colnames(data)[1]<-'tax'

data[,-1]<-lapply(data[,-1],function(x) ifelse(x>0,1,NA))
setlist<- lapply(data[ -1 ], function(i) unique(data$tax[ !is.na(i) ]))

# mycol<-c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FD8D3C')
mycol<-c( "#d9e7f2", "#eadff0")
# mycol<-c( "#ffd7d8", "#d8f2e7", "#fff2cd")
   # mycol<-c("#b4d9e5", "#91a1cf", "#716bbf","#5239a3")
library(VennDiagram)
#最多5组
# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

# 
setwd("C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/PhD project/Research-1-soil DNA extraction/4.1/rarefy")
# 
pdf('-2028.venn.pdf',height = 8,width = 8)
 display_venn(setlist,
                    lwd = 3,
                    # fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                    # Numbers
                    cex = 1.5,
                    alpha = 0.5,
                    # Set names
                    cat.cex = 2,
                    cat.fontface = "bold",
                    cat.default.pos = "outer",
                    cat.col = mycol,
                    fill = mycol,
                    col = mycol)

# 
dev.off()  # 关闭PDF文件设备，图形保存到 "venn_diagram.pdf" 文件中


print()

library(gplots)
v.table <- venn(setlist)
lengths(attributes(v.table)$intersections)
each_part<-attributes(v.table)$intersections
maxl <-max(sapply(each_part,length))
sets_table <- sapply(each_part,function(x) c(x,rep("",maxl-length(x))))
write.table(sets_table,paste0('step22.venn.pdf.txt'),sep = "\t",row.names=FALSE,quote=FALSE)

