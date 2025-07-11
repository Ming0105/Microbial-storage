df1<-read.delim('C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/R_practice/merge_tax_matadata/step2.even.5000.feature-table.txt',header = T,check.names = F)
df2<-read.delim('C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/R_practice/merge_tax_matadata/taxonomy.tsv',header = T)[,1:2]
colnames(df1)[1]<-colnames(df2)[1]<-'Feature.ID'
colnames(df2)[2]<-'taxonomy'
df<-merge(df1,df2,by='Feature.ID',all.x = T)
write.table(df,'C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/R_practice/merge_tax_matadata/step23.feature_tax_table.txt',row.names = F,sep='\t',quote = F)
df3<-read.delim('C:/Users/midu6195/OneDrive - The University of Sydney (Staff)/R_practice/merge_tax_matadata/step23.feature_tax_table.txt',header = T,check.names = F)
