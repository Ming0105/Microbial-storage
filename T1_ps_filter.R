


library(Rcpp)
library(phyloseq)
library(dada2)
library(stringr)
library(readr)
library(readxl)
library(Biostrings)
library(ggtree)
library(ggplot2)


load("D:/OneDrive - The University of Sydney (Staff)/0_Projects/8_NSW samples/Data/NSW_Transect1/T1.ps.raw.RData")
ps <- T1.ps.raw
ps

table(ps@tax_table[,1])

Bac1 <- ps %>% 
  subset_taxa(Kingdom=="Bacteria")%>%
  subset_taxa(!(Family %in% "Mitochondria") &
              !(Order %in%  "Chloroplast") )
Bac1

Depth <- data.frame(sample_sums(Bac1))

Bac1 <- prune_samples(sample_sums(Bac1)>5000,Bac1)
Bac1 <- prune_taxa(taxa_sums(Bac1)>10,Bac1 )
Bac1<- filter_taxa(Bac1,function(x)sum(x>0)>1,TRUE)
Bac1

########### Tree
library(DECIPHER)
library(phangorn)

Seqs <- getSequences(asv_data@refseq)

alignment <- AlignSeqs(DNAStringSet(Seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)

tree.NJ <- NJ(dm)
fit = pml(tree.NJ, data=phang.align)

fitGTR_0 <- update(fit, k=4, inv=0.2)
fitGTR <-optim.pml(fitGTR_0, model="GTR", optInv=TRUE, optGamma=TRUE, 
                   rearrangement = "stochastic", control = pml.control(trace = 0))

fitGTR$tree
#save(fitGTR, file="T1_Bac_fitGTR.RData")
tree <- fitGTR$tree

Bac_tree <- merge_phyloseq(Bac, phy_tree(tree))
Bac_tree

save(Bac_tree,file="Bac_tree.RData")
#save.image(file="Image_Bac_Tree.Rdata")

###
sam <- data.frame(Bac_tree@sam_data)
colnames(sam)

plot_tree(Bac_tree, "treeonly")

plot_tree(Bac_tree, nodelabf=nodeplotboot(80,0,3), color="Land.use_Field", ladderize="left")

plot_tree(Bac_tree, nodelabf=nodeplotboot(60,60,3), color="Land.use_Field", # shape="Phylum",
          ladderize="left") + 
  coord_polar(theta="y")

plot_tree(Bac_tree, size="abundance", color="Land.use_Field")

### plot

ggtree(tree)
ggtree(tree, layout="circular")

data(chiroptera, package="ape")

tree$tip.label <- Bac1@tax_table[,2]
groupInfo <- tree$tip.label
tree <- groupOTU(tree, groupInfo)
ggtree(tree, aes(color=group),layout='circular') + geom_tiplab(size=1, aes(angle=angle))



