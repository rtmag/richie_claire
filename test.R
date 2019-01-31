# spearman corr, remove normal, scale/center


options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

matrix = read.table("matrix_top38_ES.csv",sep=",",header=T,row.names=1)

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))

heatmap.2(as.matrix(matrix),col=colors,scale="row", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl")




# 
source("https://raw.githubusercontent.com/rtmag/refactor/master/R/refactor.R")

dec_matrix = read.table("raw_all.csv",sep=",",header=T,row.names=1)
all_microbes = dec_matrix[,1:36]

all_microbes_tumor = dec_matrix[,1:18]

# To tell how many microbes have no enrichmen score in any of the tumor samples run:
table(apply(all_microbes_tumor,1,var)==0)
# 144 microbes are not present in any tumor sample

# to select microbes with present in at least 1 tumor sample run:
microbes_in_tumor = all_microbes_tumor[apply(all_microbes_tumor,1,var)!=0, ]

apply(microbes_in_tumor,1,function(x) x>0)
      
rowSums(apply(microbes_in_tumor,1,function(x) x>0))
              
tumor_matrix_38 = matrix[1:18,]
              
heatmap.2(t(as.matrix(matrix)),col=colors,scale="column", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
labRow = T,xlab="", ylab="Microbes",key.title="microbe")

 pdf("heatmap_38microbe_tumor.pdf",height=10,width=7)
              tumor_matrix_a = matrix[1:18,]              
heatmap.2(t(as.matrix(tumor_matrix_a)),col=colors,scale="column", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
key.title="microbe",cexRow=.7,margins=c(8,12))
dev.off()
              
 pdf("heatmap_38microbe_normal.pdf",height=10,width=7)
tumor_matrix_b = matrix[19:36,]              
heatmap.2(t(as.matrix(tumor_matrix_b)),col=colors,scale="column", trace="none",distfun = function(x) get_dist(x,method="pearson"),srtCol=90,
ylab="Microbes",key.title="microbe",cexRow=.7,margins=c(8,12))
dev.off()
              
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))
# TOUATI's Figure
  pdf("heatmap_38microbe_all.pdf",height=10,width=7)
heatmap.2(t(as.matrix(matrix)),col=colors,scale="row", trace="none",srtCol=90,
ylab="Microbes",key.title="microbe",cexRow=.7,margins=c(8,12))
dev.off()
              
# difference
x = all_microbes_tumor[apply(all_microbes_tumor,1,var)>0.00001,]
              
              heatmap.2(t(as.matrix(x)),col=colors,scale="row", trace="none",srtCol=90,
ylab="Microbes",key.title="microbe",cexRow=.7,margins=c(8,12))
                           
# 
