# spearman corr, remove normal, scale/center


options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)

matrix = read.table("matrix_top38_ES.csv",sep=",",header=T,row.names=1)
