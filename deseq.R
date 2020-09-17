library(DESeq2)
library(tidyr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(wesanderson)
library(gridExtra)
library(VennDiagram)
library(ggtree)
library(ape)
library(pheatmap)
library(circlize)
countdata <- read.table(file="feature_count.txt",header=TRUE)
colnames(countdata)<-c("Geneid","Chr","Start","End","Strand","Length",
                       "BRA337_B1_1","BRA337_B1_2","BRA337_B1_3","BRA337_Coinoculation_1","BRA337_Coinoculation_2","BRA337_Coinoculation_3","BRA337_DAOM197198_1","BRA337_DAOM197198_2","BRA337_DAOM197198_3","CM4574_B1_1","CM4574_B1_2","CM4574_B1_3","CM4574_Coinoculation_1","CM4574_Coinoculation_2","CM4574_DAOM197198_1","CM4574_DAOM197198_2","CM4574_DAOM197198_3","COL2215_B1_1","COL2215_B1_2","COL2215_B1_3","COL2215_Coinoculation_1","COL2215_Coinoculation_2","COL2215_Coinoculation_3","COL2215_DAOM197198_1","COL2215_DAOM197198_2","COL2215_DAOM197198_3")
rownames(countdata)<-gsub("cds_","",as.character(countdata$Geneid))
head(countdata)