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
countdata <- read.table(file="feature_count_DAOM.txt",header=TRUE)
colnames(countdata)<-c("Geneid","Chr","Start","End","Strand","Length",
                       "BRA337_B1_1","BRA337_B1_2","BRA337_B1_3","BRA337_Coinoculation_1","BRA337_Coinoculation_2","BRA337_Coinoculation_3","BRA337_DAOM197198_1","BRA337_DAOM197198_2","BRA337_DAOM197198_3","CM4574_B1_1","CM4574_B1_2","CM4574_B1_3","CM4574_Coinoculation_1","CM4574_Coinoculation_2","CM4574_DAOM197198_1","CM4574_DAOM197198_2","CM4574_DAOM197198_3","COL2215_B1_1","COL2215_B1_2","COL2215_B1_3","COL2215_Coinoculation_1","COL2215_Coinoculation_2","COL2215_Coinoculation_3","COL2215_DAOM197198_1","COL2215_DAOM197198_2","COL2215_DAOM197198_3")
rownames(countdata)<-gsub("cds_","",as.character(countdata$Geneid))
head(countdata)
COL2215<-countdata[,grep("BRA337",colnames(countdata),invert=T)]
COL2215<-COL2215[,grep("CM4574",colnames(COL2215),invert=T)]
COL2215 <- COL2215[ ,7:ncol(COL2215)]
COL2215 <- as.matrix(COL2215)
conditionCOL2215 <- factor(c(rep("B1", 3), rep("Coinoculation", 3), rep("DAOM197198", 3)))
coldataCOL2215 <- data.frame(row.names=colnames(COL2215), conditionCOL2215)
##library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=COL2215,colData=coldataCOL2215,design=~conditionCOL2215)
dds
dds <- DESeq(dds)
rld <- rlogTransformation(dds)
plotPCA(rld, intgroup="conditionCOL2215")
resCOL2215 <- results(dds)
resCOL2215 <- resCOL2215[order(resCOL2215$padj),]
resdataCOL2215 <- merge(as.data.frame(resCOL2215), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdataCOL2215)[1] <- "Gene"
pvalue0.05_filtered_B1 <- resdataCOL2215[resdataCOL2215$padj<0.05,]
write.table(na.omit(pvalue0.05_filtered_B1),file="0.05_B1_rhiir.txt",sep="\t")
#### ANALYSIs FOR PAIRWISE COMPARISON (DAOM197198 and Co-inoculation)
COL2215DAOM<-COL2215[,grep("B1",colnames(COL2215),invert=T)]
COL2215DAOM <- as.matrix(COL2215DAOM)
head(COL2215DAOM)
condition <- factor(c( rep("Coinoculation", 3),rep("DAOM", 3)))
coldataDAOM <- data.frame(row.names=colnames(COL2215DAOM), condition)
dds <- DESeqDataSetFromMatrix(countData=COL2215DAOM, colData=coldataDAOM, design=~condition)
dds <- DESeq(dds)
rld <- rlogTransformation(dds)
plotPCA(rld, intgroup="condition")
resDAOM <- results(dds)
resDAOM <- resDAOM[order(resDAOM$padj),]
resdataCOL2215DAOM <- merge(as.data.frame(resDAOM), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdataCOL2215DAOM)[1] <- "Gene"
pvalue0.05_filtered_DAOM <- resdataCOL2215DAOM[resdataCOL2215DAOM$padj<0.05,]
write.table(na.omit(pvalue0.05_filtered_DAOM), file="0.05_DAOM_rhiir.txt",sep="\t")