setwd("C:/Obese_a_practise/")
dr <- read.csv("ense_2_new_2_18.csv")
#dr<- read.csv("ense_2_new_7_29.csv")
dr
head(dr)
drcoldata <- read.csv("ense_2_new_2_18.csv", row.names=1, stringsAsFactors=FALSE)
drcoldata

dim(drcoldata)

head(rownames(drcoldata))

sample<-c(rep("normalW_a",54), rep("obese_a", 47))
sample

metadata<- data.frame(sample_id = colnames(drcoldata), batch = c(5,	1,	2,	1,	2,	1,	2,	1,	1,	2,	2,	2,	1,	2,	2,	2,	5,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	3,	3,	3,	3,	3,	5,	3,	3,	3,	2,	2,	2,	2,	4,	4,	4,	5,	5,	5,	5,	5,	5,	6,	6,	6,	6,	6,	2,	1,	2,	4,	2,	1,	1,	2,	2,	1,	2,	1,	1,	2,	2,	2,	5,	2,	2,	2,	1,	1,	2,	2,	3,	3,	3,	3,	3,	5,	3,	3,	2,	2,	4,	4,	4,	4,	4,	5,	5,	5,	6,	6,	6,	6,	6))
colnames(drcoldata)
metadata
metadata$sample<- relevel(factor(sample),"normalW_a")
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=round(drcoldata), 
                             colData=metadata, 
                             design=~sample + batch)
ncol(drcoldata)
nrow(metadata)
nrow(dds)
rownames(dds)
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)
dds <- DESeq(dds)
#rld <- rlog(dds, blind=FALSE)

vsd <- vst(dds, blind = FALSE)

head(assay(vsd), 3)
sampleDists <- dist( t( assay(vsd) ) )
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <-as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( vsd$sample_id,vsd$sample, sep="-" )
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists
         , col=colors)
plotPCA(vsd, intgroup = c("sample_id","sample"))
plotPCA
library("DESeq2")
dds <- DESeq(dds)
res <- results(dds, contrast=c("sample","normalW_a","obese_a"))

#res <- results(dds)

mcols(res, use.names=TRUE)

summary(res)

res.05 <- results(dds, alpha=.05)

table(res.05$padj < .05)

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
results(dds, contrast=c("sample","normalW_a","obese_a"))
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm=TRUE)

resSig <- subset(res, padj < 0.1)



head(resSig[ order(resSig$log2FoldChange),])
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])

topGene <-rownames(res)[which.min(res$padj)]

plotCounts(dds, gene=topGene, intgroup=c("sample"))

library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res),column="SYMBOL", keytype="ENSEMBL",multiVals="first") 
res$entrez <- mapIds(org.Hs.eg.db, keys=row.names(res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
resOrdered <- res[order(res$padj),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)[1:100,]
#write.csv(resOrderedDF, file="results_2_17.csv1")
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("sample_id","sample","batch")])
pheatmap(mat, annotation_col = anno)

