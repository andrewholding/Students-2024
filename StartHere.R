library(DESeq2)
library(EnhancedVolcano)
theme_set(theme_minimal())
setwd("~/Desktop/Students-2024")
countMatrix<-readRDS("countMatrix.RDS")

counts<-countMatrix$counts
samplenames<-matrix(nrow=2,unlist(strsplit(colnames(counts),"\\.")))[1,]
colnames(counts)<-samplenames



coldata = data.frame(condition=paste0(c(rep("Fulvestrant.",8),rep("Vehicle.",8)),rep(c("Hypoxia","Normoxia"),8)),rep=gsub("[^0-9.-]", "", samplenames))
rownames(coldata)<-samplenames
coldata

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds<- DESeq(dds)

#PCA

se<-SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                         colData=colData(dds))
plotPCA( DESeqTransform( se ) )
plotPCA( DESeqTransform( se),  intgroup=c("rep")  )

######
#Rep4 looks like an outlier, removing for now.
#############


removeSamples<-grep("4",colnames(counts), value=TRUE)

filteredCounts<-counts[,colnames(counts)[!colnames(counts) %in% removeSamples]]



coldata2 = data.frame(condition=paste0(c(rep("Fulvestrant.",6),rep("Vehicle.",6)),rep(c("Hypoxia","Normoxia"),6)),rep=gsub("[^0-9.-]", "", colnames(filteredCounts)))
rownames(coldata2)<-colnames(filteredCounts)
coldata2

dds2 <- DESeqDataSetFromMatrix(countData = filteredCounts,
                               colData = coldata2,
                               design = ~ condition)
dds2 <- DESeq(dds2)


#PCA

se2<-SummarizedExperiment(log2(counts(dds2, normalized=TRUE) + 1),
                          colData=colData(dds2))
plotPCA( DESeqTransform( se2) )
plotPCA( DESeqTransform( se2),  intgroup=c("rep")  )


#Plots


resVNvFN<- results(dds2, contrast = c("condition","Vehicle.Normoxia", "Fulvestrant.Normoxia"))
resVNvFN
EnhancedVolcano(resVNvFN,
                lab = rownames(results(dds2)),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "Vehicle.Normoxia - Fulvestrant.Normoxia" 
)
