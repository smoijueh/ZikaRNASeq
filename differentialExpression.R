# Samuel Moijueh

# visualization
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# wrangle
library(dplyr)

# differential expression analysis
library(Rsubread)
library(DESeq2)


bam_files <- list.files("ZikaRNASeq_files/data/bam", pattern = ".bam$", full.names=T)
samples <- c("Mock1-1", "Mock2-1", "ZIKV1-1", "VIKV2-1")

# featureCounts
features <- featureCounts(files = bam_files,
              
              # annotation
              annot.ext = "refs/GRCh38.gtf",
              isGTFAnnotationFile = TRUE,
              GTF.featureType = "gene",
              GTF.attrType = "gene_name",
              
              # parameters specific to paired end reads
              isPairedEnd = TRUE,
              requireBothEndsMapped = TRUE,
              checkFragLength = FALSE,
              minFragLength = 50,
              maxFragLength = 600,
              
              # miscellaneous
              nthreads = 8,
              reportReads = TRUE
              )

# count matrix
cts <- features$counts
colnames(cts) <- samples

# set up experimental design
coldata <- matrix(c(rep("Mock",2), rep("Zika", 2), rep("paired-end", 4)), ncol = 2, 
                  dimnames = list(samples,c("condition", "type"))) %>% as.data.frame()

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

# pre-filter for genes that span 10 or more genes
keep <- rowSums(counts(dds)) >= 10; 
dds<- dds[keep,]

# Set the reference to be compared
dds$condition = relevel(dds$condition,"Mock")

# Run DESeq2 and format the results
dds <- DESeq(dds) 
res <- results(dds)

## MA PLOT
plotMA(res, ylim=c(-2,2), main="MA Plot")
legend("topright", 
       legend = c("Significant", "Not Significant"), 
       col = c("red", "slategray"), 
       pch = c(19,19), title = "Significance Level: p<0.1", 
       cex = 1)

# list the most statistically significant genes by the smallest p-value
resOrdered <- res[order(res$padj),]
cutoff = 5*10e-60

## Volcano Plot
resOrdered$Significant <- ifelse(resOrdered$padj < 0.05, "FDR < 0.05", "Not Sig")
l<-as.data.frame(resOrdered); l$Gene<-rownames(l)

q <- l[c(8,2,5,6,7)]; q <- q[complete.cases(q),]

q2 <- ggplot(q, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=Significant)) + ggtitle("Volcano Plot") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = 1, color = "blue", linetype = 2) + 
  geom_vline(xintercept = -1, color = "blue", linetype = 2) + 
  geom_hline(yintercept=60, linetype=2, color = "blue") +
  scale_color_manual(values=c("red", "gray"))

q2+geom_text_repel(data=filter(q, padj<cutoff), 
                   aes(label=Gene))

## PCA 
vsdB <- varianceStabilizingTransformation(dds)
plotPCA(vsdB, intgroup = "condition") + ggtitle("PCA Plot") + theme(plot.title = element_text(hjust = 0.5))

# heatmap
n = 10
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,],
                     resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )

topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','padj')]

hmcol <- brewer.pal(11,'RdBu')
nCounts <- counts(dds, normalized=TRUE)
heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(10,2), 
        main = "Top 10 Up-Regulated and Top 10 Down-Regulated Genes")

legend("topright", 
       legend = c("Up-regulated", "Down-regulated"), 
       col = c(hmcol[11], hmcol[1]), 
       pch = c(15,15), 
       cex = 1.2)

write.csv(as.data.frame(resOrdered), file="deseq2_output.csv")