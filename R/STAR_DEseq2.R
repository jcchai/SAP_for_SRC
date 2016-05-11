#  Alignment generate genome. 
STAR --runMode genomeGenerate \
--genomeDir /path/to/STAR_mm10/ \
--genomeFastaFiles ~/mm10_index/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa \
--runThreadN 16

#  Alignment to genome (or transcriptome).
STAR --genomeDir /path/to/STAR_mm10/ \
--readFilesIn read_1.fq read_2.fq \
--runThreadN 16 \
--outSAMtype BAM SortedByCoordinate 

# Sorting, munging, etc. 
samtools view -uS Aligned.out.sam -b | samtools sort -@ 8 -m 800000000 - Aligned.sorted
samtools view Aligned.sorted.bam > Aligned.sorted.sam

# raw counts generation using HTSeq-count. 
htseq-count -m union \
-r pos \
-i gene_name \
-a 10 \
--stranded=no \
Aligned.sortedByCoord.out.bam \
~/mm10_index/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf > output_basename.counts

# differential gene expression analysis using DESeq2. 
####################################################################################
# bioinformatic analysis of RNA-seq data using DESeq2
#
# combined from DESeq2 manual and vignette
# http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# and from Dave Wheeler's blog at Massey University
# http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
#
# analysis per htseq count -intersections-nonempty
# controlling for batch (different RNAseq prep libraries)
#  # notes:
#
####################################################################################

library("DESeq2")
setwd("~/path/to/working/directory/")
directory <- "/path/to/counts/directory/"

# can merge individual sample files (i.e. ctrl1.counts, ctrl2.counts, etc.)
sampleFiles <- grep("counts"list.files(directory),value=T)

# view sampleFiles
sampleFiles

# can designate different batches of samples (i.e. different sequencers,
# PE vs SE, different library preps (eg. BATCH1 vs BATCH2))
sampleBatch <- c("Batch1","Batch1","Batch1","Batch1","Batch1","Batch1",
                 "Batch2","Batch2","Batch2","Batch2")

# set sampleConditions and sampleTable for experimental conditions
sampleCondition <- c("Control, Experimental")

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition,
                          Batch = sampleBatch)

# view sampleTable
sampleTable 

ddsHTseq <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable,
                                        directory = directory,
                                        design = ~condition) 

## view ddsHTseq - should give summary of class, data, etc.
ddsHTseq

colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition,
                                      levels = c('Control','Experimental'))

# gut of DESeq2 analysis
dds <- DESeq(ddsHTseq)

res <- results(dds)
# order results by padj value (most significant to least)
res <- res[order(res$padj),]
head(res)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 

# save data results and normalized reads to csv!
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file="DATE-DESeq2-results-with-normalized.csv")

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = 'DATE_DESeq2_normalized_counts.txt', sep = '\t')

# produce DataFrame of results of statistical tests
# could way to record experimental design
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = "DATE-DESeq2-test-conditions.csv")

# replacing outlier value with estimated value as predicted by distrubution using
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates

ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.1,
             cleaned = results(ddsClean)$padj < 0.1)
addmargins(tab)
write.csv(as.data.frame(tab),file = 'DATE-DESeq2-replaceoutliers.csv')
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
write.csv(as.data.frame(resClean),file = 'DATE-DESeq2-replaceoutliers-results.csv')

####################################################################################
# Exploritory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) MA plot
# 2) rlog stabilization and variance stabiliazation
# 3) variance stabilization plot
# 4) heatmap of clustering analysis
# 5) PCA plot
#
#
####################################################################################

# MA plot of RNAseq data for entire dataset
# http://en.wikipedia.org/wiki/MA_plot
# genes with padj < 0.1 are colored Red
plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")
dev.copy(png, "DATE-DESeq2_MAplot_initial_analysis.png")
dev.off()

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld),file='DATE-DESeq2-rlog-transformed-counts.txt', sep='\t'))
write.table(as.data.frame(assay(vsd),file='DATE-DESeq2-vst-transformed-counts.txt', sep='\t'))

# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,"DATE-DESeq2_variance_stabilizing.png")
dev.off()

# clustering analysis
library("gplots")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, type, sep=" : "))
# OR
# if you want the conditions used
# rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13, 13))
dev.copy(png, "DATE-DESeq2-clustering.png")
dev.off()

# Principal components plot
# will show additional clustering of samples
# showing basic PCA function in R from DESeq2 package
# this lacks sample IDs and only broad sense of sample clustering
# its not nice - but it does the job
print(plotPCA(rld, intgroup = c("condition")))
dev.copy(png, "DATE-DESeq2_PCA_initial_analysis.png")
dev.off()

# or ggplot PCA plot
library("grDevices")
library('ggplot2')
library("genefilter")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsdMF)[select,]))

# set condition
condition <- c("condition1", 'condition2')

scores <- data.frame(sampleFiles, pca$x, condition)

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
+ geom_point(size = 5)
+ ggtitle("Principal Components")
+ scale_colour_brewer(name = " ", palette = "Set1")
+ theme(
plot.title = element_text(face = 'bold'),
legend.position = c(0,0),
legend.key = element_rect(fill = 'NA'),
legend.text = element_text(size = 10, face = "bold"),
axis.text.y = element_text(colour = "Black"),
axis.text.x = element_text(colour = "Black"),
axis.title.x = element_text(face = "bold"),
axis.title.y = element_text(face = 'bold'),
panel.grid.major.x = element_blank(),
panel.grid.major.y = element_blank(),
panel.grid.minor.x = element_blank(),
panel.grid.minor.y = element_blank(),
panel.background = element_rect(color = 'black',fill = NA)
))

ggsave(pcaplot,file="DATE-PCA_ggplot2.pdf")

# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
head(assay(rld))
plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:2],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,3:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,5:6],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,7:8],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,9:10],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,11:12],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,13:14],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,15:16],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

# heatmap of data
library("RColorBrewer")
library("gplots")
# 1000 top expressed genes with heatmap.2
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:1000]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
heatmap.2(assay(vsd)[select,], col=my_palette,
scale="row", key=T, keysize=1, symkey=T,
density.info="none", trace="none",
cexCol=0.6, labRow=F,
main="TITLE")
dev.copy(png, "DATE-DESeq2-HEATMAP.png")
dev.off()

# top 2000 genes based on row variance with heatmap3
library(heatmap3)
colsidecolors = c("darkgreen","darkgreen",
      "mediumpurple2","mediumpurple2",
      "mediumpurple2","mediumpurple2",
      "darkgreen","darkgreen",
      "mediumpurple2","mediumpurple2",
      "darkgreen","darkgreen",
      "mediumpurple2","mediumpurple2",
      "mediumpurple2","mediumpurple2")
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=T)[seq_len(min(2000,length(rv)))]
my_palette <- colorRampPalette(c("blue", "white", "red"))(1024)
heatmap3(assay(vsd)[select,],col=my_palette,
labRow = F,cexCol = 0.8,margins=c(6,6))
dev.copy(pdf, "DATE-DESeq2-heatmap3.pdf")
dev.off()

sessionInfo()

###############################################################
#
# Optional analyses
#
###############################################################

# multifactor designs
# can analysis with more than one factor influencing the counts  (e.g., different library preps, PE vs. SE)
# from manual section 1.5 for more information and rationale

# make copy of DESeq2 results
ddsMF <- dds

# change design formulate controlling for Batch
design(ddsMF) <- formula(~ Batch + condition)

# rerun DESeq analysis
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)

# order by padj values
resMF <- resMF[order(resMF$padj),]

head(resMF)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj
# save data 'resMF' to csv!
write.csv(as.data.frame(resMF),file='DATE-DESeq2_batchcontroll_initial_analysis.csv')