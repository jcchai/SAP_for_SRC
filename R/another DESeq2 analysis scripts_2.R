#
# http://uconn.edu/
# Computational Biology Core
# http://bioinformatics.uconn.edu/rnaseq-tutorial/

sickle pe \
-f Read1_R1.fastq.gz \
-r Read1_R2.fastq.gz \
-t sanger \
-o Read1_R1_trimmed.fastq.gz \
-p Read1_R2_trimmed.fastq.gz \
-q 35 \
-l 45 \
-n \
-s Read_trimmed_singles.fastq.gz
Alignment of Trimmed Reads Using STAR

After your script finishes running, we will be aligning the trimmed reads using STAR. The files that we created with sickle -o and -p are now input under “–readFilesIn.” The “–genomeDir” command must reference a directory that contains annotations in .gtf format. STAR will map the supplied reads to the genome, and create a SAM or BAM type file. SAM is human readable, while BAM is not, but is smaller due to being in binary. The BAM file will be used in HTSeq to create our raw counts.

STAR --genomeDir /your/genome/directory/ \
--readFilesIn Read1_R1_trimmed.fastq.gz Read1_R2_trimmed.fastq.gz
--runThreadN 8 \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix Your_Prefix
Convert BAM File Format to Raw Counts with HTSeq

Finally, we will use HTSeq to transform these mapped reads into counts that we can analyze with R. “-s” indicates we do not have strand specific counts. “-r” indicates the order that the reads were generated, for us it was by alignment position. “-i” indicates what feature we will be using from the GFF, here it is gene name. We identify that we are pulling in a .bam file (“-f bam”) and proceed to identify, and say where it will go.

htseq-count -s no \
-r pos \
-i gene_name \
-f bam \
Your_Prefix.Aligned.sortedByCoord.out.bam \
features.gtf > your_prefix-output_basename.counts
Analyze Counts with DESeq2

The R code below is long and slightly complicated, but I will highlight major points. This script was adapted from here and here, and much credit goes to those authors. Their code is much better commented and explained there, but my example below does come with some additional trimming of graphs I found unimportant. Some important notes:
  
  The gist of this experiment was to compare four different tissue samples from a human patient with two types of cancer. One comparison, shown below, is the fibroblast sample versus the primary tumor sample. The fibroblast sample was taken from similar, but unaffected tissue.
The most important information comes out as “-results-with-normalized.csv,” there we can see adjusted and normal p-values, as well as log2foldchange for all of the genes.
par(mar) manipulation is used to make the most appealing figures, but these values are not the same for every display or system or figure. Much documentation is available online on how to manipulate and best use par() and ggplot2 graphing parameters.
## DESeq2 Analysis of 113 Reads
## Author: Taylor Falk
## Wegrzyn Lab, University of Connecticut

library("DESeq2")
setwd("~/Documents/R/FibrovTumor/")
directory <- "~/Documents/R/FibrovTumor/"

sampleFiles<- c("113-Fibroblast-Passage-5_GTCCGC_L001_001-output_basename.counts",
                "113TuLob-primary-tumor_AGTCAA_L001_R1_001-output_basename.counts")
condition <- c("untreated","treated")
PCANames <- c("Fibroblast","Primary Tumor")
sampleTable <- data.frame(sampleName =c("Fibroblast", "Primary Tumor"), fileName = sampleFiles, condition = condition)


ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, 
                                      levels = c('untreated','treated'))

#guts
dds <- DESeq(ddsHTSeq)
res <- results(dds)

#order by padj value (most sig to least)
res <- res[order(res$padj),]
#save data results and nromalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
head(resdata)
write.csv(resdata, file = "FibrovTumorDESeq2-results-with-normalized.csv")

#send normalized count to tab delimited
mcols(res, use.names = TRUE)
write.csv(as.data.frame(mcols(res, use.name = T)),file = "FibrovTumorDESeq2-test-conditions.csv")

#produce dFrame of results of stats tests
mcols(res, use.names = TRUE)
write.csv(as.data.frame(mcols(res, use.name = T)),file = "FibrovTumorDESeq2-test-conditions.csv")

#replacing outlier val with estimated value as predicted by distribution using "trimmed mean"
#approach. rec if several replicates per treatment
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial= results(dds)$padj < 0.1, cleaned = results(ddsClean)$padj < 0.1)
addmargins(tab)
write.csv(as.data.frame(tab), file = 'FibrovTumorDESeq2-replaceoutliers.csv')
resClean <- results(ddsClean)
resclean <-resClean[order(resClean$padj),]
head(resClean)
write.csv(as.data.frame(resClean), file = "FibrovTumorDESeq2-replaceoutliers-results.csv")

##Exploritory data analysis of RNASeq data w/ DESeq2
# MA Plot, rlog stabilization and var stabilization, variance stabilization, heatmap, PCA plot

mar.old <- par('mar')
print(mar.old)

par(mar=c(5.1,4.1, 2.0, 1.0))
png("FibrovTumorDESeq2-MAplot_initial_analysis.png")
plotMA(dds, ylim = c(-8,8), main = "RNAseq Experiment")
par(mar=mar.old)
dev.off()

#normalize raw counts with rlog and variance stabilization (better for heatmaps)
rld <- rlogTransformation(dds, blind = TRUE)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
#save em!
write.table(as.data.frame(assay(rld),file = "FibrovTumorDESeq2-rlog-transformed-counts.txt", sep = '\t'))
write.table(as.data.frame(assay(vsd),file = "FibrovTumorDESeq2-vst-transformed-counts.txt", sep = '\t'))

#plot to show transformation effect
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1, length(ord), length = 50)]
last <-ord[length(ord)]
vstcol <- c('blue','black')
mar.old <- par('mar')
print(mar.old)
par(mar=c(5.1,4.1, 2.0, 1.0))
png("FibrovTumorDESeq2-variance_stabilizing.png")
matplot(px[ord], cbind(assay(vsd)[,1],log2(px))[ord, ], type = 'l', lty = 1, col = vstcol, xlab = 'n', ylab = 'f(n)', 
        main = "Variance Stabilizing Transformation")
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
par(mar=mar.old)
dev.off()

#clustering analysis
library("gplots")
library("RColorBrewer")
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition,sampleFiles , sep=" : "))
#Or if you want conditions use:
#rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
png("FibrovTumorDESeq2-clustering.png")
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))
dev.off()

#Principal components plot shows additional but rough clustering of samples
library("grDevices")
library("genefilter")

rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))
condition <- c("condition1", "condition2")
scores <- data.frame(sampleFiles, pc$x, condition)

(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(PCANames)))) 
+ geom_point(size = 5)
+ ggtitle("Principal Components")
+ scale_colour_brewer(name = " ", palette = "Set1")
+ theme(
  plot.title = element_text(face = 'bold'),
  legend.position = c(.85,.87),
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

ggsave(pcaplot,file="FibrovTumorDESeq2-PCA_ggplot2.png")

# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples
# not applicable for my data, could be for yours
# head(assay(rld))
# plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
# plot(assay(rld)[,1:2],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
# plot(assay(rld)[,3:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
# plot(assay(rld)[,5:6],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
# plot(assay(rld)[,7:8],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
# plot(assay(rld)[,9:10],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
# plot(assay(rld)[,11:12],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
# plot(assay(rld)[,13:14],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
# plot(assay(rld)[,15:16],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

# heatmap of data
library("RColorBrewer")
library("gplots")
# 1000 top expressed genes with heatmap.2
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:1000]
my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
png("FibrovTumorDESeq2-HEATMAP.png")
heatmap.2(assay(vsd)[select,], col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F,
          main="Fibroblast vs Tumor Heatmap")
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
png("FibrovTumorDESeq2-heatmap3.png")
heatmap3(assay(vsd)[select,],col=my_palette,
         labRow = F,cexCol = 0.8,margins=c(6,6))
dev.off()

mar.old <- par('mar')
print(mar.old)
par(mar=c(5.1,4.1, 2.0, 1.0))
png("FibrovTumorDESeq2-Dispersion-plot1.png")
plotDispEsts(dds, main="Dispersion Estimates")
par(mar=mar.old)
dev.off()

sessionInfo()
Merging Data and Using Biomart

To continue with analysis, we can use the .csv files we generated from the DeSEQ2 analysis and find gene ontology. This next script contains the actual biomaRt calls, and uses the .csv files to search through ensembl. If you are trying to search through other datsets, simply replace the “useMart()”  command with the dataset of your choice. Again, the biomaRt call is relatively simple, and this script is customizable in which values you want to use and retrieve.

setwd("~/Documents/R/Biomart")
library("biomaRt")
#convert .csv into .txt
list.files()
csv_files <- list.files()[grep(".csv",list.files())] 
txt_files <- c()
for (i in 1:length(csv_files)){
  write.table(read.csv(csv_files[i]), gsub(".csv",".txt",csv_files[i]))
  txt_files <- c(txt_files, gsub(".csv",".txt",csv_files[i]))
}
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
#choose the dataset you want to explore
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# input the results that include the IDs
for (i in 1:length(txt_files)){
  a <- read.table(txt_files[i], head=TRUE)
  #with GO
  b <- getBM(attributes=c('hgnc_symbol',
                          'ensembl_gene_id',
                          'ensembl_transcript_id',
                          'ensembl_peptide_id',
                          'band',
                          'chromosome_name',
                          'start_position',
                          'end_position',
                          'transcript_start',
                          'transcript_end',
                          'description',
                          'go_id',
                          'name_1006',
                          'transcript_source',
                          'status'), 
             filters=c('hgnc_symbol'), 
             values=a[1:dim(a)[1],3], 
             mart=ensembl)
  write.table(b,file=gsub(".txt","Biomart-with-GO.txt",txt_files[i]))
  
  #without GO
  b <- getBM(attributes=c('hgnc_symbol',
                          'ensembl_gene_id',
                          'ensembl_transcript_id',
                          'ensembl_peptide_id',
                          'band',
                          'chromosome_name',
                          'start_position',
                          'end_position',
                          'transcript_start',
                          'transcript_end',
                          'description',
                          'transcript_source',
                          'status'), 
             filters=c('hgnc_symbol'), 
             values=a[1:dim(a)[1],2], 
             mart=ensembl)
  write.table(b,file=gsub(".txt","Biomart-no-GO.txt",txt_files[i]))
}
After we create the .txt files with GO terms and without them, we are ready to merge everything together based on the HGNC symbol. The for loop below pulls in the .csv files and combines them with .txt files. Some data may need to be manipulated in your case, but the biomaRt documentation is well detailed in the inputs it accepts.

##Biomart DESeq2 Data merger

setwd("~/Documents/R/Biomart")
(file_list <- list.files())
csv_files <- file_list[grep(".csv",file_list)]
no_go_txts <- file_list[grep("GO.txt",file_list)]

#set up names for files
manipulation_names <- gsub(".csv","",csv_files)
file_output_names <- gsub(".csv","-MERGED.csv",csv_files)

#begin the loop
for(k in 1:length(file_output_names)){
  temp_no_go_name <- file_list[grep(paste(manipulation_names[k],"Bio",sep = ""),file_list)]
  temp_csv_name <- file_list[grep(paste(manipulation_names[k],".csv", sep = ""),file_list)]
  
  temp_no_go_file <- read.table(temp_no_go_name)
  temp_csv_file <- read.csv(temp_csv_name)
  names(temp_csv_file)[2] <- "hgnc_symbol"
  temp_merger <- merge(temp_csv_file, temp_no_go_file, by = "hgnc_symbol")
  
  #remove unnecessary cols
  temp_merger$X <- NULL
  temp_merger$transcript_source <- NULL
  temp_merger$status <- NULL
  
  #remove all but one duplicate
  for (j in names(temp_merger)[2:length(names(temp_merger))]){
    for (i in 1:length(temp_merger[[j]])){
      if (duplicated(temp_merger[[j]])[i]){
        temp_merger[[j]][i] <- ""
        as.factor(temp_merger[[j]][i])
      }
    }
  }
  write.csv(temp_merger, file = file_output_names[k])
}

