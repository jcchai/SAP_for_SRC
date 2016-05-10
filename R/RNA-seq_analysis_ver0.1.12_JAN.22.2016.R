## gene expression data analysis scripts for NCCIT
## version 0.1.9
## since: 16.SEP.2015
## last fix:11.NOV.2015

## version 0.1.6
# make drawHeatmap() function and makeDHdataset() function

## version 0.1.7
# successfully get condition information by regular expression

## version 0.1.8
# import all possible comparision function
# rename -e 's/\-/\_/g' *.count
# alsofix new function

## version 0.1.9
# SMOh request order count is wire.e

## version 0.1.12?
# for PCA plot

## data preprocessing







##-----------------------------------------------------------------------------
## load library
library("vsn")
library("DESeq2")
library(sm)
library(plotrix)
library(scatterplot3d)
library(datasets)
library(biomaRt)

# version 2.0
# draw figure for manuscipts.
library(ggplot2)
library(gplots)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(clusterProfiler)

##-----------------------------------------------------------------------------
# set working directory
setwd("./")
#setwd("~/Dropbox/Workspace/Dr.Yang")
setwd("~/Dropbox/Workspace/NCCIT_SMOh/data")

getwd()
# save work directory
save.image("./.RData")
source("./R/utils.R")
##-----------------------------------------------------------------------------
# read HTseq data by DESeq2
#directory <- getwd()
directory <- ("~/Dropbox/Workspace/BV2/for_PCA/")
directory <- ("/home/jcchai/Dropbox/Workspace/NCCIT_SMOh/data/56")

directory

## still test
makeSampleTable <- function(){
  sampleFiles <- grep("count",list.files(directory),value=TRUE)
  sampleName <- sub("*.count","\\1",sampleFiles)
  sampleCondition <- sub(".*_*.gen.count","\\2",sampleFiles)
}
rm(makeSampleTable)
## end still test


sampleFiles <- grep("count",list.files(directory),value=TRUE)
sampleFiles

#sampleName <- sub("*.gen.count","\\1",sampleFiles)
sampleName <- sub("*.count","\\1",sampleFiles)
sampleName

sampleCondition <- sub("*_.?.-.?.ref.count","\\2",sampleFiles)
sampleCondition


##-----------------------------------------------------------------------------
# sampleCondition to cell type??

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable



## set the dataset 

## end set the dataset


## test scripts ??????
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds$condition <- factor(dds$condition, levels=c("untreated","treated"))
dds$condition <- droplevels(dds$condition)
## end test


makeDDSSet <- function(sampleSet){
  "%p%"= function(x,y) paste(x,y,sep="")
  pars <- as.list(match.call()[-1])
  setType <- sub(".*_","\\1",pars[1])
  #sampleTableName <- paste("sampleTable", setType, sep = "_")
  #conditionName <- paste("condition", setType, sep = "_")
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleSet,
                                         directory = directory,
                                         design= ~ condition)
  dds <- DESeq(ddsHTSeq)
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  res <- results(dds)
  ## export expression profiles
  rld <- rlog(dds)                  # 이 데이터
  vsd <- varianceStabilizingTransformation(dds)   # 혹은 이 데이터를 사용한다. 같이 비교하는 plot이 많음.
  rlogMat <- assay(rld)
  vstMat <- assay(vsd)
  write.table(as.data.frame(assay(dds)), file=setType%p%"_exp.txt",sep="\t", quote=FALSE)
  write.table(as.data.frame(rlogMat), file=setType%p%"_exp_rlog.txt",sep="\t", quote=FALSE)
  write.table(as.data.frame(vstMat), file=setType%p%"_exp_vsd.txt",sep="\t", quote=FALSE)
  return(dds)
}

# version 0.1.12
## run

dds <- makeDDSSet(sampleTable)
res <- results(dds, contrast=c("condition","BV2-4hr-Con","BV2-4hr-LPS"))
res
rld <- rlog(dds)

# PCA plot 

plotPCA(rld, intgroup=c("condition"))
# It is also possible to customize the PCA plot using the ggplot function.
data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition)) +
  
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

## end run


##-----------------------------------------------------------------------------
# save data to text file

## make function of make sig, up, down list file
## need to fix the output name
##read dds and make sig,up,down
makeDEGslist <- function(dds, celltype, con1, con2){
  "%p%"= function(x,y) paste(x,y,sep="")
  pars <- as.list(match.call()[-1])
  res_name <- paste(celltype, con1%p%"_"%p%con2,sep = "_")
  # get sig data
  dds_res <- results(dds, contrast=c("condition",con1,con2))
  dds_res_ordered <- dds_res[order(dds_res$padj),]
  write.table(as.data.frame(dds_res_ordered), file=res_name%p%"_ordered.txt",sep="\t", quote=FALSE)
  dds_res_ordered_sig <- subset(dds_res_ordered, padj < 0.1)
  write.table(as.data.frame(dds_res_ordered_sig), file=res_name%p%"_ordered_sig.txt",sep="\t", quote=FALSE)
  # get up regulated data
  dds_res_ordered_sig_up <- subset(dds_res_ordered_sig, log2FoldChange >= 0)
  write.table(as.data.frame(dds_res_ordered_sig_up), file=res_name%p%"_ordered_sig_up.txt",sep="\t", quote=FALSE)
  # get down regulated data
  dds_res_ordered_sig_down <- subset(dds_res_ordered_sig, log2FoldChange < 0)
  write.table(as.data.frame(dds_res_ordered_sig_down), file=res_name%p%"_ordered_sig_down.txt",sep="\t", quote=FALSE)
  return(list("sig" = dds_res_ordered_sig, "up" = dds_res_ordered_sig_up, "down" = dds_res_ordered_sig_down))
}



## include in 0.1.9 version
## SMOh requeset
dds_5
dds_5_E9E_DEG <- makeDEGslist(dds_5, "human","NCCIT_EB","9_19_EB")
dds_5_E9Eet_DEG <- makeDEGslist(dds_5, "human","NCCIT_EB","9_19_EB_EtOH")
dds_5_9E9Eet_DEG <- makeDEGslist(dds_5, "human","9_19_EB","9_19_EB_EtOH")

res <- results(dds_5, contrast=c("condition","NCCIT_EB","9_19_EB"))

# get sig data

dds_res_ordered <- res[order(res$padj),]
write.table(as.data.frame(dds_res_ordered), file="_ordered_2.txt",sep="\t", quote=FALSE)
dds_res_ordered_sig <- subset(dds_res_ordered, padj < 0.1)
write.table(as.data.frame(dds_res_ordered_sig), file=res_name%p%"_ordered_sig.txt",sep="\t", quote=FALSE)

## end request


## include in 0.1.8 version
## still test 11 NOV2015
## make all possible comparision
sampleCondition
unique(sampleCondition)
str(dds_3)
unique(dds_3@colData@listData$condition)

allPC <- function(dds, celltype){
  possibleCondision <- dds@colData@listData$condition
  "%p%"= function(x,y) paste(x,y,sep="")
  pars <- as.list(match.call()[-1])
  tmp <- c()
  for (i in 1:length(unique(possibleCondision))){
    for (j in i+1:length(unique(possibleCondision))){
      if ( !is.na(unique(possibleCondision)[j])) {
        cond <- (paste(unique(possibleCondision)[i], unique(possibleCondision)[j], sep="_"))
        res_name <- paste(pars, unique(possibleCondision)[i]%p%"_"%p%unique(possibleCondision)[j],sep = "_")
        # main process part
        res_tmp <- makeDEGslist(dds, celltype, unique(possibleCondision)[i], unique(possibleCondision)[j])
        tmp <- c(tmp, res_tmp)
      }
    }
  }
  #print("finish!")
  return(tmp)
}
dds_3_PC <- allPC(dds_3,"hs")
dds_3_PC

## test version
allPC <- function(dds, celltype){
  possibleCondision <- dds@colData@listData$condition
  "%p%"= function(x,y) paste(x,y,sep="")
  pars <- as.list(match.call()[-1])
  tmp <- list()
  for (i in 1:length(unique(possibleCondision))){
    for (j in i+1:length(unique(possibleCondision))){
      if ( !is.na(unique(possibleCondision)[j])) {
        #assign(paste(unique(possibleCondision)[i], unique(possibleCondision)[j], sep="_"))
        res_name <- paste(pars, unique(possibleCondision)[i]%p%"_"%p%unique(possibleCondision)[j],sep = "_")
        res_name <- res_name[1]
        assign(pars, res_name)
        tmp <- c(tmp, pars)
        print(res_name)
        # main process part
        # makeDEGslist(dds, celltype, unique(possibleCondision)[i], unique(possibleCondision)[j])
      }
    }
  }
  #print("finish!")
  return(tmp)
}
dds_3_PC <- allPC(dds_3,"hs")
dds_3_PC
## end test version


## end 

## run 
# 3-1
unique(dds_3@colData@listData$condition)
dds_3_E10E_DEG <- makeDEGslist(dds_3, "human","NCCIT_EB","10_2_EB")
dds_3_E10Eet_DEG <- makeDEGslist(dds_3, "human","NCCIT_EB","10_2_EB_EtOH")
dds_3_10E10Eet_DEG <- makeDEGslist(dds_3, "human","10_2_EB","10_2_EB_EtOH")
# 4-1
unique(dds_4@colData@listData$condition)
dds_4_DEG <- makeDEGslist(dds_4, "human","NCCIT_EB","10_2_EB")
dds_4_DEG <- makeDEGslist(dds_4, "human","NCCIT_EB_RA","10_2_EB_RA")
dds_4_DEG <- makeDEGslist(dds_4, "human","NCCIT_EB_RA_EtOH","10_2_EB_RA_EtOH")
dds_4_DEG <- makeDEGslist(dds_4, "human","NCCIT_EB","NCCIT_EB_RA")
dds_4_DEG <- makeDEGslist(dds_4, "human","NCCIT_EB_RA","NCCIT_EB_RA_EtOH")
dds_4_DEG <- makeDEGslist(dds_4, "human","10_2_EB","10_2_EB_RA")
dds_4_DEG <- makeDEGslist(dds_4, "human","10_2_EB_RA","10_2_EB_RA_EtOH")
# 5-1
unique(dds_5@colData@listData$condition)
dds_5_E9E_DEG <- makeDEGslist(dds_5, "human","NCCIT_EB","9_19_EB")
dds_5_E9Eet_DEG <- makeDEGslist(dds_5, "human","NCCIT_EB","9_19_EB_EtOH")
dds_5_9E9Eet_DEG <- makeDEGslist(dds_5, "human","9_19_EB","9_19_EB_EtOH")
# 6-1
unique(dds_6@colData@listData$condition)
dds_6_E9E_DEG <- makeDEGslist(dds_6, "human","NCCIT_EB","9_19_EB")
dds_6_ER9ER_DEG <- makeDEGslist(dds_6, "human","NCCIT_EB_RA","9_19_EB_RA")
dds_6_ERet9ERet_DEG <- makeDEGslist(dds_6, "human","NCCIT_EB_RA_EtOH","9_19_EB_RA_EtOH")
dds_6_EER_DEG <- makeDEGslist(dds_6, "human","NCCIT_EB","NCCIT_EB_RA")
dds_6_ERERet_DEG <- makeDEGslist(dds_6, "human","NCCIT_EB_RA","NCCIT_EB_RA_EtOH")
dds_6_9E9ER_DEG <- makeDEGslist(dds_6, "human","9_19_EB","9_19_EB_RA")
dds_6_9ER9ERet_DEG <- makeDEGslist(dds_6, "human","9_19_EB_RA","9_19_EB_RA_EtOH")
## end run

## get top 25 transcripts and bottom 25 transcripts
## ?? CC request

getTB <- function(DEGs){
  pars <- as.list(match.call()[-1])
  top <- head(DEGs[[2]], 25)
  bottom <- head(DEGs[[3]], 25)
  TB <- rbind(top, bottom)
  write.table(as.data.frame(TB), file=pars$DEGs%p%"_tb.txt",sep="\t", quote=FALSE)
  TB_exp <- as.data.frame(assay(rld)[rownames(TB),])
  #TB_exp <- TB_exp[c(3:6)]
  ## check
  write.table(as.data.frame(TB_exp), file=pars$DEGs%p%"_tb_exp.txt",sep="\t", quote=FALSE)
  hmcol = colorRampPalette((brewer.pal(9, "YlOrRd")))(100)
  # cluster rows
  hc.rows <- hclust(dist(TB_exp))
  # transpose the matrix and cluster columns
  hc.cols <- hclust(dist(t(TB_exp)))
  mycl <- cutree(hc.rows, h=max(hc.rows$height)/5)
  mycolhc <- colorRampPalette(brewer.pal(9, "PRGn"))(11) # change to R colorbrowser
  mycolhc <- mycolhc[as.vector(mycl)]
  pdf(pars$DEGs%p%"_tb_heatmaps.pdf", width=20, height=30)
  heatmap.2(data.matrix(TB_exp), 
            col=hmcol, notecex=1.0, cexCol=1.0, cexRow = 1.0,
            Rowv=as.dendrogram(hc.rows), Colv=as.dendrogram(hc.cols), density.info = "none",
            srtCol=360, adjCol = c(0.5,1), 
            RowSideColors=mycolhc, trace="none", xlab= "samples", ylab= "genesymbol")
  dev.off()
}

getTB(dds_UnTX_VR_DEG)
getTB(dds_UnTX_VS_DEG)
getTB(dds_UnTX_RS_DEG)

getTB(dds_TX_VR_DEG)
getTB(dds_TX_VS_DEG)
getTB(dds_TX_RS_DEG)

getTB(dds_Vector_DEG)

getTB(dds_RuOver_DEG)

getTB(dds_shRu_DEG)

## end CC request

##-----------------------------------------------------------------------------
### heatmap scripts
# draw heatmap
# make heatmap pre-dataset
makeDHdataset <- function(DEGs){
  pars <- as.list(match.call()[-1])
  ddsData <- sub("_DEG", "\\1", pars)
  ddsData <- sub("_[91E].{2,}$","\\1",ddsData)
  rld <- rlog(get(ddsData))
  sig <- as.data.frame(assay(rld)[rownames(DEGs[[1]]),])
  up <- as.data.frame(assay(rld)[rownames(DEGs[[2]]),])
  down <- as.data.frame(assay(rld)[rownames(DEGs[[3]]),])
  dataset <- list("sig" = sig, "up" = up, "down" = down)
  return(dataset)
}

## test
tmps <- "dds_3_10E10Eet_DEG"
ddsData <- sub("_DEG", "\\1", tmps)
ddsData
ddsData2 <- sub("_[1E].{3,}$","\\1",ddsData)
ddsData2
rld <- rlog(get(ddsData))
head(assay(rld))

res_UnTX <- results(dds_UnTX)
rld_UnTX <- rlog(dds_UnTX)
head(assay(rld_UnTX))


## run
# 3-1
dds_3_E10E_DEG_eset <- makeDHdataset(dds_3_E10E_DEG)
dds_3_E10Eet_DEG_eset <- makeDHdataset(dds_3_E10Eet_DEG)
dds_3_10E10Eet_DEG_eset <- makeDHdataset(dds_3_10E10Eet_DEG)
# 4-1
dds_4_E10E_DEG_eset <- makeDHdataset(dds_4_E10E_DEG)
dds_4_ER10ER_DEG_eset <- makeDHdataset(dds_4_ER10ER_DEG)
dds_4_ERet10ERet_DEG_eset <- makeDHdataset(dds_4_ERet10ERet_DEG)
dds_4_EER_DEG_eset <- makeDHdataset(dds_4_EER_DEG)
dds_4_ERERet_DEG_eset <- makeDHdataset(dds_4_ERERet_DEG)
dds_4_10E10ER_DEG_eset <- makeDHdataset(dds_4_10E10ER_DEG)
dds_4_10ER10ERet_DEG_eset <- makeDHdataset(dds_4_10ER10ERet_DEG)
# 5-1
dds_5_E9E_DEG_eset <- makeDHdataset(dds_5_E9E_DEG)
dds_5_E9Eet_DEG_eset <- makeDHdataset(dds_5_E9Eet_DEG)
dds_5_9E9Eet_DEG_eset <- makeDHdataset(dds_5_9E9Eet_DEG)
# 6-1
dds_6_E9E_DEG_eset <- makeDHdataset(dds_6_E9E_DEG)
dds_6_ER9ER_DEG_eset <- makeDHdataset(dds_6_ER9ER_DEG)
dds_6_ERet9ERet_DEG_eset <- makeDHdataset(dds_6_ERet9ERet_DEG)
dds_6_EER_DEG_eset <- makeDHdataset(dds_6_EER_DEG)
dds_6_ERERet_DEG_eset <- makeDHdataset(dds_6_ERERet_DEG)
dds_6_9E9ER_DEG_eset <- makeDHdataset(dds_6_9E9ER_DEG)
dds_6_9ER9ERet_DEG_eset <- makeDHdataset(dds_6_9ER9ERet_DEG)
## end run



##-----------------------------------------------------------------------------
"%p%"= function(x,y) paste(x,y,sep="")

## run
# 3-1
drawHeatmap(dds_3_E10E_DEG_eset)
drawHeatmap(dds_3_E10Eet_DEG_eset)
drawHeatmap(dds_3_10E10Eet_DEG_eset)
# 4-1
drawHeatmap(dds_4_E10E_DEG_eset)
drawHeatmap(dds_4_ER10ER_DEG_eset)
drawHeatmap(dds_4_ERet10ERet_DEG_eset)
drawHeatmap(dds_4_EER_DEG_eset) 
drawHeatmap(dds_4_ERERet_DEG_eset)
drawHeatmap(dds_4_10E10ER_DEG_eset)
drawHeatmap(dds_4_10ER10ERet_DEG_eset)
# 5-1
drawHeatmap(dds_5_E9E_DEG_eset)
drawHeatmap(dds_5_E9Eet_DEG_eset)
drawHeatmap(dds_5_9E9Eet_DEG_eset)
# 6-1
drawHeatmap(dds_6_E9E_DEG_eset)
drawHeatmap(dds_6_ER9ER_DEG_eset)
drawHeatmap(dds_6_ERet9ERet_DEG_eset)
drawHeatmap(dds_6_EER_DEG_eset)
drawHeatmap(dds_6_ERERet_DEG_eset)
drawHeatmap(dds_6_9E9ER_DEG_eset)
drawHeatmap(dds_6_9ER9ERet_DEG_eset)
## end run

drawHeatmap <- function(DEG_list){
  "%p%"= function(x,y) paste(x,y,sep="")
  ## match.call return a call containing the specified arguments 
  ## and the function name also 
  ## I convert it to a list , from which I remove the first element(-1)
  ## which is the function name
  
  pars <- as.list(match.call()[-1])
  for (i in 1:length(DEG_list)){
    
    ## start option of heatmaps
    #hmcol = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)
    hmcol = colorRampPalette((brewer.pal(9, "YlOrRd")))(100)
    
    # cluster rows
    hc.rows <- hclust(dist(DEG_list[[i]]))
    #plot(hc.rows)
    # transpose the matrix and cluster columns
    hc.cols <- hclust(dist(t(DEG_list[[i]])))
    ## CC request
    #hc.cols <- hclust(dist(t(DEG_list[[i]][c(3:6)])))
    ## end CC request
    
    ## custom
    #in this case, use only full heatmap
    mycl <- cutree(hc.rows, h=max(hc.rows$height)/5)
    mycolhc <- colorRampPalette(brewer.pal(9, "PRGn"))(11) # change to R colorbrowser
    mycolhc <- mycolhc[as.vector(mycl)]
    #factor(mycl)
    
    # full heatmap
    pdfName <- paste(pars$DEG_list, names(DEG_list[i]), sep="_")
    pdf(pdfName%p%".pdf", width=20, height=30)
    heatmap.2(data.matrix(DEG_list[[i]]), 
              col=hmcol, notecex=1.0, cexCol=1.0, cexRow = 1.0, # decrease font size of row/column labels,
              Rowv=as.dendrogram(hc.rows), Colv=as.dendrogram(hc.cols), density.info = "none", #dendrogram="none",
              #srtCol = 15,adjCol = c(1,1),keysize = 0.8,
              srtCol=360, adjCol = c(0.5,1), #lmat=rbind(c(4,3,0),c(2,1,0)), lwid=c(4, 11, 1 ),
              #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 8, 2 )
              RowSideColors=mycolhc, trace="none", xlab= "samples", ylab= "genesymbol")
    dev.off()
    
    # darw heatmap by sub group.
    x <- data.matrix(DEG_list[[i]])
    for(j in 1:length(levels(factor(mycl)))){
      ysub <- x[names(mycl[mycl%in%j]),]
      if (is.null(dim(ysub))){
      }else{
        print(j)
        filename <- paste(pdfName,j%p%"_group_matrix.txt",sep="_")
        write.table(ysub,filename, sep="\t", col.names=NA, quote=F)
        if(dim(ysub)[1]>=2){
          hrsub <- hclust(dist(ysub))
          hcsub <- hclust(dist(t(ysub)))
          subpdf<-paste(pdfName,j%p%"_group_matrix.pdf",sep="_")
          pdf(subpdf, width=20, height=30)
          heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=as.dendrogram(hcsub), col=hmcol,
                    srtCol=360, adjCol = c(0.5,1),
                    density.info = "none", trace="none", xlab= "samples", ylab= "Gene symbol") #,symm = T)
          dev.off()
        }
      }
    }
  }
  ## end of heatmap scripts
  print("Finish!")
}



##-----------------------------------------------------------------------------
## Gene Ontology


library(data.table)
library(clusterProfiler)
library(biomaRt)
library(genefilter)
library(GOstats); library(GO.db); library(ath1121501.db)
library(annotate)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

## use DEGs liset dataset
hMSC_Cont_LPS_DEG
hMSC_Cont_Poly_DEG
hMSC_LPS_Poly_DEG


"%p%"= function(x,y) paste(x,y,sep="")
## Start GeneOntology 
for (i in 1:length(RA_un1_DEG)){
  #gs2id <- bitr(row.names(RA_un1_DEG[[i]]), fromType="SYMBOL", toType="ENTREZID", annoDb="org.Hs.eg.db")
  # fromType is ENTREZID, PFAM, IPI, PROSITE, ACCNUM, ALIAS, ENZYME, MAP, PATH, PMID, REFSEQ, SYMBOL, UNIGENE, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, GENENAME, UNIPROT, GO, EVIDENCE, ONTOLOGY, GOALL, EVIDENCEALL, ONTOLOGYALL, OMIM, UCSCKG.
  # gs2id <- bitr(row.names(RA_un1_DEG[[i]]), fromType="REFSEQ", toType="ENTREZID", annoDb="org.Hs.eg.db")
  gs2id <- bitr(row.names(RA_un1_DEG[[i]]), fromType="REFSEQ", toType="ENTREZID", annoDb="org.Mm.eg.db")
  gene <- (gs2id[,2])
  title <- paste("RA_un1_DEG", names(RA_un1_DEG[i]), sep = "_")
  ## end set name
  title
  for ( j in c(3,5,7,9) ){
    ggo <- groupGO(gene     = gene,
                   organism = "human",
                   ont      = "BP",
                   level    = j,
                   readable = TRUE)
    #head(summary(ggo))
    go_list <- summary(ggo)[summary(ggo)[,3]>= 1,]
    go_list_ordered <- go_list[order(go_list$Count, decreasing=T),]
    go_list_ordered$Level <- j
    if(j == 3){
      hMSC_GO_BP <- go_list_ordered
    }
    hMSC_GO_BP <- rbind(hMSC_GO_BP, go_list_ordered)
    write.table(go_list_ordered, file=title%p%"_BP_lv"%p%j%p%".txt",sep="\t", quote=FALSE)
  }
  write.table(hMSC_GO_BP, file=title%p%"_BP_total.txt",sep="\t", quote=FALSE)
  for ( k in c(3,5,7,9)){  
    ggo <- groupGO(gene     = gene,
                   organism = "human",
                   ont      = "MF",
                   level    = k,
                   readable = TRUE)
    #head(summary(ggo))
    go_list <- summary(ggo)[summary(ggo)[,3]>= 1,]
    go_list_ordered <- go_list[order(go_list$Count, decreasing=T),]
    go_list_ordered$Level <- k
    if(k == 3){
      hMSC_GO_MF <- go_list_ordered
    }
    hMSC_GO_MF <- rbind(hMSC_GO_MF, go_list_ordered)
    write.table(go_list_ordered, file=title%p%"_MF_lv"%p%k%p%".txt",sep="\t", quote=FALSE)
  }
  write.table(hMSC_GO_MF, file=title%p%"_MF_total.txt",sep="\t", quote=FALSE)
}
## End GeneOntology



## ???
bitr <- function (geneID, fromType, toType, annoDb, drop = TRUE) 
{
  idTypes <- idType(annoDb)
  msg <- paste0("should be one of ", paste(idTypes, collapse = ", "), 
                ".")
  if (!fromType %in% idTypes) {
    stop("'fromType' ", msg)
  }
  if (!all(toType %in% idTypes)) {
    stop("'toType' ", msg)
  }
  geneID %<>% as.character %>% unique
  db <- get_db_obj(annoDb)
  res <- suppressWarnings(select(db, keys = geneID, keytype = fromType, 
                                 columns = c(fromType, toType)))
  ii <- which(is.na(res[, 2]))
  if (length(ii)) {
    n <- res[ii, 1] %>% unique %>% length
    if (n) {
      warning(paste0(round(n/length(geneID) * 100, 2), 
                     "%"), " of input gene IDs are fail to map...")
    }
    if (drop) {
      res <- res[-ii, ]
    }
  }
  return(res)
}
<environment: namespace:clusterProfiler>
  
  
  
  ##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## gene ontolygy 0
library(data.table)
library(clusterProfiler)
library(biomaRt)
library( "genefilter" )
library(GOstats); library(GO.db); library(ath1121501.db)
library(annotate)
library(org.Hs.eg.db)


#  sig data
hMSC_Cont_LPS_GO_set <- list("sig" = hMSC_Cont_LPS_ordered_sig, "up" = hMSC_Cont_LPS_ordered_sig_up, "down" = hMSC_Cont_LPS_ordered_sig_down)
hMSC_Cont_LPS_GO_set

# hMSC_Cont_Poly_ordered_sig
hMSC_Cont_Poly_GO_set <- list("sig" = hMSC_Cont_Poly_ordered_sig, "up" = hMSC_Cont_Poly_ordered_sig_up, "down" = hMSC_Cont_Poly_ordered_sig_down)
hMSC_Cont_Poly_GO_set

# hMSC_LPS_Poly_ordered_sig
hMSC_LPS_Poly_GO_set <- list("sig" = hMSC_LPS_Poly_ordered_sig, "up" = hMSC_LPS_Poly_ordered_sig_up, "down" = hMSC_LPS_Poly_ordered_sig_down)
hMSC_LPS_Poly_GO_set



"%p%"= function(x,y) paste(x,y,sep="")
## Start GeneOntology 
for (i in 1:length(hMSC_Cont_Poly_GO_set)){
  #for (i in 1:length(hMSC_Cont_Poly_GO_set)){
  #  gs2id <- bitr(row.names(hMSC_Cont_Poly_GO_set[[i]]), fromType="SYMBOL", toType="ENTREZID", annoDb="org.Hs.eg.db")
  # fromType is ENTREZID, PFAM, IPI, PROSITE, ACCNUM, ALIAS, ENZYME, MAP, PATH, PMID, REFSEQ, SYMBOL, UNIGENE, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, GENENAME, UNIPROT, GO, EVIDENCE, ONTOLOGY, GOALL, EVIDENCEALL, ONTOLOGYALL, OMIM, UCSCKG.
  gs2id <- bitr(row.names(hMSC_Cont_Poly_GO_set[[i]]), fromType="REFSEQ", toType="ENTREZID", annoDb="org.Hs.eg.db")
  gene <- (gs2id[,2])
  title <- paste("hMSC_Cont_Poly_GO_set", names(hMSC_Cont_Poly_GO_set[i]), sep = "_")
  ## end set name
  title
  for ( j in c(3,5,7,9) ){
    ggo <- groupGO(gene     = gene,
                   organism = "human",
                   ont      = "BP",
                   level    = j,
                   readable = TRUE)
    #head(summary(ggo))
    go_list <- summary(ggo)[summary(ggo)[,3]>= 1,]
    go_list_ordered <- go_list[order(go_list$Count, decreasing=T),]
    go_list_ordered$Level <- j
    if(j == 3){
      hMSC_GO_BP <- go_list_ordered
    }
    hMSC_GO_BP <- rbind(hMSC_GO_BP, go_list_ordered)
    write.table(go_list_ordered, file=title%p%"_BP_lv"%p%j%p%".txt",sep="\t", quote=FALSE)
  }
  write.table(hMSC_GO_BP, file=title%p%"_BP_total.txt",sep="\t", quote=FALSE)
  for ( k in c(3,5,7,9)){  
    ggo <- groupGO(gene     = gene,
                   organism = "human",
                   ont      = "MF",
                   level    = k,
                   readable = TRUE)
    #head(summary(ggo))
    go_list <- summary(ggo)[summary(ggo)[,3]>= 1,]
    go_list_ordered <- go_list[order(go_list$Count, decreasing=T),]
    go_list_ordered$Level <- k
    if(k == 3){
      hMSC_GO_MF <- go_list_ordered
    }
    hMSC_GO_MF <- rbind(hMSC_GO_MF, go_list_ordered)
    write.table(go_list_ordered, file=title%p%"_MF_lv"%p%k%p%".txt",sep="\t", quote=FALSE)
  }
  write.table(hMSC_GO_MF, file=title%p%"_MF_total.txt",sep="\t", quote=FALSE)
}
## End GeneOntology

###############################################################################
##===
x <- res[rownames(hMSC_Cont_LPS_ordered_sig),]

hMSC_DEGs <- resMF_ordered_sig

hMSC_genelist <- unique(rownames(hMSC_Cont_LPS_ordered_sig))
length(hMSC_genelist)

hMSC_entid <- bitr(hMSC_genelist, fromType="SYMBOL", toType="ENTREZID", annoDb="org.Hs.eg.db")
head(hMSC_entid)
dim(hMSC_entid)

names(hMSC_entid)[1] <-"hgnc_symbol"
#hMSC_gene <- merge(hMSC_genelist, hMSC_entid, by = "hgnc_symbol")
#hMSC_gene

gene <- hMSC_entid[,2]
gene

# 전체 genelist에서 NA 가 아닌것 만 뽑아낸 리스트가 필요
"%p%"= function(x,y) paste(x,y,sep="")
setwd("/Users/jcchai/Dropbox/Public/prof.Chai's lab/SH Kim/hMSC_data/RNAseq_htseq")


## GO analysis each level 
## BP
for ( i in c(3,5,7,9) ){
  ggo <- groupGO(gene     = gene,
                 organism = "human",
                 ont      = "BP",
                 level    = i,
                 readable = TRUE)
  #head(summary(ggo))
  go_list <- summary(ggo)[summary(ggo)[,3]>= 1,]
  go_list_ordered <- go_list[order(go_list$Count, decreasing=T),]
  go_list_ordered$Level <- i
  if(i == 3){
    hMSC_GO_BP <- go_list_ordered
  }
  hMSC_GO_BP <- rbind(hMSC_GO_BP, go_list_ordered)
  write.table(go_list_ordered, file="hMSC_Cont_vs_LPS_DEGs_GO_BP_lv"%p%i%p%".txt",sep="\t", quote=FALSE)
}
write.table(hMSC_GO_BP, file="hMSC_Cont_vs_LPS_DEGs_GO_BP_total.txt",sep="\t", quote=FALSE)

## MF
for ( i in 1:12 ){  
  ggo <- groupGO(gene     = gene,
                 organism = "human",
                 ont      = "MF",
                 level    = i,
                 readable = TRUE)
  #head(summary(ggo))
  go_list <- summary(ggo)[summary(ggo)[,3]>= 1,]
  go_list_ordered <- go_list[order(go_list$Count, decreasing=T),]
  go_list_ordered$Level <- i
  if(i == 1){
    hMSC_GO_MF <- go_list_ordered
  }
  hMSC_GO_MF <- rbind(hMSC_GO_MF, go_list_ordered)
  write.table(go_list_ordered, file="hMSC_Cont_vs_LPS_DEGs_GO_MF_lv"%p%i%p%".txt",sep="\t", quote=FALSE)
}
write.table(hMSC_GO_MF, file="hMSC_Cont_vs_LPS_DEGs_GO_MF_total.txt",sep="\t", quote=FALSE)


d <- c(1, 100, NA, 10)
max(d, na.rm=TRUE)
# If you do want to remove all of the NAs, use this idiom instead:

d <- d[!is.na(d)]


unique((hMSC_DEGs_anno[hMSC_DEGs_anno$fit.cluster ==j])$hgnc_symbol)

#####
library(data.table)
library(clusterProfiler)
"%p%"= function(x,y) paste(x,y,sep="")
gene = rownames(hMSC_Cont_LPS_ordered_sig)
for ( i in c(3,5,7)){
  ggo <- groupGO(gene = gene, organism = "human", ont = "BP",level = i, readable = TRUE)
  #head(summary(ggo))
  go_list <- summary(ggo)[summary(ggo)[,3]>= 1,]
  go_list_ordered <- go_list[order(go_list$Count, decreasing=T),]
  go_list_ordered$Level <- i
  #?go_save <- "cluster_"%p%j%p%"_GO_BP_lv"%p%i
  ## make barplot of GO
  # file_name <- paste(i,"Cont_vs_LPS_Cluster_"%p%j%p%"_GO_BP_lv"%p%i%p%"_barplot.pdf",sep="")
  # pdf(file_name, width=20, height=60)
  # barplot(ggo, drop=TRUE, showCategory= dim(go_list_ordered)[1])
  # dev.off()
  ## merge all GO BP results
  if(i == 1){
    hMSC_GO_BP <- go_list_ordered
  }
  hMSC_GO_BP <- rbind(hMSC_GO_BP, go_list_ordered)
  write.table(go_list_ordered, file="hMSC_Cont_vs_LPS_GO_BP_lv"%p%i%p%".txt",sep="\t", quote=FALSE)
}


##-----------------------------------------------------------------------------
## ETC

## make new function
hello<-function() {
  cat(paste("Hello, ",system("whoami",T),"!\n",sep="",collapse=""))
  today<-as.list(unlist(strsplit(system("date",T)," ")))
  names(today)<-c("year","month","day","week","time","timezone")
  return(today)
}
hello()
today <- hello()
today$year
##-----------------------------------------------------------------------------
# Effects of transformations on the variance
library("vsn")
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))

library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))

heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))



##-----------------------------------------------------------------------------
## add 13OCT2015
## Gene Ontology enrichment analysis

source("http://bioconductor.org/biocLite.R")
#A set of annotation maps describing the entire Gene Ontology
biocLite("GO.db")
biocLite("topGO")
biocLite("GOstats")


#install if you haven't already
source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
biocLite("GO.db")
library(org.Hs.eg.db)
library(GO.db)

#for reference's sake
#how to get GO terms from an Entrez gene id following the manual
#http://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf

#org.Hs.egGO is an R object that provides
#mappings between entrez gene identifers and the GO
#identifers that they are directly associated with

entrez_object <- org.Hs.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(entrez_object)
# Convert to a list
entrez_to_go <- as.list(entrez_object[mapped_genes])
#http://www.ncbi.nlm.nih.gov/gene/?term=1
#entrez gene id 1 is A1BG
entrez_to_go[[1]]

#map GO terms to Entrez gene ids

go_object <- as.list(org.Hs.egGO2EG)
#test
go_object[1]
#same
go_object['GO:0000002']


#GO:0000002 -> mitochondrial genome maintenance
#entrez gene id 291 -> SLC25A4
#which if you look up has the GO term mitochondrial genome maintenance

#now let's get all the entrez gene ids for
#GO:0007411 (axon guidance)

axon_gene <- go_object['GO:0007411']
length(unlist(axon_gene, use.names=F))

length(unique(unlist(axon_gene, use.names=F)))
axon_gene <- unique(unlist(axon_gene, use.names=F))
head(axon_gene)

library("GO.db")
library("GOstats")

#as the universal list, I will use all the genes with GO terms

universe <- mapped_genes
length(axon_gene)
length(universe)

params <- new('GOHyperGParams',
              geneIds=axon_gene,
              universeGeneIds=universe,
              ontology='BP',
              pvalueCutoff=0.001,
              conditional=F,
              testDirection='over',
              annotation="org.Hs.eg.db"
)
hgOver <- hyperGTest(params)

hgOver
Gene to GO BP  test for over-representation 
4207 GO BP ids tested (997 have p < 0.001)
Selected gene set size: 319 
Gene universe size: 14844 
Annotation package: org.Hs.eg

result <- summary(hgOver)
head(result,20)

#install if necessary
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

library("biomaRt")
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

my_chr <- c(1:22, 'M', 'X', 'Y')
my_ensembl_gene <- getBM(attributes='ensembl_gene_id',
                         filters = 'chromosome_name',
                         values = my_chr,
                         mart = ensembl)

#how many entries
length(my_ensembl_gene)
[1] 52322

five_ensembl <- my_ensembl_gene[1:5,]

five_ensembl_to_entrez <- getBM(attributes=c('ensembl_gene_id', 'entrezgene'), filters = 'ensembl_gene_id', values = five_ensembl, mart = ensembl)
five_ensembl_to_entrez
ensembl_gene_id entrezgene
1 ENSG00000036549      26009
2 ENSG00000037637      54455
3 ENSG00000221126         NA
4 ENSG00000225011         NA
5 ENSG00000228176         NA

#out of interest how many entrez gene ids?
my_entrez_gene <- getBM(attributes='entrezgene',
                        filters = 'chromosome_name',
                        values = my_chr,
                        mart = ensembl)

length(my_entrez_gene[,1])


#get some more info on the entrez_gene
my_attribute <- c('entrezgene',
                  'hgnc_symbol',
                  'chromosome_name',
                  'start_position',
                  'end_position',
                  'strand')

my_entrez_gene_info  <- getBM(attributes=my_attribute,
                              filters = c('entrezgene', 'chromosome_name'),
                              values = list(entrezgene=my_entrez_gene$entrezgene, chromosome_name=my_chr),
                              mart = ensembl)

head(my_entrez_gene_info)
dim(my_entrez_gene_info)

#entrez gene id on different chromosome locations
my_entrez_gene_info[my_entrez_gene_info$entrezgene=='728369',]


