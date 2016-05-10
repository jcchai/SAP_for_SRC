## collection of analysis scripts
## version: 0.0.2
## since: 01.NOV.2015
##


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
#library(clusterProfiler)


##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
# read HTseq data by DESeq2
#directory <- getwd()
directory <- ("/home/jcchai/Dropbox/Workspace/Dr.Yang/data/2nd")
directory <- ("/home/jcchai/Dropbox/Workspace/Dr.Yang/data/2nd/UnTX")
#directory <- ("/home/jcchai/Dropbox/Workspace/SHKim/data/LPS_POLYIC/genesymbol_LP")

directory <- ("./data/2nd")
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

sampleCondition <- sub("*_*.count","\\2",sampleFiles)
sampleCondition
condition <- c("Vector_UnTX","Vector_UnTX","Vector_TX","Vector_TX","RuOver_UnTX","RuOver_UnTX","RuOver_TX","RuOver_TX","shRu_UnTX","shRu_UnTX","shRu_TX","shRu_TX")
condition_UnTX = as.factor(condition[c(1,2,5,6,9,10)])

condition_1 <- c("Vector", "Vector","Vector","Vector","RuOver","RuOver","RuOver","RuOver","shRu","shRu","shRu","shRu")
condition_2 <- c("RA", "RA","RA","RV","RV","RV","un-1")

##-----------------------------------------------------------------------------
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

##-----------------------------------------------------------------------------
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

##-----------------------------------------------------------------------------
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

##-----------------------------------------------------------------------------
### heatmap scripts
# draw heatmap
# make heatmap pre-dataset
makeDHdataset <- function(DEGs){
  sig <- as.data.frame(assay(rld)[rownames(DEGs[[1]]),])
  up <- as.data.frame(assay(rld)[rownames(DEGs[[2]]),])
  down <- as.data.frame(assay(rld)[rownames(DEGs[[3]]),])
  dataset <- list("sig" = sig, "up" = up, "down" = down)
  return(dataset)
}

##-----------------------------------------------------------------------------
## draw heatmap


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

##-----------------------------------------------------------------------------
