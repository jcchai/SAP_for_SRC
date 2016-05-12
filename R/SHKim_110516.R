## SHKim 11/5/16

dataset <- read.csv("./data/MSC_data.txt",header= T, sep ="\t")
dataset <- read.csv("./data/Poly HI_Top 50.txt",header= T, sep ="\t")
dataset <- read.csv("./data/LPS HI_Top 50.txt",header= T, sep ="\t")


dataset
dataset[1]

class(dataset)

dataset_ <- dataset[,-1]
rownames(dataset_) <- dataset[,1]
dataset <- dataset_


library(gplots)
library(RColorBrewer)

## start option of heatmaps
#hmcol = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)
hmcol = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100)

# cluster rows
hc.rows <- hclust(dist(dataset))
#plot(hc.rows)
# transpose the matrix and cluster columns
hc.cols <- hclust(dist(t(dataset)))
## CC request
#hc.cols <- hclust(dist(t(dataset[c(3:6)])))
## end CC request

## custom
#in this case, use only full heatmap
mycl <- cutree(hc.rows, h=max(hc.rows$height)/5)
mycolhc <- colorRampPalette(brewer.pal(9, "PRGn"))(11) # change to R colorbrowser
mycolhc <- mycolhc[as.vector(mycl)]
#factor(mycl)

# full heatmap
#pdfName <- paste("heatmap.pdf"", names(]), sep="_")
pdf("heatmap2.pdf", width=20, height=30)
heatmap.2(data.matrix(dataset), 
          col=hmcol, notecex=1.0, cexCol=1.0, cexRow = 1.0, # decrease font size of row/column labels,
          Rowv=NULL, #as.dendrogram(hc.rows), 
          Colv=as.dendrogram(hc.cols), 
          density.info = "none", #dendrogram="none",
          #srtCol = 15,adjCol = c(1,1),keysize = 0.8,
          srtCol=360, adjCol = c(0.5,1), #lmat=rbind(c(4,3,0),c(2,1,0)), lwid=c(4, 11, 1 ),
          #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 8, 2 )
          RowSideColors=mycolhc, trace="none", xlab= "samples", ylab= "genesymbol")
dev.off()



#__뚱__

drawHeatmap <- function(DEG_list){
  "%p%"= function(x,y) paste(x,y,sep="")
  ## match.call return a call containing the specified arguments 
  ## and the function name also 
  ## I convert it to a list , from which I remove the first element(-1)
  ## which is the function name
  
  pars <- as.list(match.call()[-1])
  for (i in 1:length(DEG_list)){
    
   
    # darw heatmap by sub group.
    x <- data.matrix(dataset)
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
