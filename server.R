# server.R
# version 0.1.4 (DEC 18, 2015)

# version 0.1.5 (FEB.11.2016)
# fix problem...

# if you want run shiny app
# runApp("~/Dropbox/Codes/R/TF/R")
library(shiny)
library(datasets)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(shiny)

library(magrittr)
library(RColorBrewer)
library(pheatmap)

## version 0.1.5

## if you want run shiny app
## runApp("~/Dropbox/Codes/R/TF/R")

# Define server logic required to analysis of gene expression profile


# Define server logic required to summarize and view the selecteds
# dataset
shinyServer(function(input, output) {
  # ./data/hMSC_Cont_VS_LPS_ordered_DEGs_sig.txt
  # upload files
  
  dataset <- reactive({
    
    # input$file1 will be NULL initially. After the user selects and uploads a
    # file, it will be a data frame with 'name', 'size', 'type', and 'datapath'
    # columns. The 'datapath' column will contain the local filenames where the
    # data can be found.
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    # read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    if (input$rowname == TRUE){
    read.csv(inFile$datapath, header=input$header, row.names=1, sep = input$sep, quote=input$quote) #quote = "\""
    }else{
      read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    }
  })
  
  # show only #n rows
  output$contents <- renderTable({
    head(dataset(), n = input$obs)
  })
  
  # Drawing heatmap
  output$heatmap <- renderPlot({
    test <- data.matrix(dataset())
    fixed <- ncol(test)
    # set color
    hmcol = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)
    # set color
    input$color
    if (input$color == 1){
      hmcol = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)
    }else if(input$color == 2){
      hmcol = colorRampPalette((brewer.pal(9, "YlOrRd")))(100)
    }else if(input$color == 3){
      hmcol = colorRampPalette(rev(brewer.pal(9, "Greys")))(100)
    }else {
      hmcol = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100)
    }
    # cluster rows
    hc.rows <- hclust(dist(test))
    #plot(hc.rows)
    # transpose the matrix and cluster columns
    hc.cols <- hclust(dist(t(test)))
    #plot(hc.cols)
    
    ## custom
    #in this case, use only full heatmap
    mycl <- cutree(hc.rows, h=max(hc.rows$height)/5) #0에 가까울수록 많이 쪼개진다.
    mycolhc <- colorRampPalette(brewer.pal(9, "RdPu"))(11) # change to R colorbrowser
    mycolhc <- mycolhc[as.vector(mycl)]
    
    #    pdf("heatmap_IRF1.pdf", width=10, height=20)
    #    pdf("test.pdf", width=10, height=20)
    if (input$fixed == TRUE){
      #test <- test[order(test[,ncol(test)]),]#,decreasing=T),]
      #test <- test[,order[test]]
      heatmap.2(test,  main=input$maintitle, col=hmcol, notecex=1.0, cexCol=1.0, cexRow = 0.85, # decrease font size of row/column labels,
                Rowv=as.dendrogram(hc.rows), Colv=FALSE, density.info = "none", #dendrogram="none",
                #srtCol = 15,adjCol = c(1,1),keysize = 0.8,
                srtCol=360, adjCol = c(0.5,1), #lmat=rbind(c(4,3,0),c(2,1,0)), lwid=c(4, 11, 1 ),
                #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 8, 2 )
                RowSideColors=mycolhc, trace="none", xlab= "samples", ylab= input$ylab) #,symm = T)
    }else{
      heatmap.2(test, main=input$maintitle, col=hmcol, notecex=1.0, cexCol=1.0, cexRow = 0.85, # decrease font size of row/column labels,
              Rowv=as.dendrogram(hc.rows), Colv=as.dendrogram(hc.cols), density.info = "none", #dendrogram="none",
              #srtCol = 15,adjCol = c(1,1),keysize = 0.8,
              srtCol=360, adjCol = c(0.5,1), #lmat=rbind(c(4,3,0),c(2,1,0)), lwid=c(4, 11, 1 ),
              #lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(1.5, 8, 2 )
              RowSideColors=mycolhc, trace="none", xlab= "samples", ylab= input$ylab) #,symm = T)
    }
  })
  
  # Download data
  output$downloadData <- downloadHandler(
    filename = function() { paste(input$file1, '.csv', sep='') },
    content = function(file) {
      pdf(filename, width=10, height=20)
      dev.off()
      #      write.csv(dataset(), file)
    }
  )
  
})


#library(shiny)
#setwd("~/Dropbox/Codes/R/TF")

#library(shiny)
#setwd("~/Dropbox/Codes/R/TF")