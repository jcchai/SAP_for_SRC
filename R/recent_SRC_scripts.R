# recent SRC scripts

library("DESeq2")
setwd("/data/works/hs_BM-MSC/hs_BM-MSC_RNA-seq_BT_Macrogen_Mar2016/star-using-UCSC-hg19-genesgtf2015version-bam-mismatch3-may112016/DESeq2")
sampleFiles=grep("tab",list.files("/data/works/hs_BM-MSC/hs_BM-MSC_RNA-seq_BT_Macrogen_Mar2016/star-using-UCSC-hg19-genesgtf2015version-bam-mismatch3-may112016/DESeq2"),value=TRUE)
sampleFiles
sampleCondition<-c("Cont","Cont","Cont","LPS-HI","LPS-HI","LPS-HI","LPS-LOW","LPS-LOW","LPS-LOW","Poly-HI","Poly-HI","Poly-HI")
sampleTable<-data.frame(sampleName=sampleFiles,fileName=sampleFiles, condition=sampleCondition)
sampleTable
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory="/data/works/hs_BM-MSC/hs_BM-MSC_RNA-seq_BT_Macrogen_Mar2016/star-using-UCSC-hg19-genesgtf2015version-bam-mismatch3-may112016/DESeq2", design=~condition)
colData(ddsHTSeq)$condition <-factor(colData(ddsHTSeq)$condition)
dds<-DESeq(ddsHTSeq)
colData(dds)
resMF<-results(dds)
head(resMF)
resMFType_BMMSC_RNAseq_BT___Cont_vs_LPSLOW <- results(dds, contrast=c("condition","LPS-LOW","Cont"))
resMFType_BMMSC_RNAseq_BT___Cont_vs_LPSHI <- results(dds, contrast=c("condition","LPS-HI","Cont"))
resMFType_BMMSC_RNAseq_BT___Cont_vs_PolyHI <- results(dds, contrast=c("condition","Poly-HI","Cont"))
resMFType_BMMSC_RNAseq_BT___LPSHI_vs_PolyHI <- results(dds, contrast=c("condition","Poly-HI","LPS-HI"))
resMFType_BMMSC_RNAseq_BT___LPSLOW_vs_LPSHI <- results(dds, contrast=c("condition","LPS-HI","LPS-LOW"))
resMFType_BMMSC_RNAseq_BT___LPSLOW_vs_PolyHI <- results(dds, contrast=c("condition","Poly-HI","LPS-LOW"))
resMFType_BMMSC_RNAseq_BT___Cont_vs_LPSLOW_ordered <- resMFType_BMMSC_RNAseq_BT___Cont_vs_LPSLOW[order(resMFType_BMMSC_RNAseq_BT___Cont_vs_LPSLOW$padj),]
resMFType_BMMSC_RNAseq_BT___Cont_vs_LPSHI_ordered <- resMFType_BMMSC_RNAseq_BT___Cont_vs_LPSHI[order(resMFType_BMMSC_RNAseq_BT___Cont_vs_LPSHI$padj),]
resMFType_BMMSC_RNAseq_BT___Cont_vs_PolyHI_ordered <- resMFType_BMMSC_RNAseq_BT___Cont_vs_PolyHI[order(resMFType_BMMSC_RNAseq_BT___Cont_vs_PolyHI$padj),]
resMFType_BMMSC_RNAseq_BT___LPSHI_vs_PolyHI_ordered <- resMFType_BMMSC_RNAseq_BT___LPSHI_vs_PolyHI[order(resMFType_BMMSC_RNAseq_BT___LPSHI_vs_PolyHI$padj),]
resMFType_BMMSC_RNAseq_BT___LPSLOW_vs_LPSHI_ordered <- resMFType_BMMSC_RNAseq_BT___LPSLOW_vs_LPSHI[order(resMFType_BMMSC_RNAseq_BT___LPSLOW_vs_LPSHI$padj),]
resMFType_BMMSC_RNAseq_BT___LPSLOW_vs_PolyHI_ordered <- resMFType_BMMSC_RNAseq_BT___LPSLOW_vs_PolyHI[order(resMFType_BMMSC_RNAseq_BT___LPSLOW_vs_PolyHI$padj),]
write.csv(as.data.frame(resMFType_BMMSC_RNAseq_BT___Cont_vs_LPSLOW_ordered), file="hs_BM-MSC_RNA-seq_BT___Cont_vs_LPS-LOW__STAR__mismatch3-genesgtf_DESeq2_May112016.csv")
write.csv(as.data.frame(resMFType_BMMSC_RNAseq_BT___Cont_vs_LPSHI_ordered), file="hs_BM-MSC_RNA-seq_BT___Cont_vs_LPS-HI__STAR__mismatch3-genesgtf_DESeq2_May112016.csv")
write.csv(as.data.frame(resMFType_BMMSC_RNAseq_BT___Cont_vs_PolyHI_ordered), file="hs_BM-MSC_RNA-seq_BT___Cont_vs_Poly-HI__STAR__mismatch3-genesgtf_DESeq2_May112016.csv")
write.csv(as.data.frame(resMFType_BMMSC_RNAseq_BT___LPSHI_vs_PolyHI_ordered), file="hs_BM-MSC_RNA-seq_BT___LPS-HI_vs_Poly-HI__STAR__mismatch3-genesgtf_DESeq2_May112016.csv")
write.csv(as.data.frame(resMFType_BMMSC_RNAseq_BT___LPSLOW_vs_LPSHI_ordered), file="hs_BM-MSC_RNA-seq_BT___LPS-LOW_vs_LPS-HI__STAR__mismatch3-genesgtf_DESeq2_May112016.csv")
write.csv(as.data.frame(resMFType_BMMSC_RNAseq_BT___LPSLOW_vs_PolyHI_ordered), file="hs_BM-MSC_RNA-seq_BT___LPS-LOW_vs_Poly-HI__STAR__mismatch3-genesgtf_DESeq2_May112016.csv")
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
head(assay(rld), 3)
write.csv(assay(rld), file="hs_BM-MSC_RNA-seq_BiolTriple___RLD__STAR__mismatch3-genesgtf_DESeq2__Macrogen_May112016.csv")
p = plotPCA (vsd, intgroup=c("condition"))
p
quit()

