# ---------------------------------------------
# Date: 01-16-2024
# Written by: Jingzhe Shang
# Summary: Fig1B
# ---------------------------------------------

# Load in libraries
suppressMessages({
    library(ggplot2)
    library(reshape2)
    library(ggpubr)
    library(corrplot)
    library(pheatmap)
})

options(stringsAsFactors = FALSE)
if(!dir.exists("../output/09.acidity")) dir.create("../output/09.acidity")
if(!dir.exists("../output/09.acidity/ph")) dir.create("../output/09.acidity/ph")

acidity_geneset <- read.table("../acidity/acidity_list_final.txt",header = TRUE,sep = "\t" )
genelist <- acidity_geneset

cor.matrix <- data.frame(ph=genelist)
    for (Cancer_Type in Cancer_Types) {
        print(Cancer_Type)
        TPMdata <- load(paste0("../output/00.RawCount/",Cancer_Type,".counts_TPM.RData"))
        Cancer_tcga_sample <- Cancer_tcga_sample[!duplicated(Cancer_tcga_sample$tcga.tcga_barcode),]
        index <- Cancer_tcga_sample$sample_T %in% "Tumor"
        Cancer_tcga_sample <- Cancer_tcga_sample[index,]
        rownames(Cancer_tcga_sample) <- Cancer_tcga_sample$tcga.tcga_barcode
        rownames(TPM) <- TPM$Genename
        TPM <- TPM[,-1]
        ## get gene tpm
        acidity_tpm <- TPM[rownames(TPM) %in% genelist,]
        ## get celltype
        celltype <- read.table(paste0("../output/08.Celltype/tumor/",Cancer_Type,".txt"), sep = "\t")
        rownames(celltype) <- gsub("\\.","-",rownames(celltype))

        index <- intersect(colnames(acidity_tpm), rownames(celltype))
        acidity_tpm <- acidity_tpm[,colnames(acidity_tpm) %in% index]
        celltype <- celltype[index,]
        ## cor
        index <-  match(colnames(acidity_tpm),rownames(celltype))
        acidity_tpm.cell <- rbind(acidity_tpm,t(celltype[index,c(24,25)]))
        acidity_tpm.cell <- acidity_tpm.cell[-nrow(acidity_tpm.cell),]
        acidity_tpm.cell.cor <- cor(t(acidity_tpm.cell),method="pearson")

        pdf(paste0("../output/09.acidity/ph/",Cancer_Type,".gene.cd8.cor.pdf"),width =8,height = 8)
        corrplot(acidity_tpm.cell.cor, type = "lower", order = "original",mar = c(4,4,4,4))
        dev.off()

        acidity_tpm.cell.cor.data <- data.frame(ph=rownames(acidity_tpm.cell.cor),CD8=acidity_tpm.cell.cor[,ncol(acidity_tpm.cell.cor)])
        cor.matrix <- merge(cor.matrix,acidity_tpm.cell.cor.data,by="ph",all=TRUE)
        colnames(cor.matrix)[ncol(cor.matrix)] <- paste0(Cancer_Type,"_CD8")

    }
    rownames(cor.matrix) <- cor.matrix$ph
    cor.matrix <- cor.matrix[rownames(cor.matrix)!="CD8_T",-1]

    pdf(paste0("../output/09.acidity/ph/allcancer.gene.cd8.cor.pdf"),width =8,height = 8)
    bk <- max(abs(cor.matrix))
    pheatmap(cor.matrix,breaks = seq(-bk,bk,length=100),cellheight = 12,
             color = colorRampPalette(c("blue","white","red"))(100) )
    dev.off()
    }
