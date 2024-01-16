# ---------------------------------------------
# Date: 01-16-2024
# Written by: Jingzhe Shang
# Summary: calculate acidity score
# ---------------------------------------------
#
# Load in libraries
suppressMessages({
    library(ggplot2)
    library(reshape2)
    library(ggpubr)
    library(cowplot)
    library(Rmisc)
    library(biomaRt)
})



acidity_score <- function(Cancer) {
    options(stringsAsFactors = FALSE)
    if(!dir.exists("../output/09.acidity")) dir.create("../output/09.acidity")
    if(!dir.exists(paste0("../output/09.acidity/",Cancer))) dir.create(paste0("../output/09.acidity/",Cancer))

    TPMdata <- load(paste0("../output/00.RawCount/",Cancer,".counts_TPM.RData"))
    Cancer_tcga_sample <- Cancer_tcga_sample[!duplicated(Cancer_tcga_sample$tcga.tcga_barcode),]
    rownames(Cancer_tcga_sample) <- Cancer_tcga_sample$tcga.tcga_barcode
    rownames(TPM) <- TPM$Genename
    TPM <- TPM[,-1]
    TPM <- TPM[order(rowMeans(TPM),decreasing = TRUE),]
    TPM <- TPM[!duplicated(rownames(TPM)),]
    #### ####
    ## geneset
    acidity_geneset <- read.table("../acidity/acidity_list_final.txt",header = TRUE,sep = "\t" )

    acidity_tpm <- TPM[rownames(TPM) %in% acidity_geneset$Genes_up,]

    ## get mean score
    gene_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    gene.score <- apply(acidity_tpm,2, gene_mean)
    score.matrix <- data.frame(Cancer=colnames(acidity_tpm),gene.score)
    score.matrix[,2] <- log2(score.matrix[,2]+1)
    write.csv(score.matrix,paste0("../output/09.acidity/",Cancer,"/acidity.matrix.csv"),quote = TRUE)

    ## boxplot
    rownames(score.matrix) <- gsub("\\.","-",rownames(score.matrix))
    index <- match(rownames(score.matrix),Cancer_tcga_sample$Cancer.sample_title)
    score.matrix.sample <- cbind(score.matrix,Cancer_tcga_sample[index,c(1:3)])

    pdf(paste0("../output/09.acidity/",Cancer,"/boxplot.pdf"),width = 8,height = 6)
    p1 <- ggboxplot(score.matrix.sample,x = "tcga.cgc_case_clinical_stage",
                    y = "score",add="jitter")
    plot(p1)
    dev.off()
}

