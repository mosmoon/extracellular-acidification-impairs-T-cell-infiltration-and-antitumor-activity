# ---------------------------------------------
# Date: 01-16-2024
# Written by: Jingzhe Shang
# Summary: preprocessing TCGA data
# ---------------------------------------------

# Load in libraries
suppressMessages({
    library(recount3)
    library(plyr)
})

rawCounts.TCGA <- function(Cancer_Type){
    options(stringsAsFactors = FALSE)
    rowdata <- load(paste0("../output/00.RawCount/",Cancer_Type,".rawdata.RData"))

    counts <- as.data.frame(assays(rse_gene_Cancer_tcga)$counts)

    Cancer_tcga_sample$sample <- substr(Cancer_tcga_sample$tcga.tcga_barcode,14,15)
    Cancer_tcga_sample$sample_T=ifelse(as.numeric(Cancer_tcga_sample$sample) %in% 1:9,"Tumor", "Normal")  # Grouping into Tumor and Normal

    ## get gene information
    Cancer_tcga_gene <- as.data.frame(rowData(rse_gene_Cancer_tcga)@listData)
    Cancer_tcga_gene <- Cancer_tcga_gene[,c(5:8,3)]

    index <- match(colnames(counts),Cancer_tcga_sample$external_id)
    colnames(counts) <- Cancer_tcga_sample$tcga.tcga_barcode[index]
    ## replace gene names
    counts$genesums <- rowSums(counts)
    index <- match(rownames(counts),Cancer_tcga_gene$gene_id)
    counts <- as.data.frame(cbind(Genename = Cancer_tcga_gene$gene_name[index],
                                  ENSEMBL = Cancer_tcga_gene$gene_id[index],
                                  counts))
    counts$Genename <- gsub("\\.\\d+$","",counts$Genename)
    ## Remove duplicate genes and zero values
    counts <- counts[order(counts$genesums,decreasing = TRUE),]
    counts <- counts[counts$genesums>0,]
    counts <- counts[!duplicated(counts$Genename),]
    ## get TPM
    TPM <- as.data.frame(assays(rse_gene_Cancer_tcga)$TPM)
    index <- match(colnames(TPM),Cancer_tcga_sample$external_id)
    colnames(TPM) <- Cancer_tcga_sample$tcga.tcga_barcode[index]
    index <- rownames(TPM) %in% rownames(counts)
    TPM <- TPM[index,]
    ## replace TPM gene names
    index <- match(rownames(TPM),Cancer_tcga_gene$gene_id)
    TPM <- as.data.frame(cbind(Genename = Cancer_tcga_gene$gene_name[index],
                               TPM))
    write.csv(counts, paste0("../output/00.RawCount/",Cancer_Type,".counts.csv"),row.names = FALSE); #save counts
    write.table(TPM, paste0("../output/00.RawCount/",Cancer_Type,".TPM.txt"),
                row.names = FALSE,sep = "\t"); #save TPM
    save(Cancer_tcga_sample,Cancer_tcga_gene,counts,TPM,
         file = paste0("../output/00.RawCount/",Cancer_Type,".counts_TPM.RData") )
}
