# ---------------------------------------------
# Date: 01-16-2024
# Written by: Jingzhe Shang
# Summary: FigS1C
# ---------------------------------------------

# Load in libraries
suppressMessages({
    library(pacman)
    library(ggplot2)
    library(reshape2)
    library(ggpubr)
    library(cowplot)
    library(Rmisc)
    library(CePa)
    library(biomaRt)
})

if(!dir.exists("../output/09.acidity")) dir.create("../output/09.acidity")


human_projects <- available_projects(organism = "human")
human_projects_tcga <- human_projects[human_projects$file_source %in% "tcga", ]
Cancer_Types <- human_projects_tcga$project

acidity_score_alltumorplot <- function(Cancer_Types) {

    Cancer_cell <- as.data.frame(matrix(NA,ncol=5,dimnames=list("", c("TCGA","CD8_T","sample_T","score","Cancer_Type"))))
    for (Cancer_Type in Cancer_Types) {
        TPMdata <- load(paste0("../output/00.RawCount/",Cancer_Type,".counts_TPM.RData"))
        Cancer_tcga_sample <- Cancer_tcga_sample[!duplicated(Cancer_tcga_sample$tcga.tcga_barcode),]
        rownames(Cancer_tcga_sample) <- Cancer_tcga_sample$tcga.tcga_barcode
        rownames(TPM) <- TPM$Genename
        TPM <- TPM[,-1]

        matrix.combined <- as.data.frame(read.csv(paste0("../output/10.acidity-16gene/",Cancer_Type,"/acidity.matrix.csv")))
        matrix <- matrix.combined[,-1]
        colnames(matrix) <- c("TCGA","score")
        ## get celltype
        celltype <- read.table(paste0("../output/08.Celltype/tumor/",Cancer_Type,".txt"),
                               sep = "\t")
        rownames(celltype) <- gsub("\\.","-",rownames(celltype))

        index <- intersect(matrix$TCGA, rownames(celltype))
        matrix <- matrix[matrix$TCGA %in% index,]
        celltype <- celltype[index,]
        ## plot cell score
        index <- match(matrix$TCGA,Cancer_tcga_sample$tcga.tcga_barcode)
        matrix.sample <- cbind(matrix,Cancer_tcga_sample[index,c(2,3,4,7)])
        matrix.sample$tcga.cgc_case_clinical_stage <- factor(matrix.sample$tcga.cgc_case_clinical_stage,
                                                             levels = sort(unique(matrix.sample$tcga.cgc_case_clinical_stage)))
        index <- match(matrix.sample$TCGA,rownames(celltype))
        matrix.sample_cell <- cbind(celltype[index,],matrix.sample)

        sub_matrix.sample_cell <- matrix.sample_cell[,c("CD8_T",
                                                        "sample_T","score"),]
        sub_matrix.sample_cell <- cbind(TCGA=rownames(sub_matrix.sample_cell),sub_matrix.sample_cell,Cancer_Type )
        Cancer_cell <- rbind(Cancer_cell,sub_matrix.sample_cell)
    }

    Cancer_cell <- Cancer_cell[-1,]
    pdf(paste0("../output/10.acidity-16gene/alltumor.cell.score.pdf"),width =13,height = 15)

    tumor <- Cancer_cell[Cancer_cell$sample_T %in% "Tumor" ,]
    Normal <- Cancer_cell[Cancer_cell$sample_T %in% "Normal" ,]
    if (nrow(tumor)>1) {
        p0 <- ggplot(tumor,aes(x=CD8_T,y=score))+
            geom_point()+stat_smooth(method = "lm",se = TRUE)+
            ylab("tumor acidity Score")+
            stat_cor(data=tumor, method = "pearson")+theme_sjz()+facet_wrap(~Cancer_Type,ncol = 5,scales = "free_x")
        plot(p0)
    }
    if (nrow(Normal)>1) {
        p1 <- ggplot(Normal,aes(x=CD8_T,y=score,))+
            geom_point()+stat_smooth(method = "lm",se = TRUE)+
            ylab("Normal acidity Score")+
            stat_cor(data=Normal, method = "pearson")+theme_sjz()+facet_wrap(~Cancer_Type,ncol = 5,scales = "free_x")
        plot(p1)
    }

    dev.off()
}

