# ---------------------------------------------
# Date: 01-16-2024
# Written by: Jingzhe Shang
# Summary: Fig2C_2D
# ---------------------------------------------

# Load in libraries
suppressMessages({
    library(gdata)
    library(ggplot2)
    library(reshape2)
    library(pheatmap)
    library(ggpubr)
    library(cowplot)
    library(common.sjz)
    library(fgsea)
})

kegg_GO_pathway.gsea <- function(setname){

    options(stringsAsFactors = FALSE)
    commondata <- load("../output/02.DEG/common_data.RData")
    ##
    geneset.remove.dup <- read.csv(paste0("../output/15.pathway/",setname,"/all.genes.counts.csv"),header = TRUE)
    geneset6 <- list(geneset.remove.dup$X)
    names(geneset6) <- setname
    ##
    pdf(paste0("../output/15.pathway/",setname,"/gsea.pathway.pdf"),width = 8,height = 2.5)
    gsea.score <-  matrix(rep(NA,9), nrow = 1,
                          dimnames = list("",c("pathway","pval","padj","log2err","ES","NES",
                                               "size","leadingEdge","Group")),
    )

    for (i in 1:ncol(gene.fc0)) {
        tempgl <- gene.fc0[,i]
        names(tempgl) <- rownames(gene.fc0)
        tempgl <- tempgl[order(tempgl,decreasing = FALSE)]
        fgseaRes <- fgsea(pathways = geneset6,
                          stats = tempgl,
                          minSize=5,
                          maxSize=3000)
        temp.score <- as.data.frame(fgseaRes[1])
        temp.score$Group <-colnames(gene.fc0)[i]
        gsea.score <- rbind(gsea.score,temp.score)

        p1 <- plotEnrichment(geneset6[[setname]],
                             tempgl) + labs(title=colnames(gene.fc0)[i])
        plot(p1)
    }
    dev.off()
}


pathway <- list(c("GO:0048870","GO","cell_motility"),
                c("GO:0050839","GO",	"GO_cell_adhesion_molecule_binding"),
                c("GO:0007229","GO","integrin-mediated_signaling_pathway"))

for (i in 1:length(pathway)) {
    kegg_GO_pathway.gsea(setname = pathway[[i]][3])
}
