# ---------------------------------------------
# Date: 01-16-2024
# Written by: Jingzhe Shang
# Summary: FigS3_2A
# ---------------------------------------------

# Load in libraries
suppressMessages({
    library(pheatmap)
    library(reshape2)
    library(ggplot2)
    library(RColorBrewer)
})


pathway <- list(c("GO:0050839","GO","GO_cell_adhesion_molecule_binding"),
                c("GO:0006935","GO","chemotaxis"),
                c("GO:2000249","GO","GO_regulation_of_actin_cytoskeleton"))
sub.fc = list(c(2:3))

for (i in 1:length(pathway)) {
    pathway.pheatmap(setname = concern_pathway[[i]][3],sub.fc = sub.fc )
}


pathway.pheatmap <- function(setname,sub.fc) {


    common.data <- load("../output/02.DEG/common_data.RData")
    #i <- 1
    for (i in 1:length(sub.fc)) {
        sub.name <- paste0(sub.fc[[i]],collapse = "")

        if (!dir.exists(paste0("../output/15.pathway/",setname,"/",sub.name))) dir.create(paste0("../output/15.pathway/",setname,"/",sub.name))
        tempname <- unlist(strsplit(colnames(deg.fc0)[sub.fc[[i]]],".VS."))
        tempname <- gsub("\\.","-",tempname)
        index <- gsub("_G\\d+$","",colnames(normcounts)) %in% tempname
        index[1:5] <- TRUE
        counts.sub <- normcounts[,index]

        deg.fc.sub <- deg.fc0[,c(sub.fc[[i]])]
        index <- apply(as.matrix(deg.fc.sub),1, function(x){
            all(x==0)
        })
        deg.fc.sub <- deg.fc.sub[!index,]
        deg.counts <- counts.sub[rownames(counts.sub) %in% rownames(deg.fc.sub),]
        deg.counts0 <- deg.counts[,-(1:5)]

        deg.counts0.log <- log2(deg.counts0+1)
        deg.counts0.log.sc <- t(scale(t(deg.counts0.log),center = TRUE,scale = FALSE))

        ###################
        geneset.remove.dup <- read.csv(paste0("../output/15.pathway/",setname,"/all.genes.counts.csv"),header = TRUE)
        rownames(geneset.remove.dup) <-geneset.remove.dup$X
        geneset.remove.dup <- geneset.remove.dup[,-1]

        ##deg.fc.sub
        index <- rownames(deg.fc.sub) %in% rownames(geneset.remove.dup)
        category.deg.fc.sub <- deg.fc.sub[index,]

        ##DEGgenes
        index <-rownames(counts.sub) %in% rownames(category.deg.fc.sub)
        category.total.deg.sub <- counts.sub[index,-c(1:5)]
        write.csv(category.total.deg.sub,paste0("../output/15.pathway/",setname,"/",sub.name,"/pathway.heatmap.csv"),
                  quote = TRUE,row.names = TRUE)
        category.total.deg.sub.log <- log2(category.total.deg.sub+1)
        category.total.deg.sub.log0 <- t(scale(t(category.total.deg.sub.log),center = TRUE,scale = FALSE))

        if (nrow(category.total.deg.sub.log)>1) {
            bk <- max(abs(category.total.deg.sub.log0))
            pdf(paste0("../output/15.pathway/",setname,"/",sub.name,"/pathway.heatmap.pdf"),width = 8,height = 8)
            pheatmap(category.total.deg.sub.log,cluster_cols = FALSE,cellwidth = 8,cellheight = 9)
            pheatmap(category.total.deg.sub.log0,cluster_cols = FALSE,cellwidth = 8,cellheight = 8,border_color = NA,
                     breaks = seq(-bk,bk,length=100),main = setname,
                     color = colorRampPalette(c("blue", "white", "red"))(100))

            dev.off()
            pdf(paste0("../output/15.pathway/",setname,"/",sub.name,"/pathway.heatmap.small.pdf"),width = 8,height = 4)

            pheatmap(category.total.deg.sub.log0,cluster_cols = FALSE,cellwidth = 8,border_color = NA,
                     breaks = seq(-bk,bk,length=100),show_rownames = FALSE,main = setname,
                     color = colorRampPalette(c("blue", "white", "red"))(100))
            dev.off()


            index <- rownames(deg.fc0) %in%  rownames(category.total.deg.sub.log)
            genes.fc <- deg.fc0[index,]

            bk <- max(abs(genes.fc))
            genes.fc[abs(genes.fc)<1] <- 0
            pdf(paste0("../output/15.pathway/",setname,"/",sub.name,"/pathway.heatmap.fc.pdf"),width = 8,height = 8)
            pheatmap(genes.fc,cluster_cols = FALSE,cellwidth = 10,cellheight = 10,
                     breaks = seq(-bk,bk,length=100),main = setname,
                     color = colorRampPalette(c("blue", "white", "red"))(100))
            dev.off()

        }
    }
}


