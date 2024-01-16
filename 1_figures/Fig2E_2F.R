# ---------------------------------------------
# Date: 01-16-2024
# Written by: Jingzhe Shang
# Summary: Fig2E_2F
# ---------------------------------------------

# Load in libraries
suppressMessages({
    library(ggplot2)
    library(reshape2)
    library(ggpubr)
    library(cowplot)
    library(Rmisc)
    library(CePa)
    library(biomaRt)
    library(ggrepel)
    library(Seurat)
    library(RColorBrewer)
    })

cluster_marker <- function(GSEdata) {
    if (!dir.exists(paste0("../output/",GSEdata,"/02.cluster"))) dir.create(paste0("../output/",GSEdata,"/02.cluster"))

    ## load data
    rawdata <- load(paste0("../output/",GSEdata,"/01.QC/",GSEdata,".QC.Rdata"))
    DefaultAssay(seuratdata) <- "RNA"

    ## specify that we will perform downstream analysis on the corrected data note that the
    ## original unmodified data still resides in the 'RNA' assay

    ## Run the standard workflow for visualization and clustering
    seuratdata <- FindVariableFeatures(seuratdata, nfeatures = 3000)

    ## Identify the 10 most highly variable genes
    top10_frist <- head(VariableFeatures(seuratdata), 10)

    ## plot variable features with and without labels
    plot1 <- VariableFeaturePlot(seuratdata)
    plot2 <- LabelPoints(plot = plot1, points = top10_frist, repel = TRUE)
    p <- plot1 + plot2
    p
    ggsave(paste0("../output/",GSEdata,"/02.cluster/0.top10.VF.pdf"),plot=p,width = 20,height =5,
           scale = 0.5)

    seuratdata <- ScaleData(seuratdata, verbose = FALSE)
    seuratdata <- RunPCA(seuratdata, verbose = FALSE)
    seuratdata <- RunUMAP(seuratdata, reduction = "pca", dims = 1:30)
    seuratdata <- FindNeighbors(seuratdata, reduction = "pca", dims = 1:30)
    seuratdata <- FindClusters(seuratdata, resolution = 0.5)
    if ("Cluster" %in% colnames(seuratdata@meta.data) ) {
        Idents(seuratdata) <- seuratdata$Cluster
    } else {
        Idents(seuratdata) <- seuratdata$seurat_clusters
    }

    # Visualization
    p1 <- DimPlot(seuratdata, reduction = "umap", group.by = "Tissue")
    p2 <- DimPlot(seuratdata, reduction = "umap", label = TRUE, repel = TRUE)
    p <- p1 + p2
    p
    ggsave(paste0("../output/",GSEdata,"/02.cluster/1.Visualization.cluster.pdf"),plot=p,width = 13,height =6)

    p3 <- DimPlot(seuratdata, reduction = "umap", split.by = "Tissue",ncol = 2)
    p3
    ggsave(paste0("../output/",GSEdata,"/02.cluster/2.Visualization.cluster.by.ident.pdf"),plot=p3,width = 6,height =6)
    save(seuratdata,file=paste0("../output/",GSEdata,"/01.QC/",GSEdata,".cluster.Rdata") )



}



celldeg <- function(GSEdata,cells="CD8",
                    genelist.dir="../genelist/genelist.txt",
                    compares) {

    theme_set(theme_cowplot())
    if (!dir.exists(paste0("../output/",GSEdata,"/06.deg"))) dir.create(paste0("../output/",GSEdata,"/06.deg"))
    cel1ldata=load(paste0("../output/",GSEdata,"/01.QC/",GSEdata,".celltype.check.Rdata") )
    ## Identify differential expressed genes across conditions
    gl.name <- gsub("\\.txt", "",basename(genelist.dir))
    genelist <- read.delim2(genelist.dir,sep = "\t",header = FALSE)
    genelist <- genelist$V1

    if (!dir.exists(paste0("../output/",GSEdata,"/06.deg/",cell))) dir.create(paste0("../output/",GSEdata,"/06.deg/",cell))
        seuratdata$singleR.cell <- Idents(seuratdata)
        ## CD8.cells <- seuratdata
        CD8.cells <- subset(seuratdata, idents = cell)
        Idents(CD8.cells) <- CD8.cells$Tissue
        ## subcell cluster
        CD8.cells <- RunPCA(CD8.cells, npcs = 30, verbose = FALSE)
        CD8.cells <- RunUMAP(CD8.cells, reduction = "pca", dims = 1:30)
        CD8.cells$celltype.stim <- paste(CD8.cells$singleR.cell,Idents(CD8.cells), sep = "_")

        avg.CD8.cells <- as.data.frame(log1p(AverageExpression(CD8.cells, verbose = FALSE)$RNA))
        avg.CD8.cells$gene <- rownames(avg.CD8.cells)

        avg.CD8.cells[avg.CD8.cells=="NaN"] <- 0
        avg.CD8.cells[avg.CD8.cells=="Inf"] <- 0

        for (i in compares) {
            Idents(CD8.cells) <- CD8.cells$Tissue
            comname <- unlist(strsplit(i," VS "))
            out_name <- paste0(comname,collapse = "_")
            comname1 <- paste(cell,comname,sep = "_")

            Idents(CD8.cells) <- "celltype.stim"
            #保存差异表格
            cd8.deg <- FindMarkers(CD8.cells, ident.1 = comname1[1], ident.2 = comname1[2],
                                   verbose = FALSE,min.pct = 0.25,min.diff.pct = -Inf,
                                   min.cells.feature =  3,min.cells.group = 3)

            cd8.deg$deg <- ifelse(cd8.deg$avg_log2FC>0,"UP", "DOWN")
            cd8.deg$deg[abs(cd8.deg$avg_log2FC)<1] <- "NA"
            cd8.deg$deg[cd8.deg$p_val_adj>0.05] <- "NA"
            cd8.deg$SYMBOL <- rownames(cd8.deg)

            cd8.deg.fc1 <- cd8.deg[abs(cd8.deg$avg_log2FC)>0.59,]
            cd8.deg.fc1 <- cd8.deg.fc1[cd8.deg.fc1$avg_log2FC!="-Inf",]
            cd8.deg.fc1 <- cd8.deg.fc1[cd8.deg.fc1$avg_log2FC!="Inf",]
            cd8.deg.fc1 <- cd8.deg.fc1[cd8.deg.fc1$p_val<0.05,]

            DEG_up <- cd8.deg.fc1[cd8.deg.fc1$avg_log2FC>0,]
            DEG_down <- cd8.deg.fc1[cd8.deg.fc1$avg_log2FC<0,]

            write.csv(DEG_up,paste0("../output/",GSEdata,"/06.deg/",cell,"/",out_name,".deg.up.csv"),
                      row.names = FALSE,quote = FALSE)
            write.csv(DEG_down,paste0("../output/",GSEdata,"/06.deg/",cell,"/",out_name,".deg.down.csv"),
                      row.names = FALSE,quote = FALSE)
            #
            avg.CD8.cells.sub <- avg.CD8.cells[avg.CD8.cells$gene %in% rownames(cd8.deg.fc1),   ]

            index <-genelist %in% rownames(cd8.deg)
            genelist.sub <- genelist[index]
            pdf(paste0("../output/",GSEdata,"/06.deg/",cell,"/",out_name,".4.genelist.DEG.featureplot.pdf"),width = 6,height = 6)
            for (j in 1:length(genelist.sub)) {
                #j <- 22
                Feature.Plot <- FeaturePlot(CD8.cells, features = genelist.sub[j],label = TRUE,# split.by = "Tissue",
                                            #blend=TRUE,
                                            #cols = c("deepskyblue1","white", "red"),
                                            #cols = brewer.pal(3,"Accent"),
                                            max.cutoff = 3
                )
                plot(Feature.Plot)
            }
            dev.off()

            if(length(unique(CD8.cells$Tissue))==3) dotcolor=c("deepskyblue1", "bisque2","pink")
            if(length(unique(CD8.cells$Tissue))==2) dotcolor=c("deepskyblue1", "pink")
            #deg vlnplot
            pdf(paste0("../output/",GSEdata,"/06.deg/",cell,"/",out_name,".5.genelist.DEG.vlnplot.pdf"),width = 6.5,height = 6)
            for (z in 1:length(genelist.sub)) {
                plots <- VlnPlot(CD8.cells, features =genelist.sub[z], split.by = "Tissue",
                                 cols = dotcolor,
                                 group.by = "singleR.cell",
                                 pt.size = 0, combine = FALSE)
                plot(plots[[1]]  )
            }
            dev.off()
          }



    }
## an example
cluster_marker(GSEdata = "GSE108989")

celldeg(GSEdata = "GSE108989",cells = "CD8",
        genelist.dir = "../genelist/genelist.txt",
        compares = c("TTC VS PTC"))
