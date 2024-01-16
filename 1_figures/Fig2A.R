# ---------------------------------------------
# Date: 01-16-2024
# Written by: Jingzhe Shang
# Summary: Fig2A
# ---------------------------------------------

# Load in libraries
suppressMessages({
    library(ggvenn)
    library(Vennerable)
    library(gridGraphics)
})


venn.analysis <- function(compares,FC,vennsub){

    normCounts <- read.csv("../output/02.DEG/2.normCounts.csv",header = T)
    colnames(normCounts) <- gsub("\\.","-",colnames(normCounts))
    mrnacounts <- normCounts[normCounts$gene_biotype=="protein_coding",]
    if (!dir.exists("../output/04.mRNA_Venn")) dir.create("../output/04.mRNA_Venn")

    for (i in 1:length(compares)) {

        assign(compares[i],
               read.csv(paste0("../output/02.DEG/",gsub(" ","_",compares[i]),"/",gsub(" ","_",compares[i]),"_DEG.csv"),header = T))
    }

    list.deg <- list(get(compares[1])$GeneName)

    for (i in 2:length(compares)) {
        compare.name=get(compares[i])$GeneName
        list.deg <- c(list.deg,list(compare.name))
    }
    names(list.deg) <- compares
    venn.list <- list.deg
    ## up and down
    get.FC.deg.list <- function(FC,compares) {
        if (FC > 0) {
            index <- get(compares[1])$log2fdc > FC
        }
        if (FC < 0) {
            index <- get(compares[1])$log2fdc< FC
        }
        list.deg.reg <- list(get(compares[1])$GeneName[index])
        for (i in 2:length(compares)) {
            if (FC > 0) {
                index <- get(compares[i])$log2fdc > FC
            }
            if (FC < 0) {
                index <- get(compares[i])$log2fdc< FC
            }
            compare.name=get(compares[i])$GeneName[index]
            list.deg.reg <- c(list.deg.reg,list(compare.name))
        }
        names(list.deg.reg) <- compares
        list.deg.reg
    }

    list.deg.down <- get.FC.deg.list(-FC,compares)
    list.deg.up <- get.FC.deg.list(FC,compares)

    ## plot
    plot.venn <- function(venn.list,reg,nc) {

        p1 <- ggvenn(venn.list,names(venn.list)[nc],text_size=6,set_name_size = 2,
                     stroke_color = NA,show_percentage = FALSE)
        p1 <- p1+labs(title=reg)
        plot(p1)
        Vstem <- Venn(venn.list[nc])
        venn_partitions <- Vstem@IntersectionSets[-1]

        r <- max(unlist(lapply(venn_partitions,length)))
        c <- length(venn_partitions)
        venn.matrix <- matrix("", nrow =r , ncol = c)
        for (i in 1:length(venn_partitions)) {
            gene <- unlist(venn_partitions[i])
            l <- length(gene)
            if (l>0) venn.matrix[1:l,i] <- gene
        }
        colnames(venn.matrix) <- names(venn_partitions)
        write.csv(venn.matrix,row.names = FALSE,
                  paste0("../output/04.mRNA_Venn/",paste(nc,collapse = ""),reg,"matrix.csv"))

        venn.name <- lapply(names(venn_partitions),function(x) unlist(strsplit(x,split = "")))
        venn.name <- do.call(rbind,venn.name)
        colnames(venn.name) <- names(venn.list)[nc]
        venn_partitions <- cbind(venn.name,unlist(lapply(venn_partitions,function(x) {
            paste(x,collapse = ", ")
        })))
        write.csv(venn_partitions,row.names = FALSE,
                  paste0("../output/04.mRNA_Venn/",paste(nc,collapse = ""),reg,".csv"))
    }

    plot.sub <- function(vennsub){
        sub.n <- paste0(vennsub,collapse = "")
        pdf(paste0("../output/04.mRNA_Venn/Venn",sub.n,".pdf"),width = 6,height = 4)
        plot.venn(venn.list,"",vennsub)
        plot.venn(list.deg.up,"up",vennsub)
        plot.venn(list.deg.down,"down",vennsub)
        dev.off()
    }
    lapply(vennsub, plot.sub)

}

compares <- c("5LA VS CTRL","10LA VS CTRL","6HCL VS CTRL","8HCL VS CTRL")
venn.analysis(compares = compares,FC = 1,vennsub = list(c(2,3),1:4) )




