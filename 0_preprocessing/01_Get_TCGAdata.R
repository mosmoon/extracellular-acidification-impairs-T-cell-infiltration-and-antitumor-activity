# ---------------------------------------------
# Date: 01-16-2024
# Written by: Jingzhe Shang
# Summary: Get Tumor rawdata
# ---------------------------------------------

# Load in libraries
suppressMessages({
    library(recount3)
    library(plyr)
})

# get Tumor data function
recounts_TCGA <- function(Cancer_Type) {
    #recount3_cache_files()
    ## Find all available human projects
    human_projects <- available_projects()
    ## Find the project you are interested in,
    Cancer_tcga <- subset(
        human_projects,
        project == Cancer_Type & project_type == "data_sources"

    )
    ## Create a RangedSummarizedExperiment (RSE) object at the gene level
    rse_gene_Cancer_tcga <- create_rse(project_info = Cancer_tcga)
    assay(rse_gene_Cancer_tcga, "counts") <- transform_counts(rse_gene_Cancer_tcga)
    assay(rse_gene_Cancer_tcga, "TPM") <- recount::getTPM(rse_gene_Cancer_tcga,length_var = "score")

    ## Sample metadata
    Cancer_tcga_cols <- rownames(colData(rse_gene_Cancer_tcga))
    Cancer_tcga_colDatas <- colData(rse_gene_Cancer_tcga)@listData
    index <- c("external_id","tcga.tcga_barcode","tcga.cgc_case_clinical_stage")
    Cancer_tcga_sample <- as.data.frame(Cancer_tcga_colDatas[index])

    ## Sample grouping（tumor，normal）
    Cancer_tcga_sample$sample <- substr(Cancer_tcga_sample$tcga.tcga_barcode,14,15)

    ## Filter paraffin embedding
    Cancer_tcga_sample$vial <- substr(Cancer_tcga_sample$tcga.tcga_barcode,16,16)
    index <- Cancer_tcga_sample$vial %in% "B"
    Cancer_tcga_sample <- Cancer_tcga_sample[!index,]

    ## Filter non-RNA data
    Cancer_tcga_sample$Analyte <- substr(Cancer_tcga_sample$tcga.tcga_barcode,20,20)
    index <- Cancer_tcga_sample$Analyte %in% c("R","H","T")
    Cancer_tcga_sample <- Cancer_tcga_sample[index,]
    Cancer_tcga_sample$sample_T <- ifelse(Cancer_tcga_sample$sample %in% 1:9,"Tumor", "Normal")  #ifelse实现分组

    index <-   rse_gene_Cancer_tcga$tcga.tcga_barcode %in% Cancer_tcga_sample$tcga.tcga_barcode
    rse_gene_Cancer_tcga <- rse_gene_Cancer_tcga[,index]

    save(rse_gene_Cancer_tcga,Cancer_tcga_sample,
         file = paste0("../output/00.RawCount/",Cancer_Type,".rawdata.RData"))
}

## get Tumor data

human_projects <- available_projects(organism = "human")
human_projects_tcga <- human_projects[human_projects$file_source %in% "tcga", ]
lapply(human_projects_tcga$project,recounts_TCGA)








