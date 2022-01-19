###############################################################################
## Initialise renv

if (!requireNamespace("remotes")){
  install.packages("remotes")
}

remotes::install_github("rstudio/renv")

if (!file.exists("renv.lock")){
    renv::init()
} else {
    renv::restore(prompt=FALSE)
}


###############################################################################
## Prepare subset for sub-clustering                                         ##
library(Seurat)
## load Seurat object with basedata
FN <- "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/471_scRNAseq_newSamples_Epi6_Epi10_SC20193_SC21201/workdir/A4A6A7A8.Seurat.Robj"
load(FN)

## Data selection for subsetting:
sampleID: primaryTumor


# Remove clusters: 4, 8, 11, 13.
OsC_sel1 <- subset(x = OsC, subset = seurat_clusters %in% c(0:3, 5:7, 9:10, 12) )

sampleIDs <- unique(OsC_sel1$sampleID)

for (i in 1:length(sampleIDs)){
    FNout <- paste0(
      "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/471_scRNAseq_newSamples_Epi6_Epi10_SC20193_SC21201/basedata/",
      "input_",
      sampleIDs[i],
      "_C012356791012.txt"
    )
    
    OsC_temp <- subset(x= OsC_sel1, subset = sampleID == sampleIDs[i])
    dfMatrix <- OsC_temp[["RNA"]]@counts
    dfMatrix <- dfMatrix[!(Matrix::rowSums(dfMatrix) == 0),]
    print(paste0(sampleIDs[i], ": ", ncol(dfMatrix)))
    if (nrow(dfMatrix) > 10){
        write.table(dfMatrix, FNout, sep = "\t")
    }
}

OsC@meta.data[["meta_Subclustering_A4A6A7A8"]] <- "Rest"
OsC@meta.data[OsC@meta.data$seurat_clusters %in% c(0:3, 5:7, 9:10, 12), "meta_Subclustering_A4A6A7A8"] <- "Selected"

save(
  OsC, 
  file = FN
)


renv::snapshot()
