###############################################################################
## Initialise renv

if (!requireNamespace("remotes")){
  install.packages("remotes")
}

remotes::install_github("rstudio/renv")

if (!file.exists("renv.lock")){
    renv::init()
} else {
    renv::restore()
}


###############################################################################
## Prepare subset for sub-clustering                                         ##
library(Seurat)
## load Seurat object with basedata
FN <- "/camp/stp/babs/working/boeings/Projects/lovellbadger/emily.frost/358B_RLL_EF_single_cell_p27kip1_in_ovarian_granulosa_cells_SC19235/workdir/rll358B.Seurat.Robj"
load(FN)

## Data selection for subsetting:
sampleID: primaryTumor



OsC_sel1 <- subset(x = OsC, subset = seurat_clusters %in% c(1,2,3,8,14) )

sampleIDs <- unique(OsC_sel1$sampleID)

for (i in 1:length(sampleIDs)){
    FNout <- paste0(
      "/camp/stp/babs/working/boeings/Projects/lovellbadger/emily.frost/358B_RLL_EF_single_cell_p27kip1_in_ovarian_granulosa_cells_SC19235/basedata/",
      "input_",
      sampleIDs[i],
      "_C123814.txt"
    )
    
    OsC_temp <- subset(x= OsC_sel1, subset = sampleID == sampleIDs[i])
    dfMatrix <- OsC_temp[["RNA"]]@counts
    dfMatrix <- dfMatrix[!(Matrix::rowSums(dfMatrix) == 0),]
    print(paste0(sampleIDs[i], ": ", ncol(dfMatrix)))
    if (nrow(dfMatrix) > 10){
        write.table(dfMatrix, FNout, sep = "\t")
    }
}

OsC@meta.data[["meta_Subclustering_rll358B"]] <- "Rest"
OsC@meta.data[OsC@meta.data$seurat_clusters %in% c(1,2,3,8,14), "meta_Subclustering_rll358B"] <- "Selected"

save(
  OsC, 
  file = FN
)


renv::snapshot()
