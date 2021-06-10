###############################################################################
## Prepare subset for sub-clustering                                         ##
library(Seurat)
## load Seurat object with basedata
FN <- "/camp/stp/babs/working/boeings/Projects/lovellbadger/emily.frost/358B_RLL_EF_single_cell_p27kip1_in_ovarian_granulosa_cells_SC19235/workdir/rll358B.Seurat.Robj"
load(FN)

## Data selection for subsetting:
sampleID: primaryTumor



OsC_sel1 <- subset(x = OsC, subset = seurat_clusters %in% c(1,2,7,9,15) )

sampleIDs <- unique(OsC_sel1$sampleID)

for (i in 1:length(sampleIDs)){
    FNout <- paste0(
      "/camp/stp/babs/working/boeings/Projects/lovellbadger/emily.frost/358B_RLL_EF_single_cell_p27kip1_in_ovarian_granulosa_cells_SC19235/basedata/",
      "input_",
      sampleIDs[i],
      "_C127915.txt"
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
OsC@meta.data[OsC@meta.data$seurat_clusters %in% c(1,2,7,9,15), "meta_Subclustering_rll358B"] <- "Selected"

save(
  OsC, 
  file = FN
)



library(Seurat)
library(ggplot2)
library(tidyverse)
library(knitr)

VersionPdfExt <- paste0(".V", gsub("-", "", Sys.Date()), ".pdf")

if (dir.exists("/Volumes/babs/working/boeings/")){
  hpc.mount <- "/Volumes/babs/working/boeings/"
} else if (dir.exists("Y:/working/boeings/")){
  hpc.mount <- "Y:/working/boeings/"
} else if (dir.exists("/camp/stp/babs/working/boeings/")){
  hpc.mount <- "/camp/stp/babs/working/boeings/"
} else {
  hpc.mount <- ""
}


FN <- paste0(hpc.mount, "Projects/reference_data/documentation/BC.parameters.txt")
dbTable <- read.delim(
  FN, 
  sep = "\t",
  stringsAsFactors = F
)

db.pwd <- as.vector(dbTable[1,1])

figureCount <- 1

source("assets/R/SBwebtools.pckg.r")

if (length(.libPaths()) > 2){
  .libPaths(.libPaths()[2:3])
}
## Create biologic Object for visualization ##

ObioFN <- paste0("../", list.files("..")[grep(".bioLOGIC.Robj", list.files(".."))])

load(ObioFN)

checkFile = paste0(
  Obio@parameterList$project_id,
  ".bioLOGIC.Robj"
)


Obio <- setMountingPoint(Obio)
Obio <- setAnalysisPaths(Obio)
Obio <- setCrickGenomeAndGeneNameTable(Obio)
Obio <- createAnalysisFolders(
  Obio
)
Obio <- setDataBaseParameters(Obio)
## Upload new metadata table ##
#```{r child = 'src/modules/db_tools/upload.meta.data.table.to.DB.Rmd', eval=TRUE}
#
#```