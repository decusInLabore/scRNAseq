###############################################################################
# module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;
###############################################################################

###############################################################################
## Initialise renv

# if (!requireNamespace("remotes")){
#   install.packages("remotes")
# }

# remotes::install_github("rstudio/renv")

if (!file.exists("renv.lock")){
    renv::init()
} else {
    renv::restore(prompt=FALSE)
}


###############################################################################
## Prepare subset for sub-clustering                                         ##
library(Seurat)
## load Seurat object with basedata
FN <- "/camp/stp/babs/working/boeings/Projects/hillc/danielle.park/450_singleCell_scPDAC_PRJCA001063full/workdir/scPDAC_Wu_all.Seurat.Robj"
load(FN)

## Expected name of the Seurat object: OsC ##

## Data selection for subsetting:
#sampleID: primaryTumor

# selectClusters for subclustering
selClust <- c(4,5)
samplestring <- paste0("C", paste0(selClust , collapse = ""))
subClusterID <- "Fib"

# Select only cluster 1
OsC_sel1 <- subset(x = OsC, subset = seurat_clusters %in% selClust)

sampleIDs <- unique(OsC_sel1$sampleID)

## Remove skin sample ##
#sampleIDs <- sampleIDs[sampleIDs != "Skin01sub"]

## make samplestring


for (i in 1:length(sampleIDs)){
    FNout <- paste0(
      "/camp/stp/babs/working/boeings/Projects/hillc/danielle.park/450_singleCell_scPDAC_PRJCA001063full/basedata/",
      "input_",
      sampleIDs[i],
      "_",
      samplestring,
      subClusterID,
      ".txt"
    )
    
    OsC_temp <- subset(x= OsC_sel1, subset = sampleID == sampleIDs[i])
    dfMatrix <- OsC_temp[["RNA"]]@counts
    dfMatrix <- dfMatrix[!(Matrix::rowSums(dfMatrix) == 0),]
    print(paste0(sampleIDs[i], ": ", ncol(dfMatrix)))
    if (nrow(dfMatrix) > 10){
        write.table(dfMatrix, FNout, sep = "\t")
    }
}

OsC@meta.data[[paste0("meta_Subclustering_", subClusterID)]] <- "Rest"
OsC@meta.data[OsC@meta.data$seurat_clusters %in% selClust, paste0("meta_Subclustering_",subClusterID)] <- "Selected"

save(
  OsC, 
  file = FN
)


#renv::snapshot()
