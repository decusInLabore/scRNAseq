###############################################################################
# module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;
###############################################################################

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
FN <- "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449E_subclustering_449C_SC20193_C1_2_3_4_5/workdir/pbl449EsubCluster.Seurat.Robj"
load(FN)

## Data selection for subsetting:
#sampleID: primaryTumor


# Select only cluster 1
OsC_sel1 <- subset(x = OsC, subset = seurat_clusters %in% c(1))

sampleIDs <- unique(OsC_sel1$sampleID)

## Remove skin sample ##
#sampleIDs <- sampleIDs[sampleIDs != "Skin01sub"]

for (i in 1:length(sampleIDs)){
    FNout <- paste0(
      "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449E_subclustering_449C_SC20193_C1_2_3_4_5/basedata/",
      "input_",
      sampleIDs[i],
      "_C1only.txt"
    )
    
    OsC_temp <- subset(x= OsC_sel1, subset = sampleID == sampleIDs[i])
    dfMatrix <- OsC_temp[["RNA"]]@counts
    dfMatrix <- dfMatrix[!(Matrix::rowSums(dfMatrix) == 0),]
    print(paste0(sampleIDs[i], ": ", ncol(dfMatrix)))
    if (nrow(dfMatrix) > 10){
        write.table(dfMatrix, FNout, sep = "\t")
    }
}

OsC@meta.data[["meta_Subclustering_C1only"]] <- "Rest"
OsC@meta.data[OsC@meta.data$seurat_clusters %in% c(1), "meta_Subclustering_C1only"] <- "Selected"

save(
  OsC, 
  file = FN
)


renv::snapshot()
