###############################################################################
## Prepare subset for sub-clustering                                         ##

## load Seurat object with basedata
FN <- "/camp/stp/babs/working/boeings/Projects/hillc/danielle.park/437_fibroblast_GSE154778/workdir/scPDAC.Seurat.Robj"
load(FN)

## Data selection for subsetting:
sampleID: primaryTumor


library(Seurat)
OsC_sel1 <- subset(x = OsC, subset = meta_Region == "Primary_Tumor" & clusterName == "C4")

sampleIDs <- unique(OsC_sel1$sampleID)

for (i in 1:length(sampleIDs)){
    FNout <- paste0(
      "/camp/stp/babs/working/boeings/Projects/hillc/danielle.park/437subC_fibroblast_GSE154778_primTumor_C4/basedata/",
      "input_",
      sampleIDs[i],
      "_C4.txt"
    )
    
    OsC_temp <- subset(x= OsC_sel1, subset = sampleID == sampleIDs[i])
    dfMatrix <- OsC_temp[["RNA"]]@counts
    dfMatrix <- dfMatrix[!(Matrix::rowSums(dfMatrix) == 0),]
    print(paste0(sampleIDs[i], ": ", ncol(dfMatrix)))
    if (nrow(dfMatrix) > 10){
        write.table(dfMatrix, FNout, sep = "\t")
    }
}

