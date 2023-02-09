###############################################################################
################################################################################

## The subset analysis needs to run from the main analysis module, to make use of the 
# renv lock file there.


###############################################################################
## Initialise renv

# if (!requireNamespace("remotes")){
#   install.packages("remotes")
# }
# 
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
FN <- "/camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/496_scRNAseq_myc_timecourse/workdir/MYC.Seurat.Robj"

load(FN)

## Data selection for subsetting:
#sampleID: primaryTumor

selectedClusters <- c(12, 14, 17)

# Select only cluster 1
OsC_sel1 <- subset(x = OsC, subset = seurat_clusters %in% selectedClusters)


allSampleIDs <- unique(OsC_sel1$sampleID)

## aggregate by sample
## remove everything after the last underscore in sampleIDs
sampleIDs <- unique(gsub('(.*)_\\w+', '\\1', allSampleIDs))

## Remove skin sample ##
#sampleIDs <- sampleIDs[sampleIDs != "Skin01sub"]

smallSubsets <- as.vector(NULL, mode = "character" )

for (i in 1:length(sampleIDs)){
    FNout <- paste0(
      "/camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/496_scRNAseq_myc_timecourse/basedata/"
      ,
      "input_",
      sampleIDs[i],
      "_496V3C",
      paste0(selectedClusters, collapse = "C"),
      ".txt"
    )
    
    sampleGroup <- allSampleIDs[grep(paste0(sampleIDs[i], "_[A-Za-z]"), allSampleIDs)]
    print(paste0("Sample in group ", sampleIDs[i], ":", paste0(sampleGroup, collapse = ", ")))
    
    OsC_temp <- subset(x= OsC_sel1, subset = sampleID %in% sampleGroup)
    dfMatrix <- OsC_temp[["RNA"]]@counts
    dfMatrix <- dfMatrix[!(Matrix::rowSums(dfMatrix) == 0),]
    print(paste0(sampleIDs[i], ": ", ncol(dfMatrix)))
    if (ncol(dfMatrix) > 200){
        write.table(dfMatrix, FNout, sep = "\t")
        print("Processed.")
    } else {
        smallSubsets <- c(
            smallSubsets, 
            sampleIDs[i]
        )
        print(smallSubsets)
    }
}

if (length(smallSubsets) > 0){
  FNout <- paste0(
    "/camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/496_scRNAseq_myc_timecourse/basedata/",
    "input_combined_small_subsets_",
    paste0(smallSubsets, collapse = "_"),
    "_496V3C",
    paste0(selectedClusters, collapse = "C"),
    ".txt"
  )
  
  OsC_temp <- subset(x= OsC_sel1, subset = sampleID %in% smallSubsets)
  dfMatrix <- OsC_temp[["RNA"]]@counts
  dfMatrix <- dfMatrix[!(Matrix::rowSums(dfMatrix) == 0),]
  print(paste0(sampleIDs[i], ": ", ncol(dfMatrix)))
  if (nrow(dfMatrix) > 200){
    write.table(dfMatrix, FNout, sep = "\t")
  }
}

# save(
#   OsC, 
#   file = FN
# )



