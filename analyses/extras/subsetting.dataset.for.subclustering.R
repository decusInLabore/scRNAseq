###############################################################################
# module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;
###############################################################################

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


# Select only cluster 1
OsC_sel1 <- subset(x = OsC, subset = seurat_clusters %in% c(5,8))


sampleIDs <- unique(OsC_sel1$sampleID)

## Remove skin sample ##
#sampleIDs <- sampleIDs[sampleIDs != "Skin01sub"]

for (i in 1:length(sampleIDs)){
    FNout <- paste0(
      "/camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/496_scRNAseq_myc_timecourse/basedata/"
      ,
      "input_",
      sampleIDs[i],
      "_C5C8.txt"
    )
    
    OsC_temp <- subset(x= OsC_sel1, subset = sampleID == sampleIDs[i])
    dfMatrix <- OsC_temp[["RNA"]]@counts
    dfMatrix <- dfMatrix[!(Matrix::rowSums(dfMatrix) == 0),]
    print(paste0(sampleIDs[i], ": ", ncol(dfMatrix)))
    if (nrow(dfMatrix) > 10){
        write.table(dfMatrix, FNout, sep = "\t")
    }
}


save(
  OsC, 
  file = FN
)



