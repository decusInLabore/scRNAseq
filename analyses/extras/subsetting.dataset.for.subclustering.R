
## Run this module in /scripts/scRNAseq/analyses/Main_Analysis

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

baseDir <- "/camp/stp/babs/working/boeings/Projects/sahaie/rebecca.lee/506_mouse_liver_single_cell_SC20154_plus_SC22147/"
FN <- paste0(baseDir, "workdir/SC22147.Seurat.Robj")

load(FN)

## Data selection for subsetting:
#sampleID: primaryTumor

# 3.	To subcluster Neutrophil (19), Hepatocyte (8) and B_cell (5) clusters in samples 1-5 

## Review table
unique(OsC@meta.data[,c("clusterName", "seurat_clusters")])


# Select only cluster 1
selClusters <- c(5,8, 19)
OsC_sel1 <- subset(x = OsC, subset = seurat_clusters %in% selClusters)


sampleIDs <- unique(OsC_sel1$sampleID)

## Remove skin sample ##
#sampleIDs <- sampleIDs[sampleIDs != "Skin01sub"]

outDataDir <- paste0(baseDir, "basedata/")

if (!dir.exists(outDataDir)){
    dir.create(outDataDir)
}


sampleVec <- as.vector(NULL, mode = "character")
pathVec   <- as.vector(NULL, mode = "character")

for (i in 1:length(sampleIDs)){
    FNout <- paste0(
      outDataDir,
      "input_",
      sampleIDs[i],
      "_C",
      paste0(selClusters, collapse = "C"),
      ".txt"
    )
    
    
    
    OsC_temp <- subset(x= OsC_sel1, subset = sampleID == sampleIDs[i])
    dfMatrix <- OsC_temp[["RNA"]]@counts
    dfMatrix <- dfMatrix[!(Matrix::rowSums(dfMatrix) == 0),]
    print(paste0(sampleIDs[i], ": ", ncol(dfMatrix)))
    if (ncol(dfMatrix) > 10){
        write.table(dfMatrix, FNout, sep = "\t")
        print(paste0(sampleIDs[i], " written to file."))
        
        pathVec <- c(
            pathVec,
            FNout
        )
        
        sampleVec <- c(
            sampleVec,
            sampleIDs[i]
        )
    }
    
    
    
}


# save(
#   OsC, 
#   file = FN
# )

## Create text template with file paths
dfTemplate <- data.frame(t(data.frame(sample=sampleVec, path=pathVec)))

FnT <- "../../../scRNAseq/design/pathTemplate.txt"
write.table(
    dfTemplate,
    FnT, 
    sep ="\t"
)



