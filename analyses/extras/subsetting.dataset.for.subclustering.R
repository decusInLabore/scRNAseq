
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

# To subcluster melanoma cells (current clusters Melanoma_1 and renamed Melanoma_2) 
# in samples 2,4, and 5

# BackgroundLiverM3	  1
# LiverMetM3	        2
# LiverNormalM1	      3
# SubcutM3	          4
# SubcutNoLiverMetM2	5

# and to compare proportions of different sub clusters. Deliverables: Subclustering of Melanoma_1 and Melanoma_2 subclusters in samples 2, 3, and 5 in analysis SC22147, summarized in a standard single-cell analysis report. 

## Review table
unique(OsC@meta.data[,c("clusterName", "seurat_clusters")])

# Melanoma               3
#  Melanoma_2           13

# Select only cluster 1
selClusters <- c(3, 13)
OsC_sel1 <- subset(x = OsC, subset = seurat_clusters %in% selClusters)

# select samples 2, 4 and 5
selSamples <- c(
  "LiverMetM3R1", "LiverMetM3R2",
  "SubcutM3R1", "SubcutM3R2",
  "SubcutNoLiverMetM2R1", "SubcutNoLiverMetM2R2"
)

OsC_sel1 <- subset(x = OsC_sel1, subset = sampleName %in% selSamples)

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


# LiverMetM3R2 features only 10 cells in this selection - excluded. 

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



