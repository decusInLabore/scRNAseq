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
FN <- "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/484_scRNAseq_Characterisation_EpCAM_CD90_thymic_epithelial_cells_SC22096/workdir/SC22096.Seurat.Robj"
load(FN)

## Data selection for subsetting:
#sampleID: primaryTumor

## Assume that all unassigned cells are human
OsC@meta.data[OsC@meta.data$meta_inferred_species == 0, "meta_inferred_species"] <- "human"

# Select only cluster 1
OsC_sel1 <- subset(x = OsC, subset = meta_inferred_species %in% c("human"))
OsC_sel1@meta.data[["sampleID"]] <- paste0(OsC_sel1@meta.data$sampleID, "_", OsC_sel1@meta.data$meta_Hash)
OsC_sel1@meta.data$sampleID <- gsub("_0", "", OsC_sel1@meta.data$sampleID)
OsC_sel1@meta.data$sampleID <- gsub("EpCAMpos", "", OsC_sel1@meta.data$sampleID)
OsC_sel1@meta.data$sampleID <- gsub("EpCAMneg", "", OsC_sel1@meta.data$sampleID)



sampleIDs <- unique(OsC_sel1$sampleID)

## Remove skin sample ##
#sampleIDs <- sampleIDs[sampleIDs != "Skin01sub"]

for (i in 1:length(sampleIDs)){
    FNout <- paste0(
      "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/484_scRNAseq_Characterisation_EpCAM_CD90_thymic_epithelial_cells_SC22096/basedata/"
      ,
      "input_",
      sampleIDs[i],
      "_hs.txt"
    )
    
    OsC_temp <- subset(x= OsC_sel1, subset = sampleID == sampleIDs[i])
    dfMatrix <- OsC_temp[["RNA"]]@counts
    dfMatrix <- dfMatrix[!(Matrix::rowSums(dfMatrix) == 0),]
    print(paste0(sampleIDs[i], ": ", ncol(dfMatrix)))
    if (nrow(dfMatrix) > 10){
        write.table(dfMatrix, FNout, sep = "\t")
    }
}

OsC@meta.data[["meta_Subclustering_T_cells"]] <- "Rest"
OsC@meta.data[OsC@meta.data$seurat_clusters %in% c(1), "meta_Subclustering_T_cells"] <- "Selected"

save(
  OsC, 
  file = FN
)


renv::snapshot()
