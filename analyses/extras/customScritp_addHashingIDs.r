###############################################################################
## Assign human/mouse cells                                                  ##

# The Goal of this script is to add hashing tags to existing single-cells
# in a single-cell object. 



###############################################################################
## Processing of the cellranger count output and assignment to cellIDs in the##
## integrated single-cell experiment                                         ##

# This is the R-version we're using
# module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a


sourceList <- list(
    "SC22096multiA9" = "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/484A_scRNAseq_Characterisation_EpCAM_CD90_thymic_epithelial_cells_SC22096_hash/FASTQ_files/SC22096multiA9/outs/multi/multiplexing_analysis/tag_calls_per_cell.csv",
    "SC22096multiA10" = "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/484A_scRNAseq_Characterisation_EpCAM_CD90_thymic_epithelial_cells_SC22096_hash/FASTQ_files/SC22096multiA10/outs/multi/multiplexing_analysis/tag_calls_per_cell.csv",
    "SC22096multiA11" = "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/484A_scRNAseq_Characterisation_EpCAM_CD90_thymic_epithelial_cells_SC22096_hash/FASTQ_files/SC22096multiA11/outs/multi/multiplexing_analysis/tag_calls_per_cell.csv",
    "SC22096multiA12" = "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/484A_scRNAseq_Characterisation_EpCAM_CD90_thymic_epithelial_cells_SC22096_hash/FASTQ_files/SC22096multiA12/outs/multi/multiplexing_analysis/tag_calls_per_cell.csv"
)


for (i in 1:length(sourceList)){
    df <- read.delim(sourceList[[i]], sep=",", stringsAsFactors = F)
    
    dfTest <- df[grep("\\|", df$feature_call),]
    for (j in 1:nrow(dfTest)){
        feature_callVec <- unlist(strsplit(dfTest[j, "feature_call"], "\\|"))
        num_umisVec <- as.numeric(unlist(strsplit(dfTest[j, "num_umis"], "\\|")))
        pos <- which.max(num_umisVec)
        dfTest[j, "feature_call"] <- feature_callVec[pos]
        dfTest[j, "num_umis"] <- num_umisVec[pos]
        
    }
    
    rmPos <- dfTest$cell_barcode
    df <- df[!(df$cell_barcode %in% rmPos),]
    
    df <- rbind(
        df,
        dfTest
    )
    
    df[["sampleID"]] <- names(sourceList)[i]
    
    if (i==1){
        dfRes <- df
    } else {
        dfRes <- rbind(
            dfRes,
            df
        )
    }
}

dfRes[["limsID"]] <- dfRes$sampleID

##  Done making cellcode list                                                ##
###############################################################################

###############################################################################
## Adding cell codes to the integrated object                                ##

## This is the cellranger aggregate suffix
dfRes[["suffix"]] <- ""

dfRes[dfRes$sampleID == "SC22096multiA9", "suffix"] <- "-1"
dfRes[dfRes$sampleID == "SC22096multiA10", "suffix"] <- "-2"
dfRes[dfRes$sampleID == "SC22096multiA11", "suffix"] <- "-1"
dfRes[dfRes$sampleID == "SC22096multiA12", "suffix"] <- "-2"

dfRes$cellID <- gsub("-1", "", dfRes$cellID)
dfRes$cellID <- paste0(dfRes$cellID, dfRes$suffix)

## Switch lims IDs to merged IDS
dfRes[dfRes$sampleID == "SC22096multiA9", "sampleID"] <- "A1A2MERGED"
dfRes[dfRes$sampleID == "SC22096multiA10", "sampleID"] <- "A1A2MERGED"
dfRes[dfRes$sampleID == "SC22096multiA11", "sampleID"] <- "A3A4MERGED"
dfRes[dfRes$sampleID == "SC22096multiA12", "sampleID"] <- "A3A4MERGED"

## Remove original suffix
dfRes[["cellID"]] <- sapply(dfRes$cell_barcode, function(x) unlist(strsplit(x, "-"))[1])
dfRes[["cellID"]] <- paste0(dfRes$cellID, dfRes$suffix)
dfRes$suffix <- NULL
# Done in accounting for cellranger aggreate
## Done                                                                      ##
###############################################################################

###############################################################################
## Match with IDs in integrated dataset                                      ##

load("../SC22096.Seurat.Robj")

## This suffix is from the integration in Seurat
dfID <- unique(OsC@meta.data[,c("cellID", "sampleID")])
dfID[["suffix"]] <- sapply(dfID$cellID, function(x) unlist(strsplit(x, "_"))[2])

dfID <- unique(dfID[,c("suffix", "sampleID")])


sampleAnnotation <- list(
    "THY135" = "A1A2MERGED",
    "THY136" = "A3A4MERGED"
)

for (i in 1:length(sampleAnnotation)){
    dfRes[dfRes$sampleID == sampleAnnotation[[i]], "sampleID"] <- names(sampleAnnotation)[i]
}



for (i in 1:length(sampleAnnotation)){
    dfRes[dfRes$sampleID == names(sampleAnnotation)[i], "cellID"] <- paste0(dfRes[dfRes$sampleID == names(sampleAnnotation)[i], "cellID"], "_", as.vector(dfID[dfID$sampleID == names(sampleAnnotation)[i], "suffix"]))
}


## Done                                                                      ##
###############################################################################

###############################################################################
## Add Meta columns for CMO tag                                              ##

dfRes[["meta_Hash"]] <- ""
dfRes[dfRes$feature_call == "CMO305", "meta_Hash"] <- "CD49f"
dfRes[dfRes$feature_call == "CMO306", "meta_Hash"] <- "CD90hiEpCAMneg"
dfRes[dfRes$feature_call == "CMO307", "meta_Hash"] <- "CD90negEpCAMpos"
dfRes[dfRes$feature_call == "CMO308", "meta_Hash"] <- "CD49f"
dfRes[dfRes$feature_call == "CMO309", "meta_Hash"] <- "CD90hiEpCAMneg"
dfRes[dfRes$feature_call == "CMO310", "meta_Hash"] <- "CD90negEpCAMpos"

names(dfRes) <- gsub("feature_call", "meta_hash_feature_call", names(dfRes))
names(dfRes) <- gsub("num_umis", "meta_hash_num_umis", names(dfRes))

dfRes$cell_barcode <- NULL
dfRes$num_features <- NULL
dfRes$sampleID <- NULL
dfRes$limsID <- NULL
dfRes$suffix <- NULL
## Done                                                                      ##
###############################################################################

row.names(dfRes) <- dfRes$cellID
dfRes$cellID <- NULL



OsC <- biologicToolsSC::addDf2seuratMetaData(
    obj = OsC, 
    dfAdd = dfRes
)

OsC@meta.data[is.na(OsC@meta.data)] <- 0



save(OsC,
     file = paste0(
         Obio@parameterList$localWorkDir,
         Obio@parameterList$project_id,
         ".Seurat.Robj"
     )
)


## Path in which the cellranger count outputs recide:
pathToCountDir <- c("/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/484A_scRNAseq_Characterisation_EpCAM_CD90_thymic_epithelial_cells_SC22096_hash/FASTQ_files/")

## Names of the cellranger count ids from the alignment to the human/mouse transcriptomes
sampleIDvec <- c(
    "RAG4782A1mmhs",
    "RAG4782A2mmhs",
    "RAG4782A3mmhs",
    "RAG4782A4mmhs",
    "RAG4782A5mmhs",
    "RAG4782A6mmhs",                         
    "RAG4782A7mmhs",
    "RAG4782A8mmhs"
)

for (i in 1:length(sampleIDvec)){
        sampleID <- sampleIDvec[i]
        dataDir <- paste0(pathToCountDir,sampleID,"/outs/filtered_feature_bc_matrix")
        
        
        
        assign(
            "fullMat", #names(obj@parameterList[[obj@parameterList$inputMode]])[i],
            Seurat::Read10X(data.dir = dataDir, gene.column = 2)
        )
        
        mmMat <- fullMat[grep("mm10", row.names(fullMat)), ]
        hsMat <- fullMat[grep("GRCh38", row.names(fullMat)), ]
        
        Omm = Seurat::CreateSeuratObject(
            counts = mmMat,
            project = sampleID,
            min.cells = 0,
            min.features = 0
        )
        
        Ohs = Seurat::CreateSeuratObject(
            counts = hsMat ,
            project = sampleID,
            min.cells = 0,
            min.features = 0
        )
        
        NmmGenes <- nrow(Omm)
        NhsGenes <- nrow(Ohs)
        
        dfMetaHS <- Ohs@meta.data
        dfMetaHS$orig.ident <- NULL
        dfMetaHS[["Perc_of_all_genes_in_cell"]] <- round(100*(dfMetaHS$nFeature_RNA / NhsGenes),2)
        names(dfMetaHS) <- paste0("hs_", names(dfMetaHS))
        dfMetaHS[["cellID"]] <- row.names(dfMetaHS)
        
        dfMetaMM <- Omm@meta.data
        dfMetaMM$orig.ident <- NULL
        dfMetaMM[["Perc_of_all_genes_in_cell"]] <- round(100*(dfMetaMM$nFeature_RNA / NmmGenes),2)
        names(dfMetaMM) <- paste0("mm_", names(dfMetaMM))
        dfMetaMM[["cellID"]] <- row.names(dfMetaMM)
        
        dfMeta <- merge(
            dfMetaHS, 
            dfMetaMM, 
            by.x = "cellID",
            by.y = "cellID"
        )
        
        dfMeta[["inferred_species"]] <- ifelse(dfMeta$mm_Perc_of_all_genes_in_cell > dfMeta$hs_Perc_of_all_genes_in_cell, "mouse", "human")
        dfMeta[["limsID"]] <- sampleID
        
        if (i ==1){
            dfRes <- dfMeta
        } else {
            dfRes <- rbind(
                dfRes,
                dfMeta
            )
        }
}

## Get sample annotation from design file ##

## Load Seurat object to annotate mouse cells in ##
## Load Seurat object

library(biologicSeqTools2)

ObioFN <- paste0("../", list.files("..")[grep(".bioLOGIC.Robj", list.files(".."))])

if (file.exists(ObioFN)){
    load(paste0(ObioFN))
    print(paste0("Obio object ", Obio@parameterList$localWorkDir, ObioFN, " exists and is loaded."))
} else {
    exit()
}

## Reset paths to local environment
Obio <- setMountingPoint(Obio)
Obio <- setAnalysisPaths(Obio)
Obio <- setCrickGenomeAndGeneNameTable(Obio)
Obio <- createAnalysisFolders(
    Obio
)
Obio <- setDataBaseParameters(Obio)




SeuratFN <- paste0(Obio@parameterList$localWorkDir,list.files(Obio@parameterList$localWorkDir)[grep(".Seurat.Robj", list.files(Obio@parameterList$localWorkDir))])


if (file.exists(SeuratFN)){
    load(SeuratFN)
    print(paste0("Obio object ", Obio@parameterList$localWorkDir,SeuratFN, " exists and is loaded."))
    
} else {
    exit()
}

sampleAnnotation <- list()

for (i in 1:length(Obio@sampleDetailList)){
    sampleAnnotation[[names(Obio@sampleDetailList)[i]]] <- Obio@sampleDetailList[[i]][["limsId"]]
}
 
## Challenge here: deal with cell ranger aggregate ##
## In dfRes, cellIDs are -1_1 for RAG4782A1 for A1A2MERGED
##           and         -2_1 for RAG4782A2 for A1A2MERGED
## In dfRes, cellIDs are -1_1 for RAG4782A3 for A3A4MERGED
##           and         -2_1 for RAG4782A4 for A3A4MERGED
dfRes[["suffix"]] <- ""

dfRes$limsID <- gsub("mmhs", "", dfRes$limsID)
dfRes[dfRes$limsID == "RAG4782A1", "suffix"] <- "-1"
dfRes[dfRes$limsID == "RAG4782A2", "suffix"] <- "-2"
dfRes[dfRes$limsID == "RAG4782A3", "suffix"] <- "-1"
dfRes[dfRes$limsID == "RAG4782A4", "suffix"] <- "-2"

dfRes$cellID <- gsub("-1", "", dfRes$cellID)
dfRes$cellID <- paste0(dfRes$cellID, dfRes$suffix)

## Switch lims IDs to merged IDS
dfRes[dfRes$limsID == "RAG4782A1", "limsID"] <- "A1A2MERGED"
dfRes[dfRes$limsID == "RAG4782A2", "limsID"] <- "A1A2MERGED"
dfRes[dfRes$limsID == "RAG4782A3", "limsID"] <- "A3A4MERGED"
dfRes[dfRes$limsID == "RAG4782A4", "limsID"] <- "A3A4MERGED"

for (i in 1:length(sampleAnnotation)){
    dfRes[dfRes$limsID == sampleAnnotation[[i]], "sampleID"] <- names(sampleAnnotation)[i]
}




dfID <- unique(OsC@meta.data[,c("cellID", "sampleID")])
dfID[["suffix"]] <- sapply(dfID$cellID, function(x) unlist(strsplit(x, "_"))[2])

dfID <- unique(dfID[,c("suffix", "sampleID")])

for (i in 1:length(sampleAnnotation)){
    dfRes[dfRes$sampleID == names(sampleAnnotation)[i], "cellID"] <- paste0(dfRes[dfRes$sampleID == names(sampleAnnotation)[i], "cellID"], "_", as.vector(dfID[dfID$sampleID == names(sampleAnnotation)[i], "suffix"]))
}


#names(dfRes) <- gsub("hs_Perc_of_all_genes_in_cell", "meta_hs_Perc_of_all_genes_in_cell",names(dfRes))
#names(dfRes) <- gsub("mm_Perc_of_all_genes_in_cell", "meta_mm_Perc_of_all_genes_in_cell",names(dfRes))
names(dfRes) <- gsub("meta_meta_", "meta_",names(dfRes))



row.names(dfRes) <- dfRes$cellID
dfRes$cellID <- NULL
dfRes$sampleID <- NULL
dfRes$limsID <- NULL


OsC <- biologicToolsSC::addDf2seuratMetaData(
    obj = OsC, 
    dfAdd = dfRes
)

OsC@meta.data[is.na(OsC@meta.data)] <- 0



save(OsC,
     file = paste0(
         Obio@parameterList$localWorkDir,
         Obio@parameterList$project_id,
         ".Seurat.Robj"
     )
)

# devtools::install_github("decusInLabore/biologicSeqTools")
# 
# 
# ObioFN <- "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449B_scRNAseq_human_thymus_cultures_SC20193_skin_thymus_custom_anchors/workdir/SC20193customA.bioLOGIC.Robj"
# load(ObioFN)
# 
# dfdbTable <- OsC@meta.data
# columnDbCatList <- biologicSeqTools::inferDBcategories(dfdbTable)
# 
# if (dir.exists("/Volumes/babs/working/boeings/")){
#     hpc.mount <- "/Volumes/babs/working/boeings/"
# } else if (dir.exists("Y:/working/boeings/")){
#     hpc.mount <- "Y:/working/boeings/"
# } else if (dir.exists("/camp/stp/babs/working/boeings/")){
#     hpc.mount <- "/camp/stp/babs/working/boeings/"
# } else {
#     hpc.mount <- ""
# }
# #Create the environment and load a suitable version of R, e.g. so:
# FN <- paste0(hpc.mount, "Projects/reference_data/documentation/BC.parameters.txt")
# dbTable <- read.delim(
#     FN, 
#     sep = "\t",
#     stringsAsFactors = F
# )
# db.pwd <- as.vector(dbTable[1,1])
# 
# 
# biologicSeqTools::upload.datatable.to.database(
#     host = Obio@dbDetailList$host,
#     user = Obio@dbDetailList$db.user,
#     password = db.pwd,
#     prim.data.db = Obio@dbDetailList$primDataDB,
#     dbTableName = Obio@parameterList$PCAdbTableName,
#     df.data = dfdbTable,
#     db.col.parameter.list = columnDbCatList,
#     new.table = TRUE
# )
# biologicSeqTools::killDbConnections()