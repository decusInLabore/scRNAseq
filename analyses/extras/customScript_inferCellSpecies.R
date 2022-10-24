###############################################################################
## Assign human/mouse cells                                                  ##

# The goal of this script is to create a table that infers the species of 
# each cell in a single-cell experiment

# In step 1 we run cellranger count using the combined human/mouse transcriptome
# on each individual sample. 

# Then we transform that into a table with the correct cellIDs in the integrated
# single cell experiment. 

###############################################################################
## cellranger count commands                                                 ##

# ml CellRanger/5.0.0-bcl2fastq-2.20.0
# 
# sbatch --time=72:00:00 --wrap "cellranger count --id=RAG4782A1mmhs --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-and-mm10-2020-A --fastqs=/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220426_A01366_0185_BHCH3GDMXY/fastq,/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220505_A01366_0193_BHKC3MDSX3/fastq --sample=RAG4782A1" --job-name="CR.A1" -c 16 --mem-per-cpu=7000 -o CR.A1.slurm >> commands.txt
# sbatch --time=72:00:00 --wrap "cellranger count --id=RAG4782A2mmhs --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-and-mm10-2020-A --fastqs=/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220426_A01366_0185_BHCH3GDMXY/fastq,/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220505_A01366_0193_BHKC3MDSX3/fastq --sample=RAG4782A2" --job-name="CR.A2" -c 16 --mem-per-cpu=7000 -o CR.A2.slurm >> commands.txt
# sbatch --time=72:00:00 --wrap "cellranger count --id=RAG4782A3mmhs --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-and-mm10-2020-A --fastqs=/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220426_A01366_0185_BHCH3GDMXY/fastq,/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220505_A01366_0193_BHKC3MDSX3/fastq --sample=RAG4782A3" --job-name="CR.A3" -c 16 --mem-per-cpu=7000 -o CR.A3.slurm >> commands.txt
# sbatch --time=72:00:00 --wrap "cellranger count --id=RAG4782A4mmhs --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-and-mm10-2020-A --fastqs=/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220426_A01366_0185_BHCH3GDMXY/fastq,/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220505_A01366_0193_BHKC3MDSX3/fastq --sample=RAG4782A4" --job-name="CR.A4" -c 16 --mem-per-cpu=7000 -o CR.A4.slurm >> commands.txt
# sbatch --time=72:00:00 --wrap "cellranger count --id=RAG4782A5mmhs --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-and-mm10-2020-A --fastqs=/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220426_A01366_0185_BHCH3GDMXY/fastq,/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220505_A01366_0193_BHKC3MDSX3/fastq --sample=RAG4782A5" --job-name="CR.A5" -c 16 --mem-per-cpu=7000 -o CR.A5.slurm >> commands.txt
# sbatch --time=72:00:00 --wrap "cellranger count --id=RAG4782A6mmhs --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-and-mm10-2020-A --fastqs=/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220426_A01366_0185_BHCH3GDMXY/fastq,/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220505_A01366_0193_BHKC3MDSX3/fastq --sample=RAG4782A6" --job-name="CR.A6" -c 16 --mem-per-cpu=7000 -o CR.A6.slurm >> commands.txt
# sbatch --time=72:00:00 --wrap "cellranger count --id=RAG4782A7mmhs --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-and-mm10-2020-A --fastqs=/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220426_A01366_0185_BHCH3GDMXY/fastq,/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220505_A01366_0193_BHKC3MDSX3/fastq --sample=RAG4782A7" --job-name="CR.A7" -c 16 --mem-per-cpu=7000 -o CR.A7.slurm >> commands.txt
# sbatch --time=72:00:00 --wrap "cellranger count --id=RAG4782A8mmhs --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-and-mm10-2020-A --fastqs=/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220426_A01366_0185_BHCH3GDMXY/fastq,/camp/stp/sequencing/outputs/babs/data/bonfantip/roberta.ragazzini/SC22096/primary_data/220505_A01366_0193_BHKC3MDSX3/fastq --sample=RAG4782A8" --job-name="CR.A8" -c 16 --mem-per-cpu=7000 -o CR.A8.slurm >> commands.txt

## Run cellranger count on the command line                                  ##
###############################################################################

###############################################################################
## Processing of the cellranger count output and assignment to cellIDs in the##
## integrated single-cell experiment                                         ##

# This is the R-version we're using
# module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a

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

## Infering the species for each cell in each individual sample ##

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
        
        
        ## Whatever species has a higher percentage of all possible species genes detected will be called as such
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

## Done infering species for each individual sample                          ##
###############################################################################

###############################################################################
## Add cellIDs of integrated experiment                                      ##

## Load integrated single cell data object ##
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

## Samples that had to be merged with cellranger aggr will feature a -1, -2, 
## -n suffix. Here true for samples A1, A2, A3, A4.

dfRes$limsID <- gsub("mmhs", "", dfRes$limsID)
dfRes[dfRes$limsID == "RAG4782A1", "suffix"] <- "-1"
dfRes[dfRes$limsID == "RAG4782A2", "suffix"] <- "-2"
dfRes[dfRes$limsID == "RAG4782A3", "suffix"] <- "-1"
dfRes[dfRes$limsID == "RAG4782A4", "suffix"] <- "-2"

dfRes$cellID <- gsub("-1", "", dfRes$cellID)
dfRes$cellID <- paste0(dfRes$cellID, dfRes$suffix)

## After making the cellIDs for A1/A2 and A3/A4 unique, as cellranger aggr does,
## we can change the sampleID to the merged sample ID
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


names(dfRes) <- gsub("hs_nFeature_RNA", "meta_hs_nFeature_RNA",names(dfRes))
names(dfRes) <- gsub("mm_nFeature_RNA", "meta_mm_nFeature_RNA",names(dfRes))

names(dfRes) <- gsub("inferred_species", "meta_inferred_species",names(dfRes))
names(dfRes) <- gsub("inferred_species", "meta_inferred_species",names(dfRes))

#names(dfRes) <- gsub("hs_Perc_of_all_genes_in_cell", "meta_hs_Perc_of_all_genes_in_cell",names(dfRes))
#names(dfRes) <- gsub("mm_Perc_of_all_genes_in_cell", "meta_mm_Perc_of_all_genes_in_cell",names(dfRes))
names(dfRes) <- gsub("meta_meta_", "meta_",names(dfRes))


## Done creating assignment table                                            ##
###############################################################################

###############################################################################
## Save table to be integrated into biologicViewerSC                         ##

dfRes$sampleID <- NULL
dfRes$limsID <- NULL
dfRes$suffix <- NULL

FNout <- "../temp/species.inference.txt"
write.table(
    dfRes, 
    FNout, 
    sep="\t"
)

## The above table will be picked up and integarted by the sc_PartC_ module
###############################################################################