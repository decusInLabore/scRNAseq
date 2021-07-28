################################

## CellRanger commands ##
## ml CellRanger/3.0.2-bcl2fastq-2.20.0

##sbatch --time=71:00:00 --wrap "cellranger count --id=RAG1452A1 --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-and-mm10-3.1.0 --fastqs=/camp/stp/babs/inputs/sequencing/fastq/201127_K00102_0540_AHK7CNBBXY/fastq/SC20193,/camp/stp/babs/inputs/sequencing/fastq/201210_K00102_0544_AHKJWFBBXY/fastq/SC20193 --sample=RAG1452A1" --job-name="CRRAG1452A1" -c 16 --mem-per-cpu=7000 -o cellranger.RAG1452A1.slurm >> commands.txt
##sbatch --time=71:00:00 --wrap "cellranger count --id=RAG1452A2 --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-and-mm10-3.1.0 --fastqs=/camp/stp/babs/inputs/sequencing/fastq/201127_K00102_0540_AHK7CNBBXY/fastq/SC20193 --sample=RAG1452A2" --job-name="CRRAG1452A2" -c 16 --mem-per-cpu=7000 -o cellranger.RAG1452A2.slurm >> commands.txt
##sbatch --time=71:00:00 --wrap "cellranger count --id=RAG1452A3 --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-and-mm10-3.1.0 --fastqs=/camp/stp/babs/inputs/sequencing/fastq/201127_K00102_0540_AHK7CNBBXY/fastq/SC20193,/camp/stp/babs/inputs/sequencing/fastq/201210_K00102_0544_AHKJWFBBXY/fastq/SC20193 --sample=RAG1452A3" --job-name="CRRAG1452A3" -c 16 --mem-per-cpu=7000 -o cellranger.RAG1452A3.slurm >> commands.txt
##sbatch --time=71:00:00 --wrap "cellranger count --id=RAG1452A4 --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-and-mm10-3.1.0 --fastqs=/camp/stp/babs/inputs/sequencing/fastq/201127_K00102_0540_AHK7CNBBXY/fastq/SC20193,/camp/stp/babs/inputs/sequencing/fastq/201210_K00102_0544_AHKJWFBBXY/fastq/SC20193 --sample=RAG1452A4" --job-name="CRRAG1452A4" -c 16 --mem-per-cpu=7000 -o cellranger.RAG1452A4.slurm >> commands.txt
##sbatch --time=71:00:00 --wrap "cellranger count --id=RAG1452A5 --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-and-mm10-3.1.0 --fastqs=/camp/stp/babs/inputs/sequencing/fastq/201127_K00102_0540_AHK7CNBBXY/fastq/SC20193,/camp/stp/babs/inputs/sequencing/fastq/201210_K00102_0544_AHKJWFBBXY/fastq/SC20193 --sample=RAG1452A5" --job-name="CRRAG1452A5" -c 16 --mem-per-cpu=7000 -o cellranger.RAG1452A5.slurm >> commands.txt
##sbatch --time=71:00:00 --wrap "cellranger count --id=RAG1452A6 --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-and-mm10-3.1.0 --fastqs=/camp/stp/babs/inputs/sequencing/fastq/201127_K00102_0540_AHK7CNBBXY/fastq/SC20193,/camp/stp/babs/inputs/sequencing/fastq/201210_K00102_0544_AHKJWFBBXY/fastq/SC20193 --sample=RAG1452A6" --job-name="CRRAG1452A6" -c 16 --mem-per-cpu=7000 -o cellranger.RAG1452A6.slurm >> commands.txt
##sbatch --time=71:00:00 --wrap "cellranger count --id=RAG1452A7 --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-and-mm10-3.1.0 --fastqs=/camp/stp/babs/inputs/sequencing/fastq/201127_K00102_0540_AHK7CNBBXY/fastq/SC20193,/camp/stp/babs/inputs/sequencing/fastq/201210_K00102_0544_AHKJWFBBXY/fastq/SC20193 --sample=RAG1452A7" --job-name="CRRAG1452A7" -c 16 --mem-per-cpu=7000 -o cellranger.RAG1452A7.slurm >> commands.txt
##sbatch --time=71:00:00 --wrap "cellranger count --id=RAG1452A8 --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-and-mm10-3.1.0 --fastqs=/camp/stp/babs/inputs/sequencing/fastq/201127_K00102_0540_AHK7CNBBXY/fastq/SC20193,/camp/stp/babs/inputs/sequencing/fastq/201210_K00102_0544_AHKJWFBBXY/fastq/SC20193 --sample=RAG1452A8" --job-name="CRRAG1452A8" -c 16 --mem-per-cpu=7000 -o cellranger.RAG1452A8.slurm >> commands.txt


sampleIDvec <- c(
    "RAG1452A1",
    "RAG1452A2",
    "RAG1452A3",
    "RAG1452A4",
    "RAG1452A5",
    "RAG1452A6",                         
    "RAG1452A7",
    "RAG1452A8"
)

for (i in 1:length(sampleIDvec)){
        sampleID <- sampleIDvec[i]
        dataDir <- paste0("/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/413_scRNAseq_human_thymus_cultures_SC20193/FASTQ_files_hs_mm/",sampleID,"/outs/filtered_feature_bc_matrix")
        
        
        
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
 
sampleAnnotation <- list(
    "Skin01" = "RAG1452A1",
    "Oesophagus20" =  "RAG1452A2",
    "ThyP89" =  "RAG1452A3",
    "ThyP94" =  "RAG1452A4",
    "MedullaThy89" = "RAG1452A5",
    "MedullaThy94" =  "RAG1452A6",
    "CortexThy89" =  "RAG1452A7",
    "CortexThy94" =  "RAG1452A8"
)

dfRes[["sampleID"]] <- ""

for (i in 1:length(sampleAnnotation)){
    dfRes[dfRes$limsID == sampleAnnotation[[i]], "sampleID"] <- names(sampleAnnotation)[i]
}


OFN <- "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449B_scRNAseq_human_thymus_cultures_SC20193_skin_thymus_custom_anchors/workdir/SC20193customA.Seurat.Robj"
load(OFN)

dfID <- unique(OsC@meta.data[,c("cellID", "sampleID")])
dfID[["suffix"]] <- sapply(dfID$cellID, function(x) unlist(strsplit(x, "_"))[2])

dfID <- unique(dfID[,c("suffix", "sampleID")])

for (i in 1:length(sampleAnnotation)){
    dfRes[dfRes$sampleID == names(sampleAnnotation)[i], "cellID"] <- gsub("-1", paste0("_", as.vector(dfID[dfID$sampleID == names(sampleAnnotation)[i], "suffix"])),dfRes[dfRes$sampleID == names(sampleAnnotation)[i], "cellID"])
}

row.names(dfRes) <- dfRes$cellID
dfRes$cellID <- NULL
dfRes$sampleID <- NULL
dfRes$limsID <- NULL

devtools::install_github("decusInLabore/biologicToolsSC")

OsC <- biologicToolsSC::addDf2seuratMetaData(
    obj = OsC, 
    dfAdd = dfRes
)

OsC@meta.data[is.na(OsC@meta.data)] <- 0

names(OsC@meta.data) <- gsub("inferred_species", "meta_inferred_species",names(OsC@meta.data))
names(OsC@meta.data) <- gsub("hs_Perc_of_all_genes_in_cell", "meta_hs_Perc_of_all_genes_in_cell",names(OsC@meta.data))
names(OsC@meta.data) <- gsub("mm_Perc_of_all_genes_in_cell", "meta_mm_Perc_of_all_genes_in_cell",names(OsC@meta.data))
names(OsC@meta.data) <- gsub("meta_meta_", "meta_",names(OsC@meta.data))
OsC@meta.data <- OsC@meta.data[,!duplicated(names(OsC@meta.data))]


save(
    OsC,
    file = OFN
)

devtools::install_github("decusInLabore/biologicSeqTools")


ObioFN <- "/camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449B_scRNAseq_human_thymus_cultures_SC20193_skin_thymus_custom_anchors/workdir/SC20193customA.bioLOGIC.Robj"
load(ObioFN)

dfdbTable <- OsC@meta.data
columnDbCatList <- biologicSeqTools::inferDBcategories(dfdbTable)

if (dir.exists("/Volumes/babs/working/boeings/")){
    hpc.mount <- "/Volumes/babs/working/boeings/"
} else if (dir.exists("Y:/working/boeings/")){
    hpc.mount <- "Y:/working/boeings/"
} else if (dir.exists("/camp/stp/babs/working/boeings/")){
    hpc.mount <- "/camp/stp/babs/working/boeings/"
} else {
    hpc.mount <- ""
}
#Create the environment and load a suitable version of R, e.g. so:
FN <- paste0(hpc.mount, "Projects/reference_data/documentation/BC.parameters.txt")
dbTable <- read.delim(
    FN, 
    sep = "\t",
    stringsAsFactors = F
)
db.pwd <- as.vector(dbTable[1,1])


biologicSeqTools::upload.datatable.to.database(
    host = Obio@dbDetailList$host,
    user = Obio@dbDetailList$db.user,
    password = db.pwd,
    prim.data.db = Obio@dbDetailList$primDataDB,
    dbTableName = Obio@parameterList$PCAdbTableName,
    df.data = dfdbTable,
    db.col.parameter.list = columnDbCatList,
    new.table = TRUE
)
biologicSeqTools::killDbConnections()