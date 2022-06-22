################################################################################
## R-function to make cellranger count commands                               ##

cmdVec <- "ml CellRanger/5.0.0-bcl2fastq-2.20.0"
transcriptome <- "/camp/svc/reference/Genomics/10x/10x_transcriptomes/Gallus_gallus-6.0-release-97"

folderVec <- c(
    "/camp/stp/sequencing/outputs/babs/data/lovell-badger/gunes.taylor/SC18208/primary_data/181102_D00446_0276_ACCHYNANXX/fastq",
    "/camp/stp/sequencing/outputs/babs/data/lovell-badger/gunes.taylor/SC18208/primary_data/181221_K00102_0289_AH32Y7BBXY/fastq",  
    
    "/camp/stp/sequencing/outputs/babs/data/lovell-badger/gunes.taylor/SC18208/primary_data/190328_K00102_0324_AH5KG3BBXY/fastq",
    "/camp/stp/sequencing/outputs/babs/data/lovell-badger/gunes.taylor/SC18208/primary_data/190705_K00102_0362_BHC53JBBXY/fastq",
    "/camp/stp/sequencing/outputs/babs/data/lovell-badger/gunes.taylor/SC18208/primary_data/190823_K00102_0384_BHCLJKBBXY/fastq"
)


fileList <- list()

fileList <- sapply(folderVec, function(x) list.files(x))

fileList <- lapply(fileList, function(x) sapply(x, function(y) unlist(strsplit(y, "_"))[1] ))

fileList <- lapply(fileList, function(x) unique(as.vector(x)))

## Make a folder list for each sample ##

## First, get all sample IDs


allIDs <- sort(unique(unlist(fileList)))

## remove IDs that have already been successfully processed:
processedIDs <- list.files("/camp/stp/babs/working/boeings/Projects/lovellbadger/gunes.taylor/490_scRNAseq_SC18208/FASTQ_files/")
processedIDs <- processedIDs[processedIDs %in% allIDs]

allIDs <- allIDs[!(allIDs %in% processedIDs)]

for (i in 1:length(allIDs)){
    sampleID <- allIDs[i]
    pos <- unlist(lapply(fileList, function(x)  grep(paste0("^", sampleID,"$"), x)))
    
    if (length(pos) > 0){
       folderString <- paste0(names(pos), collapse = ",")
    }
    
    cmd <- paste0(
        'sbatch --time=72:00:00 --wrap "cellranger count --id=',
        sampleID,
        ' --transcriptome=',
        transcriptome,
        ' --fastqs=',
        folderString,
        ' --sample=',
        sampleID, 
        '" --job-name="'
        ,sampleID,
        '" -c 16 --mem-per-cpu=7000 -o CR.',sampleID,'.slurm >> commands.txt'
    )
    
    cmdVec <- c(
        cmdVec, 
        cmd
    )
    
}

## write to file

sink("../../FASTQ_files/cellranger.commands.part2.sh")

for (i in 1:length(cmdVec)){
    cat(cmdVec[i]); cat("\n")
}

sink()