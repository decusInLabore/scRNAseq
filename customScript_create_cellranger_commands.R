# module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;R

################################################################################
## R-function to make cellranger count commands                               ##

cmdVec <- "ml CellRanger/5.0.0-bcl2fastq-2.20.0"
transcriptome <- "/camp/svc/reference/Genomics/10x/10x_transcriptomes/Gallus_gallus-6.0-release-97"

folderVec <- c(
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

sink("../../FASTQ_files/cellranger.commands.sh")

for (i in 1:length(cmdVec)){
    cat(cmdVec[i]); cat("\n")
}

sink()
