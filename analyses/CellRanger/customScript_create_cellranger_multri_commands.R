# module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;R

################################################################################
## R-function to make cellranger count commands                               ##

folder <- "/camp/stp/babs/working/boeings/Projects/hillc/foteini.papaleonidopoulou/505_scRNAseq_nodal_signaling_in_mouse_ES_cells_SC22177/"

cmdVec <- "module load CellRanger/6.0.1-bcl2fastq-2.20.0;ml use /camp/apps/eb/dev/modules/all;ml --ignore-cache bcl2fastq2/2.20.0-foss-2018b"
transcriptome <- "/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-mm10-2020-A"

folderVec <- c(
    "/camp/stp/sequencing/outputs/babs/data/hillc/foteini.papaleonidopoulou/SC22177/primary_data/220726_A01366_0242_BHN3GKDSX3/fastq",
    "/camp/stp/sequencing/outputs/babs/data/hillc/foteini.papaleonidopoulou/SC22177/primary_data/220727_A01366_0244_BHCJJCDMXY/fastq",
    "/camp/stp/sequencing/outputs/babs/data/hillc/foteini.papaleonidopoulou/SC22177/primary_data/220811_A01366_0259_AHCGJVDMXY/fastq"
)

## Create a list to determine which fastq folders to quote
fileList <- list()
fileList <- sapply(folderVec, function(x) list.files(x))
fileList <- lapply(fileList, function(x) sapply(x, function(y) unlist(strsplit(y, "_"))[1] ))
fileList <- lapply(fileList, function(x) unique(as.vector(x)))


featureTypeList <- list(
    "gene expression" = c("PAP4934A1","PAP4934A2","PAP4934A3","PAP4934A4"),
    "Multiplexing Capture" = c("PAP4934A5","PAP4934A6","PAP4934A7","PAP4934A8")
)

samplesList <- list(
    "PAP4934A1" = list(
        "PAP4934A1CMO301" = c("CMO301", "d0_WT"),
        "PAP4934A1CMO302" = c("CMO302", "d0_KO_S2"),
        "PAP4934A1CMO303" = c("CMO303", "d0_KO_S3"),
        "PAP4934A1CMO304" = c("CMO304", "d0_KO_S4")
    ),
    "PAP4934A2" = list(
        "PAP4934A2CMO305" = c("CMO305", "d2_WT"),
        "PAP4934A2CMO302" = c("CMO306", "d2_KO_S2"),
        "PAP4934A2CMO303" = c("CMO307", "d2_KO_S3"),
        "PAP4934A2CMO304" = c("CMO308", "d2_KO_S4")
    ),
    "PAP4934A3" = list(
        "PAP4934A3CMO309" = c("CMO309", "d4_WT"),
        "PAP4934A3CMO310" = c("CMO310", "d4_KO_S2"),
        "PAP4934A3CMO311" = c("CMO311", "d4_KO_S3"),
        "PAP4934A3CMO312" = c("CMO312", "d4_KO_S4")
    ),
    "PAP4934A4" = list(
        "PAP4934A4CMO305" = c("CMO305", "d6_WT"),
        "PAP4934A4CMO306" = c("CMO306", "d6_KO_S2"),
        "PAP4934A4CMO307" = c("CMO307", "d6_KO_S3"),
        "PAP4934A4CMO308" = c("CMO308", "d6_KO_S4")
    )
)

## Create on cellranger multi config file per list entry

cmdList <- purrr::map( 1:length( samplesList ), function( i ) {
    sampleID <- names(samplesList)[i]
    ## Create header
    dfConfig <- rbind(
        c("[gene-expression]","",	"# Section describing configuration of Gene Expression libraries only"),
        c("reference", transcriptome, "# Path of folder containing 10x-compatible reference. Required for gene expression and Feature Barcode libraries."),
        c("no-bam",	"TRUE",	"# Disables BAM(.bai) file generation"),
        c("", "", "")
    )

    dfConfig <- data.frame(cbind(
        dfConfig,
        rep("", nrow(dfConfig))
    ))

    ## Add libraries section
    dfLib <- data.frame(rbind(
        c("[libraries]", "", "# Table describing the input sequencing libraries. Must contain the header and 1 or more rows. Required.", ""),
        c("fastq_id", "fastqs", "lanes", "feature_types")
    ))

    dfTemp <- purrr::map( 1:length( featureTypeList ), function( i ) {
            dfTemp2 <- purrr::map( 1:length( featureTypeList[[i]] ), function( j ) {
                ## Get all relevant fastq file folder
                sampleID <- featureTypeList[[i]][j]
                pos <- unlist(lapply(fileList, function(x)  grep(paste0("^", sampleID,"$"), x)))
                folderVec <- names(pos)

                for (m in 1:length(folderVec)){
                    newRow <- c(sampleID, folderVec[m], "any", names(featureTypeList)[i])

                    if (m ==1){
                        dfTempRes <- newRow
                    } else {
                        dfTempRes <- data.frame(rbind(
                            dfTempRes,
                            newRow
                        ))
                    }
                }

                row.names(dfTempRes) <- NULL

                return(dfTempRes)
            })

        dfTemp3 <- do.call(rbind, dfTemp2)
        return(dfTemp3)
    }
    )

    dfTemp5 <- do.call(rbind, dfTemp)

    dfLib <- rbind(
        dfLib,
        dfTemp5
    )




    ## Create samples section
    tempRes <- purrr::map( 1:length( samplesList[[i]] ), function( j ) {
        newRow <- c(names( samplesList[[i]][j]), unlist(samplesList[[i]][[j]]))
        return(newRow)
    }
    )
    res2 <- do.call(rbind, tempRes)

    dfSampleSection <- rbind(
        c("[samples]","",""),
        c("sample_id",	"cmo_ids",	"description"),
        res2
    )

    dfSampleSection <- data.frame(cbind(
        dfSampleSection,
        rep("", nrow(dfSampleSection))
    ))

    names(dfSampleSection) <- names(dfLib)
    names(dfConfig) <- names(dfLib)

    dfSettings <- rbind(
        dfConfig,
        dfLib,
        dfSampleSection
    )

    FNconf <- paste0(folder, "FASTQ_files/", sampleID, ".cellranger.multi.config.csv")
    readr::write_csv(dfSettings, FNconf)

    cmd <- paste0(
        'sbatch --time=72:00:00 --wrap "cellranger multi --id ',
        sampleID,
        ' --csv ',
        FNconf,
        '" -c 16 --mem-per-cpu=7000 -o CR.',sampleID,'.slurm >> commands.txt'
    )

    return(cmd)

}
)



cmdVec <- c(
    "module load CellRanger/6.0.1-bcl2fastq-2.20.0;ml use /camp/apps/eb/dev/modules/all;ml --ignore-cache bcl2fastq2/2.20.0-foss-2018b",
    unlist(cmdList)
)



## write to file
FN <- paste0(folder, "FASTQ_files/cellranger.multi.commands.sh")
sink(FN)

for (i in 1:length(cmdVec)){
    cat(cmdVec[i]); cat("\n")
}

sink()
