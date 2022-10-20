if (dir.exists("/Volumes/babs/working/boeings/")){
    hpc.mount <- "/Volumes/babs/working/boeings/"
} else if (dir.exists("Y:/working/boeings/")){
    hpc.mount <- "Y:/working/boeings/"
} else if (dir.exists("/camp/stp/babs/working/boeings/")){
    hpc.mount <- "/camp/stp/babs/working/boeings/"
} else {
    hpc.mount <- ""
}

FN <- paste0(hpc.mount, "Projects/reference_data/documentation/BC.parameters.txt")
dbTable <- read.delim(
    FN,
    sep = "\t",
    stringsAsFactors = F
)

db.pwd <- as.vector(dbTable[1,1])


library(bioLOGIC)

###############################################################################
## Create animation file                                                     ##


createAnimationFrame <- function(
    cellID = "cellID",
    UMAP_1 = "UMAP_1",
    UMAP_2 = "UMAP_2",
    clusterName = "clusterName",
    clusterColor = "clusterColor",
    state = "stage name 1",
    dfMeta = OsC@meta.data
){
    
    dfTemp <- dfMeta
    dfTemp <- dfMeta[,c(cellID,UMAP_1, UMAP_2, clusterName, clusterColor)]
    dfTemp[["state"]] <- state
    
    dfTemp <- dfTemp[!is.na(dfTemp[,UMAP_1]), ]
    dfTemp <- dfTemp[!is.na(dfTemp[,UMAP_2]), ]
    dfTemp <- dfTemp[dfTemp[,UMAP_1] != 0 & dfTemp[,UMAP_2] != 0, ]
    
    if (is.null(clusterColor)){
        colOpts <- sort(unique(dfTemp[,clusterName]))
        cols <- scales::hue_pal()(length(colOpts))
        dfCol <- data.frame(colOpts, cols)
        names(dfCol) <- c(clusterName, "clusterColor")
        
        dfTemp <- dplyr::full_join(
            dfTemp,
            dfCol,
            by = clusterName
        )
        
    }
    
    names(dfTemp) <- gsub(cellID, "cellID", names(dfTemp))
    names(dfTemp) <- gsub(UMAP_1, "UMAP_1", names(dfTemp))
    names(dfTemp) <- gsub(UMAP_2, "UMAP_2", names(dfTemp))
    names(dfTemp) <- gsub(clusterName, "clusterName", names(dfTemp))
    
    if (!is.null(clusterColor)){
        names(dfTemp) <- gsub(clusterColor, "clusterColor", names(dfTemp))
    }
    
    dfAdd <- dfTemp[, c("cellID","UMAP_1", "UMAP_2", "clusterName", "clusterColor", "state")]
    return(dfAdd)
}



## Done                                                                      ##
###############################################################################


###############################################################################
## First frame                                                               ##

stateName <- "Experiment A4A6A7A8 Before Subclustering"
orderVec <- c(stateName)

dfMeta <- import.db.table.from.db(
    dbtable = "A4A6A7A8_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
)

dfAdd <- createAnimationFrame(
    cellID = "cellID",
    UMAP_1 = "UMAP_1",
    UMAP_2 = "UMAP_2",
    clusterName = "clusterName",
    clusterColor = NULL,
    state = stateName,
    dfMeta = dfMeta
)

## Get cluster colors ##
FNcol <- paste0(hpc.mount, "projects/bonfantip/roberta.ragazzini/471_scRNAseq_newSamples_Epi6_Epi10_SC20193_SC21201/workdir/scRNAseq/design/clusterAnnotation.txt")

dfCol <- read.delim(
    FNcol, 
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

dfCol <- unique(dfCol[,c("clusterName", "clusterColor")])

dfAdd$clusterColor <- NULL

dfAdd <- dplyr::full_join(
    dfAdd,
    dfCol,
    by = clusterName
)



dfReel <- dfAdd

## Done                                                                      ##
###############################################################################

###############################################################################
## Second Frame                                                              ##
## Change color for clusters to be removed to black
## Clusters marked for removal: C4, C8, C11, C13

stateName <- "Contaminating cells to be removed are marked black"
orderVec <- c(orderVec, stateName)

dfMeta <- import.db.table.from.db(
    dbtable = "A4A6A7A8_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
)

dfAdd <- createAnimationFrame(
    cellID = "cellID",
    UMAP_1 = "UMAP_1",
    UMAP_2 = "UMAP_2",
    clusterName = "clusterName",
    clusterColor = NULL,
    state = stateName,
    dfMeta = dfMeta
)

## Get cluster colors ##
FNcol <- paste0(hpc.mount, "projects/bonfantip/roberta.ragazzini/471_scRNAseq_newSamples_Epi6_Epi10_SC20193_SC21201/workdir/scRNAseq/design/clusterAnnotation.txt")

dfCol <- read.delim(
    FNcol, 
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

dfCol <- unique(dfCol[,c("clusterName", "clusterColor")])

dfCol[dfCol$clusterName %in% c("C4", "C8", "C11", "C13"),"clusterColor"] <- "#000000"

dfAdd$clusterColor <- NULL

dfAdd <- dplyr::full_join(
    dfAdd,
    dfCol,
    by = clusterName
)

dfReel <- rbind(
    dfReel, 
    dfAdd
)

## Done                                                                      ##
###############################################################################

###############################################################################
## Third Frame                                                               ##
## Change color for clusters to be removed to black
## Clusters marked for removal: C4, C8, C11, C13
stateName <- "Contaminating cells removed from A4A6A7A8"
orderVec <- c(orderVec, stateName)

dfMeta <- import.db.table.from.db(
    dbtable = "A4A6A7A8_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
)

dfAdd <- createAnimationFrame(
    cellID = "cellID",
    UMAP_1 = "UMAP_1",
    UMAP_2 = "UMAP_2",
    clusterName = "clusterName",
    clusterColor = NULL,
    state = stateName,
    dfMeta = dfMeta
)

## Get cluster colors ##
FNcol <- paste0(hpc.mount, "projects/bonfantip/roberta.ragazzini/471_scRNAseq_newSamples_Epi6_Epi10_SC20193_SC21201/workdir/scRNAseq/design/clusterAnnotation.txt")

dfCol <- read.delim(
    FNcol, 
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

dfCol <- unique(dfCol[,c("clusterName", "clusterColor")])

dfCol[dfCol$clusterName %in% c("C4", "C8", "C11", "C13"),"clusterColor"] <- "#000000"

dfAdd$clusterColor <- NULL

dfAdd <- dplyr::full_join(
    dfAdd,
    dfCol,
    by = clusterName
)

dfAdd <- dfAdd[!(dfAdd$clusterName %in% c("C4", "C8", "C11", "C13")), ]

dfReel <- rbind(
    dfReel, 
    dfAdd
)

## Done                                                                      ##
###############################################################################

###############################################################################
## Fortht Frame                                                              ##
## Change color for clusters to be removed to black
## Clusters marked for removal: C4, C8, C11, C13

stateName <- "New UMAP for Cell Selection A4A6A7A8sub"
orderVec <- c(orderVec, stateName)

dfMeta <- import.db.table.from.db(
    dbtable = "A4A6A7A8sub_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
)

dfAdd <- createAnimationFrame(
    cellID = "cellID",
    UMAP_1 = "UMAP_1",
    UMAP_2 = "UMAP_2",
    clusterName = "clusterName",
    clusterColor = NULL,
    state = stateName,
    dfMeta = dfMeta
)

## Get cluster colors ##
FNcol <- paste0(hpc.mount, "projects/bonfantip/roberta.ragazzini/471_subset_scRNAseq_newSamples_Epi6_Epi10_SC20193_SC21201/workdir/scRNAseq/design/clusterAnnotation.txt")

dfCol <- read.delim(
    FNcol, 
    header = T,
    sep = "\t",
    stringsAsFactors = F
)

dfCol <- dfCol[,c("clusterName", "clusterColor")]

dfAdd$clusterColor <- NULL

dfAdd <- dplyr::full_join(
    dfAdd,
    dfCol,
    by = clusterName
)

dfReel <- rbind(
    dfReel, 
    dfAdd
)

## Done                                                                      ##
###############################################################################

###############################################################################
## Fifth Frame                                                               ##
## Change color for clusters to be removed to black
## Clusters marked for removal: C4, C8, C11, C13

stateName <- "Cell Cycle Phases A4A6A7A8sub - No regression"
orderVec <- c(orderVec, stateName)

dfMeta <- import.db.table.from.db(
    dbtable = "A4A6A7A8sub_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
)

dfAdd <- createAnimationFrame(
    cellID = "cellID",
    UMAP_1 = "UMAP_1_Without_Regression",
    UMAP_2 = "UMAP_2_Without_Regression",
    clusterName = "Phase",
    clusterColor = NULL,
    state = stateName,
    dfMeta = dfMeta
)



dfReel <- rbind(
    dfReel, 
    dfAdd
)

## Done                                                                      ##
###############################################################################

###############################################################################
## Sixtht Frame                                                              ##
## Change color for clusters to be removed to black
## Clusters marked for removal: C4, C8, C11, C13

stateName <- "Cell Cycle Phases A4A6A7A8sub - Full Cell Cycle Regression"
orderVec <- c(orderVec, stateName)

dfMeta <- import.db.table.from.db(
    dbtable = "A4A6A7A8sub_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
)

dfAdd <- createAnimationFrame(
    cellID = "cellID",
    UMAP_1 = "UMAP_1_CellCycle_Regression",
    UMAP_2 = "UMAP_2_CellCycle_Regression",
    clusterName = "Phase",
    clusterColor = NULL,
    state = stateName,
    dfMeta = dfMeta
)



dfReel <- rbind(
    dfReel, 
    dfAdd
)

## Done                                                                      ##
###############################################################################


dfReel$state <- factor(dfReel$state, levels = orderVec)

#> Loading required package: ggplot2

# We'll start with a static plot

maxX <- 1.1 * max(dfReel$UMAP_1)
minX <- 1.1 * min(dfReel$UMAP_1)
maxY <- 1.1 * max(dfReel$UMAP_2)
minY <- 1.1 * min(dfReel$UMAP_2)

## Set dotsize on first frame
avgSize <- round(nrow(dfReel) / length(orderVec))
dotsize  = 1
if (avgSize  > 10000){
    dotsize  = 0.75
} else if (avgSize  > 20000){
    dotsize = 0.5
} else if (avgSize  > 50000){
    dotsize = 0.25
}

p <- ggplot(dfReel, aes(x = UMAP_1, y = UMAP_2, color = clusterColor)
    ) + geom_point(size = as.numeric(dotsize)) + ggplot2::theme_bw(
    ) +  ggplot2::theme(
        axis.text.y   = ggplot2::element_text(size=8),
        axis.text.x   = ggplot2::element_text(size=8),
        axis.title.y  = ggplot2::element_text(size=8),
        axis.title.x  = ggplot2::element_text(size=8),
        axis.line = ggplot2::element_line(colour = "black"),
        panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
        legend.title = ggplot2::element_blank()
    ) + scale_colour_identity() + labs(title = 'Stage: {closest_state}'
    ) + ggplot2::coord_fixed(ratio=1
    ) + xlim(c(minX, maxX)) + ylim(c(minY, maxY))

## Check
p2 <- p + facet_wrap( ~ state,  ncol = 2)

###############################################################################
## Add cluster labels in plot                                                ##
# library(dplyr)
# dfTextPos <- dfReel %>% group_by(clusterName, clusterColor) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)
# p <- p + ggrepel::geom_label_repel(data = dfTextPos,aes(label = clusterName))
##
###############################################################################


anim <- p + gganimate::transition_states(
    state,
    transition_length = 1,
    state_length = 5,
    wrap = FALSE
)

animate(anim, fps=5)


###############################################################################
## Save                                                                      ##
dir.create("test_anim")
animate(anim, nframes = 24, device = "png",
        renderer = file_renderer("test_anim", prefix = "A4A5_UMAP", overwrite = TRUE))


#library(av)
png1 <- paste0("test_anim/", list.files("test_anim"))

videoDir <- "/Volumes/babs/www/boeings/bioLOGIC_external/data/A4A6A7A8sub/html/report_videos/"

if (!dir.exists(videoDir)){
    dir.create(videoDir)
}

FNtsVideo <- paste0(videoDir, "processAnime.mp4")

png_files <- c(png1) #sprintf("input%03d.png", 1:10)
av::av_encode_video(png_files, FNtsVideo , framerate = 1)



## Done                                                                      ##
###############################################################################

