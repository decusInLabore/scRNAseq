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

## Load first 

dfMeta <- import.db.table.from.db(
    dbtable = "A4A6A7A8_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
)

dfMeta[["seuratSampleID"]] <- sapply(dfMeta$cellID, function(x) unlist(strsplit(x, "_"))[2])
dfMeta[["mergeID"]] <- sapply(dfMeta$cellID, function(x) unlist(strsplit(x, "_"))[1])

dfMeta[["merge"]] <- paste0(dfMeta$sampleID, "_", dfMeta$mergeID)

#dfMeta <- dfMeta[,c("cellID","UMAP_1", "UMAP_2", "sampleName", "sampleColor", "clusterName", "clusterColor")]

dfMeta <- dfMeta[,c("cellID","UMAP_1", "UMAP_2", "clusterName", "clusterColor")]
dfMeta[["state"]] <- "Before Subclustering"

## Load second
dfMeta2 <- import.db.table.from.db(
    dbtable = "A4A6A7A8sub_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
)

dfMeta2 <- dfMeta2[,c("cellID","UMAP_1", "UMAP_2", "sampleName", "sampleColor", "clusterName", "clusterColor")]
dfMeta2 <- dfMeta2[,c("cellID","UMAP_1", "UMAP_2", "clusterName", "clusterColor")]
dfMeta2[["state"]] <- "New UMAP After Subclustering A4A6A7A8sub"





## Add transition state
dfMetaT1 <- dfMeta
dfMetaT1[["state"]] <- "Cells Marked for Removal in A4A6A7A8"
dfMetaT1$clusterName <- gsub("C4", "Remove",dfMetaT1$clusterName)
dfMetaT1$clusterName <- gsub("C8", "Remove",dfMetaT1$clusterName)
dfMetaT1$clusterName <- gsub("C11", "Remove",dfMetaT1$clusterName)
dfMetaT1$clusterName <- gsub("C13", "Remove",dfMetaT1$clusterName)
dfMetaT1[dfMetaT1$clusterName == "Remove", "clusterColor"] <- "#000000"



dfMetaT2 <- dfMetaT1
dfMetaT2[["state"]] <- "UMAP A4A6A7A8 After Cell Removal"
dfMetaT2 <- dfMetaT2[-grep("Remove", dfMetaT2$clusterName), ]

###############################################################################
## Add cell cycle states                                                     ##
dfMetaCS <- import.db.table.from.db(
    dbtable = "A4A6A7A8sub_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
)

dfWO <- dfMetaCS[dfMetaCS$UMAP_1_Without_Regression != 0 & dfMetaCS$UMAP_2_Without_Regression != 0 ,c("cellID", "UMAP_1_Without_Regression", "UMAP_2_Without_Regression", "ClusterNames_CellCycle_Regression",  "clusterColor")]

## Done                                                                      ##
###############################################################################


dfReel <- rbind(
    dfMeta, 
    dfMetaT1,
    dfMetaT2,
    dfMeta2
)

dfReel$state <- factor(dfReel$state, levels = c(
    "Before Subclustering", 
    "Cells Marked for Removal in A4A6A7A8", 
    "UMAP A4A6A7A8 After Cell Removal", 
    "New UMAP After Subclustering A4A6A7A8sub"
    )
)

library(gganimate)
#> Loading required package: ggplot2

# We'll start with a static plot

maxX <- 1.1 * max(dfReel$UMAP_1)
minX <- 1.1 * min(dfReel$UMAP_1)
maxY <- 1.1 * max(dfReel$UMAP_2)
minY <- 1.1 * min(dfReel$UMAP_2)

p <- ggplot(dfReel, aes(x = UMAP_1, y = UMAP_2, color = clusterColor)) +
    geom_point() + ggplot2::theme_bw(
    )  +  ggplot2::theme(
        axis.text.y   = ggplot2::element_text(size=8),
        axis.text.x   = ggplot2::element_text(size=8),
        axis.title.y  = ggplot2::element_text(size=8),
        axis.title.x  = ggplot2::element_text(size=8),
        axis.line = ggplot2::element_line(colour = "black"),
        panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
        legend.title = ggplot2::element_blank()
    ) + scale_colour_identity() + labs(title = 'Stage: {closest_state}') + ggplot2::coord_fixed(ratio=1
    )  + xlim(c(minX, maxX)) + ylim(c(minY, maxY))

## Check
p2 <- p + facet_wrap( ~ state,  ncol = 2)

###############################################################################
## Add cluster labels in plot                                                ##
library(dplyr)
dfTextPos <- dfReel %>% group_by(clusterName, clusterColor) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)
p <- p + ggrepel::geom_label_repel(data = dfTextPos,aes(label = clusterName))
##
###############################################################################


library(gganimate)
anim <- p +
    gganimate::transition_states(state,
                      transition_length = 1,
                      state_length = 5,
                      wrap = FALSE)

animate(anim, fps=5)


###############################################################################
## Save                                                                      ##
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

