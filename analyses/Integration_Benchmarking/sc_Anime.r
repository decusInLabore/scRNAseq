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

if (!requireNamespace("remotes")){
    install.packages("remotes")
}


if (!requireNamespace("renv")){
    remotes::install_github("rstudio/renv")
}

if (!file.exists("renv.lock")){
    renv::init()
} else {
    renv::restore(prompt = FALSE)
}


## Load data
FN <- "../../../../workdir/temp/integration.coordinates.txt"

df <- read.delim(
    FN, 
    sep = "\t",
    stringsAsFactors = F
)

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


dfMeta <- biologicSeqTools2::import.db.table.from.db(
    dbtable = "integration01_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
)
stateName <- "Before Integration"

#stateName <- "Contaminating cells to be removed are marked black"
orderVec <-  stateName

dfMeta$cellID <- gsub('(.*)_\\w+', '\\1', dfMeta$cellID)
dfMeta$cellID <- paste0(dfMeta$cellID, "_",dfMeta$orig_ident)

#dfMeta <- df[df$integration_method == "cca", ]

dfClust <- dfMeta[,c("cellID", "original_clusterName")]

df <- merge(
    df, 
    dfClust, 
    by.x = "cellID",
    by.y = "cellID"
)

## Get original colors for clusters
library(dplyr)

dfCol1 <- biologicSeqTools2::import.db.table.from.db(
    dbtable = "pbl359subCluster07_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
) %>% dplyr::select(clusterName, clusterColor) %>% dplyr::distinct()

dfCol2 <- biologicSeqTools2::import.db.table.from.db(
    dbtable = "SC22096subset_PCA",
    dbname = "pbl_data",
    user     = "boeingS",
    password = db.pwd,
    host     = "10.27.241.82"
) %>% dplyr::select(clusterName, clusterColor) %>% dplyr::distinct()

dfCol <- rbind(
    dfCol1, 
    dfCol2
)

dfCol$clusterName <- paste0("org_cluster_", dfCol$clusterName)

df <- dplyr::full_join(
    df, 
    dfCol, 
    by = c("original_clusterName" = "clusterName")
)

dfAdd <- createAnimationFrame(
    cellID = "cellID",
    UMAP_1 = "UMAP_1_original",
    UMAP_2 = "UMAP_2_original",
    clusterName = "original_clusterName",
    clusterColor = "clusterColor",
    state = stateName,
    dfMeta = dfMeta
)




## Get cluster colors ##
# FNcol <- paste0(hpc.mount, "projects/bonfantip/roberta.ragazzini/471_scRNAseq_newSamples_Epi6_Epi10_SC20193_SC21201/workdir/scRNAseq/design/clusterAnnotation.txt")
# 
# dfCol <- read.delim(
#     FNcol, 
#     header = T,
#     sep = "\t",
#     stringsAsFactors = F
# )
# 
# dfCol <- unique(dfCol[,c("clusterName", "clusterColor")])
# 
# dfAdd$clusterColor <- NULL
# 
# dfAdd <- dplyr::full_join(
#     dfAdd,
#     dfCol,
#     by = clusterName
# )



dfReel <- dfAdd

## Done                                                                      ##
###############################################################################



###############################################################################
## Second frame                                                               ##


# dfMeta <- import.db.table.from.db(
#     dbtable = "A4A6A7A8_PCA",
#     dbname = "pbl_data",
#     user     = "boeingS",
#     password = db.pwd,
#     host     = "10.27.241.82"
# )
stateName <- "CCA Integration"

#stateName <- "Contaminating cells to be removed are marked black"
orderVec <-  c(orderVec, stateName)

dfMeta <- df[df$integration_method == "cca", ]

dfAdd <- createAnimationFrame(
    cellID = "cellID",
    UMAP_1 = "UMAP_1",
    UMAP_2 = "UMAP_2",
    clusterName = "original_clusterName",
    clusterColor = "clusterColor",
    state = stateName,
    dfMeta = dfMeta
)




## Get cluster colors ##
# FNcol <- paste0(hpc.mount, "projects/bonfantip/roberta.ragazzini/471_scRNAseq_newSamples_Epi6_Epi10_SC20193_SC21201/workdir/scRNAseq/design/clusterAnnotation.txt")
# 
# dfCol <- read.delim(
#     FNcol, 
#     header = T,
#     sep = "\t",
#     stringsAsFactors = F
# )
# 
# dfCol <- unique(dfCol[,c("clusterName", "clusterColor")])
# 
# dfAdd$clusterColor <- NULL
# 
# dfAdd <- dplyr::full_join(
#     dfAdd,
#     dfCol,
#     by = clusterName
# )



dfReel <- rbind(
    dfReel, 
    dfAdd
)

## Done                                                                      ##
###############################################################################

###############################################################################
## Second Frame                                                              ##
## Change color for clusters to be removed to black
## Clusters marked for removal: C4, C8, C11, C13

# stateName <- "Contaminating cells to be removed are marked black"
# orderVec <- c(orderVec, stateName)

# dfMeta <- import.db.table.from.db(
#     dbtable = "A4A6A7A8_PCA",
#     dbname = "pbl_data",
#     user     = "boeingS",
#     password = db.pwd,
#     host     = "10.27.241.82"
# )


stateName <- "RPCA Integration"
orderVec <- c(orderVec, stateName)


dfMeta <- df[df$integration_method == "rpca", ]

dfAdd <- createAnimationFrame(
    cellID = "cellID",
    UMAP_1 = "UMAP_1",
    UMAP_2 = "UMAP_2",
    clusterName = "original_clusterName",
    clusterColor = "clusterColor",
    state = stateName,
    dfMeta = dfMeta
)

## Get cluster colors ##
# FNcol <- paste0(hpc.mount, "projects/bonfantip/roberta.ragazzini/471_scRNAseq_newSamples_Epi6_Epi10_SC20193_SC21201/workdir/scRNAseq/design/clusterAnnotation.txt")
# 
# dfCol <- read.delim(
#     FNcol, 
#     header = T,
#     sep = "\t",
#     stringsAsFactors = F
# )
# 
# dfCol <- unique(dfCol[,c("clusterName", "clusterColor")])
# 
# dfCol[dfCol$clusterName %in% c("C4", "C8", "C11", "C13"),"clusterColor"] <- "#000000"
# 
# dfAdd$clusterColor <- NULL
# 
# dfAdd <- dplyr::full_join(
#     dfAdd,
#     dfCol,
#     by = clusterName
# )

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

# stateName <- "Contaminating cells to be removed are marked black"
# orderVec <- c(orderVec, stateName)

# dfMeta <- import.db.table.from.db(
#     dbtable = "A4A6A7A8_PCA",
#     dbname = "pbl_data",
#     user     = "boeingS",
#     password = db.pwd,
#     host     = "10.27.241.82"
# )


stateName <- "STACAS Integration"
orderVec <- c(orderVec, stateName)


dfMeta <- df[df$integration_method == "stacas", ]

dfAdd <- createAnimationFrame(
    cellID = "cellID",
    UMAP_1 = "UMAP_1",
    UMAP_2 = "UMAP_2",
    clusterName = "original_clusterName",
    clusterColor = "clusterColor",
    state = stateName,
    dfMeta = dfMeta
)

## Get cluster colors ##
# FNcol <- paste0(hpc.mount, "projects/bonfantip/roberta.ragazzini/471_scRNAseq_newSamples_Epi6_Epi10_SC20193_SC21201/workdir/scRNAseq/design/clusterAnnotation.txt")
# 
# dfCol <- read.delim(
#     FNcol, 
#     header = T,
#     sep = "\t",
#     stringsAsFactors = F
# )
# 
# dfCol <- unique(dfCol[,c("clusterName", "clusterColor")])
# 
# dfCol[dfCol$clusterName %in% c("C4", "C8", "C11", "C13"),"clusterColor"] <- "#000000"
# 
# dfAdd$clusterColor <- NULL
# 
# dfAdd <- dplyr::full_join(
#     dfAdd,
#     dfCol,
#     by = clusterName
# )

dfReel <- rbind(
    dfReel, 
    dfAdd
)

## Done                                                                      ##
###############################################################################

dfReel$clusterColor <- NULL
dfReel <- merge(
    dfReel, 
    dfCol, 
    by.x = "clusterName",
    by.y = "clusterName"
)

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
    dotsize  = 0.25
} else if (avgSize  > 20000){
    dotsize = 0.1
} else if (avgSize  > 50000){
    dotsize = 0.05
}

library(ggplot2)

colVec <- dfCol$clusterColor
names(colVec) <- dfCol$clusterName
p <- ggplot(dfReel, aes(x = UMAP_1, y = UMAP_2, color = clusterName)
    ) + geom_point(size = as.numeric(dotsize)) + ggplot2::theme_bw(
    ) +  ggplot2::theme(
        axis.text.y   = ggplot2::element_text(size=8),
        axis.text.x   = ggplot2::element_text(size=8),
        axis.title.y  = ggplot2::element_text(size=8),
        axis.title.x  = ggplot2::element_text(size=8),
        axis.line = ggplot2::element_line(colour = "black"),
        panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
        legend.title = ggplot2::element_blank(),
        legend.position = "none"
    ) + scale_color_manual(values=colVec) + labs(title = '{closest_state}'
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

gganimate::animate(anim, fps=5)


###############################################################################
## Save                                                                      ##

animDir <- "../../../../workdir/test_anim/"
dir.create(animDir)

gganimate::animate(anim, nframes = 24, device = "png",
        renderer = gganimate::file_renderer(animDir, prefix = "Integrationtest_", overwrite = TRUE))


#library(av)
png1 <- paste0(animDir, list.files(animDir))

videoDir <- "/Volumes/babs/www/boeings/bioLOGIC_external/data/integration01/html/report_videos/"

if (!dir.exists(videoDir)){
    dir.create(videoDir)
}

FNtsVideo <- paste0(videoDir, "processAnime.mp4")

png_files <- c(png1) #sprintf("input%03d.png", 1:10)
av::av_encode_video(png_files, FNtsVideo , framerate = 1)



## Done                                                                      ##
###############################################################################

