###############################################################################
## Make Pseudotime Heatmaps                                                  ##

createLineageHeatmap <- function(
    OsC,  
    lineageSelection = NULL,
    pseudotimeLineageName = "Pseudotime_MC_Lineage_1",
    heatmapGeneVec = NULL,
    clusterRows = TRUE,
    plotTopNgenes = 100,
    plotAbsoluteValues = FALSE,  # otherwise: relative
    highColor = "#c4463a",
    midColor = "#fffbbc",
    lowColor = "#3060cf"
){
    
    if (!is.null(lineageSelection)){
        OsC_branch1 <- subset(x = OsC, subset = seurat_clusters %in% lineageSelection )
    } else {
        OsC_branch1 <- OsC
    }
    
    
    if (is.null(heatmapGeneVec)){
        Y <- data.matrix( OsC_branch1@assays[["RNA"]]@data)
        
        Y <-  Y[OsC_branch1@assays$integrated@var.features,]
        
        t <- OsC_branch1@meta.data[,pseudotimeLineageName]
        
        # Fit GAM for each gene using pseudotime as independent variable.
        
        gam.pval <- apply(Y, 1, function(z){
            d <- data.frame(z=z, t=t)
            tmp <- gam::gam(z ~ gam::lo(t), data=d)
            p <- summary(tmp)[4][[1]][1,5]
            p
        })
        topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:plotTopNgenes]
    } else {
        topgenes <- heatmapGeneVec
    }
    
    topgenes <- topgenes[topgenes %in% row.names(OsC_branch1)]
    
    # Identify genes with the most significant time-dependent model fit.
    
    
    
    
    # Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.
    #require(clusterExperiment)
    
    
    
    ## Get cells in pseudotime order ##
    dfOrder <- data.frame(OsC_branch1@meta.data)
    dfOrder <- dfOrder[,c("cellID", pseudotimeLineageName, "clusterName")]
    dfOrder <- dfOrder[order(dfOrder[,pseudotimeLineageName], decreasing = F),]
    cellOrder <- dfOrder$cellID
    cellCluster <- dfOrder$clusterName
    
    
    if (plotAbsoluteValues){
        dfDat <- OsC_branch1@assays$RNA@data
    } else {
        dfDat <- OsC_branch1@assays$RNA@scale.data    
    }
    
    dfDat <- dfDat[topgenes,cellOrder]
    
    PTcol <- names(dfOrder)[grep("Pseudotime", names(dfOrder))]
    dfOrder <- dfOrder[,c("cellID", PTcol, "clusterName")]
    dfOrder <- dfOrder[order(dfOrder[,PTcol], decreasing = F),]
    
    
    cellOrder <- as.character(dfOrder$cellID)
    cellCluster <- as.character(dfOrder$clusterName)
    PTorder <- dfOrder[,PTcol]
    names(PTorder) <- dfOrder$cellID
    names(dfOrder) <- gsub(PTcol, "Pseudotime", names(dfOrder))
    
    ## Make heatmap
    
    ## Average Pseudotime data into 100 bins ##
    mHmBase <- data.matrix(dfDat)
    
    if ( nrow(mHmBase) < 110){
        showRowNames <- TRUE
    } else {
        showRowNames <- FALSE
    }
    
    ## Create heatmap plot ##
    #library(ComplexHeatmap)
    #library(circlize)
    f1 = circlize::colorRamp2(seq(-4, 4, length = 3), c(lowColor, midColor, highColor))    
    
    if (plotAbsoluteValues){
        f1 = circlize::colorRamp2(seq(0, 6, length = 2), c(lowColor, highColor)) 
    } else {
        f1 = circlize::colorRamp2(seq(-4, 4, length = 3), c(lowColor, midColor, highColor))    
    }
    
    ## Create top annotation and colorbars ##
    # from https://www.biostars.org/p/368265/
    # from https://bioconductor.statistik.tu-dortmund.de/packages/3.1/bioc/vignettes/ComplexHeatmap/inst/doc/ComplexHeatmap.html
    
    
    anno <- as.data.frame(colnames(mHmBase))
    colnames(anno) <- "cellID"
    anno$Group <- cellCluster
    
    ## Color sample groups in line with the designated sample group color ##
    
    #library(scales)
    #hue_pal()(2)
    dfSel <- OsC_branch1@meta.data
    row.names(dfSel) <- dfSel$cellID
    dfSel <- dfSel[cellOrder, ]
    df <- unique(data.frame(OsC_branch1@meta.data[,c("clusterName", "clusterColor")]))
    
    GroupVec <- as.vector(df$clusterColor)
    names(GroupVec) <- as.vector(df$clusterName)
    
    
    df2 <- unique(dfSel[,c("cellID","clusterName", "clusterColor")])
    df2 <- na.omit(df2)
    row.names(df2) <- df2$cellID
    df2 <- df2[cellOrder, ]
    ClusterVec <- as.vector(df2[,c("clusterName")])
    df2 <- data.frame(df2[,c("clusterName")])
    names(df2) <- "Cluster"
    
    # df2 <- unique(data.frame(OsC_branch1@meta.data[,c("cellID","clusterName", "clusterColor")]))
    # row.names(df2) <- df2$cellID
    # df2 <- df2[cellOrder, ]
    # df2 <- data.frame(df2[,c("clusterName")])
    # names(df2) <- "Group"
    
    col_fun = circlize::colorRamp2(c(0,  ceiling(max(PTorder))), c("lightgrey", "black"))
    
    # ha = ComplexHeatmap::HeatmapAnnotation(
    #     df = df2, col = list(Group = GroupVec)
    # )
    
    ha = ComplexHeatmap::HeatmapAnnotation(
        Pseudotime = PTorder, 
        Cluster = ClusterVec,
        col = list(
            Pseudotime = col_fun,
            Cluster = GroupVec
        )
    )
    
    ComplexHeatmap::ht_opt(
        legend_border = "black",
        heatmap_border = TRUE,
        annotation_border = TRUE
    )
    
    h1 = ComplexHeatmap::Heatmap(
        mHmBase,
        column_title = gsub(
            "_", 
            " ", 
            paste0("Heatmap_", pseudotimeLineageName)
        ),
        name = pseudotimeLineageName, 
        #row_km = 5,
        col = f1,
        
        show_column_names = F,
        show_row_names = showRowNames,
        border = TRUE,
        cluster_columns = F,
        cluster_rows = clusterRows,
        
        #Dendrogram configurations: columns
        clustering_distance_columns="euclidean",
        clustering_method_columns="complete",
        column_dend_height=unit(10,"mm"),
        
        #Dendrogram configurations: rows
        clustering_distance_rows="euclidean",
        clustering_method_rows="complete",
        row_dend_width=grid::unit(10,"mm"),
        top_annotation = ha,
        show_heatmap_legend = TRUE,
        #row_title = NULL,
        #show_row_dend = FALSE,
        row_names_gp = grid::gpar(fontsize = 5)
    ) 
    
    ComplexHeatmap::ht_opt(RESET = TRUE)
    return(h1) 
}      

##  End of function                                                          ##
###############################################################################