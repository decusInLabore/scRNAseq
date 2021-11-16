---
  
output: 
    html_document:
        highlight: default
        theme: paper
        code_folding: hide
        df_print: tibble
        toc: true
        toc_depth: 3
        toc_float: true
        css: /camp/stp/babs/working/boeings/Stefan/protocol_files/github/boeings/templates/style/style.css

---
  
  
  
```{r setup_mon, include=FALSE}
knitr::opts_chunk$set(
  tidy = TRUE,
  tidy.opts = list(width.cutoff = 120),
  message = FALSE,
  warning = FALSE
)
```


```{r hpc_notes_mon, include=FALSE}

## Get interactive session ##
#  srun --time=08:00:00 --mem=40G -p int --pty bash

# module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;R;

# sbatch --time=12:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;Rscript runA.r" --job-name="rA"  --mem=100G -o rA.slurm >> commands.txt

# sbatch --time=24:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;Rscript runDGE.r" --job-name="rDGE" -p hmem --mem=300G -o rDGE.slurm >> commands.txt

# --mem-per-cpu=14G -p hmem --pty bash

```

```{r populate_meta_data_database_Mon, eval=TRUE, echo=F, results=F}
#Create the environment and load a suitable version of R, e.g. so:
VersionPdfExt <- paste0(".V", gsub("-", "", Sys.Date()), ".pdf")

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

figureCount <- 1



devtools::install_github("decusInLabore/biologicSeqTools")
devtools::install_github("decusInLabore/biologicToolsSC")
library(Seurat)
library(biologicSeqTools)
library(biologicToolsSC)


# source("assets/R/SBwebtools.pckg.r")
# source("assets/R/scTools.r")
# 
# if (length(.libPaths()) > 2){
#     .libPaths(.libPaths()[2:3])
# }
## Create biologic Object for visualization ##

ObioFN <- paste0("../", list.files("..")[grep(".bioLOGIC.Robj", list.files(".."))])

load(ObioFN)

checkFile = paste0(
  Obio@parameterList$project_id,
  ".bioLOGIC.Robj"
)


Obio <- setMountingPoint(Obio)
Obio <- setAnalysisPaths(Obio)
Obio <- setCrickGenomeAndGeneNameTable(Obio)
Obio <- createAnalysisFolders(
  Obio
)
Obio <- setDataBaseParameters(Obio)

## Upload metadata table > p315_PCA
# Obio@parameterList$host <- "10.152.22.193"
# Obio@parameterList$db.user <- "boeingS"
# db.pwd <- "5+3f4nB040420"


ObioFN <- paste0("../", list.files("..")[grep(".Seurat.Robj", list.files(".."))])

load(ObioFN)

## Create url string
if (Obio@parameterList$host == "10.27.241.234"){
  urlString <- "biologic.thecrick.org"
} else {
  urlString <- "biologic.crick.ac.uk"
}

legendDotSize <- 5


shinyURL <- paste0(
  "https://shiny-bioinformatics.crick.ac.uk/shiny/boeings/",
  Obio@parameterList$project_id,
  "_app/"
)            


## Set file paths ##
baseFN <- paste0(
  Obio@parameterList$project_id, 
  ".trajectory.table.xlsx"
)


outPutFN <- paste0(
  Obio@parameterList$reportTableDir,
  baseFN
)


FNrel <- paste0("report_tables/", baseFN)


## Create table link string ##

tableLink <- paste0('<a href="https://',urlString,'/mdata/',Obio@parameterList$project_id, '/html/', FNrel,' target="_blank">here</a>')  

tableString <- paste0('An Excel table with the DGE results can be downloaded ',
                      tableLink
)
```

# ### Re-dimension reduction for 3D rendering
# 
# if (Dim = "3D") {
#   
#   print ("Running UMAP 3D")
#   
#   seurat <- RunUMAP(object = seurat, reduction = "pca", dims = 1:20, n.components = 3)
#   
#   print("Clustering 3D")
#   
#   seurat <- FindNeighbors(object=seurat, dims=1:20)
#   seurat <- FindClusters(object=seurat, resolution=0.5)
#   seurat[[sprintf("ClusterNames_%.1f_%dPC", 0.5, 20)]] <- Idents(object = seurat)
#   
# }
# 

### Building the necessary parts for a basic cds

# part one, gene annotations



seurat <- seuratobj


library(Seurat)
library(monocle3)
library(htmlwidgets)

seurat <- readRDS("WHEREVER/YOU/KEEP/YOURS")

gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information

cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]], row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix

New_matrix <- seurat@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix


### Construct the basic cds object

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


### Construct and assign the made up partition 
###### I DO NOT ADVISE THIS

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


### Assign the cluster info

list_cluster <- seurat@meta.data[[sprintf("seurat_clusters")]]
#list_cluster <- seurat@meta.data[[sprintf("ClusterNames_%s_%sPC", 0.5, 20)]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster


### Could be a space-holder, but essentially fills out louvain parameters

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


### Assign UMAP coordinate

cds_from_seurat@reducedDims@listData[["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings


### Assign feature loading for downstream module analysis

cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings


### Learn graph, this step usually takes a significant period of time for larger samples

print("Learning graph, which can take a while depends on the sample")

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)

plot_cells(cds_from_seurat)

####Here I chose to save the gene_metadata, rds, etc.`