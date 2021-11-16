

processing file: sc_PartB_Analysis.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |.                                                                     |   2%
  ordinary text without R code

  |                                                                              |...                                                                   |   4%
label: setup (with options) 
List of 1
 $ include: logi FALSE

  |                                                                              |....                                                                  |   6%
  ordinary text without R code

  |                                                                              |.....                                                                 |   8%
label: unnamed-chunk-1 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B.1.set.parameters.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B.1.set.parameters.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..............                                                        |  20%
  ordinary text without R code

  |                                                                              |............................                                          |  40%
label: set_parameters_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..........................................                            |  60%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: B_set_directories (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |.......                                                               |  10%
  ordinary text without R code

  |                                                                              |........                                                              |  12%
label: unnamed-chunk-2 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B.2.retrieve.reference.gene.sets.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B.2.retrieve.reference.gene.sets.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..............                                                        |  20%
  ordinary text without R code

  |                                                                              |............................                                          |  40%
label: partB_retrieve_reference_gene_sets_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..........................................                            |  60%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: partB_create_sample_list (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

Loading required package: DBI
  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |..........                                                            |  14%
  ordinary text without R code

  |                                                                              |...........                                                           |  16%
label: unnamed-chunk-3 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B.3.create.integrated.seurat.object.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B.3.create.integrated.seurat.object.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_create_integrated_Seurat_object_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: B3_partB_create_sample_list (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
Computing 2000 integration features
Scaling features for provided objects
Finding all pairwise anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 3718 anchors
Filtering anchors
	Retained 2321 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 3231 anchors
Filtering anchors
	Retained 1299 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 4620 anchors
Filtering anchors
	Retained 1491 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 4672 anchors
Filtering anchors
	Retained 1123 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7627 anchors
Filtering anchors
	Retained 1204 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 5035 anchors
Filtering anchors
	Retained 2016 anchors
Extracting within-dataset neighbors
Merging dataset 1 into 2
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 3 into 4
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 4 3 into 2 1
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Warning: Adding a command log without an assay associated with it
Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
  |                                                                              |..........................................................            |  83%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: partB_create_PCA_plots (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
22:27:33 UMAP embedding parameters a = 0.9922 b = 1.112
22:27:33 Read 8807 rows and found 30 numeric columns
22:27:33 Using Annoy for neighbor search, n_neighbors = 30
22:27:33 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
22:27:35 Writing NN index file to temp file /tmp/Rtmpk2IyZx/file11283161e5779
22:27:35 Searching Annoy index using 1 thread, search_k = 3000
22:27:38 Annoy recall = 100%
22:27:40 Commencing smooth kNN distance calibration using 1 thread
22:27:45 Initializing from normalized Laplacian + noise
22:27:46 Commencing optimization for 500 epochs, with 367252 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
22:28:22 Optimization finished
Computing nearest neighbor graph
Computing SNN
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|

Attaching package: 'scales'

The following object is masked from 'package:purrr':

    discard

The following object is masked from 'package:readr':

    col_factor

Calculating cluster Polykeratins
Calculating cluster Ionocytes
Calculating cluster NEDC_I_NEUROD1
Calculating cluster mTEC_Myo
Calculating cluster mTEC_Diff
Calculating cluster NEDC_III_SOX2
Calculating cluster HB_region
Calculating cluster NEDC_II_NEUROD1
Calculating cluster Polykeratin_prolif
Calculating cluster cTEC_II
Calculating cluster cTEC_Neurons_act
Calculating cluster T_cells_cont
Calculating cluster cTEC_diff
Calculating cluster NEDC_IV_SOX2
Calculating cluster cTEC_I
Calculating cluster cTEC_III

  |                                                                              |............                                                          |  18%
  ordinary text without R code

  |                                                                              |..............                                                        |  20%
label: unnamed-chunk-4 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B.4.create.variation.plots.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B.4.create.variation.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..........                                                            |  14%
  ordinary text without R code

  |                                                                              |....................                                                  |  29%
label: partB_create_variation_plots (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..............................                                        |  43%
  ordinary text without R code

  |                                                                              |........................................                              |  57%
label: partB_create-individual-varFeatures (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Finished averaging RNA for cluster all
Finished averaging integrated for cluster all
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Finished averaging RNA for cluster Epip6sub
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |=========                                                             |  12%  |                                                                              |============                                                          |  17%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  33%  |                                                                              |==========================                                            |  38%  |                                                                              |=============================                                         |  42%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |=========================================                             |  58%  |                                                                              |============================================                          |  62%  |                                                                              |===============================================                       |  67%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Finished averaging RNA for cluster Epip7sub
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |=====                                                                 |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  15%  |                                                                              |=============                                                         |  19%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |======================                                                |  31%  |                                                                              |========================                                              |  35%  |                                                                              |===========================                                           |  38%  |                                                                              |==============================                                        |  42%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  58%  |                                                                              |===========================================                           |  62%  |                                                                              |==============================================                        |  65%  |                                                                              |================================================                      |  69%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |=========================================================             |  81%  |                                                                              |===========================================================           |  85%  |                                                                              |==============================================================        |  88%  |                                                                              |=================================================================     |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Finished averaging RNA for cluster Epip10sub
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   5%  |                                                                              |======                                                                |   9%  |                                                                              |==========                                                            |  14%  |                                                                              |=============                                                         |  18%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |=============================                                         |  41%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  59%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |=========================================================             |  82%  |                                                                              |============================================================          |  86%  |                                                                              |================================================================      |  91%  |                                                                              |===================================================================   |  95%  |                                                                              |======================================================================| 100%
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Finished averaging RNA for cluster Epip14sub
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   5%  |                                                                              |======                                                                |   9%  |                                                                              |==========                                                            |  14%  |                                                                              |=============                                                         |  18%  |                                                                              |================                                                      |  23%  |                                                                              |===================                                                   |  27%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |=============================                                         |  41%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |=========================================                             |  59%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |===================================================                   |  73%  |                                                                              |======================================================                |  77%  |                                                                              |=========================================================             |  82%  |                                                                              |============================================================          |  86%  |                                                                              |================================================================      |  91%  |                                                                              |===================================================================   |  95%  |                                                                              |======================================================================| 100%
  |                                                                              |..................................................                    |  71%
   inline R code fragments

  |                                                                              |............................................................          |  86%
label: create-var-feature-plot-II (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |...............                                                       |  22%
  ordinary text without R code

  |                                                                              |................                                                      |  24%
label: unnamed-chunk-5 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B.6.create.pca.plots.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B.6.create.pca.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |.........                                                             |  12%
  ordinary text without R code

  |                                                                              |..................                                                    |  25%
label: partB_pcaB6_plots (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..........................                                            |  38%
  ordinary text without R code

  |                                                                              |...................................                                   |  50%
label: partB_upload_markers (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

  |                                                                              |............................................                          |  62%
  ordinary text without R code

  |                                                                              |....................................................                  |  75%
label: data_prep_sample_UMAP_rendering (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: logi FALSE

Registered S3 method overwritten by 'treeio':
  method     from
  root.phylo ape 
ggtree v1.16.6  For help: https://yulab-smu.github.io/treedata-book/

If you use ggtree in published research, please cite the most appropriate paper(s):

[36m-[39m Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods for mapping and visualizing associated data on phylogeny using ggtree. Molecular Biology and Evolution 2018, 35(12):3041-3043. doi: 10.1093/molbev/msy194
[36m-[39m Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam. ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution 2017, 8(1):28-36, doi:10.1111/2041-210X.12628



Attaching package: 'ggtree'

The following object is masked from 'package:tidyr':

    expand

Finished averaging RNA for cluster Epip6sub
Finished averaging RNA for cluster Epip7sub
Finished averaging RNA for cluster Epip10sub
Finished averaging RNA for cluster Epip14sub
  |                                                                              |.............................................................         |  88%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: Plot_tsne_data_plotting (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |..................                                                    |  25%
  ordinary text without R code

  |                                                                              |...................                                                   |  27%
label: unnamed-chunk-6 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B7.calculate.average.expression.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B7.calculate.average.expression.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..................                                                    |  25%
  ordinary text without R code

  |                                                                              |...................................                                   |  50%
label: partB_calculate_avg_expr (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |....................................................                  |  75%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: Identify differential expressed genes across conditions (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

Finished averaging RNA for cluster all
Finished averaging integrated for cluster all
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
Finished averaging RNA for cluster Polykeratins_Epip6sub
Finished averaging RNA for cluster Ionocytes_Epip6sub
Finished averaging RNA for cluster NEDC_I_NEUROD1_Epip6sub
Finished averaging RNA for cluster mTEC_Myo_Epip6sub
Finished averaging RNA for cluster mTEC_Diff_Epip6sub
Finished averaging RNA for cluster NEDC_III_SOX2_Epip6sub
Finished averaging RNA for cluster HB_region_Epip6sub
Finished averaging RNA for cluster NEDC_II_NEUROD1_Epip6sub
Finished averaging RNA for cluster Polykeratin_prolif_Epip6sub
Finished averaging RNA for cluster cTEC_II_Epip6sub
Finished averaging RNA for cluster cTEC_Neurons_act_Epip6sub
Finished averaging RNA for cluster T_cells_cont_Epip6sub
Finished averaging RNA for cluster cTEC_diff_Epip6sub
Finished averaging RNA for cluster NEDC_III_SOX2_Epip7sub
Finished averaging RNA for cluster NEDC_IV_SOX2_Epip7sub
Finished averaging RNA for cluster NEDC_II_NEUROD1_Epip7sub
Finished averaging RNA for cluster NEDC_I_NEUROD1_Epip7sub
Finished averaging RNA for cluster cTEC_Neurons_act_Epip7sub
Finished averaging RNA for cluster Polykeratins_Epip7sub
Finished averaging RNA for cluster Ionocytes_Epip7sub
Finished averaging RNA for cluster mTEC_Diff_Epip7sub
Finished averaging RNA for cluster mTEC_Myo_Epip7sub
Finished averaging RNA for cluster HB_region_Epip7sub
Finished averaging RNA for cluster cTEC_diff_Epip7sub
Finished averaging RNA for cluster Polykeratin_prolif_Epip7sub
Finished averaging RNA for cluster cTEC_I_Epip7sub
Finished averaging RNA for cluster cTEC_II_Epip7sub
Finished averaging RNA for cluster T_cells_cont_Epip7sub
Finished averaging RNA for cluster cTEC_III_Epip7sub
Finished averaging RNA for cluster cTEC_III_Epip10sub
Finished averaging RNA for cluster cTEC_I_Epip10sub
Finished averaging RNA for cluster Polykeratins_Epip10sub
Finished averaging RNA for cluster cTEC_Neurons_act_Epip10sub
Finished averaging RNA for cluster cTEC_II_Epip10sub
Finished averaging RNA for cluster cTEC_diff_Epip10sub
Finished averaging RNA for cluster T_cells_cont_Epip10sub
Finished averaging RNA for cluster NEDC_II_NEUROD1_Epip10sub
Finished averaging RNA for cluster NEDC_III_SOX2_Epip10sub
Finished averaging RNA for cluster Polykeratin_prolif_Epip10sub
Finished averaging RNA for cluster mTEC_Diff_Epip10sub
Finished averaging RNA for cluster HB_region_Epip10sub
Finished averaging RNA for cluster Ionocytes_Epip10sub
Finished averaging RNA for cluster NEDC_I_NEUROD1_Epip10sub
Finished averaging RNA for cluster NEDC_IV_SOX2_Epip10sub
Finished averaging RNA for cluster mTEC_Myo_Epip10sub
Finished averaging RNA for cluster cTEC_I_Epip14sub
Finished averaging RNA for cluster cTEC_II_Epip14sub
Finished averaging RNA for cluster Polykeratins_Epip14sub
Finished averaging RNA for cluster cTEC_Neurons_act_Epip14sub
Finished averaging RNA for cluster cTEC_diff_Epip14sub
Finished averaging RNA for cluster T_cells_cont_Epip14sub
Finished averaging RNA for cluster cTEC_III_Epip14sub
Finished averaging RNA for cluster NEDC_III_SOX2_Epip14sub
Finished averaging RNA for cluster Polykeratin_prolif_Epip14sub
Finished averaging RNA for cluster NEDC_II_NEUROD1_Epip14sub
Finished averaging RNA for cluster mTEC_Diff_Epip14sub
Finished averaging RNA for cluster HB_region_Epip14sub
Finished averaging RNA for cluster NEDC_I_NEUROD1_Epip14sub
Finished averaging RNA for cluster Ionocytes_Epip14sub
Finished averaging integrated for cluster Polykeratins_Epip6sub
Finished averaging integrated for cluster Ionocytes_Epip6sub
Finished averaging integrated for cluster NEDC_I_NEUROD1_Epip6sub
Finished averaging integrated for cluster mTEC_Myo_Epip6sub
Finished averaging integrated for cluster mTEC_Diff_Epip6sub
Finished averaging integrated for cluster NEDC_III_SOX2_Epip6sub
Finished averaging integrated for cluster HB_region_Epip6sub
Finished averaging integrated for cluster NEDC_II_NEUROD1_Epip6sub
Finished averaging integrated for cluster Polykeratin_prolif_Epip6sub
Finished averaging integrated for cluster cTEC_II_Epip6sub
Finished averaging integrated for cluster cTEC_Neurons_act_Epip6sub
Finished averaging integrated for cluster T_cells_cont_Epip6sub
Finished averaging integrated for cluster cTEC_diff_Epip6sub
Finished averaging integrated for cluster NEDC_III_SOX2_Epip7sub
Finished averaging integrated for cluster NEDC_IV_SOX2_Epip7sub
Finished averaging integrated for cluster NEDC_II_NEUROD1_Epip7sub
Finished averaging integrated for cluster NEDC_I_NEUROD1_Epip7sub
Finished averaging integrated for cluster cTEC_Neurons_act_Epip7sub
Finished averaging integrated for cluster Polykeratins_Epip7sub
Finished averaging integrated for cluster Ionocytes_Epip7sub
Finished averaging integrated for cluster mTEC_Diff_Epip7sub
Finished averaging integrated for cluster mTEC_Myo_Epip7sub
Finished averaging integrated for cluster HB_region_Epip7sub
Finished averaging integrated for cluster cTEC_diff_Epip7sub
Finished averaging integrated for cluster Polykeratin_prolif_Epip7sub
Finished averaging integrated for cluster cTEC_I_Epip7sub
Finished averaging integrated for cluster cTEC_II_Epip7sub
Finished averaging integrated for cluster T_cells_cont_Epip7sub
Finished averaging integrated for cluster cTEC_III_Epip7sub
Finished averaging integrated for cluster cTEC_III_Epip10sub
Finished averaging integrated for cluster cTEC_I_Epip10sub
Finished averaging integrated for cluster Polykeratins_Epip10sub
Finished averaging integrated for cluster cTEC_Neurons_act_Epip10sub
Finished averaging integrated for cluster cTEC_II_Epip10sub
Finished averaging integrated for cluster cTEC_diff_Epip10sub
Finished averaging integrated for cluster T_cells_cont_Epip10sub
Finished averaging integrated for cluster NEDC_II_NEUROD1_Epip10sub
Finished averaging integrated for cluster NEDC_III_SOX2_Epip10sub
Finished averaging integrated for cluster Polykeratin_prolif_Epip10sub
Finished averaging integrated for cluster mTEC_Diff_Epip10sub
Finished averaging integrated for cluster HB_region_Epip10sub
Finished averaging integrated for cluster Ionocytes_Epip10sub
Finished averaging integrated for cluster NEDC_I_NEUROD1_Epip10sub
Finished averaging integrated for cluster NEDC_IV_SOX2_Epip10sub
Finished averaging integrated for cluster mTEC_Myo_Epip10sub
Finished averaging integrated for cluster cTEC_I_Epip14sub
Finished averaging integrated for cluster cTEC_II_Epip14sub
Finished averaging integrated for cluster Polykeratins_Epip14sub
Finished averaging integrated for cluster cTEC_Neurons_act_Epip14sub
Finished averaging integrated for cluster cTEC_diff_Epip14sub
Finished averaging integrated for cluster T_cells_cont_Epip14sub
Finished averaging integrated for cluster cTEC_III_Epip14sub
Finished averaging integrated for cluster NEDC_III_SOX2_Epip14sub
Finished averaging integrated for cluster Polykeratin_prolif_Epip14sub
Finished averaging integrated for cluster NEDC_II_NEUROD1_Epip14sub
Finished averaging integrated for cluster mTEC_Diff_Epip14sub
Finished averaging integrated for cluster HB_region_Epip14sub
Finished averaging integrated for cluster NEDC_I_NEUROD1_Epip14sub
Finished averaging integrated for cluster Ionocytes_Epip14sub
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
Finished averaging RNA for cluster Polykeratins
Finished averaging RNA for cluster Ionocytes
Finished averaging RNA for cluster NEDC_I_NEUROD1
Finished averaging RNA for cluster mTEC_Myo
Finished averaging RNA for cluster mTEC_Diff
Finished averaging RNA for cluster NEDC_III_SOX2
Finished averaging RNA for cluster HB_region
Finished averaging RNA for cluster NEDC_II_NEUROD1
Finished averaging RNA for cluster Polykeratin_prolif
Finished averaging RNA for cluster cTEC_II
Finished averaging RNA for cluster cTEC_Neurons_act
Finished averaging RNA for cluster T_cells_cont
Finished averaging RNA for cluster cTEC_diff
Finished averaging RNA for cluster NEDC_IV_SOX2
Finished averaging RNA for cluster cTEC_I
Finished averaging RNA for cluster cTEC_III
Finished averaging integrated for cluster Polykeratins
Finished averaging integrated for cluster Ionocytes
Finished averaging integrated for cluster NEDC_I_NEUROD1
Finished averaging integrated for cluster mTEC_Myo
Finished averaging integrated for cluster mTEC_Diff
Finished averaging integrated for cluster NEDC_III_SOX2
Finished averaging integrated for cluster HB_region
Finished averaging integrated for cluster NEDC_II_NEUROD1
Finished averaging integrated for cluster Polykeratin_prolif
Finished averaging integrated for cluster cTEC_II
Finished averaging integrated for cluster cTEC_Neurons_act
Finished averaging integrated for cluster T_cells_cont
Finished averaging integrated for cluster cTEC_diff
Finished averaging integrated for cluster NEDC_IV_SOX2
Finished averaging integrated for cluster cTEC_I
Finished averaging integrated for cluster cTEC_III
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
Finished averaging RNA for cluster Epip6sub
Finished averaging RNA for cluster Epip7sub
Finished averaging RNA for cluster Epip10sub
Finished averaging RNA for cluster Epip14sub
Finished averaging integrated for cluster Epip6sub
Finished averaging integrated for cluster Epip7sub
Finished averaging integrated for cluster Epip10sub
Finished averaging integrated for cluster Epip14sub
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
Finished averaging RNA for cluster Epip6sub
Finished averaging RNA for cluster Epip7sub
Finished averaging RNA for cluster Epip10sub
Finished averaging RNA for cluster Epip14sub

  |                                                                              |.....................                                                 |  29%
  ordinary text without R code

  |                                                                              |......................                                                |  31%
label: unnamed-chunk-7 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B8.plot.pca.and.tsne.results.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B8.plot.pca.and.tsne.results.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_pca_tsne_plots (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: Plot_dimred_by_cluster_data_prep (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: logi FALSE

Finished averaging RNA for cluster Epip6sub
Finished averaging RNA for cluster Epip7sub
Finished averaging RNA for cluster Epip10sub
Finished averaging RNA for cluster Epip14sub
Finished averaging RNA for cluster Polykeratins
Finished averaging RNA for cluster Ionocytes
Finished averaging RNA for cluster NEDC_I_NEUROD1
Finished averaging RNA for cluster mTEC_Myo
Finished averaging RNA for cluster mTEC_Diff
Finished averaging RNA for cluster NEDC_III_SOX2
Finished averaging RNA for cluster HB_region
Finished averaging RNA for cluster NEDC_II_NEUROD1
Finished averaging RNA for cluster Polykeratin_prolif
Finished averaging RNA for cluster cTEC_II
Finished averaging RNA for cluster cTEC_Neurons_act
Finished averaging RNA for cluster T_cells_cont
Finished averaging RNA for cluster cTEC_diff
Finished averaging RNA for cluster NEDC_IV_SOX2
Finished averaging RNA for cluster cTEC_I
Finished averaging RNA for cluster cTEC_III
  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: Plot_dim_red_by_cluster (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |.......................                                               |  33%
  ordinary text without R code

  |                                                                              |.........................                                             |  35%
label: unnamed-chunk-8 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B9.euclidian.distances.figure.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B9.euclidian.distances.figure.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..............                                                        |  20%
  ordinary text without R code

  |                                                                              |............................                                          |  40%
label: partB_euclidian_distances (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..........................................                            |  60%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: Plot_umap_sleepwalk (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

Estimating 'maxdist' for feature matrix 1
  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |..........................                                            |  37%
  ordinary text without R code

  |                                                                              |...........................                                           |  39%
label: unnamed-chunk-9 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B10.doublet.estimation.plots.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B10.doublet.estimation.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_doublet_estimates (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: doublet-data (with options) 
List of 5
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: logi FALSE
 $ message: logi FALSE

  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: doublet-plot (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |.............................                                         |  41%
  ordinary text without R code

  |                                                                              |..............................                                        |  43%
label: unnamed-chunk-10 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B11.cell.cycle.phase.plots.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B11.cell.cycle.phase.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_cell_cycle_phase (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: Plot_dimred_by_cell_cycle_phase_data_prep (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: logi FALSE

Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
Finished averaging RNA for cluster Epip6sub
Finished averaging RNA for cluster Epip7sub
Finished averaging RNA for cluster Epip10sub
Finished averaging RNA for cluster Epip14sub
Finished averaging RNA for cluster Polykeratins
Finished averaging RNA for cluster Ionocytes
Finished averaging RNA for cluster NEDC_I_NEUROD1
Finished averaging RNA for cluster mTEC_Myo
Finished averaging RNA for cluster mTEC_Diff
Finished averaging RNA for cluster NEDC_III_SOX2
Finished averaging RNA for cluster HB_region
Finished averaging RNA for cluster NEDC_II_NEUROD1
Finished averaging RNA for cluster Polykeratin_prolif
Finished averaging RNA for cluster cTEC_II
Finished averaging RNA for cluster cTEC_Neurons_act
Finished averaging RNA for cluster T_cells_cont
Finished averaging RNA for cluster cTEC_diff
Finished averaging RNA for cluster NEDC_IV_SOX2
Finished averaging RNA for cluster cTEC_I
Finished averaging RNA for cluster cTEC_III
  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: Plot_dim_red_by_cell_cycle_phase (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |................................                                      |  45%
  ordinary text without R code

  |                                                                              |.................................                                     |  47%
label: unnamed-chunk-11 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B12.Ncell.barchart.plots.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B12.Ncell.barchart.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_N_cell_barchart (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: make_percent_N_cells_barchart (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: plot_barchart (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |..................................                                    |  49%
  ordinary text without R code

  |                                                                              |....................................                                  |  51%
label: unnamed-chunk-12 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B13.AUC.cat.enrichments.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B13.AUC.cat.enrichments.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_AUC_cat_enrichment (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: AUC_prep_from_file (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

Quantiles for the number of genes detected by cell: 
(Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC).
  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: enrichment_plots (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |.....................................                                 |  53%
  ordinary text without R code

  |                                                                              |......................................                                |  55%
label: unnamed-chunk-13 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B14.linear.regression.data.prep.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B14.linear.regression.data.prep.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_linear_regression (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: data_pre_linear_regression (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: plot_linear_regression (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |........................................                              |  57%
  ordinary text without R code

  |                                                                              |.........................................                             |  59%
label: unnamed-chunk-14 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B15.heatmaps.and.dotplots.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B15.heatmaps.and.dotplots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_heatmaps_and_dotplots (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: Cluster_overview_data_prep (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
Scale for 'fill' is already present. Adding another scale for 'fill', which
will replace the existing scale.
  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: Cluster_overview (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |...........................................                           |  61%
  ordinary text without R code

  |                                                                              |............................................                          |  63%
label: unnamed-chunk-15 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B16.make.cluster.table.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B16.make.cluster.table.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_make_cluster_table (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: ClusterTable_data_prep (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


Attaching package: 'DT'

The following object is masked from 'package:Seurat':

    JS

  |                                                                              |..........................................................            |  83%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: render_ClusterTable (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |.............................................                         |  65%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: unnamed-chunk-16 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B17.cell.type.estimate.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B17.cell.type.estimate.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_cell_type_estimate (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: CellTypePlot_dataprep (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F


Registered S3 method overwritten by 'enrichplot':
  method               from
  fortify.enrichResult DOSE
clusterProfiler v3.12.0  For help: https://guangchuangyu.github.io/software/clusterProfiler

If you use clusterProfiler in published research, please cite:
Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.

Attaching package: 'clusterProfiler'

The following object is masked from 'package:purrr':

    simplify

Loading required package: GSEABase
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: 'BiocGenerics'

The following objects are masked from 'package:parallel':

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from 'package:dplyr':

    combine, intersect, setdiff, union

The following objects are masked from 'package:stats':

    IQR, mad, sd, var, xtabs

The following objects are masked from 'package:base':

    anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    union, unique, unsplit, which, which.max, which.min

Loading required package: Biobase
Welcome to Bioconductor

    Vignettes contain introductory material; view with
    'browseVignettes()'. To cite Bioconductor, see
    'citation("Biobase")', and for packages 'citation("pkgname")'.

Loading required package: annotate
Loading required package: AnnotationDbi
Loading required package: stats4
Loading required package: IRanges
Loading required package: S4Vectors

Attaching package: 'S4Vectors'

The following object is masked from 'package:ggtree':

    expand

The following objects are masked from 'package:dplyr':

    first, rename

The following object is masked from 'package:tidyr':

    expand

The following object is masked from 'package:base':

    expand.grid


Attaching package: 'IRanges'

The following object is masked from 'package:ggtree':

    collapse

The following objects are masked from 'package:dplyr':

    collapse, desc, slice

The following object is masked from 'package:purrr':

    reduce


Attaching package: 'AnnotationDbi'

The following object is masked from 'package:dplyr':

    select

Loading required package: XML
Loading required package: graph

Attaching package: 'graph'

The following object is masked from 'package:XML':

    addNode

The following object is masked from 'package:stringr':

    boundary

--> No gene can be mapped....
--> Expected input gene ID: CCL22,TNFRSF4,ZBTB32,TREM2,MS4A1,SLC12A1
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: FAM198B,DAPK2,HK3,P2RY13,P2RY13,SIGLEC1
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: DPEP2,GZMB,FRK,STAP1,P2RY13,BIRC3
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: ALOX15,CCL14,P2RX5,STAP1,MMP12,RGS1
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: RGS13,ETV3,BLK,CD79B,FFAR2,BANK1
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: TNFSF14,VPREB3,CD19,ZNF222,CCL4,TMEM255A
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: IGHM,CD22,CCL23,TNFRSF17,CFP,PSG2
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: IL1A,PLA2G7,MYB,IGHD,IL18RAP,CPA3
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: IGKC,CCL13,CCL5,NMBR,ALOX15,IGSF6
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: PDCD1LG2,IGSF6,NOD2,TNFAIP6,NCF2,PDE6C
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: CD1C,CCL22,GPR18,FLVCR2,CCL8,SLC2A6
--> return NULL...
  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: CellTypePlot (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |................................................                      |  69%
  ordinary text without R code

  |                                                                              |.................................................                     |  71%
label: unnamed-chunk-17 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B18.gene.set.heatmaps.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B18.gene.set.heatmaps.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_gene_set_heatmaps (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: MakeHM (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |..........................................................            |  83%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: calculate-markdown-heatmap-chunks-dynamically (with options) 
List of 4
 $ echo   : symbol F
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |...................................................                   |  73%
  ordinary text without R code

  |                                                                              |....................................................                  |  75%
label: unnamed-chunk-18 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B19.example.feature.plots.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B19.example.feature.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_example_feature_plots (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: ExampleFeaturePlot_dataprep (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |..........................................................            |  83%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: ExampleFeaturePlot (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |......................................................                |  76%
  ordinary text without R code

  |                                                                              |.......................................................               |  78%
label: unnamed-chunk-19 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B20.pca.cell.and.loading.plots.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B20.pca.cell.and.loading.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |....                                                                  |   6%
  ordinary text without R code

  |                                                                              |.........                                                             |  12%
label: partB_pca_loadings (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |.............                                                         |  19%
  ordinary text without R code

  |                                                                              |..................                                                    |  25%
label: elbow_plot (with options) 
List of 3
 $ fig.height: num 6
 $ fig.width : num 10
 $ fig.cap   : language paste0("**Figure, ", figureCount, ": Elbowplot** Variance explained per PCA dimension")

  |                                                                              |......................                                                |  31%
  ordinary text without R code

  |                                                                              |..........................                                            |  38%
label: PCA_dataprep (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F


Attaching package: 'ggpubr'

The following object is masked from 'package:ggtree':

    rotate


Attaching package: 'gridExtra'

The following object is masked from 'package:Biobase':

    combine

The following object is masked from 'package:BiocGenerics':

    combine

The following object is masked from 'package:dplyr':

    combine

  |                                                                              |...............................                                       |  44%
   inline R code fragments

  |                                                                              |...................................                                   |  50%
label: plot_PCADimPlots (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |.......................................                               |  56%
  ordinary text without R code

  |                                                                              |............................................                          |  62%
label: pca_boxplot_dataprep (with options) 
List of 4
 $ echo   : symbol T
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

  |                                                                              |................................................                      |  69%
   inline R code fragments

  |                                                                              |....................................................                  |  75%
label: pca_boxplot (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |.........................................................             |  81%
  ordinary text without R code

  |                                                                              |.............................................................         |  88%
label: create-pca-enrichment-data (with options) 
List of 4
 $ echo   : symbol T
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

--> No gene can be mapped....
--> Expected input gene ID: CD72,TREM2,CCL17,IL3,CD33,ADAMDEC1
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: BCL7A,RAB27B,CXCL3,LRMP,CD37,CCL22
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: FCRL2,RGS13,DHRS11,IL2RA,FCER2,P2RY14
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: P2RY14,ADAMDEC1,FCGR2B,ABCB4,CRYBB1,CTSG
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: ZBTB32,CD38,FCER1A,SIGLEC1,EBI3,CCL18
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: CXCL3,BRAF,IL18R1,AQP9,ATP8B4,ADAMDEC1
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: CCL13,BCL7A,AIM2,CD86,RGS13,MAK
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: CCL8,RGS13,P2RX1,C3AR1,STAP1,CD38
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: MYB,MS4A3,AIF1,CLC,CLC,UGT1A8
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: PTGIR,CD80,IGKC,HPGDS,MEP1A,IFI44L
--> return NULL...
  |                                                                              |..................................................................    |  94%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: create-pca-enrichment-plot (with options) 
List of 4
 $ echo   : symbol T
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |........................................................              |  80%
  ordinary text without R code

  |                                                                              |..........................................................            |  82%
label: unnamed-chunk-20 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B21.pseudotime.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B21.pseudotime.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: partB_pseudotime (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: pseudotime-enrichment-data (with options) 
List of 4
 $ echo   : symbol T
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: pseudotime-plot (with options) 
List of 4
 $ echo   : symbol T
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |...........................................................           |  84%
  ordinary text without R code

  |                                                                              |............................................................          |  86%
label: unnamed-chunk-21 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B22.add.extra.tables.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B22.add.extra.tables.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..................                                                    |  25%
  ordinary text without R code

  |                                                                              |...................................                                   |  50%
label: partB_add_extra_tables (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |....................................................                  |  75%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: Additional tables (with options) 
List of 4
 $ echo   : symbol F
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F


  |                                                                              |..............................................................        |  88%
  ordinary text without R code

  |                                                                              |...............................................................       |  90%
label: saveobject (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol T
 $ results: symbol F

  |                                                                              |.................................................................     |  92%
  ordinary text without R code

  |                                                                              |..................................................................    |  94%
label: create_report_params (with options) 
List of 2
 $ eval   : symbol T
 $ results: chr "asis"

  |                                                                              |...................................................................   |  96%
  ordinary text without R code

  |                                                                              |..................................................................... |  98%
label: documentation (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol T
 $ results: symbol T

  |                                                                              |......................................................................| 100%
   inline R code fragments


output file: sc_PartB_Analysis.knit.md

/camp/apps/eb/dev/software/Anaconda2/2019.10/bin/pandoc +RTS -K512m -RTS sc_PartB_Analysis.utf8.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/446B_scRNAseq_subClustering_pbl359B_12clusters/workdir/html_local/sc_PartB_Analysis.html --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 5 --variable toc_float=1 --variable toc_selectors=h1,h2,h3,h4,h5 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --css src/assets/style/style.css --variable theme=paper --include-in-header /tmp/Rtmpk2IyZx/rmarkdown-str11283286e70ab.html --mathjax --variable 'mathjax-url:https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/446B_scRNAseq_subClustering_pbl359B_12clusters/workdir/html_local/sc_PartB_Analysis.html
There were 15 warnings (use warnings() to see them)
