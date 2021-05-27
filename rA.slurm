

processing file: sc_PartA_QC.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |...                                                                   |   4%
  ordinary text without R code

  |                                                                              |......                                                                |   9%
label: setup (with options) 
List of 1
 $ include: logi FALSE

  |                                                                              |.........                                                             |  13%
  ordinary text without R code

  |                                                                              |............                                                          |  17%
label: unnamed-chunk-1 (with options) 
List of 2
 $ child: chr "src/modules/settings/set.parameters.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/settings/set.parameters.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |.......                                                               |  10%
  ordinary text without R code

  |                                                                              |..............                                                        |  20%
label: set_parameters_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |.....................                                                 |  30%
  ordinary text without R code

  |                                                                              |............................                                          |  40%
label: set_parameters_load_design (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |..........................................                            |  60%
label: set_parameters_create_Obio_object (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |.................................................                     |  70%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: set_parameters_create_design_table (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |...............................................................       |  90%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: render_sample_Table (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |...............                                                       |  22%
  ordinary text without R code

  |                                                                              |..................                                                    |  26%
label: unnamed-chunk-2 (with options) 
List of 2
 $ child: chr "src/modules/qc/create.10X.QC.table.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/qc/create.10X.QC.table.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |.....                                                                 |   8%
  ordinary text without R code

  |                                                                              |...........                                                           |  15%
label: set_QC_parameters_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |................                                                      |  23%
  ordinary text without R code

  |                                                                              |......................                                                |  31%
label: create_10X_QC_table (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |...........................                                           |  38%
   inline R code fragments

  |                                                                              |................................                                      |  46%
label: cellRangerCommands (with options) 
List of 4
 $ echo   : symbol F
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |......................................                                |  54%
   inline R code fragments

  |                                                                              |...........................................                           |  62%
label: render_QCTable (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |................................................                      |  69%
  ordinary text without R code

  |                                                                              |......................................................                |  77%
label: QC_raw_data (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

  |                                                                              |...........................................................           |  85%
   inline R code fragments

  |                                                                              |.................................................................     |  92%
label: Plot_tsne_data_plotting (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |.....................                                                 |  30%
  ordinary text without R code

  |                                                                              |........................                                              |  35%
label: unnamed-chunk-3 (with options) 
List of 2
 $ child: chr "src/modules/qc/A.create.sample.list.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/qc/A.create.sample.list.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..............                                                        |  20%
  ordinary text without R code

  |                                                                              |............................                                          |  40%
label: set_create_sample_list_parameters_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..........................................                            |  60%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: create_sample_list (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

Loading required package: DBI
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
15:30:21 UMAP embedding parameters a = 0.9922 b = 1.112
15:30:21 Read 2965 rows and found 30 numeric columns
15:30:21 Using Annoy for neighbor search, n_neighbors = 30
15:30:21 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:30:21 Writing NN index file to temp file /tmp/slurm_25471725/RtmpWs41h1/file1a416c75434f
15:30:21 Searching Annoy index using 1 thread, search_k = 3000
15:30:22 Annoy recall = 100%
15:30:23 Commencing smooth kNN distance calibration using 1 thread
15:30:24 Initializing from normalized Laplacian + noise
15:30:24 Commencing optimization for 500 epochs, with 118306 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:30:34 Optimization finished
Computing nearest neighbor graph
Computing SNN
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
15:30:53 UMAP embedding parameters a = 0.9922 b = 1.112
15:30:53 Read 2147 rows and found 30 numeric columns
15:30:53 Using Annoy for neighbor search, n_neighbors = 30
15:30:53 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:30:53 Writing NN index file to temp file /tmp/slurm_25471725/RtmpWs41h1/file1a4175f010b
15:30:53 Searching Annoy index using 1 thread, search_k = 3000
15:30:54 Annoy recall = 100%
15:30:55 Commencing smooth kNN distance calibration using 1 thread
15:30:56 Initializing from normalized Laplacian + noise
15:30:56 Commencing optimization for 500 epochs, with 84802 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:31:03 Optimization finished
Computing nearest neighbor graph
Computing SNN
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
15:31:35 UMAP embedding parameters a = 0.9922 b = 1.112
15:31:35 Read 3418 rows and found 30 numeric columns
15:31:35 Using Annoy for neighbor search, n_neighbors = 30
15:31:35 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:31:35 Writing NN index file to temp file /tmp/slurm_25471725/RtmpWs41h1/file1a4130e1972e
15:31:35 Searching Annoy index using 1 thread, search_k = 3000
15:31:36 Annoy recall = 100%
15:31:37 Commencing smooth kNN distance calibration using 1 thread
15:31:38 Initializing from normalized Laplacian + noise
15:31:38 Commencing optimization for 500 epochs, with 138060 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:31:49 Optimization finished
Computing nearest neighbor graph
Computing SNN
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
15:32:15 UMAP embedding parameters a = 0.9922 b = 1.112
15:32:15 Read 2730 rows and found 30 numeric columns
15:32:15 Using Annoy for neighbor search, n_neighbors = 30
15:32:15 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:32:15 Writing NN index file to temp file /tmp/slurm_25471725/RtmpWs41h1/file1a4155d0224d
15:32:15 Searching Annoy index using 1 thread, search_k = 3000
15:32:16 Annoy recall = 100%
15:32:16 Commencing smooth kNN distance calibration using 1 thread
15:32:18 Initializing from normalized Laplacian + noise
15:32:18 Commencing optimization for 500 epochs, with 104842 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:32:27 Optimization finished
Computing nearest neighbor graph
Computing SNN
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
15:32:40 UMAP embedding parameters a = 0.9922 b = 1.112
15:32:40 Read 1431 rows and found 30 numeric columns
15:32:40 Using Annoy for neighbor search, n_neighbors = 30
15:32:40 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:32:40 Writing NN index file to temp file /tmp/slurm_25471725/RtmpWs41h1/file1a4174a81ea8
15:32:40 Searching Annoy index using 1 thread, search_k = 3000
15:32:40 Annoy recall = 100%
15:32:41 Commencing smooth kNN distance calibration using 1 thread
15:32:42 Initializing from normalized Laplacian + noise
15:32:42 Commencing optimization for 500 epochs, with 54038 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:32:47 Optimization finished
Computing nearest neighbor graph
Computing SNN
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
15:33:11 UMAP embedding parameters a = 0.9922 b = 1.112
15:33:11 Read 2574 rows and found 30 numeric columns
15:33:11 Using Annoy for neighbor search, n_neighbors = 30
15:33:11 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:33:11 Writing NN index file to temp file /tmp/slurm_25471725/RtmpWs41h1/file1a4114e91ac0
15:33:11 Searching Annoy index using 1 thread, search_k = 3000
15:33:12 Annoy recall = 100%
15:33:12 Commencing smooth kNN distance calibration using 1 thread
15:33:14 Initializing from normalized Laplacian + noise
15:33:14 Commencing optimization for 500 epochs, with 100416 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:33:23 Optimization finished
Computing nearest neighbor graph
Computing SNN
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
15:33:26 UMAP embedding parameters a = 0.9922 b = 1.112
15:33:26 Read 345 rows and found 30 numeric columns
15:33:26 Using Annoy for neighbor search, n_neighbors = 30
15:33:26 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:33:26 Writing NN index file to temp file /tmp/slurm_25471725/RtmpWs41h1/file1a411c8b4090
15:33:26 Searching Annoy index using 1 thread, search_k = 3000
15:33:26 Annoy recall = 100%
15:33:27 Commencing smooth kNN distance calibration using 1 thread
15:33:28 Initializing from normalized Laplacian + noise
15:33:28 Commencing optimization for 500 epochs, with 13946 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:33:30 Optimization finished
Computing nearest neighbor graph
Computing SNN
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
15:34:09 UMAP embedding parameters a = 0.9922 b = 1.112
15:34:09 Read 4174 rows and found 30 numeric columns
15:34:09 Using Annoy for neighbor search, n_neighbors = 30
15:34:09 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:34:10 Writing NN index file to temp file /tmp/slurm_25471725/RtmpWs41h1/file1a412d0f34c9
15:34:10 Searching Annoy index using 1 thread, search_k = 3000
15:34:11 Annoy recall = 100%
15:34:11 Commencing smooth kNN distance calibration using 1 thread
15:34:12 Initializing from normalized Laplacian + noise
15:34:12 Commencing optimization for 500 epochs, with 169554 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:34:26 Optimization finished
Computing nearest neighbor graph
Computing SNN
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |...........................                                           |  39%
  ordinary text without R code

  |                                                                              |..............................                                        |  43%
label: unnamed-chunk-4 (with options) 
List of 2
 $ child: chr "src/modules/qc/A.create.sample.level.QC.plots.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/qc/A.create.sample.level.QC.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |.....                                                                 |   7%
  ordinary text without R code

  |                                                                              |.........                                                             |  13%
label: sample_level_QC_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..............                                                        |  20%
  ordinary text without R code

  |                                                                              |...................                                                   |  27%
label: percent_mt (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

mixtools package, version 1.2.0, Released 2020-02-05
This package is based upon work supported by the National Science Foundation under Grant No. SES-0518772.

  |                                                                              |.......................                                               |  33%
   inline R code fragments

  |                                                                              |............................                                          |  40%
label: Plot_percentMT_hist (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |.................................                                     |  47%
  ordinary text without R code

  |                                                                              |.....................................                                 |  53%
label: nFeature_data (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

  |                                                                              |..........................................                            |  60%
   inline R code fragments

  |                                                                              |...............................................                       |  67%
label: Plot_nFeatureRNA_hist (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |...................................................                   |  73%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: UMAP_data (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

  |                                                                              |.............................................................         |  87%
   inline R code fragments

  |                                                                              |.................................................................     |  93%
label: Plot_UMAP_SL (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |.................................                                     |  48%
  ordinary text without R code

  |                                                                              |.....................................                                 |  52%
label: unnamed-chunk-5 (with options) 
List of 2
 $ child: chr "src/modules/qc/A.create.sample.level.cell.cycle.QC.plots.Rmd"
 $ eval : logi FALSE

  |                                                                              |........................................                              |  57%
  ordinary text without R code

  |                                                                              |...........................................                           |  61%
label: unnamed-chunk-6 (with options) 
List of 2
 $ child: chr "src/modules/qc/A.perform.cell.cycle.regression.Rmd"
 $ eval : logi FALSE

  |                                                                              |..............................................                        |  65%
  ordinary text without R code

  |                                                                              |.................................................                     |  70%
label: unnamed-chunk-7 (with options) 
List of 2
 $ child: chr "src/modules/qc/A.create.integrated.cat.feature.plots.Rmd"
 $ eval : logi FALSE

  |                                                                              |....................................................                  |  74%
  ordinary text without R code

  |                                                                              |.......................................................               |  78%
label: saveobject (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol T
 $ results: symbol F

  |                                                                              |..........................................................            |  83%
  ordinary text without R code

  |                                                                              |.............................................................         |  87%
label: documentation (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol T
 $ results: symbol T

  |                                                                              |................................................................      |  91%
  ordinary text without R code

  |                                                                              |...................................................................   |  96%
label: create_report_params (with options) 
List of 2
 $ eval   : symbol T
 $ results: chr "asis"

Project Database could not be reached or has no entry in Obio@parameterList$lims.id for this analysis.
  |                                                                              |......................................................................| 100%
   inline R code fragments


output file: sc_PartA_QC.knit.md

/camp/apps/eb/dev/software/Anaconda2/2019.10/bin/pandoc +RTS -K512m -RTS sc_PartA_QC.utf8.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/hillc/danielle.park/439_scRNAseq_KPC_PDAC_GSE129455/workdir/html_local/sc_PartA_QC.html --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 5 --variable toc_float=1 --variable toc_selectors=h1,h2,h3,h4,h5 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --css src/assets/style/style.css --variable theme=paper --include-in-header /tmp/slurm_25471725/RtmpWs41h1/rmarkdown-str1a415531e0b.html --mathjax --variable 'mathjax-url:https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/hillc/danielle.park/439_scRNAseq_KPC_PDAC_GSE129455/workdir/html_local/sc_PartA_QC.html
There were 35 warnings (use warnings() to see them)
