

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
19:11:28 UMAP embedding parameters a = 0.9922 b = 1.112
19:11:28 Read 1226 rows and found 15 numeric columns
19:11:28 Using Annoy for neighbor search, n_neighbors = 30
19:11:28 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
19:11:28 Writing NN index file to temp file /tmp/Rtmp13dtpq/filedbcf23c6c257
19:11:28 Searching Annoy index using 1 thread, search_k = 3000
19:11:28 Annoy recall = 100%
19:11:29 Commencing smooth kNN distance calibration using 1 thread
19:11:31 Initializing from normalized Laplacian + noise
19:11:31 Commencing optimization for 500 epochs, with 47194 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
19:11:37 Optimization finished
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
19:12:05 UMAP embedding parameters a = 0.9922 b = 1.112
19:12:05 Read 2176 rows and found 15 numeric columns
19:12:05 Using Annoy for neighbor search, n_neighbors = 30
19:12:05 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
19:12:05 Writing NN index file to temp file /tmp/Rtmp13dtpq/filedbcf4bf40a64
19:12:05 Searching Annoy index using 1 thread, search_k = 3000
19:12:06 Annoy recall = 100%
19:12:06 Commencing smooth kNN distance calibration using 1 thread
19:12:08 Initializing from normalized Laplacian + noise
19:12:08 Commencing optimization for 500 epochs, with 85796 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
19:12:17 Optimization finished
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
19:12:32 UMAP embedding parameters a = 0.9922 b = 1.112
19:12:32 Read 1289 rows and found 15 numeric columns
19:12:32 Using Annoy for neighbor search, n_neighbors = 30
19:12:32 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
19:12:32 Writing NN index file to temp file /tmp/Rtmp13dtpq/filedbcf6e8d4518
19:12:32 Searching Annoy index using 1 thread, search_k = 3000
19:12:32 Annoy recall = 100%
19:12:33 Commencing smooth kNN distance calibration using 1 thread
19:12:35 Initializing from normalized Laplacian + noise
19:12:35 Commencing optimization for 500 epochs, with 54322 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
19:12:40 Optimization finished
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
19:13:08 UMAP embedding parameters a = 0.9922 b = 1.112
19:13:08 Read 2420 rows and found 15 numeric columns
19:13:08 Using Annoy for neighbor search, n_neighbors = 30
19:13:08 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
19:13:09 Writing NN index file to temp file /tmp/Rtmp13dtpq/filedbcf31e076d9
19:13:09 Searching Annoy index using 1 thread, search_k = 3000
19:13:10 Annoy recall = 100%
19:13:11 Commencing smooth kNN distance calibration using 1 thread
19:13:12 Initializing from normalized Laplacian + noise
19:13:12 Commencing optimization for 500 epochs, with 102908 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
19:13:22 Optimization finished
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

/camp/apps/eb/dev/software/Anaconda2/2019.10/bin/pandoc +RTS -K512m -RTS sc_PartA_QC.utf8.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/459_scRNAseq_pbl359sub_TcellclusterRemoval/workdir/html_local/sc_PartA_QC.html --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 5 --variable toc_float=1 --variable toc_selectors=h1,h2,h3,h4,h5 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --css src/assets/style/style.css --variable theme=paper --include-in-header /tmp/Rtmp13dtpq/rmarkdown-strdbcf6cb473f7.html --mathjax --variable 'mathjax-url:https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/459_scRNAseq_pbl359sub_TcellclusterRemoval/workdir/html_local/sc_PartA_QC.html
There were 35 warnings (use warnings() to see them)
