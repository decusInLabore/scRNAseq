

processing file: sc_PartA_QC.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |...                                                                   |   5%
  ordinary text without R code

  |                                                                              |.......                                                               |  10%
label: setup (with options) 
List of 1
 $ include: logi FALSE

  |                                                                              |..........                                                            |  14%
  ordinary text without R code

  |                                                                              |.............                                                         |  19%
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


  |                                                                              |.................                                                     |  24%
  ordinary text without R code

  |                                                                              |....................                                                  |  29%
label: unnamed-chunk-2 (with options) 
List of 2
 $ child: chr "src/modules/qc/create.10X.QC.table.Rmd"
 $ eval : logi FALSE

  |                                                                              |.......................                                               |  33%
  ordinary text without R code

  |                                                                              |...........................                                           |  38%
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


Attaching package: 'Seurat'

The following object is masked from 'package:SummarizedExperiment':

    Assays

Warning: The 'show_progress' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)
Calculating cell attributes from input UMI matrix: log_umi
Variance stabilizing transformation of count matrix of size 18889 by 10185
Model formula is y ~ log_umi
Get Negative Binomial regression parameters per gene
Using 2000 genes, 10185 cells
There are 1 estimated thetas smaller than 1e-07 - will be set to 1e-07
Found 87 outliers - those will be ignored in fitting/regularization step

Second step: Get residuals using fitted parameters for 18889 genes
Computing corrected count matrix for 18889 genes
Calculating gene attributes
Wall clock passed: Time difference of 4.809521 mins
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Loading required package: DBI
Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
15:02:01 UMAP embedding parameters a = 0.9922 b = 1.112
15:02:01 Read 10185 rows and found 20 numeric columns
15:02:01 Using Annoy for neighbor search, n_neighbors = 30
15:02:01 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:02:04 Writing NN index file to temp file /tmp/slurm_25072569/RtmpubRsc8/file1091a6931b8
15:02:04 Searching Annoy index using 1 thread, search_k = 3000
15:02:08 Annoy recall = 100%
15:02:09 Commencing smooth kNN distance calibration using 1 thread
15:02:11 Initializing from normalized Laplacian + noise
15:02:11 Commencing optimization for 200 epochs, with 429476 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:02:28 Optimization finished
Computing nearest neighbor graph
Computing SNN
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |..............................                                        |  43%
  ordinary text without R code

  |                                                                              |.................................                                     |  48%
label: unnamed-chunk-4 (with options) 
List of 2
 $ child: chr "src/modules/qc/A.create.sample.level.QC.plots.Rmd"
 $ eval : logi FALSE

  |                                                                              |.....................................                                 |  52%
  ordinary text without R code

  |                                                                              |........................................                              |  57%
label: unnamed-chunk-5 (with options) 
List of 2
 $ child: chr "src/modules/qc/A.perform.cell.cycle.regression.Rmd"
 $ eval : logi FALSE

  |                                                                              |...........................................                           |  62%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: unnamed-chunk-6 (with options) 
List of 2
 $ child: chr "src/modules/qc/A.create.integrated.cat.feature.plots.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/qc/A.create.integrated.cat.feature.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |.........                                                             |  12%
  ordinary text without R code

  |                                                                              |..................                                                    |  25%
label: A_create_integrated_cat_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..........................                                            |  38%
  ordinary text without R code

  |                                                                              |...................................                                   |  50%
label: A_createIntegrated (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol T
 $ results: symbol F

Warning: The 'show_progress' argument is deprecated as of v0.3. Use 'verbosity' instead. (in sctransform::vst)
Calculating cell attributes from input UMI matrix: log_umi
Variance stabilizing transformation of count matrix of size 18766 by 9200
Model formula is y ~ log_umi
Get Negative Binomial regression parameters per gene
Using 2000 genes, 9200 cells
Found 74 outliers - those will be ignored in fitting/regularization step

Second step: Get residuals using fitted parameters for 18766 genes
Computing corrected count matrix for 18766 genes
Calculating gene attributes
Wall clock passed: Time difference of 4.21929 mins
Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:07:47 UMAP embedding parameters a = 0.9922 b = 1.112
15:07:47 Read 9200 rows and found 20 numeric columns
15:07:47 Using Annoy for neighbor search, n_neighbors = 30
15:07:47 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:07:49 Writing NN index file to temp file /tmp/slurm_25072569/RtmpubRsc8/file1091323ca63c
15:07:49 Searching Annoy index using 1 thread, search_k = 3000
15:07:53 Annoy recall = 100%
15:07:54 Commencing smooth kNN distance calibration using 1 thread
15:07:56 Initializing from normalized Laplacian + noise
15:07:56 Commencing optimization for 500 epochs, with 396282 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
15:08:33 Optimization finished
Computing nearest neighbor graph
Computing SNN
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
  |                                                                              |............................................                          |  62%
  ordinary text without R code

  |                                                                              |....................................................                  |  75%
label: AUC_prep_from_file (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

Quantiles for the number of genes detected by cell: 
(Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC).
  |                                                                              |.............................................................         |  88%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: A_enrichment_plots (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |..................................................                    |  71%
  ordinary text without R code

  |                                                                              |.....................................................                 |  76%
label: saveobject (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol T
 $ results: symbol F

  |                                                                              |.........................................................             |  81%
  ordinary text without R code

  |                                                                              |............................................................          |  86%
label: documentation (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol T
 $ results: symbol T

  |                                                                              |...............................................................       |  90%
  ordinary text without R code

  |                                                                              |...................................................................   |  95%
label: create_report_params (with options) 
List of 2
 $ eval   : symbol T
 $ results: chr "asis"

Project Database could not be reached or has no entry in Obio@parameterList$lims.id for this analysis.
  |                                                                              |......................................................................| 100%
   inline R code fragments


output file: sc_PartA_QC.knit.md

/camp/apps/eb/dev/software/Anaconda2/2019.10/bin/pandoc +RTS -K512m -RTS sc_PartA_QC.utf8.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/boeings/stefan.boeing/426_scRNAseq_PBMC_10k/workdir/html_local/sc_PartA_QC.html --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 5 --variable toc_float=1 --variable toc_selectors=h1,h2,h3,h4,h5 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --css src/assets/style/style.css --variable theme=paper --include-in-header /tmp/slurm_25072569/RtmpubRsc8/rmarkdown-str10914673a5a5.html --mathjax --variable 'mathjax-url:https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/boeings/stefan.boeing/426_scRNAseq_PBMC_10k/workdir/html_local/sc_PartA_QC.html
There were 50 or more warnings (use warnings() to see the first 50)
