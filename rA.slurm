

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
20:59:37 UMAP embedding parameters a = 0.9922 b = 1.112
20:59:37 Read 50 rows and found 10 numeric columns
20:59:37 Using Annoy for neighbor search, n_neighbors = 30
20:59:37 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
20:59:38 Writing NN index file to temp file /tmp/slurm_25076982/RtmpImOPLj/file2dcf48f323
20:59:38 Searching Annoy index using 1 thread, search_k = 3000
20:59:38 Annoy recall = 100%
20:59:38 Commencing smooth kNN distance calibration using 1 thread
20:59:40 Initializing from normalized Laplacian + noise
20:59:40 Commencing optimization for 500 epochs, with 1534 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
20:59:41 Optimization finished
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
20:59:42 UMAP embedding parameters a = 0.9922 b = 1.112
20:59:42 Read 114 rows and found 10 numeric columns
20:59:42 Using Annoy for neighbor search, n_neighbors = 30
20:59:42 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
20:59:42 Writing NN index file to temp file /tmp/slurm_25076982/RtmpImOPLj/file2dcf3d1a16db
20:59:42 Searching Annoy index using 1 thread, search_k = 3000
20:59:43 Annoy recall = 100%
20:59:43 Commencing smooth kNN distance calibration using 1 thread
20:59:45 Initializing from normalized Laplacian + noise
20:59:45 Commencing optimization for 500 epochs, with 4348 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
20:59:46 Optimization finished
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
20:59:47 UMAP embedding parameters a = 0.9922 b = 1.112
20:59:47 Read 45 rows and found 10 numeric columns
20:59:47 Using Annoy for neighbor search, n_neighbors = 30
20:59:47 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
20:59:47 Writing NN index file to temp file /tmp/slurm_25076982/RtmpImOPLj/file2dcf7e637b16
20:59:47 Searching Annoy index using 1 thread, search_k = 3000
20:59:47 Annoy recall = 100%
20:59:47 Commencing smooth kNN distance calibration using 1 thread
20:59:49 Initializing from normalized Laplacian + noise
20:59:49 Commencing optimization for 500 epochs, with 1304 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
20:59:50 Optimization finished
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
20:59:55 UMAP embedding parameters a = 0.9922 b = 1.112
20:59:55 Read 260 rows and found 30 numeric columns
20:59:55 Using Annoy for neighbor search, n_neighbors = 30
20:59:55 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
20:59:55 Writing NN index file to temp file /tmp/slurm_25076982/RtmpImOPLj/file2dcf3b6d44aa
20:59:55 Searching Annoy index using 1 thread, search_k = 3000
20:59:55 Annoy recall = 100%
20:59:56 Commencing smooth kNN distance calibration using 1 thread
20:59:57 Initializing from normalized Laplacian + noise
20:59:57 Commencing optimization for 500 epochs, with 10600 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
20:59:59 Optimization finished
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
21:00:01 UMAP embedding parameters a = 0.9922 b = 1.112
21:00:01 Read 147 rows and found 10 numeric columns
21:00:01 Using Annoy for neighbor search, n_neighbors = 30
21:00:01 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:00:01 Writing NN index file to temp file /tmp/slurm_25076982/RtmpImOPLj/file2dcf5e5c1bb3
21:00:01 Searching Annoy index using 1 thread, search_k = 3000
21:00:01 Annoy recall = 100%
21:00:02 Commencing smooth kNN distance calibration using 1 thread
21:00:03 Initializing from normalized Laplacian + noise
21:00:03 Commencing optimization for 500 epochs, with 5600 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:00:04 Optimization finished
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
21:00:07 UMAP embedding parameters a = 0.9922 b = 1.112
21:00:07 Read 171 rows and found 10 numeric columns
21:00:07 Using Annoy for neighbor search, n_neighbors = 30
21:00:07 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:00:07 Writing NN index file to temp file /tmp/slurm_25076982/RtmpImOPLj/file2dcfe776404
21:00:07 Searching Annoy index using 1 thread, search_k = 3000
21:00:07 Annoy recall = 100%
21:00:08 Commencing smooth kNN distance calibration using 1 thread
21:00:09 Initializing from normalized Laplacian + noise
21:00:09 Commencing optimization for 500 epochs, with 6674 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:00:10 Optimization finished
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
21:00:13 UMAP embedding parameters a = 0.9922 b = 1.112
21:00:13 Read 217 rows and found 10 numeric columns
21:00:13 Using Annoy for neighbor search, n_neighbors = 30
21:00:13 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:00:13 Writing NN index file to temp file /tmp/slurm_25076982/RtmpImOPLj/file2dcf63d10e78
21:00:13 Searching Annoy index using 1 thread, search_k = 3000
21:00:13 Annoy recall = 100%
21:00:14 Commencing smooth kNN distance calibration using 1 thread
21:00:15 Initializing from normalized Laplacian + noise
21:00:15 Commencing optimization for 500 epochs, with 8308 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:00:17 Optimization finished
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
21:00:19 UMAP embedding parameters a = 0.9922 b = 1.112
21:00:19 Read 130 rows and found 10 numeric columns
21:00:19 Using Annoy for neighbor search, n_neighbors = 30
21:00:19 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:00:19 Writing NN index file to temp file /tmp/slurm_25076982/RtmpImOPLj/file2dcf49ecc99e
21:00:19 Searching Annoy index using 1 thread, search_k = 3000
21:00:19 Annoy recall = 100%
21:00:20 Commencing smooth kNN distance calibration using 1 thread
21:00:21 Initializing from normalized Laplacian + noise
21:00:21 Commencing optimization for 500 epochs, with 4432 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:00:22 Optimization finished
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

Computing 2000 integration features
Scaling features for provided objects
Finding all pairwise anchors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 250 anchors
Filtering anchors
	Retained 163 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 205 anchors
Filtering anchors
	Retained 203 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 220 anchors
Filtering anchors
	Retained 220 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 250 anchors
Filtering anchors
	Retained 100 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 563 anchors
Filtering anchors
	Retained 212 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 220 anchors
Filtering anchors
	Retained 101 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 250 anchors
Filtering anchors
	Retained 137 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 524 anchors
Filtering anchors
	Retained 274 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 220 anchors
Filtering anchors
	Retained 125 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 704 anchors
Filtering anchors
	Retained 344 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 250 anchors
Filtering anchors
	Retained 139 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 539 anchors
Filtering anchors
	Retained 259 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 220 anchors
Filtering anchors
	Retained 133 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 806 anchors
Filtering anchors
	Retained 398 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 643 anchors
Filtering anchors
	Retained 335 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 250 anchors
Filtering anchors
	Retained 111 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 548 anchors
Filtering anchors
	Retained 221 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 220 anchors
Filtering anchors
	Retained 104 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 923 anchors
Filtering anchors
	Retained 366 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 678 anchors
Filtering anchors
	Retained 261 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 749 anchors
Filtering anchors
	Retained 303 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 250 anchors
Filtering anchors
	Retained 148 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 492 anchors
Filtering anchors
	Retained 281 anchors
Extracting within-dataset neighbors
Running CCA
Warning in irlba(A = mat3, nv = num.cc) :
  You're computing too large a percentage of total singular values, use a standard svd instead.
Merging objects
Finding neighborhoods
Finding anchors
	Found 220 anchors
Filtering anchors
	Retained 141 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 637 anchors
Filtering anchors
	Retained 354 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 561 anchors
Filtering anchors
	Retained 318 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 582 anchors
Filtering anchors
	Retained 332 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 624 anchors
Filtering anchors
	Retained 362 anchors
Extracting within-dataset neighbors
Merging dataset 3 into 2
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 1 into 2 3
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 8 into 7
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 5 into 2 3 1
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 6 into 4
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 2 3 1 5 into 4 6
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 7 8 into 4 6 2 3 1 5
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Warning: Adding a command log without an assay associated with it
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:01:46 UMAP embedding parameters a = 0.9922 b = 1.112
21:01:46 Read 1133 rows and found 20 numeric columns
21:01:46 Using Annoy for neighbor search, n_neighbors = 30
21:01:46 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:01:47 Writing NN index file to temp file /tmp/slurm_25076982/RtmpImOPLj/file2dcf46ba8428
21:01:47 Searching Annoy index using 1 thread, search_k = 3000
21:01:47 Annoy recall = 100%
21:01:48 Commencing smooth kNN distance calibration using 1 thread
21:01:49 Initializing from normalized Laplacian + noise
21:01:49 Commencing optimization for 500 epochs, with 50826 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
21:01:54 Optimization finished
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

/camp/apps/eb/dev/software/Anaconda2/2019.10/bin/pandoc +RTS -K512m -RTS sc_PartA_QC.utf8.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/hillc/danielle.park/437subC_fibroblast_GSE154778_primTumor_C4/workdir/html_local/sc_PartA_QC.html --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 5 --variable toc_float=1 --variable toc_selectors=h1,h2,h3,h4,h5 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --css src/assets/style/style.css --variable theme=paper --include-in-header /tmp/slurm_25076982/RtmpImOPLj/rmarkdown-str2dcf72a91db7.html --mathjax --variable 'mathjax-url:https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/hillc/danielle.park/437subC_fibroblast_GSE154778_primTumor_C4/workdir/html_local/sc_PartA_QC.html
There were 28 warnings (use warnings() to see them)
