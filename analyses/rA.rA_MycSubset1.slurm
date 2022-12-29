
The following have been reloaded with a version change:
  1) FFTW/3.3.4-gompi-2016b => FFTW/3.3.8-gompi-2020a
  2) GCC/5.4.0-2.26 => GCC/9.3.0
  3) GCCcore/5.4.0 => GCCcore/9.3.0
  4) OpenBLAS/0.2.18-GCC-5.4.0-2.26-LAPACK-3.6.1 => OpenBLAS/0.3.9-GCC-9.3.0
  5) OpenMPI/1.10.3-GCC-5.4.0-2.26 => OpenMPI/4.0.3-GCC-9.3.0
  6) ScaLAPACK/2.0.2-gompi-2016b-OpenBLAS-0.2.18-LAPACK-3.6.1 => ScaLAPACK/2.1.0-gompi-2020a
  7) binutils/2.26-GCCcore-5.4.0 => binutils/2.34-GCCcore-9.3.0
  8) foss/2016b => foss/2020a
  9) gompi/2016b => gompi/2020a
 10) hwloc/1.11.3-GCC-5.4.0-2.26 => hwloc/2.2.0-GCCcore-9.3.0
 11) numactl/2.0.11-GCC-5.4.0-2.26 => numactl/2.0.13-GCCcore-9.3.0



processing file: QC.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |...                                                                   |   4%
  ordinary text without R code

  |                                                                              |......                                                                |   8%
label: setup (with options) 
List of 1
 $ include: logi FALSE

  |                                                                              |........                                                              |  12%
  ordinary text without R code

  |                                                                              |...........                                                           |  16%
label: unnamed-chunk-1 (with options) 
List of 2
 $ child: chr "01-set.parameters.Rmd"
 $ eval : logi TRUE



processing file: ./01-set.parameters.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |......                                                                |   8%
  ordinary text without R code

  |                                                                              |............                                                          |  17%
label: set_parameters_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..................                                                    |  25%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: set_parameters_load_design (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |.............................                                         |  42%
  ordinary text without R code

  |                                                                              |...................................                                   |  50%
label: set_parameters_create_Obio_object (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |.........................................                             |  58%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: set_parameters_create_design_table (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |....................................................                  |  75%
  ordinary text without R code

  |                                                                              |..........................................................            |  83%
label: render_sample_Table (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |................................................................      |  92%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: parameter_table (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |..............                                                        |  20%
  ordinary text without R code

  |                                                                              |.................                                                     |  24%
label: unnamed-chunk-2 (with options) 
List of 2
 $ child: chr "02-retrieve.reference.gene.sets.Rmd"
 $ eval : logi TRUE



processing file: ./02-retrieve.reference.gene.sets.Rmd
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
 $ results: logi FALSE

Loading required package: DBI
  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |....................                                                  |  28%
  ordinary text without R code

  |                                                                              |......................                                                |  32%
label: unnamed-chunk-3 (with options) 
List of 2
 $ child: chr "03-create.10X.QC.table.Rmd"
 $ eval : logi TRUE



processing file: ./03-create.10X.QC.table.Rmd
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
 $ echo   : logi FALSE
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
 $ results: logi FALSE

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


  |                                                                              |.........................                                             |  36%
  ordinary text without R code

  |                                                                              |............................                                          |  40%
label: unnamed-chunk-4 (with options) 
List of 2
 $ child: chr "04-create.sample.list.Rmd"
 $ eval : logi TRUE



processing file: ./04-create.sample.list.Rmd
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

Registered S3 method overwritten by 'spatstat.geom':
  method     from
  print.boxx cli 
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
This message will be shown once per session
10:11:27 UMAP embedding parameters a = 0.9922 b = 1.112
10:11:27 Read 3296 rows and found 30 numeric columns
10:11:27 Using Annoy for neighbor search, n_neighbors = 30
10:11:27 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:11:28 Writing NN index file to temp file /tmp/slurm_47112575/RtmpxgvIpe/file5cb373184efe
10:11:28 Searching Annoy index using 1 thread, search_k = 3000
10:11:29 Annoy recall = 100%
10:11:30 Commencing smooth kNN distance calibration using 1 thread
10:11:32 Initializing from normalized Laplacian + noise
10:11:32 Commencing optimization for 500 epochs, with 144966 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:11:45 Optimization finished
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
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
10:11:59 UMAP embedding parameters a = 0.9922 b = 1.112
10:11:59 Read 642 rows and found 30 numeric columns
10:11:59 Using Annoy for neighbor search, n_neighbors = 30
10:11:59 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:11:59 Writing NN index file to temp file /tmp/slurm_47112575/RtmpxgvIpe/file5cb35c219475
10:11:59 Searching Annoy index using 1 thread, search_k = 3000
10:11:59 Annoy recall = 100%
10:12:00 Commencing smooth kNN distance calibration using 1 thread
10:12:02 Initializing from normalized Laplacian + noise
10:12:02 Commencing optimization for 500 epochs, with 24924 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:12:05 Optimization finished
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
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
10:12:18 UMAP embedding parameters a = 0.9922 b = 1.112
10:12:18 Read 655 rows and found 30 numeric columns
10:12:18 Using Annoy for neighbor search, n_neighbors = 30
10:12:18 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:12:18 Writing NN index file to temp file /tmp/slurm_47112575/RtmpxgvIpe/file5cb32f4f6988
10:12:18 Searching Annoy index using 1 thread, search_k = 3000
10:12:18 Annoy recall = 100%
10:12:19 Commencing smooth kNN distance calibration using 1 thread
10:12:21 Initializing from normalized Laplacian + noise
10:12:21 Commencing optimization for 500 epochs, with 26044 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:12:24 Optimization finished
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
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
10:13:45 UMAP embedding parameters a = 0.9922 b = 1.112
10:13:48 Read 4915 rows and found 30 numeric columns
10:13:48 Using Annoy for neighbor search, n_neighbors = 30
10:13:48 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:13:48 Writing NN index file to temp file /tmp/slurm_47112575/RtmpxgvIpe/file5cb31be38544
10:13:48 Searching Annoy index using 1 thread, search_k = 3000
10:13:50 Annoy recall = 100%
10:13:51 Commencing smooth kNN distance calibration using 1 thread
10:13:53 Initializing from normalized Laplacian + noise
10:13:53 Commencing optimization for 500 epochs, with 211618 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:14:12 Optimization finished
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
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
10:15:15 UMAP embedding parameters a = 0.9922 b = 1.112
10:15:15 Read 3172 rows and found 30 numeric columns
10:15:15 Using Annoy for neighbor search, n_neighbors = 30
10:15:15 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:15:16 Writing NN index file to temp file /tmp/slurm_47112575/RtmpxgvIpe/file5cb32c8ff213
10:15:16 Searching Annoy index using 1 thread, search_k = 3000
10:15:17 Annoy recall = 100%
10:15:18 Commencing smooth kNN distance calibration using 1 thread
10:15:20 Initializing from normalized Laplacian + noise
10:15:20 Commencing optimization for 500 epochs, with 137878 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:15:32 Optimization finished
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
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
10:16:13 UMAP embedding parameters a = 0.9922 b = 1.112
10:16:13 Read 2639 rows and found 30 numeric columns
10:16:13 Using Annoy for neighbor search, n_neighbors = 30
10:16:13 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:16:13 Writing NN index file to temp file /tmp/slurm_47112575/RtmpxgvIpe/file5cb36401d4c7
10:16:13 Searching Annoy index using 1 thread, search_k = 3000
10:16:14 Annoy recall = 100%
10:16:15 Commencing smooth kNN distance calibration using 1 thread
10:16:17 Initializing from normalized Laplacian + noise
10:16:17 Commencing optimization for 500 epochs, with 114598 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:16:28 Optimization finished
Computing nearest neighbor graph
Computing SNN
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |...............................                                       |  44%
  ordinary text without R code

  |                                                                              |..................................                                    |  48%
label: unnamed-chunk-5 (with options) 
List of 2
 $ child: chr "05-create.sample.level.QC.plots.Rmd"
 $ eval : logi TRUE



processing file: ./05-create.sample.level.QC.plots.Rmd
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
 $ results: logi FALSE

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
 $ results: logi FALSE

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
 $ results: logi FALSE

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


  |                                                                              |....................................                                  |  52%
  ordinary text without R code

  |                                                                              |.......................................                               |  56%
label: unnamed-chunk-6 (with options) 
List of 2
 $ child: chr "06-create.sample.level.cell.cycle.QC.plots.Rmd"
 $ eval : logi TRUE



processing file: ./06-create.sample.level.cell.cycle.QC.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..........                                                            |  14%
  ordinary text without R code

  |                                                                              |....................                                                  |  29%
label: sample_level_QC__CC_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..............................                                        |  43%
  ordinary text without R code

  |                                                                              |........................................                              |  57%
label: UMAP_data_CC (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: logi FALSE

Attaching SeuratObject
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
  |                                                                              |..................................................                    |  71%
   inline R code fragments

  |                                                                              |............................................................          |  86%
label: Plot_UMAP_CC (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |..........................................                            |  60%
  ordinary text without R code

  |                                                                              |.............................................                         |  64%
label: unnamed-chunk-7 (with options) 
List of 2
 $ child: chr "07-run.doublet.finder.Rmd"
 $ eval : logi TRUE



processing file: ./07-run.doublet.finder.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |.....                                                                 |   7%
  ordinary text without R code

  |                                                                              |.........                                                             |  13%
label: doublet_finder_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..............                                                        |  20%
  ordinary text without R code

  |                                                                              |...................                                                   |  27%
label: DF_data (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: logi FALSE

Loading required package: fields
Loading required package: spam
Loading required package: dotCall64
Loading required package: grid

Attaching package: 'grid'

The following object is masked from 'package:mixtools':

    depth

Spam version 2.5-1 (2019-12-12) is loaded.
Type 'help( Spam)' or 'demo( spam)' for a short introduction 
and overview of this package.
Help for individual functions is also obtained by adding the
suffix '.spam' to the function name, e.g. 'help( chol.spam)'.

Attaching package: 'spam'

The following objects are masked from 'package:base':

    backsolve, forwardsolve

Loading required package: maps
See https://github.com/NCAR/Fields for
 an extensive vignette, other supplements and source code 
Loading required package: KernSmooth
KernSmooth 2.23 loaded
Copyright M. P. Wand 1997-2009
Loading required package: sctransform
Calculating cell attributes from input UMI matrix: log_umi
Variance stabilizing transformation of count matrix of size 16834 by 4395
Model formula is y ~ log_umi
Get Negative Binomial regression parameters per gene
Using 2000 genes, 4395 cells
There are 1 estimated thetas smaller than 1e-07 - will be set to 1e-07
Found 113 outliers - those will be ignored in fitting/regularization step

Second step: Get residuals using fitted parameters for 16834 genes
Computing corrected count matrix for 16834 genes
Calculating gene attributes
Wall clock passed: Time difference of 2.071194 mins
Determine variable features
Place corrected count matrix in counts slot
Centering data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
Set default assay to SCT
PC_ 1 
Positive:  Ccl5, S100a4, Cd52, Lgals1, Nkg7, Ms4a4b, Tmsb10, Trbc1, Il7r, Trbc2 
	   Cd3g, Rgs1, Vps37b, Ptprc, Cd8b1, Cxcr6, Vim, Junb, Coro1a, AW112010 
	   Srgn, Cd3d, Arhgdib, Ltb, Zfp36l2, Rac2, Icos, Hcst, S100a10, Tnfaip3 
Negative:  Sftpa1, Sftpb, Sftpc, Sftpd, Wfdc2, Cxcl15, Cd74, Lyz2, Npc2, Scd1 
	   Slc34a2, Napsa, Chil1, Ctsh, Hp, Mgst1, Atp1b1, Cbr2, Lcn2, Cldn18 
	   Cldn3, S100g, Krt18, Retnla, Aqp5, Ager, Alcam, Hc, Lpcat1, Sfta2 
PC_ 2 
Positive:  S100a4, S100a6, Ccl5, Lgals1, Sftpb, Wfdc2, Sftpd, Cd52, Nkg7, Lgals3 
	   Sftpa1, Cxcr6, Cxcl15, Uba52, Ctsh, Cldn3, Gm42418, Tmsb4x, Areg, Rgs1 
	   Aqp5, Krt18, Napsa, Alcam, Mgst1, Krt8, Tmsb10, Crip1, Cyp2f2, Zfp36l2 
Negative:  Cldn5, Ramp2, Hpgd, Ptprb, Tmem100, Calcrl, Cdh5, Egfl7, Pltp, Clec14a 
	   Gpihbp1, Cavin2, Ly6a, Acvrl1, Tspan7, Tm4sf1, Aqp1, Ifitm3, Cd93, BC028528 
	   Clec1a, Ier3, Ace, Icam2, Ehd4, Pecam1, Cav1, Cxcl12, Sema3c, Glp1r 
PC_ 3 
Positive:  Clu, Scgb1a1, Krt19, Scgb3a1, S100a6, Krt18, Krt7, Ltf, Krt8, Anxa1 
	   Cxcl17, Cyp2f2, Cp, Tspan1, Areg, Epcam, Rbp4, Gsto1, Scgb3a2, Basp1 
	   F3, Cd24a, Ffar4, Fxyd3, Cd9, Wfdc2, Selenbp1, Lgals3, Gabrp, Adam28 
Negative:  Sftpc, Lyz2, Cxcl15, Slc34a2, Scd1, Napsa, Ccl5, Chil1, Npc2, Hc 
	   Lamp3, Elovl1, Lpcat1, Lcn2, Sftpa1, S100g, Cd74, Rnase4, Sfta2, Il33 
	   Lgi3, Ppp1r14c, Fasn, Chia1, Ly6c2, Abca3, Cox6a2, Ager, Egfl6, Fabp5 
PC_ 4 
Positive:  Ccl5, Nkg7, Ly6c2, Ms4a4b, Cd8b1, Klf2, Klrd1, Scgb1a1, Gzma, Cd8a 
	   Gimap4, Gramd3, Klre1, Epsti1, Rps24, Scgb3a2, AW112010, Xcl1, Rpl13, Ms4a6b 
	   Gimap7, Rps16, Gm2682, Itga4, Ctsw, Sell, Bcl2, Klrc1, Fcer1g, Prf1 
Negative:  S100a4, Gadd45b, Junb, Lmo4, Ramp1, Pxdc1, Cxcr6, Maf, Furin, Ccr2 
	   Ltb4r1, Nfkbia, Il17a, Serpinb1a, Tmem176a, Tcrg-C1, Icos, Bcl2a1b, Blk, Ramp3 
	   Pim1, Cd163l1, Lgals3, Hk2, Ckb, Nr4a1, Il2ra, Il7r, Bcl2a1d, Actn2 
PC_ 5 
Positive:  Ccl5, S100a6, Nkg7, Cd8b1, Lgals1, S100a4, Cldn18, Ctsh, Areg, Egr1 
	   Cd8a, Aqp5, Jun, Cd52, Vim, Tspan8, Actb, Phlda1, Krt18, S100a14 
	   Napsa, Ager, Lgals3, Atf3, Lmo7, Crip1, Fos, Ms4a4b, Krt8, Lpcat1 
Negative:  Scgb1a1, Retnla, Scgb3a2, Cbr2, Scgb3a1, Hp, Trf, Cldn10, Cp, Lyz2 
	   Selenbp1, Cyp2f2, B430010I23Rik, Pigr, Ldhb, mt-Co3, Ftl1, Pon1, Nupr1, mt-Atp6 
	   Fabp5, Lypd2, mt-Nd2, mt-Nd4, Fth1, Npc2, Bpifb1, Tff2, Gsta3, Inmt 
Calculating cell attributes from input UMI matrix: log_umi
Variance stabilizing transformation of count matrix of size 14870 by 856
Model formula is y ~ log_umi
Get Negative Binomial regression parameters per gene
Using 2000 genes, 856 cells
There are 1 estimated thetas smaller than 1e-07 - will be set to 1e-07
Found 81 outliers - those will be ignored in fitting/regularization step

Second step: Get residuals using fitted parameters for 14870 genes
Computing corrected count matrix for 14870 genes
Calculating gene attributes
Wall clock passed: Time difference of 25.88137 secs
Determine variable features
Place corrected count matrix in counts slot
Centering data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
Set default assay to SCT
PC_ 1 
Positive:  Sftpc, Sftpb, Sftpa1, Sftpd, Cxcl15, Npc2, Scd1, Lyz2, Cd74, Wfdc2 
	   Napsa, Slc34a2, Chil1, Atp1b1, S100g, Cldn3, Sfta2, Cbr2, Cldn18, Ctsh 
	   Ager, Lcn2, Lpcat1, Aqp5, Hp, Mgst1, Alcam, Lamp3, Dram1, Muc1 
Negative:  Cd52, Cd3g, S100a4, Ptprc, Trbc2, Lgals1, Il7r, Coro1a, Srgn, Tmsb10 
	   Cd3d, Arhgdib, Ms4a4b, Ltb, Vps37b, Rac2, Cd3e, Ptprcap, Laptm5, Tnfaip3 
	   Ets1, AW112010, Thy1, Lsp1, Cxcr6, Rgs2, H2-Q7, Fxyd5, Selplg, S100a10 
PC_ 2 
Positive:  Egfl7, Calcrl, Ramp2, Pecam1, Cdh5, BC028528, Cavin2, Ptprb, Ly6a, Ifitm3 
	   Ace, Aqp1, Tspan7, Acvrl1, Pltp, Cldn5, Podxl, Ly6c1, Flt1, Cav1 
	   Tmem100, Cd93, Cyyr1, Foxf1, Ecscr, Adgrl4, Eng, Icam2, Ier3, Gpihbp1 
Negative:  Cbr2, Lyz2, Cxcl15, Npc2, Trf, Hp, Selenbp1, S100g, Scgb1a1, Rps24 
	   Lcn2, Chchd10, Mgst1, Retnla, AW112010, Ldhb, Cyp2f2, Wfdc2, Nupr1, Ly6c2 
	   Scp2, Chil1, Cox6a2, Ccl5, Cd52, Nkg7, Lgals1, Me1, Cxcl17, H2afj 
PC_ 3 
Positive:  Maf, Igf1r, Lmo4, Fosl2, Ccr2, mt-Co1, Junb, Il7r, Nr4a1, Ikzf3 
	   Zfp36l1, Tcrg-C1, Pxdc1, Il2ra, Samsn1, Furin, mt-Co3, Serpinb1a, Hk2, Blk 
	   Gadd45b, mt-Atp6, Emb, mt-Co2, Rora, Ramp1, Pim1, Nfkbid, Uhrf2, Tmem64 
Negative:  Ms4a4b, Ccl5, Nkg7, AW112010, Cd8b1, Klrd1, Klf2, Ms4a6b, Cd8a, Ctsw 
	   Gimap4, Ly6c2, Ctla2a, Gimap6, Epsti1, Gramd3, Ly6a, B2m, Cd2, Trac 
	   Psmb9, Rac2, Cd52, Klrc1, Cd3d, Gm8369, Prf1, Rps24, Gimap7, Arhgdib 
PC_ 4 
Positive:  Ms4a4b, Lamp3, Lpcat1, mt-Co1, Slc34a2, Ptprc, Ccl5, Napsa, Scd1, Hc 
	   Ager, Cldn18, Abca3, mt-Atp6, mt-Co3, Itga4, Fasn, Egfl6, Xist, mt-Co2 
	   Scd2, Elovl1, Soat1, Gramd3, Nkg7, Chil1, Mme, Malt1, Rnase4, Dennd4a 
Negative:  Scgb1a1, S100a6, Tmem176a, Tmem176b, Cd24a, Cxcr6, Ifitm3, S100a4, Capg, Tcrg-C1 
	   Scgb3a1, Pxdc1, Actn2, Ramp1, Ltb4r1, Gsto1, Igf1r, Serpinb1a, Lgals3, Rbp4 
	   Ckb, Krt83, Blk, Cp, Gpx2, Aldh1a1, Rorc, Anxa1, Muc5b, Ly6g5b 
PC_ 5 
Positive:  S100a4, Lgals1, Cd52, Cd3g, Capg, S100a10, Vim, Sftpc, Cxcr6, Napsa 
	   Sh3bgrl3, Rnase4, Cd74, Ramp1, B2m, Npc2, Lyz2, Lamp3, Ptprcap, Ager 
	   Chil1, S100a6, Il33, Scd1, Actb, Lsp1, H2-Aa, H2-Ab1, Elovl1, Sfta2 
Negative:  Scgb1a1, mt-Co2, mt-Co3, mt-Atp6, mt-Co1, mt-Cytb, Gm26917, Cd24a, mt-Nd1, Cyp2f2 
	   Neat1, Cp, mt-Nd4, mt-Nd2, Malat1, B430010I23Rik, Rbp4, Aldh1a1, Ffar4, mt-Nd5 
	   Dennd4a, Gsto1, Ltf, Tspan1, Clu, Vps37b, Krt19, Scgb3a1, Pigr, F3 
Calculating cell attributes from input UMI matrix: log_umi
Variance stabilizing transformation of count matrix of size 14712 by 873
Model formula is y ~ log_umi
Get Negative Binomial regression parameters per gene
Using 2000 genes, 873 cells
Found 83 outliers - those will be ignored in fitting/regularization step

Second step: Get residuals using fitted parameters for 14712 genes
Computing corrected count matrix for 14712 genes
Calculating gene attributes
Wall clock passed: Time difference of 26.24273 secs
Determine variable features
Place corrected count matrix in counts slot
Centering data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
Set default assay to SCT
PC_ 1 
Positive:  Sftpc, Sftpb, Sftpa1, Sftpd, Cd74, Npc2, Cxcl15, Wfdc2, Scd1, Lyz2 
	   Napsa, Slc34a2, Atp1b1, Chil1, Sfta2, S100g, Cbr2, Ctsh, Cldn3, Cldn18 
	   Alcam, Lcn2, Ager, Hp, Lpcat1, Mgst1, H2-Ab1, Muc1, H2-Aa, Dram1 
Negative:  Cd52, Cd3g, Lgals1, Ms4a4b, Ptprc, S100a4, Trbc2, Coro1a, Arhgdib, Srgn 
	   Tmsb10, Rac2, Cd3d, Il7r, AW112010, Nkg7, Fxyd5, Cd3e, H2-Q7, Laptm5 
	   Vim, Ccl5, Ptprcap, S100a10, Trac, Vps37b, Ets1, Lsp1, Cxcr6, Thy1 
PC_ 2 
Positive:  Ramp2, Egfl7, Cdh5, Calcrl, Acvrl1, Cav1, Cavin2, Tmem100, Ace, Cldn5 
	   Pecam1, Ecscr, Icam2, Klf4, Cd93, Aqp1, Tm4sf1, Tspan7, Gpihbp1, Pltp 
	   Flt1, Ptprb, BC028528, Eng, Clec14a, Ly6a, Thbd, Plvap, Cyyr1, Slco2a1 
Negative:  Nkg7, Lgals1, Ccl5, Cd52, S100a6, Cd8b1, Ms4a4b, Cd3g, Cd8a, Arhgdib 
	   Hcst, AW112010, Sh3bgrl3, Nupr1, Trac, Rac2, S100a4, Trbc2, Sftpa1, Cd3d 
	   Thy1, Emp3, Klrd1, Ctsw, Dusp2, Rps19, Cxcl15, Retnla, Cbr2, Npc2 
PC_ 3 
Positive:  Ccl5, Nkg7, Ms4a4b, Cd8b1, AW112010, Klf2, Cd8a, Ms4a6b, Vim, Lgals1 
	   Ctla2a, Epsti1, Klrd1, Ctsw, Cd52, B2m, Arhgdib, Cd3d, Gimap4, Ly6a 
	   Rac2, Hcst, Itga1, Sh3bgrl3, Cd3g, Tuba1a, Trac, Crip1, Ly6c2, Tmsb10 
Negative:  mt-Co3, mt-Atp6, mt-Co2, mt-Co1, mt-Cytb, mt-Nd1, Il7r, Lmo4, Ccr2, mt-Nd2 
	   Furin, Maf, mt-Nd4, Samsn1, Igf1r, Il2ra, Emb, Gadd45b, Nrip1, Zfp36l1 
	   Ramp3, Pim1, Ikzf3, Ltb4r1, Icos, Fosl2, Nfkbid, Uhrf2, St6galnac3, Tmem64 
PC_ 4 
Positive:  Slc34a2, Lamp3, Napsa, Lpcat1, Egfl6, Rnase4, Hc, Scd1, Fasn, Abca3 
	   Scd2, Cldn18, Lyz2, Mme, Cd74, Ager, Chil1, H2-Ab1, Soat1, H2-Aa 
	   Tgoln1, Elovl1, Alcam, Ptprf, Ppp1r14c, Lgi3, Muc1, Sftpc, Slco4c1, Sfta2 
Negative:  Scgb1a1, Scgb3a2, Cd24a, Ifitm3, B430010I23Rik, Cyp2f2, Cp, Cbr2, Gsto1, Cldn10 
	   Gsta4, Aldh1a1, Gabrp, Pigr, Cckar, Plpp3, Tmem176a, Gpx2, Hes1, Pam 
	   Hspb1, Rbp4, Clu, Gstm1, Rassf9, Gstm2, Ldhb, Krt19, Ifitm2, Qsox1 
PC_ 5 
Positive:  Gramd3, Gm2682, Txk, Ly6c2, Slamf6, Klf3, Pde2a, Lef1, Ms4a4b, Ripor2 
	   Klf2, Sidt1, S1pr1, Bcl2, mt-Co3, mt-Atp6, Vps37b, Atp11b, Itga4, mt-Co2 
	   Ugcg, Cd27, Sell, Gimap4, Macf1, Arhgap15, mt-Co1, Dennd4a, Cmah, Scgb1a1 
Negative:  S100a4, S100a6, Cxcr6, Lgals1, Cd52, Vim, Capg, S100a10, Cd3g, Fgl2 
	   Sh3bgrl3, Lgals3, Jaml, Ramp1, Gng2, Thy1, Itm2b, Krt83, Ptprcap, Fxyd5 
	   Srgn, Ctla2a, Emp3, Id2, Ltb4r1, Bcl2a1b, Rgs1, Ltb, Itgae, Lsp1 
Calculating cell attributes from input UMI matrix: log_umi
Variance stabilizing transformation of count matrix of size 17086 by 6553
Model formula is y ~ log_umi
Get Negative Binomial regression parameters per gene
Using 2000 genes, 5000 cells
Found 112 outliers - those will be ignored in fitting/regularization step

Second step: Get residuals using fitted parameters for 17086 genes
Computing corrected count matrix for 17086 genes
Calculating gene attributes
Wall clock passed: Time difference of 3.22684 mins
Determine variable features
Place corrected count matrix in counts slot
Centering data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
Set default assay to SCT
PC_ 1 
Positive:  Sftpc, Sftpa1, Lyz2, Scgb1a1, Sftpb, Cxcl15, Sftpd, Scd1, Retnla, Wfdc2 
	   Slc34a2, Npc2, Hp, Napsa, Chil1, Lcn2, Trf, Cd74, S100g, Ctsh 
	   Clu, Cbr2, Cyp2f2, Atp1b1, Bsg, Hc, Cldn18, Selenbp1, Nupr1, Krt18 
Negative:  Ccl5, S100a4, Cd52, Ms4a4b, Cd3g, Trbc2, Tmsb10, Nkg7, Il7r, Lgals1 
	   Cd3d, Ltb, Trbc1, AW112010, Coro1a, Arhgdib, Cxcr6, Ptprcap, Rac2, Ptprc 
	   Hcst, S100a10, H2-Q7, Icos, Trac, Thy1, Rgs1, Ms4a6b, Cd3e, Cd8b1 
PC_ 2 
Positive:  Scgb1a1, Cldn5, Hpgd, Scgb3a1, Ly6a, Ramp2, Calcrl, Cdh5, Tmem100, Ifitm3 
	   Ptprb, Ier3, Cavin2, Egfl7, Sema3c, Clu, Tm4sf1, Cxcl12, Icam2, Pltp 
	   Aqp1, Cd93, Ltf, Gpihbp1, Clec14a, Scgb3a2, Thbd, BC028528, Tspan7, Acvrl1 
Negative:  Sftpc, Cxcl15, Sftpa1, Slc34a2, Scd1, Ccl5, Napsa, Sftpb, Chil1, Sftpd 
	   Lcn2, Lyz2, S100a4, Cd52, Nkg7, Npc2, Lamp3, Hc, Elovl1, Lpcat1 
	   S100g, Rnase4, Ms4a4b, Cd3g, Lgals1, Trbc2, Uba52, Fasn, Sfta2, Ly6c2 
PC_ 3 
Positive:  Sftpc, Lyz2, Cxcl15, Slc34a2, Scd1, Cldn5, Sftpa1, Hpgd, Napsa, Chil3 
	   Ramp2, Calcrl, Elovl1, Tmem100, Fabp5, Chil1, Cdh5, Ptprb, Egfl7, Lamp3 
	   Ccl6, Rnase4, Cd36, Cxcl12, Sema3c, Cavin2, Ly6c1, Icam2, Cd93, Hc 
Negative:  Scgb1a1, Clu, Ltf, Scgb3a1, Wfdc2, Scgb3a2, S100a6, Hp, Cp, Krt18 
	   Areg, Bpifb1, Anxa1, Cyp2f2, Rbp4, Krt19, Krt8, Trf, Muc5b, Tspan1 
	   Agr2, Cxcl17, Krt7, Ifitm1, Ffar4, Tff2, Nupr1, Basp1, Epcam, Fxyd3 
PC_ 4 
Positive:  Sftpb, Hpgd, Cldn5, Scd1, Ptprb, Cdh5, Cxcl15, Calcrl, Napsa, Sema3c 
	   Slc34a2, Tmem100, Sftpd, Cxcl12, Egfl7, Sftpa1, Ly6a, Ramp2, Cavin2, Gpihbp1 
	   Cd93, Clec14a, Tspan7, Cd36, Ifitm3, Aqp1, Itga1, Slco2a1, Scn7a, Ace 
Negative:  Chil3, Ccl6, Lyz2, Fth1, Tyrobp, Ftl1, Ctsd, Ctss, Cybb, Fcer1g 
	   Plet1, Mpeg1, Ear2, Lpl, Fabp5, Psap, Lgals3, Actb, Mrc1, Cxcl2 
	   Atp6v0d2, Cebpb, Abcg1, Il18, Alox5ap, Fabp4, Msrb1, Marco, Sirpa, Cd68 
PC_ 5 
Positive:  Ccl5, Nkg7, Ms4a4b, Ly6c2, Cd8b1, Klrd1, Gzma, Scgb1a1, Ms4a6b, Epsti1 
	   Gramd3, Gimap4, Klre1, Gimap7, Klf2, Klrk1, Cd8a, Bcl2, Gm2682, Ccl4 
	   Ctsw, Xcl1, Itga4, Dusp2, Gzmb, Txk, AW112010, Eomes, Retnla, Klrc2 
Negative:  S100a4, Cxcr6, Il17a, Trdc, Krt83, Maf, Gadd45b, Icos, Ramp1, Tcrg-C1 
	   Cd163l1, S100a6, Actn2, Ltb4r1, Tmem176a, Capg, Junb, Il7r, Ccr2, Serpinb1a 
	   Blk, Pxdc1, Lmo4, Il2ra, Ly6g5b, Lgals1, Tmem176b, Jaml, Nfkbia, Ltb 
Calculating cell attributes from input UMI matrix: log_umi
Variance stabilizing transformation of count matrix of size 16179 by 4229
Model formula is y ~ log_umi
Get Negative Binomial regression parameters per gene
Using 2000 genes, 4229 cells
Found 93 outliers - those will be ignored in fitting/regularization step

Second step: Get residuals using fitted parameters for 16179 genes
Computing corrected count matrix for 16179 genes
Calculating gene attributes
Wall clock passed: Time difference of 2.15837 mins
Determine variable features
Place corrected count matrix in counts slot
Centering data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
Set default assay to SCT
PC_ 1 
Positive:  Sftpc, Sftpa1, Lyz2, Sftpb, Cxcl15, Sftpd, Slc34a2, Scd1, Npc2, Retnla 
	   Napsa, Wfdc2, Chil1, Cd74, Lcn2, Hp, S100g, Scgb1a1, Hc, Cbr2 
	   Ctsh, Ager, Elovl1, Lpcat1, Cldn18, Atp1b1, Lamp3, Rnase4, Chia1, Mgst1 
Negative:  Ccl5, Cd52, S100a4, Cd3g, Lgals1, Tmsb10, Ms4a4b, Nkg7, Trbc2, Il7r 
	   AW112010, Coro1a, Ptprc, Ltb, Cd3d, S100a6, Vim, Arhgdib, Cxcr6, Rac2 
	   Rgs1, Rgs2, Srgn, Ptprcap, Icos, S100a10, Trac, Junb, Fxyd5, Trbc1 
PC_ 2 
Positive:  Sftpc, Lyz2, Sftpa1, Cxcl15, Slc34a2, Sftpb, Scd1, S100a4, Napsa, Cd52 
	   Sftpd, Chil1, Lgals1, Lcn2, Cd3g, Nkg7, Npc2, Ccl5, Cxcr6, Hc 
	   Rgs1, Icos, Elovl1, Trbc2, Rgs2, Lpcat1, S100g, Rnase4, Tmsb4x, Lamp3 
Negative:  Cldn5, Ramp2, Ier3, Ly6a, Calcrl, Cxcl12, Ptprb, Cdh5, Hpgd, Ifitm3 
	   Gpihbp1, Egfl7, Cavin2, Tmem100, Scgb1a1, Pltp, Icam2, Tm4sf1, Clec14a, Ace 
	   Aqp1, Acvrl1, Scn7a, Cav1, Sema3c, Pecam1, Tspan7, Klf4, Ly6c1, Cd93 
PC_ 3 
Positive:  Scgb1a1, Scgb3a2, Scgb3a1, Cyp2f2, Hp, Clu, Ltf, Cp, Wfdc2, Trf 
	   Krt19, Krt18, Krt8, Cbr2, Rbp4, Retnla, S100a6, Anxa1, Gsto1, Pigr 
	   Fxyd3, Tspan1, Bpifb1, Krt7, Selenbp1, Cxcl17, Ly6d, Muc5b, Gpx2, Basp1 
Negative:  Sftpc, Lyz2, Cxcl15, Sftpa1, Slc34a2, Scd1, Napsa, Chil1, Lcn2, Elovl1 
	   Lamp3, Lrg1, Cldn5, Cd74, Hc, Rnase4, Fabp5, Sftpb, Lpcat1, Npc2 
	   Ramp2, Sftpd, Cxcl12, Cdh5, Calcrl, Egfl7, Ptprb, Gpihbp1, Tmem100, Cd36 
PC_ 4 
Positive:  S100a4, Gadd45b, Junb, S100a6, Icos, Cxcr6, Maf, Nfkbia, Samsn1, Areg 
	   Ccr2, Furin, Ltb4r1, Ramp1, Serpinb1a, Capg, Nr4a1, Il17a, Ikzf3, Actn2 
	   Pxdc1, Lgals3, Rgs2, Uhrf2, Krt83, Tmem176a, Lgals1, Rora, Rgs1, Tcrg-C1 
Negative:  Ccl5, Ms4a4b, Nkg7, Ly6c2, Klrd1, Cd8b1, Ncr1, Gramd3, Xcl1, Gimap4 
	   Bcl2, Ctsw, Klre1, Sell, Rps24, Txk, Gimap6, Gm2682, Epsti1, AW112010 
	   Retnla, Rps16, Rpl13, Cd7, Ms4a6b, Gimap7, Cd27, Klf2, Pde2a, Rpsa 
PC_ 5 
Positive:  S100a6, Ccl5, Clu, Nkg7, Cd8b1, Cd8a, Krt8, Anxa1, Krt18, Sftpa1 
	   Krt19, Areg, Fos, Krt7, Sftpd, Atf3, Egr1, Basp1, Wfdc2, Ly6d 
	   Sftpb, Ltf, Lgals1, AW112010, Ctsw, Gsto1, Jun, Lgals3, Cd52, Aqp5 
Negative:  Scgb1a1, Retnla, Scgb3a2, Il7r, Hp, Maf, Gm42418, Samsn1, Cbr2, Icos 
	   Tnfrsf4, Ikzf2, Ramp3, Il2ra, Pim1, Il1rl1, B430010I23Rik, Serpinb1a, Gata3, Trbc1 
	   Cldn10, Ccr2, Gm26740, Ltb4r1, Gadd45b, Slamf6, Rora, Nfkbid, Trdc, Ifi27l2a 
Calculating cell attributes from input UMI matrix: log_umi
Variance stabilizing transformation of count matrix of size 15767 by 3519
Model formula is y ~ log_umi
Get Negative Binomial regression parameters per gene
Using 2000 genes, 3519 cells
Found 155 outliers - those will be ignored in fitting/regularization step

Second step: Get residuals using fitted parameters for 15767 genes
Computing corrected count matrix for 15767 genes
Calculating gene attributes
Wall clock passed: Time difference of 1.683764 mins
Determine variable features
Place corrected count matrix in counts slot
Centering data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
Set default assay to SCT
PC_ 1 
Positive:  Sftpa1, Sftpc, Sftpb, Sftpd, Wfdc2, Cxcl15, Npc2, Lyz2, Cd74, Scd1 
	   Slc34a2, Napsa, Chil1, Cbr2, Cldn3, Ctsh, Cldn18, Atp1b1, S100g, Hp 
	   Mgst1, Alcam, Sfta2, Retnla, Cd63, Lcn2, Dram1, Lpcat1, Aqp5, Krt18 
Negative:  Ccl5, S100a4, Cd52, Nkg7, Tmsb10, Ms4a4b, Lgals1, Trbc1, Trbc2, Cd3g 
	   Il7r, Ly6c2, AW112010, Cd8b1, Cxcr6, Ltb, Cd3d, Junb, Ptprc, Hcst 
	   Coro1a, Rgs2, Vps37b, Trac, Ifi27l2a, Ms4a6b, Rps24, S100a10, Icos, Ptprcap 
PC_ 2 
Positive:  Ccl5, S100a6, S100a4, Lgals1, Nkg7, Tmsb4x, Cd52, Cxcr6, Ly6c2, Crip1 
	   Rgs1, Lgals3, Tmsb10, S100a10, Uba52, Zfp36l2, Klrk1, S100a11, AW112010, Bcl2a1b 
	   Capg, Cd3g, Lsp1, Rgs2, Thy1, Icos, Ptprcap, Ccr2, Wfdc2, Ms4a4b 
Negative:  Hpgd, Ramp2, Cldn5, Ptprb, Tmem100, Calcrl, Pltp, Clec14a, Ifitm3, Cdh5 
	   Tm4sf1, Ier3, Gpihbp1, Cavin2, Ly6a, Tspan7, Jun, Sema3c, Ehd4, Aqp1 
	   Egfl7, Cd93, Plvap, Itga1, Acvrl1, Eng, Clec1a, Cd36, Bmpr2, Cxcl12 
PC_ 3 
Positive:  Ccl5, Ly6c2, Nkg7, Ms4a4b, Cd8b1, Klrd1, Gzma, Klf2, Rpl13, Fcer1g 
	   Gramd3, Rps24, Rps16, Bcl2, Rps15a, Rpl23, Epsti1, Gimap4, Klre1, Rpsa 
	   Tyrobp, Rps7, Rps10, Ctsw, Rpl12, Sell, Txk, Rpl32, Rps20, Ms4a6b 
Negative:  S100a4, Gadd45b, S100a6, Cxcr6, Tmem176a, Ramp1, Lmo4, Tmem176b, Icos, Capg 
	   Junb, Pxdc1, Maf, Serpinb1a, Ccr2, Il17a, Bcl2a1b, Ltb4r1, Rora, Furin 
	   Actn2, Ckb, Rgcc, Nfkbia, Lgals3, Areg, Lgals1, Ramp3, Nr4a1, Bcl2a1d 
PC_ 4 
Positive:  Sftpc, Lyz2, Cxcl15, Slc34a2, Napsa, Lamp3, Rnase4, Chil1, Elovl1, Lcn2 
	   Sftpa1, Scd1, Hc, Lpcat1, Lgi3, Abca3, Fasn, Fabp5, Cd74, Egfl6 
	   Sfta2, Ppp1r14c, Il33, Soat1, H2-Aa, Mid1ip1, Lrp2, Itih4, Dmkn, Lrg1 
Negative:  Scgb1a1, Cyp2f2, Clu, Krt18, Krt19, Rbp4, Krt8, Cxcl17, Scgb3a1, Cp 
	   Krt7, Tspan1, Anxa1, Hp, Fxyd3, Selenbp1, Ltf, Ffar4, Pigr, Cd24a 
	   Scgb3a2, F3, Gsta4, Mgst1, Cd9, S100a6, Gsto1, Muc5b, B430010I23Rik, Aldh1a1 
PC_ 5 
Positive:  Scgb1a1, mt-Co3, Scgb3a2, Retnla, mt-Atp6, Hp, Cbr2, Trf, mt-Co2, mt-Nd4 
	   Sftpc, mt-Nd2, mt-Cytb, Lyz2, Scgb3a1, Nupr1, mt-Co1, mt-Nd1, Cyp2f2, Inmt 
	   Cldn10, Sftpa1, Selenbp1, Pon1, Aldh1a1, Cp, Tff2, Pigr, Cxcl15, Fmo2 
Negative:  Ccl5, S100a6, Lgals1, S100a4, Nkg7, Cd8b1, Vim, Lgals3, Actb, Tmsb4x 
	   Crip1, Gzmb, Cd8a, Cd52, Gzma, S100a10, Rgs1, H2-Ab1, Actg1, Jun 
	   Tspan8, Cxcr6, Egr1, Cldn18, Ms4a4b, Fgl2, Sparc, H2-Aa, Klrk1, H2-Eb1 
  |                                                                              |.......................                                               |  33%
   inline R code fragments

  |                                                                              |............................                                          |  40%
label: Plot_DF (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |.................................                                     |  47%
  ordinary text without R code

  |                                                                              |.....................................                                 |  53%
label: UMAP_percentage (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: logi FALSE

  |                                                                              |..........................................                            |  60%
   inline R code fragments

  |                                                                              |...............................................                       |  67%
label: Plot_percent_MT (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |...................................................                   |  73%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: UMAP_nFeatures (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: logi FALSE

  |                                                                              |.............................................................         |  87%
   inline R code fragments

  |                                                                              |.................................................................     |  93%
label: Plot_nFeatureRNA (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |................................................                      |  68%
  ordinary text without R code

  |                                                                              |..................................................                    |  72%
label: unnamed-chunk-8 (with options) 
List of 2
 $ child: chr "08-perform.cell.cycle.regression.Rmd"
 $ eval : logi TRUE



processing file: ./08-perform.cell.cycle.regression.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |........                                                              |  11%
  ordinary text without R code

  |                                                                              |................                                                      |  22%
label: create_cell_cycle_reg_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |.......................                                               |  33%
  ordinary text without R code

  |                                                                              |...............................                                       |  44%
label: createIntegrated (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi TRUE
 $ results: logi FALSE

Computing 2000 integration features
Scaling features for provided objects
Finding all pairwise anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2714 anchors
Filtering anchors
	Retained 2490 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2766 anchors
Filtering anchors
	Retained 2483 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 1932 anchors
Filtering anchors
	Retained 1759 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7491 anchors
Filtering anchors
	Retained 5099 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2758 anchors
Filtering anchors
	Retained 1872 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2799 anchors
Filtering anchors
	Retained 1751 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 6855 anchors
Filtering anchors
	Retained 4909 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2712 anchors
Filtering anchors
	Retained 1916 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2714 anchors
Filtering anchors
	Retained 1818 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7996 anchors
Filtering anchors
	Retained 5255 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 6484 anchors
Filtering anchors
	Retained 4947 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2466 anchors
Filtering anchors
	Retained 1819 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2616 anchors
Filtering anchors
	Retained 1902 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 6676 anchors
Filtering anchors
	Retained 4513 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 5951 anchors
Filtering anchors
	Retained 4081 anchors
Merging dataset 2 into 1
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 3 into 6
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 6 3 into 1 2
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 5 into 4
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 1 2 6 3 into 4 5
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:36:36 UMAP embedding parameters a = 0.9922 b = 1.112
Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
Also defined by 'spam'
10:36:36 Read 15319 rows and found 30 numeric columns
10:36:36 Using Annoy for neighbor search, n_neighbors = 30
Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
Also defined by 'spam'
10:36:36 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:36:38 Writing NN index file to temp file /tmp/slurm_47112575/RtmpxgvIpe/file5cb3f33c8e3
10:36:38 Searching Annoy index using 1 thread, search_k = 3000
10:36:43 Annoy recall = 100%
10:36:45 Commencing smooth kNN distance calibration using 1 thread
10:36:48 Initializing from normalized Laplacian + noise
10:36:49 Commencing optimization for 200 epochs, with 696986 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:37:13 Optimization finished
Computing nearest neighbor graph
Computing SNN
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Computing 2000 integration features
Scaling features for provided objects
Finding all pairwise anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2714 anchors
Filtering anchors
	Retained 2490 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2766 anchors
Filtering anchors
	Retained 2483 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 1932 anchors
Filtering anchors
	Retained 1759 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7491 anchors
Filtering anchors
	Retained 5099 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2758 anchors
Filtering anchors
	Retained 1872 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2799 anchors
Filtering anchors
	Retained 1751 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 6855 anchors
Filtering anchors
	Retained 4909 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2712 anchors
Filtering anchors
	Retained 1916 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2714 anchors
Filtering anchors
	Retained 1818 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7996 anchors
Filtering anchors
	Retained 5255 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 6484 anchors
Filtering anchors
	Retained 4947 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2466 anchors
Filtering anchors
	Retained 1819 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2616 anchors
Filtering anchors
	Retained 1902 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 6676 anchors
Filtering anchors
	Retained 4513 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 5951 anchors
Filtering anchors
	Retained 4081 anchors
Merging dataset 2 into 1
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 3 into 6
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 6 3 into 1 2
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 5 into 4
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 1 2 6 3 into 4 5
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:52:54 UMAP embedding parameters a = 0.9922 b = 1.112
Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
Also defined by 'spam'
10:52:54 Read 15319 rows and found 30 numeric columns
10:52:54 Using Annoy for neighbor search, n_neighbors = 30
Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
Also defined by 'spam'
10:52:54 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:52:56 Writing NN index file to temp file /tmp/slurm_47112575/RtmpxgvIpe/file5cb372f9cc25
10:52:56 Searching Annoy index using 1 thread, search_k = 3000
10:53:01 Annoy recall = 100%
10:53:03 Commencing smooth kNN distance calibration using 1 thread
10:53:06 Initializing from normalized Laplacian + noise
10:53:06 Commencing optimization for 200 epochs, with 701834 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
10:53:30 Optimization finished
Computing nearest neighbor graph
Computing SNN
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Warning: The following features are not present in the object: Fam64a, Hn1, not searching for symbol synonyms
Computing 2000 integration features
Scaling features for provided objects
Finding all pairwise anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2714 anchors
Filtering anchors
	Retained 2490 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2766 anchors
Filtering anchors
	Retained 2483 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 1932 anchors
Filtering anchors
	Retained 1759 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7491 anchors
Filtering anchors
	Retained 5099 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2758 anchors
Filtering anchors
	Retained 1872 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2799 anchors
Filtering anchors
	Retained 1751 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 6855 anchors
Filtering anchors
	Retained 4909 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2712 anchors
Filtering anchors
	Retained 1916 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2714 anchors
Filtering anchors
	Retained 1818 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7996 anchors
Filtering anchors
	Retained 5255 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 6484 anchors
Filtering anchors
	Retained 4947 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2466 anchors
Filtering anchors
	Retained 1819 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 2616 anchors
Filtering anchors
	Retained 1902 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 6676 anchors
Filtering anchors
	Retained 4513 anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 5951 anchors
Filtering anchors
	Retained 4081 anchors
Merging dataset 2 into 1
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 3 into 6
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 6 3 into 1 2
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 5 into 4
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 1 2 6 3 into 4 5
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Calculating gene variances
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Calculating feature variances of standardized and clipped values
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
11:11:19 UMAP embedding parameters a = 0.9922 b = 1.112
Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
Also defined by 'spam'
11:11:19 Read 15319 rows and found 30 numeric columns
11:11:19 Using Annoy for neighbor search, n_neighbors = 30
Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
Also defined by 'spam'
11:11:20 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
11:11:22 Writing NN index file to temp file /tmp/slurm_47112575/RtmpxgvIpe/file5cb3f95cce0
11:11:22 Searching Annoy index using 1 thread, search_k = 3000
11:11:27 Annoy recall = 100%
11:11:29 Commencing smooth kNN distance calibration using 1 thread
11:11:32 Initializing from normalized Laplacian + noise
11:11:33 Commencing optimization for 200 epochs, with 697768 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
11:11:57 Optimization finished
Computing nearest neighbor graph
Computing SNN
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
  |                                                                              |.......................................                               |  56%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: unnamed-chunk-9 (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: logi FALSE


Attaching package: 'dplyr'

The following objects are masked from 'package:stats':

    filter, lag

The following objects are masked from 'package:base':

    intersect, setdiff, setequal, union

  |                                                                              |......................................................                |  78%
   inline R code fragments

  |                                                                              |..............................................................        |  89%
label: unnamed-chunk-10 (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |.....................................................                 |  76%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: saveobject (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi TRUE
 $ results: logi FALSE

  |                                                                              |...........................................................           |  84%
  ordinary text without R code

  |                                                                              |..............................................................        |  88%
label: documentation (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi TRUE
 $ results: logi TRUE

  |                                                                              |................................................................      |  92%
  ordinary text without R code

  |                                                                              |...................................................................   |  96%
label: create_report_params (with options) 
List of 2
 $ eval   : logi TRUE
 $ results: chr "asis"

  |                                                                              |......................................................................| 100%
   inline R code fragments


output file: QC.knit.md

/camp/apps/eb/software/pandoc/2.2.3.2-foss-2016b/bin/pandoc +RTS -K512m -RTS QC.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/528_scRNAseq_subclustering_of_496_C7_C12_subset1/html_local/QC.html --lua-filter /camp/stp/babs/working/boeings/Projects/reference_data/renv/cache/v5/R-4.0/x86_64-pc-linux-gnu/rmarkdown/2.13/ac78f4d2e0289d4cba73b88af567b8b1/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/stp/babs/working/boeings/Projects/reference_data/renv/cache/v5/R-4.0/x86_64-pc-linux-gnu/rmarkdown/2.13/ac78f4d2e0289d4cba73b88af567b8b1/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 5 --variable toc_float=1 --variable toc_selectors=h1,h2,h3,h4,h5 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/stp/babs/working/boeings/Projects/reference_data/renv/cache/v5/R-4.0/x86_64-pc-linux-gnu/rmarkdown/2.13/ac78f4d2e0289d4cba73b88af567b8b1/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --variable theme=paper --css src/style/style.css --mathjax --variable 'mathjax-url=https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --include-in-header /tmp/slurm_47112575/RtmpxgvIpe/rmarkdown-str5cb33f3f1299.html --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/528_scRNAseq_subclustering_of_496_C7_C12_subset1/html_local/QC.html
There were 50 or more warnings (use warnings() to see the first 50)
