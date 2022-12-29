
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
03:00:26 UMAP embedding parameters a = 0.9922 b = 1.112
03:00:26 Read 3296 rows and found 30 numeric columns
03:00:26 Using Annoy for neighbor search, n_neighbors = 30
03:00:26 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:00:27 Writing NN index file to temp file /tmp/slurm_47105226/RtmpCkULvu/file12f8b323768e4
03:00:27 Searching Annoy index using 1 thread, search_k = 3000
03:00:28 Annoy recall = 100%
03:00:29 Commencing smooth kNN distance calibration using 1 thread
03:00:31 Initializing from normalized Laplacian + noise
03:00:31 Commencing optimization for 500 epochs, with 144966 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:00:44 Optimization finished
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
03:00:58 UMAP embedding parameters a = 0.9922 b = 1.112
03:00:58 Read 642 rows and found 30 numeric columns
03:00:58 Using Annoy for neighbor search, n_neighbors = 30
03:00:58 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:00:58 Writing NN index file to temp file /tmp/slurm_47105226/RtmpCkULvu/file12f8b4bc136ae
03:00:58 Searching Annoy index using 1 thread, search_k = 3000
03:00:58 Annoy recall = 100%
03:00:59 Commencing smooth kNN distance calibration using 1 thread
03:01:01 Initializing from normalized Laplacian + noise
03:01:01 Commencing optimization for 500 epochs, with 24924 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:01:04 Optimization finished
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
03:01:17 UMAP embedding parameters a = 0.9922 b = 1.112
03:01:17 Read 655 rows and found 30 numeric columns
03:01:17 Using Annoy for neighbor search, n_neighbors = 30
03:01:17 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:01:17 Writing NN index file to temp file /tmp/slurm_47105226/RtmpCkULvu/file12f8b350979eb
03:01:17 Searching Annoy index using 1 thread, search_k = 3000
03:01:17 Annoy recall = 100%
03:01:18 Commencing smooth kNN distance calibration using 1 thread
03:01:20 Initializing from normalized Laplacian + noise
03:01:20 Commencing optimization for 500 epochs, with 26044 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:01:23 Optimization finished
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
03:02:45 UMAP embedding parameters a = 0.9922 b = 1.112
03:02:45 Read 4915 rows and found 30 numeric columns
03:02:45 Using Annoy for neighbor search, n_neighbors = 30
03:02:45 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:02:45 Writing NN index file to temp file /tmp/slurm_47105226/RtmpCkULvu/file12f8b2fe95a2e
03:02:45 Searching Annoy index using 1 thread, search_k = 3000
03:02:47 Annoy recall = 100%
03:02:48 Commencing smooth kNN distance calibration using 1 thread
03:02:50 Initializing from normalized Laplacian + noise
03:02:50 Commencing optimization for 500 epochs, with 211618 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:03:09 Optimization finished
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
03:03:58 UMAP embedding parameters a = 0.9922 b = 1.112
03:03:58 Read 3172 rows and found 30 numeric columns
03:03:58 Using Annoy for neighbor search, n_neighbors = 30
03:03:58 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:03:58 Writing NN index file to temp file /tmp/slurm_47105226/RtmpCkULvu/file12f8b320a5eef
03:03:58 Searching Annoy index using 1 thread, search_k = 3000
03:03:59 Annoy recall = 100%
03:04:00 Commencing smooth kNN distance calibration using 1 thread
03:04:02 Initializing from normalized Laplacian + noise
03:04:02 Commencing optimization for 500 epochs, with 137878 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:04:15 Optimization finished
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
03:04:55 UMAP embedding parameters a = 0.9922 b = 1.112
03:04:55 Read 2639 rows and found 30 numeric columns
03:04:55 Using Annoy for neighbor search, n_neighbors = 30
03:04:55 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:04:56 Writing NN index file to temp file /tmp/slurm_47105226/RtmpCkULvu/file12f8b7f076a79
03:04:56 Searching Annoy index using 1 thread, search_k = 3000
03:04:56 Annoy recall = 100%
03:04:57 Commencing smooth kNN distance calibration using 1 thread
03:04:59 Initializing from normalized Laplacian + noise
03:05:00 Commencing optimization for 500 epochs, with 114598 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
03:05:10 Optimization finished
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
