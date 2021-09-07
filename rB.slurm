

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

Loading required package: biologicSeqTools

Attaching package: 'biologicSeqTools'

The following objects are masked _by_ '.GlobalEnv':

    add2labCatSelectionDBtable, add2vec,
    add.category.to.lab.reference.table.hs, addCoVar,
    addDGEcomparisons2DesignFile, addGeneAnnotation,
    addProject2ProjectTable, addTPMorFPKMtable,
    completeDesignBasedOnSampleID, createAlignmentShellScript,
    createAnalysisFolders, createAndFormatExcelOutputFiles,
    createbulkRNASeqAnalysisBashScripts,
    createConcatenateFASTQfilesShellScript, createDdsObject,
    createDesignFileCrickASFsamples,
    create.gmt.file.from.ref.data.table, create.gsea.rnk.files,
    create.GSEA.table, createNewProject, createPCAloadingPlots,
    createRMDscript, createSampleDendrogram, createSettingsFile,
    createSeuratSettingsFile, createSRAdownloadScript,
    createSRRdownloadScript, create.timecourse.cat.lines,
    create.tpm.and.fpkm.tables, createTrimGaloreShellScript,
    create.website.parameters, datatable.to.website.ptm, dbcat2gmt,
    doDGEanalysis, doLinearFittings, doLRTanalysis, doPCA, DotPlotSB,
    excel, hpcClusterSubmisison, import.db.table.from.db,
    killDbConnections, list.db.table.col.names, list.db.tables.in.db,
    listExistingProjects, make.hm, msigdb.gmt2refDB,
    organizeFastqFiles, prepareRSEMcountMatrix, profileCluster,
    readAndPrepareCountMatrix, retrieve.gene.category.from.db,
    selectHeatmapGenes, setAnalysisPaths,
    setCrickGenomeAndGeneNameTable, setDataBaseParameters,
    setMountingPoint, upload.datatable.to.database,
    upload.pca.table.to.db

  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |.......                                                               |   9%
  ordinary text without R code

  |                                                                              |........                                                              |  11%
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


  |                                                                              |.........                                                             |  13%
  ordinary text without R code

  |                                                                              |...........                                                           |  15%
label: unnamed-chunk-3 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B.3.create.integrated.seurat.object.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B.3.create.integrated.seurat.object.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..............                                                        |  20%
  ordinary text without R code

  |                                                                              |............................                                          |  40%
label: partB_create_integrated_Seurat_object_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..........................................                            |  60%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: B3_partB_create_sample_list (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: symbol F

Computing 2000 integration features
Scaling features for provided objects
Finding all pairwise anchors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 13347 anchors
Filtering anchors
	Retained 4579 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 13318 anchors
Filtering anchors
	Retained 3693 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 11130 anchors
Filtering anchors
	Retained 7031 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 13347 anchors
Filtering anchors
	Retained 4579 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 12361 anchors
Filtering anchors
	Retained 9708 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 11135 anchors
Filtering anchors
	Retained 6452 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 13318 anchors
Filtering anchors
	Retained 3693 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 11130 anchors
Filtering anchors
	Retained 7031 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 11826 anchors
Filtering anchors
	Retained 9372 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 11130 anchors
Filtering anchors
	Retained 7031 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 8005 anchors
Filtering anchors
	Retained 4719 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7150 anchors
Filtering anchors
	Retained 5147 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7285 anchors
Filtering anchors
	Retained 4522 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7150 anchors
Filtering anchors
	Retained 5147 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 7285 anchors
Filtering anchors
	Retained 4522 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 11741 anchors
Filtering anchors
	Retained 3719 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 10223 anchors
Filtering anchors
	Retained 5112 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 10170 anchors
Filtering anchors
	Retained 4969 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 10223 anchors
Filtering anchors
	Retained 5112 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 10170 anchors
Filtering anchors
	Retained 4969 anchors
Extracting within-dataset neighbors
Running CCA
Merging objects
Finding neighborhoods
Finding anchors
	Found 6800 anchors
Filtering anchors
	Retained 2952 anchors
Extracting within-dataset neighbors
Merging dataset 5 into 3
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 4 into 2
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 6 into 2 4
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 3 5 into 2 4 6
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 7 into 2 4 6 3 5
Extracting anchors for merged samples
Finding integration vectors
Finding integration vector weights
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Integrating data
Merging dataset 1 into 2 4 6 3 5 7
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
  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.............                                                         |  19%
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
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
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
Finished averaging RNA for cluster Skin01sub
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
Finished averaging RNA for cluster ThyP89sub
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
Finished averaging RNA for cluster ThyP94sub
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
Finished averaging RNA for cluster MedullaThy89sub
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
Finished averaging RNA for cluster MedullaThy94sub
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
Finished averaging RNA for cluster CortexThy89sub
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   5%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  14%  |                                                                              |=============                                                         |  19%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  29%  |                                                                              |=======================                                               |  33%  |                                                                              |===========================                                           |  38%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  48%  |                                                                              |=====================================                                 |  52%  |                                                                              |========================================                              |  57%  |                                                                              |===========================================                           |  62%  |                                                                              |===============================================                       |  67%  |                                                                              |==================================================                    |  71%  |                                                                              |=====================================================                 |  76%  |                                                                              |=========================================================             |  81%  |                                                                              |============================================================          |  86%  |                                                                              |===============================================================       |  90%  |                                                                              |===================================================================   |  95%  |                                                                              |======================================================================| 100%
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
Finished averaging RNA for cluster CortexThy94sub
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%  |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
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


  |                                                                              |...............                                                       |  21%
  ordinary text without R code

  |                                                                              |................                                                      |  23%
label: unnamed-chunk-5 (with options) 
List of 2
 $ child: chr "src/modules/section_B/B.5.cluster.DGE.Rmd"
 $ eval : logi TRUE



processing file: ./src/modules/section_B/B.5.cluster.DGE.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..................                                                    |  25%
  ordinary text without R code

  |                                                                              |...................................                                   |  50%
label: partB_pca_plots (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |....................................................                  |  75%
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
16:41:04 UMAP embedding parameters a = 0.9922 b = 1.112
16:41:04 Read 23764 rows and found 30 numeric columns
16:41:04 Using Annoy for neighbor search, n_neighbors = 30
16:41:04 Building Annoy index with metric = cosine, n_trees = 50
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
16:41:09 Writing NN index file to temp file /tmp/RtmpgJKUOq/file84c94bbffd1d
16:41:09 Searching Annoy index using 1 thread, search_k = 3000
16:41:20 Annoy recall = 100%
16:41:21 Commencing smooth kNN distance calibration using 1 thread
16:41:24 Initializing from normalized Laplacian + noise
16:41:26 Commencing optimization for 200 epochs, with 1048306 positive edges
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
16:42:08 Optimization finished
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

Calculating cluster C1
Calculating cluster C5
Calculating cluster C9
Calculating cluster C3
Calculating cluster C2
Calculating cluster C0
Calculating cluster C8
Calculating cluster C6
Calculating cluster C4
Calculating cluster C10
Calculating cluster C7

  |                                                                              |.................                                                     |  25%
  ordinary text without R code

  |                                                                              |..................                                                    |  26%
label: unnamed-chunk-6 (with options) 
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

The following object is masked from 'package:IRanges':

    collapse

The following object is masked from 'package:S4Vectors':

    expand

The following object is masked from 'package:tidyr':

    expand

Finished averaging RNA for cluster Skin01sub
Finished averaging RNA for cluster ThyP89sub
Finished averaging RNA for cluster ThyP94sub
Finished averaging RNA for cluster MedullaThy89sub
Finished averaging RNA for cluster MedullaThy94sub
Finished averaging RNA for cluster CortexThy89sub
Finished averaging RNA for cluster CortexThy94sub
  |                                                                              |.............................................................         |  88%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: Plot_tsne_data_plotting (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |....................                                                  |  28%
  ordinary text without R code

  |                                                                              |.....................                                                 |  30%
label: unnamed-chunk-7 (with options) 
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
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
Finished averaging RNA for cluster C1_Skin01sub
Finished averaging RNA for cluster C5_Skin01sub
Finished averaging RNA for cluster C9_Skin01sub
Finished averaging RNA for cluster C3_Skin01sub
Finished averaging RNA for cluster C2_Skin01sub
Finished averaging RNA for cluster C0_Skin01sub
Finished averaging RNA for cluster C8_Skin01sub
Finished averaging RNA for cluster C6_Skin01sub
Finished averaging RNA for cluster C4_Skin01sub
Finished averaging RNA for cluster C10_Skin01sub
Finished averaging RNA for cluster C7_Skin01sub
Finished averaging RNA for cluster C6_ThyP89sub
Finished averaging RNA for cluster C5_ThyP89sub
Finished averaging RNA for cluster C1_ThyP89sub
Finished averaging RNA for cluster C4_ThyP89sub
Finished averaging RNA for cluster C0_ThyP89sub
Finished averaging RNA for cluster C7_ThyP89sub
Finished averaging RNA for cluster C3_ThyP89sub
Finished averaging RNA for cluster C2_ThyP89sub
Finished averaging RNA for cluster C9_ThyP89sub
Finished averaging RNA for cluster C8_ThyP89sub
Finished averaging RNA for cluster C10_ThyP89sub
Finished averaging RNA for cluster C0_ThyP94sub
Finished averaging RNA for cluster C2_ThyP94sub
Finished averaging RNA for cluster C6_ThyP94sub
Finished averaging RNA for cluster C1_ThyP94sub
Finished averaging RNA for cluster C5_ThyP94sub
Finished averaging RNA for cluster C3_ThyP94sub
Finished averaging RNA for cluster C4_ThyP94sub
Finished averaging RNA for cluster C7_ThyP94sub
Finished averaging RNA for cluster C9_ThyP94sub
Finished averaging RNA for cluster C6_MedullaThy89sub
Finished averaging RNA for cluster C5_MedullaThy89sub
Finished averaging RNA for cluster C1_MedullaThy89sub
Finished averaging RNA for cluster C4_MedullaThy89sub
Finished averaging RNA for cluster C0_MedullaThy89sub
Finished averaging RNA for cluster C7_MedullaThy89sub
Finished averaging RNA for cluster C3_MedullaThy89sub
Finished averaging RNA for cluster C2_MedullaThy89sub
Finished averaging RNA for cluster C9_MedullaThy89sub
Finished averaging RNA for cluster C8_MedullaThy89sub
Finished averaging RNA for cluster C10_MedullaThy89sub
Finished averaging RNA for cluster C0_MedullaThy94sub
Finished averaging RNA for cluster C2_MedullaThy94sub
Finished averaging RNA for cluster C6_MedullaThy94sub
Finished averaging RNA for cluster C1_MedullaThy94sub
Finished averaging RNA for cluster C5_MedullaThy94sub
Finished averaging RNA for cluster C3_MedullaThy94sub
Finished averaging RNA for cluster C4_MedullaThy94sub
Finished averaging RNA for cluster C7_MedullaThy94sub
Finished averaging RNA for cluster C9_MedullaThy94sub
Finished averaging RNA for cluster C4_CortexThy89sub
Finished averaging RNA for cluster C8_CortexThy89sub
Finished averaging RNA for cluster C7_CortexThy89sub
Finished averaging RNA for cluster C0_CortexThy89sub
Finished averaging RNA for cluster C2_CortexThy89sub
Finished averaging RNA for cluster C3_CortexThy89sub
Finished averaging RNA for cluster C9_CortexThy89sub
Finished averaging RNA for cluster C1_CortexThy89sub
Finished averaging RNA for cluster C5_CortexThy89sub
Finished averaging RNA for cluster C6_CortexThy89sub
Finished averaging RNA for cluster C8_CortexThy94sub
Finished averaging RNA for cluster C2_CortexThy94sub
Finished averaging RNA for cluster C4_CortexThy94sub
Finished averaging RNA for cluster C3_CortexThy94sub
Finished averaging RNA for cluster C0_CortexThy94sub
Finished averaging RNA for cluster C1_CortexThy94sub
Finished averaging RNA for cluster C5_CortexThy94sub
Finished averaging RNA for cluster C6_CortexThy94sub
Finished averaging RNA for cluster C7_CortexThy94sub
Finished averaging RNA for cluster C9_CortexThy94sub
Finished averaging RNA for cluster C10_CortexThy94sub
Finished averaging integrated for cluster C1_Skin01sub
Finished averaging integrated for cluster C5_Skin01sub
Finished averaging integrated for cluster C9_Skin01sub
Finished averaging integrated for cluster C3_Skin01sub
Finished averaging integrated for cluster C2_Skin01sub
Finished averaging integrated for cluster C0_Skin01sub
Finished averaging integrated for cluster C8_Skin01sub
Finished averaging integrated for cluster C6_Skin01sub
Finished averaging integrated for cluster C4_Skin01sub
Finished averaging integrated for cluster C10_Skin01sub
Finished averaging integrated for cluster C7_Skin01sub
Finished averaging integrated for cluster C6_ThyP89sub
Finished averaging integrated for cluster C5_ThyP89sub
Finished averaging integrated for cluster C1_ThyP89sub
Finished averaging integrated for cluster C4_ThyP89sub
Finished averaging integrated for cluster C0_ThyP89sub
Finished averaging integrated for cluster C7_ThyP89sub
Finished averaging integrated for cluster C3_ThyP89sub
Finished averaging integrated for cluster C2_ThyP89sub
Finished averaging integrated for cluster C9_ThyP89sub
Finished averaging integrated for cluster C8_ThyP89sub
Finished averaging integrated for cluster C10_ThyP89sub
Finished averaging integrated for cluster C0_ThyP94sub
Finished averaging integrated for cluster C2_ThyP94sub
Finished averaging integrated for cluster C6_ThyP94sub
Finished averaging integrated for cluster C1_ThyP94sub
Finished averaging integrated for cluster C5_ThyP94sub
Finished averaging integrated for cluster C3_ThyP94sub
Finished averaging integrated for cluster C4_ThyP94sub
Finished averaging integrated for cluster C7_ThyP94sub
Finished averaging integrated for cluster C9_ThyP94sub
Finished averaging integrated for cluster C6_MedullaThy89sub
Finished averaging integrated for cluster C5_MedullaThy89sub
Finished averaging integrated for cluster C1_MedullaThy89sub
Finished averaging integrated for cluster C4_MedullaThy89sub
Finished averaging integrated for cluster C0_MedullaThy89sub
Finished averaging integrated for cluster C7_MedullaThy89sub
Finished averaging integrated for cluster C3_MedullaThy89sub
Finished averaging integrated for cluster C2_MedullaThy89sub
Finished averaging integrated for cluster C9_MedullaThy89sub
Finished averaging integrated for cluster C8_MedullaThy89sub
Finished averaging integrated for cluster C10_MedullaThy89sub
Finished averaging integrated for cluster C0_MedullaThy94sub
Finished averaging integrated for cluster C2_MedullaThy94sub
Finished averaging integrated for cluster C6_MedullaThy94sub
Finished averaging integrated for cluster C1_MedullaThy94sub
Finished averaging integrated for cluster C5_MedullaThy94sub
Finished averaging integrated for cluster C3_MedullaThy94sub
Finished averaging integrated for cluster C4_MedullaThy94sub
Finished averaging integrated for cluster C7_MedullaThy94sub
Finished averaging integrated for cluster C9_MedullaThy94sub
Finished averaging integrated for cluster C4_CortexThy89sub
Finished averaging integrated for cluster C8_CortexThy89sub
Finished averaging integrated for cluster C7_CortexThy89sub
Finished averaging integrated for cluster C0_CortexThy89sub
Finished averaging integrated for cluster C2_CortexThy89sub
Finished averaging integrated for cluster C3_CortexThy89sub
Finished averaging integrated for cluster C9_CortexThy89sub
Finished averaging integrated for cluster C1_CortexThy89sub
Finished averaging integrated for cluster C5_CortexThy89sub
Finished averaging integrated for cluster C6_CortexThy89sub
Finished averaging integrated for cluster C8_CortexThy94sub
Finished averaging integrated for cluster C2_CortexThy94sub
Finished averaging integrated for cluster C4_CortexThy94sub
Finished averaging integrated for cluster C3_CortexThy94sub
Finished averaging integrated for cluster C0_CortexThy94sub
Finished averaging integrated for cluster C1_CortexThy94sub
Finished averaging integrated for cluster C5_CortexThy94sub
Finished averaging integrated for cluster C6_CortexThy94sub
Finished averaging integrated for cluster C7_CortexThy94sub
Finished averaging integrated for cluster C9_CortexThy94sub
Finished averaging integrated for cluster C10_CortexThy94sub
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
Finished averaging RNA for cluster C1
Finished averaging RNA for cluster C5
Finished averaging RNA for cluster C9
Finished averaging RNA for cluster C3
Finished averaging RNA for cluster C2
Finished averaging RNA for cluster C0
Finished averaging RNA for cluster C8
Finished averaging RNA for cluster C6
Finished averaging RNA for cluster C4
Finished averaging RNA for cluster C10
Finished averaging RNA for cluster C7
Finished averaging integrated for cluster C1
Finished averaging integrated for cluster C5
Finished averaging integrated for cluster C9
Finished averaging integrated for cluster C3
Finished averaging integrated for cluster C2
Finished averaging integrated for cluster C0
Finished averaging integrated for cluster C8
Finished averaging integrated for cluster C6
Finished averaging integrated for cluster C4
Finished averaging integrated for cluster C10
Finished averaging integrated for cluster C7
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
Finished averaging RNA for cluster Skin01sub
Finished averaging RNA for cluster ThyP89sub
Finished averaging RNA for cluster ThyP94sub
Finished averaging RNA for cluster MedullaThy89sub
Finished averaging RNA for cluster MedullaThy94sub
Finished averaging RNA for cluster CortexThy89sub
Finished averaging RNA for cluster CortexThy94sub
Finished averaging integrated for cluster Skin01sub
Finished averaging integrated for cluster ThyP89sub
Finished averaging integrated for cluster ThyP94sub
Finished averaging integrated for cluster MedullaThy89sub
Finished averaging integrated for cluster MedullaThy94sub
Finished averaging integrated for cluster CortexThy89sub
Finished averaging integrated for cluster CortexThy94sub
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
Finished averaging RNA for cluster Skin01sub
Finished averaging RNA for cluster ThyP89sub
Finished averaging RNA for cluster ThyP94sub
Finished averaging RNA for cluster MedullaThy89sub
Finished averaging RNA for cluster MedullaThy94sub
Finished averaging RNA for cluster CortexThy89sub
Finished averaging RNA for cluster CortexThy94sub

  |                                                                              |......................                                                |  32%
  ordinary text without R code

  |                                                                              |........................                                              |  34%
label: unnamed-chunk-8 (with options) 
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

Finished averaging RNA for cluster Skin01sub
Finished averaging RNA for cluster ThyP89sub
Finished averaging RNA for cluster ThyP94sub
Finished averaging RNA for cluster MedullaThy89sub
Finished averaging RNA for cluster MedullaThy94sub
Finished averaging RNA for cluster CortexThy89sub
Finished averaging RNA for cluster CortexThy94sub
Finished averaging RNA for cluster C1
Finished averaging RNA for cluster C5
Finished averaging RNA for cluster C9
Finished averaging RNA for cluster C3
Finished averaging RNA for cluster C2
Finished averaging RNA for cluster C0
Finished averaging RNA for cluster C8
Finished averaging RNA for cluster C6
Finished averaging RNA for cluster C4
Finished averaging RNA for cluster C10
Finished averaging RNA for cluster C7
  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: Plot_dim_red_by_cluster (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |.........................                                             |  36%
  ordinary text without R code

  |                                                                              |..........................                                            |  38%
label: unnamed-chunk-9 (with options) 
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


  |                                                                              |............................                                          |  40%
  ordinary text without R code

  |                                                                              |.............................                                         |  42%
label: unnamed-chunk-10 (with options) 
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


  |                                                                              |..............................                                        |  43%
  ordinary text without R code

  |                                                                              |................................                                      |  45%
label: unnamed-chunk-11 (with options) 
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
Finished averaging RNA for cluster Skin01sub
Finished averaging RNA for cluster ThyP89sub
Finished averaging RNA for cluster ThyP94sub
Finished averaging RNA for cluster MedullaThy89sub
Finished averaging RNA for cluster MedullaThy94sub
Finished averaging RNA for cluster CortexThy89sub
Finished averaging RNA for cluster CortexThy94sub
Finished averaging RNA for cluster C1
Finished averaging RNA for cluster C5
Finished averaging RNA for cluster C9
Finished averaging RNA for cluster C3
Finished averaging RNA for cluster C2
Finished averaging RNA for cluster C0
Finished averaging RNA for cluster C8
Finished averaging RNA for cluster C6
Finished averaging RNA for cluster C4
Finished averaging RNA for cluster C10
Finished averaging RNA for cluster C7
  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: Plot_dim_red_by_cell_cycle_phase (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |.................................                                     |  47%
  ordinary text without R code

  |                                                                              |..................................                                    |  49%
label: unnamed-chunk-12 (with options) 
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


  |                                                                              |....................................                                  |  51%
  ordinary text without R code

  |                                                                              |.....................................                                 |  53%
label: unnamed-chunk-13 (with options) 
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


  |                                                                              |......................................                                |  55%
  ordinary text without R code

  |                                                                              |........................................                              |  57%
label: unnamed-chunk-14 (with options) 
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


  |                                                                              |.........................................                             |  58%
  ordinary text without R code

  |                                                                              |..........................................                            |  60%
label: unnamed-chunk-15 (with options) 
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


  |                                                                              |............................................                          |  62%
  ordinary text without R code

  |                                                                              |.............................................                         |  64%
label: unnamed-chunk-16 (with options) 
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


  |                                                                              |..............................................                        |  66%
  ordinary text without R code

  |                                                                              |................................................                      |  68%
label: unnamed-chunk-17 (with options) 
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

The following object is masked from 'package:DelayedArray':

    simplify

The following object is masked from 'package:purrr':

    simplify

Loading required package: GSEABase
Loading required package: annotate
Loading required package: AnnotationDbi

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
--> Expected input gene ID: AQP9,TNS2,NOTCH4,FMO3,TRAC,MT1F
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: NTRK1,LTB,FAIM3,SLC15A3,CD27,PLA1A
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: CD22,RASSF4,CCL19,RAB27B,SLC12A8,ADAM28
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: MMP12,CHI3L1,FAM174B,KIAA0226L,DAPK2,CCL19
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: CPA3,AIF1,DENND5B,GPC4,CCL18,CD209
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: TNFRSF1B,AGT,CACYBP,CLU,F8,GZMA
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: CD180,DENND5B,KIAA0226L,TREM2,CCL8,NIPSNAP3B
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: CTSC,CASP8,PDIA2,GPLD1,TTR,CARF
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: SEPT8,ST3GAL6,FCGR2B,CD244,CD22,EMR2
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: RSAD2,GPR65,IL5,NCF2,SLC12A8,CSF2
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: PSG2,FAM198B,TRPM6,P2RX1,CCL20,CHI3L1
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


  |                                                                              |.................................................                     |  70%
  ordinary text without R code

  |                                                                              |..................................................                    |  72%
label: unnamed-chunk-18 (with options) 
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


  |                                                                              |....................................................                  |  74%
  ordinary text without R code

  |                                                                              |.....................................................                 |  75%
label: unnamed-chunk-19 (with options) 
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


  |                                                                              |......................................................                |  77%
  ordinary text without R code

  |                                                                              |.......................................................               |  79%
label: unnamed-chunk-20 (with options) 
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
--> Expected input gene ID: CHI3L1,CCL5,NCF2,SELL,MMP12,KIAA0226L
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: CLC,SLC12A8,CD33,CRISP3,FFAR2,FAM198B
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: TREM2,CCL5,IGSF6,GGT5,IFI44L,HDC
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: FLVCR2,PLEKHG3,HHEX,IL1B,KYNU,BACH2
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: LILRA3,PDE6C,TNFRSF17,PRG2,NOX3,TRPM6
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: P2RX1,STAP1,ACHE,CD38,CLIC2,BPI
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: CD19,CD1C,PPFIBP1,SELL,TMEM255A,KIAA0754
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: FRK,PLEKHG3,RGS1,FAM174B,RRP12,ADRB2
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: KIAA0754,AQP9,C3AR1,GPC4,BEND5,CCL19
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: STAP1,FES,TREM2,EGR2,FZD2,HPGDS
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: ZNF286A,CD72,LTB,MYB,FAIM3,SLC2A6
--> return NULL...
--> No gene can be mapped....
--> Expected input gene ID: EPN2,IL18RAP,IGLL3P,BACH2,PTPRCAP,TREM2
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


  |                                                                              |.........................................................             |  81%
  ordinary text without R code

  |                                                                              |..........................................................            |  83%
label: unnamed-chunk-21 (with options) 
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


  |                                                                              |...........................................................           |  85%
  ordinary text without R code

  |                                                                              |.............................................................         |  87%
label: unnamed-chunk-22 (with options) 
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


  |                                                                              |..............................................................        |  89%
  ordinary text without R code

  |                                                                              |...............................................................       |  91%
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

Project Database could not be reached or has no entry in Obio@parameterList$lims.id for this analysis.
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

/camp/apps/eb/dev/software/Anaconda2/2019.10/bin/pandoc +RTS -K512m -RTS sc_PartB_Analysis.utf8.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449D_subclustering_449C_SC20193_C1_2_3_4_5/workdir/html_local/sc_PartB_Analysis.html --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 5 --variable toc_float=1 --variable toc_selectors=h1,h2,h3,h4,h5 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --css src/assets/style/style.css --variable theme=paper --include-in-header /tmp/RtmpgJKUOq/rmarkdown-str84c93007f54a.html --mathjax --variable 'mathjax-url:https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449D_subclustering_449C_SC20193_C1_2_3_4_5/workdir/html_local/sc_PartB_Analysis.html
There were 23 warnings (use warnings() to see them)
