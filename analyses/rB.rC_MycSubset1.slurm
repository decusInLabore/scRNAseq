
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



processing file: Main_Analysis.Rmd
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
 $ child: chr "01-set.parameters.Rmd"
 $ eval : logi TRUE



processing file: ./01-set.parameters.Rmd
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


  |                                                                              |.......                                                               |   9%
  ordinary text without R code

  |                                                                              |........                                                              |  11%
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


  |                                                                              |.........                                                             |  13%
  ordinary text without R code

  |                                                                              |...........                                                           |  15%
label: unnamed-chunk-3 (with options) 
List of 2
 $ child: chr "03-create.integrated.seurat.object.Rmd"
 $ eval : logi FALSE

  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.............                                                         |  19%
label: unnamed-chunk-4 (with options) 
List of 2
 $ eval: logi FALSE
 $ echo: logi FALSE

  |                                                                              |...............                                                       |  21%
  ordinary text without R code

  |                                                                              |................                                                      |  23%
label: unnamed-chunk-5 (with options) 
List of 2
 $ child: chr "04-create.variation.plots.Rmd"
 $ eval : logi FALSE

  |                                                                              |.................                                                     |  25%
  ordinary text without R code

  |                                                                              |..................                                                    |  26%
label: unnamed-chunk-6 (with options) 
List of 2
 $ child: chr "05-create.pca.plots.Rmd"
 $ eval : logi FALSE

  |                                                                              |....................                                                  |  28%
  ordinary text without R code

  |                                                                              |.....................                                                 |  30%
label: unnamed-chunk-7 (with options) 
List of 2
 $ eval: logi TRUE
 $ echo: logi TRUE

  |                                                                              |......................                                                |  32%
  ordinary text without R code

  |                                                                              |........................                                              |  34%
label: unnamed-chunk-8 (with options) 
List of 2
 $ child: chr "06-calculate.average.expression.Rmd"
 $ eval : logi TRUE



processing file: ./06-calculate.average.expression.Rmd
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

The following grouping variables have 1 value and will be ignored: ident
All grouping variables have 1 value only. Computing across all cells.
Warning: The following arguments are not used: row.names
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   9%  |                                                                              |=========                                                             |  13%  |                                                                              |============                                                          |  17%  |                                                                              |===============                                                       |  22%  |                                                                              |==================                                                    |  26%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  48%  |                                                                              |=====================================                                 |  52%  |                                                                              |========================================                              |  57%  |                                                                              |===========================================                           |  61%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  78%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  87%  |                                                                              |================================================================      |  91%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
Warning: The following arguments are not used: row.names
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   9%  |                                                                              |=========                                                             |  13%  |                                                                              |============                                                          |  17%  |                                                                              |===============                                                       |  22%  |                                                                              |==================                                                    |  26%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  48%  |                                                                              |=====================================                                 |  52%  |                                                                              |========================================                              |  57%  |                                                                              |===========================================                           |  61%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  78%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  87%  |                                                                              |================================================================      |  91%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
Warning: The following arguments are not used: row.names
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   9%  |                                                                              |=========                                                             |  13%  |                                                                              |============                                                          |  17%  |                                                                              |===============                                                       |  22%  |                                                                              |==================                                                    |  26%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  48%  |                                                                              |=====================================                                 |  52%  |                                                                              |========================================                              |  57%  |                                                                              |===========================================                           |  61%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  78%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  87%  |                                                                              |================================================================      |  91%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
Warning: The following arguments are not used: row.names
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   9%  |                                                                              |=========                                                             |  13%  |                                                                              |============                                                          |  17%  |                                                                              |===============                                                       |  22%  |                                                                              |==================                                                    |  26%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |===========================                                           |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |=================================                                     |  48%  |                                                                              |=====================================                                 |  52%  |                                                                              |========================================                              |  57%  |                                                                              |===========================================                           |  61%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  74%  |                                                                              |=======================================================               |  78%  |                                                                              |==========================================================            |  83%  |                                                                              |=============================================================         |  87%  |                                                                              |================================================================      |  91%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%

  |                                                                              |.........................                                             |  36%
  ordinary text without R code

  |                                                                              |..........................                                            |  38%
label: unnamed-chunk-9 (with options) 
List of 2
 $ child: chr "07-plot.pca.and.tsne.results.Rmd"
 $ eval : logi TRUE



processing file: ./07-plot.pca.and.tsne.results.Rmd
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


Attaching package: 'scales'

The following object is masked from 'package:purrr':

    discard

The following object is masked from 'package:readr':

    col_factor

ggtree v3.5.0.900 For help: https://yulab-smu.top/treedata-book/

If you use the ggtree package suite in published research, please cite
the appropriate paper(s):

Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
ggtree: an R package for visualization and annotation of phylogenetic
trees with their covariates and other associated data. Methods in
Ecology and Evolution. 2017, 8(1):28-36. doi:10.1111/2041-210X.12628

S Xu, Z Dai, P Guo, X Fu, S Liu, L Zhou, W Tang, T Feng, M Chen, L
Zhan, T Wu, E Hu, Y Jiang, X Bo, G Yu. ggtreeExtra: Compact
visualization of richly annotated phylogenetic data. Molecular Biology
and Evolution. 2021, 38(9):4039-4042. doi: 10.1093/molbev/msab166

Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods
for mapping and visualizing associated data on phylogeny using ggtree.
Molecular Biology and Evolution. 2018, 35(12):3041-3043.
doi:10.1093/molbev/msy194

Attaching package: 'ggtree'

The following object is masked from 'package:tidyr':

    expand

  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: Plot_dim_red_by_cluster (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |............................                                          |  40%
  ordinary text without R code

  |                                                                              |.............................                                         |  42%
label: unnamed-chunk-10 (with options) 
List of 2
 $ child: chr "08-euclidian.distances.figure.Rmd"
 $ eval : logi TRUE



processing file: ./08-euclidian.distances.figure.Rmd
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


  |                                                                              |..............................                                        |  43%
  ordinary text without R code

  |                                                                              |................................                                      |  45%
label: unnamed-chunk-11 (with options) 
List of 2
 $ child: chr "09-doublet.estimation.plots.Rmd"
 $ eval : logi TRUE



processing file: ./09-doublet.estimation.plots.Rmd
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


  |                                                                              |.................................                                     |  47%
  ordinary text without R code

  |                                                                              |..................................                                    |  49%
label: unnamed-chunk-12 (with options) 
List of 2
 $ child: chr "10-cell.cycle.phase.plots.Rmd"
 $ eval : logi TRUE



processing file: ./10-cell.cycle.phase.plots.Rmd
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

Warning: The following features are not present in the object: MCM5, PCNA, TYMS, FEN1, MCM2, MCM4, RRM1, UNG, GINS2, MCM6, CDCA7, DTL, PRIM1, UHRF1, MLF1IP, HELLS, RFC2, RPA2, NASP, RAD51AP1, GMNN, WDR76, SLBP, CCNE2, UBR7, POLD3, MSH2, ATAD2, RAD51, RRM2, CDC45, CDC6, EXO1, TIPIN, DSCC1, BLM, CASP8AP2, USP1, CLSPN, POLA1, CHAF1B, BRIP1, E2F8, not searching for symbol synonyms
Warning: The following features are not present in the object: HMGB2, CDK1, NUSAP1, UBE2C, BIRC5, TPX2, TOP2A, NDC80, CKS2, NUF2, CKS1B, MKI67, TMPO, CENPF, TACC3, FAM64A, SMC4, CCNB2, CKAP2L, CKAP2, AURKB, BUB1, KIF11, ANP32E, TUBB4B, GTSE1, KIF20B, HJURP, CDCA3, HN1, CDC20, TTK, CDC25C, KIF2C, RANGAP1, NCAPD2, DLGAP5, CDCA2, CDCA8, ECT2, KIF23, HMMR, AURKA, PSRC1, ANLN, LBR, CKAP5, CENPE, CTCF, NEK2, G2E3, GAS2L3, CBX5, CENPA, not searching for symbol synonyms
  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: Plot_dim_red_by_cell_cycle_phase (with options) 
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
 $ child: chr "11-Ncell.barchart.plots.Rmd"
 $ eval : logi TRUE



processing file: ./11-Ncell.barchart.plots.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |............                                                          |  17%
  ordinary text without R code

  |                                                                              |.......................                                               |  33%
label: unnamed-chunk-27 (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |...................................                                   |  50%
  ordinary text without R code

  |                                                                              |...............................................                       |  67%
label: unnamed-chunk-28 (with options) 
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


  |                                                                              |......................................                                |  55%
  ordinary text without R code

  |                                                                              |........................................                              |  57%
label: unnamed-chunk-14 (with options) 
List of 2
 $ child: chr "B133.timeseries.area.plot.Rmd"
 $ eval : logi TRUE



processing file: ./B133.timeseries.area.plot.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |..........                                                            |  14%
  ordinary text without R code

  |                                                                              |....................                                                  |  29%
label: create_sankey_init (with options) 
List of 3
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE

  |                                                                              |..............................                                        |  43%
  ordinary text without R code

  |                                                                              |........................................                              |  57%
label: unnamed-chunk-29 (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: logi FALSE

`summarise()` has grouped output by 'timepoint'. You can override using the `.groups` argument.
  |                                                                              |..................................................                    |  71%
   inline R code fragments

  |                                                                              |............................................................          |  86%
label: unnamed-chunk-30 (with options) 
List of 4
 $ echo   : logi TRUE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |......................................................................| 100%
  ordinary text without R code


  |                                                                              |.........................................                             |  58%
  ordinary text without R code

  |                                                                              |..........................................                            |  60%
label: unnamed-chunk-15 (with options) 
List of 2
 $ child: chr "12-AUC.cat.enrichments.Rmd"
 $ eval : logi TRUE



processing file: ./12-AUC.cat.enrichments.Rmd
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
 $ results: logi FALSE

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


  |                                                                              |............................................                          |  62%
  ordinary text without R code

  |                                                                              |.............................................                         |  64%
label: unnamed-chunk-16 (with options) 
List of 2
 $ child: chr "13-linear.regression.data.prep.Rmd"
 $ eval : logi TRUE



processing file: ./13-linear.regression.data.prep.Rmd
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
 $ results: logi FALSE

  |                                                                              |..........................................................            |  83%
   inline R code fragments

  |                                                                              |......................................................................| 100%
label: plot_linear_regression (with options) 
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
 $ child: chr "14-heatmaps.and.dotplots.Rmd"
 $ eval : logi TRUE



processing file: ./14-heatmaps.and.dotplots.Rmd
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
 $ results: logi FALSE

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


  |                                                                              |.................................................                     |  70%
  ordinary text without R code

  |                                                                              |..................................................                    |  72%
label: unnamed-chunk-18 (with options) 
List of 2
 $ child: chr "15-make.cluster.table.Rmd"
 $ eval : logi FALSE

  |                                                                              |....................................................                  |  74%
  ordinary text without R code

  |                                                                              |.....................................................                 |  75%
label: unnamed-chunk-19 (with options) 
List of 2
 $ child: chr "16-gene.set.heatmaps.Rmd"
 $ eval : logi TRUE



processing file: ./16-gene.set.heatmaps.Rmd
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
 $ echo   : logi FALSE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"


  |                                                                              |......................................................                |  77%
  ordinary text without R code

  |                                                                              |.......................................................               |  79%
label: unnamed-chunk-20 (with options) 
List of 2
 $ child: chr "17-example.feature.plots.Rmd"
 $ eval : logi TRUE



processing file: ./17-example.feature.plots.Rmd
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


  |                                                                              |.........................................................             |  81%
  ordinary text without R code

  |                                                                              |..........................................................            |  83%
label: unnamed-chunk-21 (with options) 
List of 2
 $ child: chr "18-pca.cell.and.loading.plots.Rmd"
 $ eval : logi FALSE

  |                                                                              |...........................................................           |  85%
  ordinary text without R code

  |                                                                              |.............................................................         |  87%
label: unnamed-chunk-22 (with options) 
List of 2
 $ child: chr "19-add.extra.tables.Rmd"
 $ eval : logi FALSE

  |                                                                              |..............................................................        |  89%
  ordinary text without R code

  |                                                                              |...............................................................       |  91%
label: saveobject (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi TRUE
 $ results: logi FALSE

  |                                                                              |.................................................................     |  92%
  ordinary text without R code

  |                                                                              |..................................................................    |  94%
label: create_report_params (with options) 
List of 2
 $ eval   : logi TRUE
 $ results: chr "asis"

  |                                                                              |...................................................................   |  96%
  ordinary text without R code

  |                                                                              |..................................................................... |  98%
label: documentation (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi TRUE
 $ results: logi TRUE

  |                                                                              |......................................................................| 100%
   inline R code fragments


output file: Main_Analysis.knit.md

/camp/apps/eb/software/pandoc/2.2.3.2-foss-2016b/bin/pandoc +RTS -K512m -RTS Main_Analysis.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/528_scRNAseq_subclustering_of_496_C7_C12_subset1/html_local/Main_Analysis.html --lua-filter /camp/stp/babs/working/boeings/Projects/home_extra/R/x86_64-pc-linux-gnu-library/4.0/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/stp/babs/working/boeings/Projects/home_extra/R/x86_64-pc-linux-gnu-library/4.0/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 5 --variable toc_float=1 --variable toc_selectors=h1,h2,h3,h4,h5 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/stp/babs/working/boeings/Projects/home_extra/R/x86_64-pc-linux-gnu-library/4.0/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --variable theme=paper --css src/style/style.css --mathjax --variable 'mathjax-url=https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --include-in-header /tmp/slurm_47148333/Rtmp6fiaqH/rmarkdown-str13844243b6078.html --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/528_scRNAseq_subclustering_of_496_C7_C12_subset1/html_local/Main_Analysis.html
There were 15 warnings (use warnings() to see them)
