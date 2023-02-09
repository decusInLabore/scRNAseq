
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

# Bootstrapping renv 0.15.4-11 -----------------------------------------------
* Downloading renv 0.15.4-11 from GitHub ... OK
* Installing renv 0.15.4-11 ... Done!
* Successfully installed and loaded renv 0.15.4-11.
Installing BiocManager [1.30.19] ...
	OK [linked cache]
The following package(s) will be updated:

# Bioconductor =======================
- AUCell                 [* -> 1.12.0]
- AnnotationDbi          [* -> 1.52.0]
- Biobase                [* -> 2.50.0]
- BiocGenerics           [* -> 0.36.1]
- BiocParallel           [* -> 1.24.1]
- BiocVersion            [* -> 3.12.0]
- DESeq2                 [* -> 1.30.1]
- DO.db                  [* -> 2.9]
- DOSE                   [* -> 3.16.0]
- DelayedArray           [* -> 0.16.3]
- GO.db                  [* -> 3.12.1]
- GOSemSim               [* -> 2.16.1]
- GSEABase               [* -> 1.52.1]
- GenomeInfoDb           [* -> 1.26.7]
- GenomeInfoDbData       [* -> 1.2.4]
- GenomicRanges          [* -> 1.42.0]
- IRanges                [* -> 2.24.1]
- MatrixGenerics         [* -> 1.2.1]
- S4Vectors              [* -> 0.28.1]
- SummarizedExperiment   [* -> 1.20.0]
- XVector                [* -> 0.30.0]
- annotate               [* -> 1.68.0]
- clusterProfiler        [* -> 3.18.1]
- enrichplot             [* -> 1.10.2]
- fgsea                  [* -> 1.16.0]
- genefilter             [* -> 1.72.1]
- geneplotter            [* -> 1.68.0]
- graph                  [* -> 1.68.0]
- qvalue                 [* -> 2.22.0]
- zlibbioc               [* -> 1.36.0]

# CRAN ===============================
- BiocManager            [1.30.19 -> 1.30.16]
- DBI                    [1.1.0 -> 1.1.1]
- MASS                   [7.3-51.6 -> 7.3-54]
- Matrix                 [1.2-18 -> 1.3-4]
- R.methodsS3            [1.8.0 -> 1.8.1]
- R.oo                   [1.23.0 -> 1.24.0]
- R.utils                [2.9.2 -> 2.11.0]
- R6                     [2.4.1 -> 2.5.1]
- RCurl                  [1.98-1.2 -> 1.98-1.5]
- Rcpp                   [1.0.4.6 -> 1.0.9]
- ape                    [5.3 -> 5.5]
- callr                  [3.4.3 -> 3.7.0]
- cli                    [2.0.2 -> 3.1.0]
- clipr                  [0.7.0 -> 0.7.1]
- codetools              [0.2-16 -> 0.2-18]
- colorspace             [1.4-1 -> 2.0-2]
- cowplot                [1.0.0 -> 1.1.1]
- crayon                 [1.3.4 -> 1.4.2]
- curl                   [4.3 -> 4.3.2]
- deldir                 [0.1-25 -> 1.0-6]
- digest                 [0.6.25 -> 0.6.28]
- dplyr                  [0.8.5 -> 1.0.7]
- ellipsis               [0.3.0 -> 0.3.2]
- evaluate               [0.14 -> 0.15]
- fansi                  [0.4.1 -> 0.5.0]
- farver                 [2.0.3 -> 2.1.0]
- fastmap                [1.0.1 -> 1.1.0]
- fastmatch              [1.1-0 -> 1.1-3]
- formatR                [1.7 -> 1.11]
- fs                     [1.4.1 -> 1.5.0]
- future                 [1.17.0 -> 1.19.1]
- generics               [0.0.2 -> 0.1.1]
- ggforce                [0.3.1 -> 0.3.3]
- ggplot2                [3.3.0 -> 3.3.5]
- ggraph                 [2.0.2 -> 2.0.5]
- gh                     [1.1.0 -> 1.3.0]
- globals                [0.12.5 -> 0.14.0]
- glue                   [1.4.0 -> 1.5.0]
- graphlayouts           [0.7.0 -> 0.7.1]
- highr                  [0.8 -> 0.9]
- hms                    [0.5.3 -> 1.1.1]
- htmltools              [0.4.0 -> 0.5.2]
- httpuv                 [1.5.2 -> 1.6.3]
- httr                   [1.4.1 -> 1.4.2]
- isoband                [0.2.1 -> 0.2.5]
- jsonlite               [1.6.1 -> 1.7.2]
- knitr                  [1.28 -> 1.39]
- labeling               [0.3 -> 0.4.2]
- later                  [1.0.0 -> 1.3.0]
- lattice                [0.20-41 -> 0.20-45]
- lifecycle              [0.2.0 -> 1.0.1]
- magrittr               [1.5 -> 2.0.1]
- matrixStats            [0.56.0 -> 0.61.0]
- mgcv                   [1.8-31 -> 1.8-38]
- mime                   [0.9 -> 0.12]
- nlme                   [3.1-147 -> 3.1-153]
- openssl                [1.4.1 -> 1.4.5]
- pillar                 [1.4.3 -> 1.6.4]
- processx               [3.4.2 -> 3.5.2]
- promises               [1.1.0 -> 1.2.0.1]
- ps                     [1.3.2 -> 1.6.0]
- rappdirs               [0.3.1 -> 0.3.3]
- rcmdcheck              [1.3.3 -> 1.4.0]
- remotes                [2.1.1 -> 2.4.1]
- rlang                  [0.4.5 -> 0.4.12]
- rmarkdown              [2.1 -> 2.13]
- roxygen2               [7.1.0 -> 7.1.2]
- rprojroot              [1.3-2 -> 2.0.2]
- rstudioapi             [0.11 -> 0.13]
- scales                 [1.1.0 -> 1.1.1]
- sctransform            [0.2.1 -> 0.3.2]
- sessioninfo            [1.1.1 -> 1.2.1]
- shiny                  [1.4.0.2 -> 1.7.1]
- snow                   [0.4-3 -> 0.4-4]
- spatstat.data          [1.4-3 -> 2.1-0]
- spatstat.utils         [1.17-0 -> 2.2-0]
- stringi                [1.4.6 -> 1.7.5]
- sys                    [3.3 -> 3.4]
- tibble                 [3.0.1 -> 3.1.6]
- tidygraph              [1.1.2 -> 1.2.0]
- tidyr                  [1.0.2 -> 1.1.4]
- tidyselect             [1.0.0 -> 1.1.1]
- tinytex                [0.22 -> 0.34]
- tweenr                 [1.0.1 -> 1.0.2]
- usethis                [1.6.1 -> 2.1.3]
- utf8                   [1.1.4 -> 1.2.2]
- uuid                   [0.1-4 -> 1.0-3]
- vctrs                  [0.2.4 -> 0.3.8]
- viridisLite            [0.3.0 -> 0.4.0]
- withr                  [2.2.0 -> 2.4.2]
- xfun                   [0.13 -> 0.30]
- zip                    [2.0.4 -> 2.2.0]
- RMySQL                 [* -> 0.10.22]
- RSpectra               [* -> 0.16-0]
- RcppAnnoy              [* -> 0.0.19]
- Seurat                 [* -> 4.0.5]
- SeuratObject           [* -> 4.0.3]
- aplot                  [* -> 0.1.1]
- attempt                [* -> 0.3.1]
- bslib                  [* -> 0.3.1]
- cachem                 [* -> 1.0.6]
- config                 [* -> 0.3.1]
- cpp11                  [* -> 0.4.1]
- credentials            [* -> 1.3.1]
- dockerfiler            [* -> 0.1.4]
- dqrng                  [* -> 0.3.0]
- fontawesome            [* -> 0.2.2]
- gdtools                [* -> 0.2.3]
- gert                   [* -> 1.4.3]
- ggfun                  [* -> 0.0.6]
- ggplotify              [* -> 0.1.0]
- gitcreds               [* -> 0.1.1]
- golem                  [* -> 0.3.1]
- gridGraphics           [* -> 0.5-1]
- here                   [* -> 1.0.1]
- jquerylib              [* -> 0.1.4]
- jrc                    [* -> 0.5.0]
- officer                [* -> 0.4.1]
- rvg                    [* -> 0.2.5]
- sass                   [* -> 0.4.0]
- scattermore            [* -> 0.7]
- scatterpie             [* -> 0.1.7]
- shadowtext             [* -> 0.0.9]
- sitmo                  [* -> 2.0.2]
- sleepwalk              [* -> 0.3.2]
- spatstat.core          [* -> 2.3-1]
- spatstat.geom          [* -> 2.3-0]
- spatstat.sparse        [* -> 2.0-0]
- systemfonts            [* -> 1.0.3]
- thinkr                 [* -> 0.15]
- tidytree               [* -> 0.3.9]
- uwot                   [* -> 0.1.10]
- yulab.utils            [* -> 0.0.4]

# GitHub =============================
- renv                   [0.15.4-11 -> rstudio/renv@main]
- biologicSeqTools2      [* -> decusinlabore/biologicSeqTools2@main]
- biologicToolsSC        [* -> decusInLabore/biologicToolsSC@main]
- biologicViewerSC       [* -> decusInLabore/biologicViewerSC@main]
- ggsankey               [* -> davidsjoberg/ggsankey@main]
- ggtree                 [* -> YuLab-SMU/ggtree]
- rvcheck                [* -> GuangchuangYu/rvcheck]
- treeio                 [* -> GuangchuangYu/treeio]

Installing BiocManager [1.30.16] ...
	OK [linked cache]
Installing renv [0.15.4-11] ...
	OK [linked cache]
Installing lattice [0.20-45] ...
	OK [linked cache]
Installing nlme [3.1-153] ...
	OK [linked cache]
Installing Rcpp [1.0.9] ...
	OK [linked cache]
Installing ape [5.5] ...
	OK [linked cache]
Installing ps [1.6.0] ...
	OK [linked cache]
Installing R6 [2.5.1] ...
	OK [linked cache]
Installing processx [3.5.2] ...
	OK [linked cache]
Installing callr [3.7.0] ...
	OK [linked cache]
Installing glue [1.5.0] ...
	OK [linked cache]
Installing cli [3.1.0] ...
	OK [linked cache]
Installing clipr [0.7.1] ...
	OK [linked cache]
Installing codetools [0.2-18] ...
	OK [linked cache]
Installing colorspace [2.0-2] ...
	OK [linked cache]
Installing digest [0.6.28] ...
	OK [linked cache]
Installing isoband [0.2.5] ...
	OK [linked cache]
Installing MASS [7.3-54] ...
	OK [linked cache]
Installing Matrix [1.3-4] ...
	OK [linked cache]
Installing mgcv [1.8-38] ...
	OK [linked cache]
Installing rlang [0.4.12] ...
	OK [linked cache]
Installing farver [2.1.0] ...
	OK [linked cache]
Installing labeling [0.4.2] ...
	OK [linked cache]
Installing lifecycle [1.0.1] ...
	OK [linked cache]
Installing viridisLite [0.4.0] ...
	OK [linked cache]
Installing scales [1.1.1] ...
	OK [linked cache]
Installing ellipsis [0.3.2] ...
	OK [linked cache]
Installing fansi [0.5.0] ...
	OK [linked cache]
Installing magrittr [2.0.1] ...
	OK [linked cache]
Installing crayon [1.4.2] ...
	OK [linked cache]
Installing utf8 [1.2.2] ...
	OK [linked cache]
Installing vctrs [0.3.8] ...
	OK [linked cache]
Installing pillar [1.6.4] ...
	OK [linked cache]
Installing tibble [3.1.6] ...
	OK [linked cache]
Installing withr [2.4.2] ...
	OK [linked cache]
Installing ggplot2 [3.3.5] ...
	OK [linked cache]
Installing cowplot [1.1.1] ...
	OK [linked cache]
Installing curl [4.3.2] ...
	OK [linked cache]
Installing DBI [1.1.1] ...
	OK [linked cache]
Installing deldir [1.0-6] ...
	OK [linked cache]
Installing generics [0.1.1] ...
	OK [linked cache]
Installing tidyselect [1.1.1] ...
	OK [linked cache]
Installing dplyr [1.0.7] ...
	OK [linked cache]
Installing evaluate [0.15] ...
	OK [linked cache]
Installing fastmap [1.1.0] ...
	OK [linked cache]
Installing fastmatch [1.1-3] ...
	OK [linked cache]
Installing formatR [1.11] ...
	OK [linked cache]
Installing fs [1.5.0] ...
	OK [linked cache]
Installing globals [0.14.0] ...
	OK [linked cache]
Installing future [1.19.1] ...
	OK [linked cache]
Installing tweenr [1.0.2] ...
	OK [linked cache]
Installing ggforce [0.3.3] ...
	OK [linked cache]
Installing cpp11 [0.4.1] ...
	OK [linked cache]
Installing tidyr [1.1.4] ...
	OK [linked cache]
Installing tidygraph [1.2.0] ...
	OK [linked cache]
Installing graphlayouts [0.7.1] ...
	OK [linked cache]
Installing ggraph [2.0.5] ...
	OK [linked cache]
Installing gitcreds [0.1.1] ...
	OK [linked cache]
Installing jsonlite [1.7.2] ...
	OK [linked cache]
Installing mime [0.12] ...
	OK [linked cache]
Installing sys [3.4] ...
	OK [linked cache]
Installing openssl [1.4.5] ...
	OK [linked cache]
Installing httr [1.4.2] ...
	OK [linked cache]
Installing gh [1.3.0] ...
	OK [linked cache]
Installing xfun [0.30] ...
	OK [linked cache]
Installing highr [0.9] ...
	OK [linked cache]
Installing hms [1.1.1] ...
	OK [linked cache]
Installing htmltools [0.5.2] ...
	OK [linked cache]
Installing later [1.3.0] ...
	OK [linked cache]
Installing promises [1.2.0.1] ...
	OK [linked cache]
Installing httpuv [1.6.3] ...
	OK [linked cache]
Installing stringi [1.7.5] ...
	OK [linked cache]
Installing knitr [1.39] ...
	OK [linked cache]
Installing matrixStats [0.61.0] ...
	OK [linked cache]
Installing R.methodsS3 [1.8.1] ...
	OK [linked cache]
Installing R.oo [1.24.0] ...
	OK [linked cache]
Installing R.utils [2.11.0] ...
	OK [linked cache]
Installing rappdirs [0.3.3] ...
	OK [linked cache]
Installing rprojroot [2.0.2] ...
	OK [linked cache]
Installing sessioninfo [1.2.1] ...
	OK [linked cache]
Installing rcmdcheck [1.4.0] ...
	OK [linked cache]
Installing RCurl [1.98-1.5] ...
	OK [linked cache]
Installing remotes [2.4.1] ...
	OK [linked cache]
Installing sass [0.4.0] ...
	OK [linked cache]
Installing jquerylib [0.1.4] ...
	OK [linked cache]
Installing bslib [0.3.1] ...
	OK [linked cache]
Installing tinytex [0.34] ...
	OK [linked cache]
Installing rmarkdown [2.13] ...
	OK [linked cache]
Installing rstudioapi [0.13] ...
	OK [linked cache]
Installing roxygen2 [7.1.2] ...
	OK [linked cache]
Installing sctransform [0.3.2] ...
	OK [linked cache]
Installing fontawesome [0.2.2] ...
	OK [linked cache]
Installing cachem [1.0.6] ...
	OK [linked cache]
Installing shiny [1.7.1] ...
	OK [linked cache]
Installing snow [0.4-4] ...
	OK [linked cache]
Installing spatstat.utils [2.2-0] ...
	OK [linked cache]
Installing spatstat.data [2.1-0] ...
	OK [linked cache]
Installing credentials [1.3.1] ...
	OK [linked cache]
Installing zip [2.2.0] ...
	OK [linked cache]
Installing gert [1.4.3] ...
	OK [linked cache]
Installing usethis [2.1.3] ...
	OK [linked cache]
Installing uuid [1.0-3] ...
	OK [linked cache]
Installing BiocGenerics [0.36.1] ...
	OK [linked cache]
Installing S4Vectors [0.28.1] ...
	OK [linked cache]
Installing Biobase [2.50.0] ...
	OK [linked cache]
Installing IRanges [2.24.1] ...
	OK [linked cache]
Installing AnnotationDbi [1.52.0] ...
	OK [linked cache]
Installing annotate [1.68.0] ...
	OK [linked cache]
Installing graph [1.68.0] ...
	OK [linked cache]
Installing GSEABase [1.52.1] ...
	OK [linked cache]
Installing GenomeInfoDbData [1.2.4] ...
	OK [linked cache]
Installing GenomeInfoDb [1.26.7] ...
	OK [linked cache]
Installing MatrixGenerics [1.2.1] ...
	OK [linked cache]
Installing DelayedArray [0.16.3] ...
	OK [linked cache]
Installing zlibbioc [1.36.0] ...
	OK [linked cache]
Installing XVector [0.30.0] ...
	OK [linked cache]
Installing GenomicRanges [1.42.0] ...
	OK [linked cache]
Installing SummarizedExperiment [1.20.0] ...
	OK [linked cache]
Installing AUCell [1.12.0] ...
	OK [linked cache]
Installing BiocParallel [1.24.1] ...
	OK [linked cache]
Installing BiocVersion [3.12.0] ...
	OK [linked cache]
Installing genefilter [1.72.1] ...
	OK [linked cache]
Installing geneplotter [1.68.0] ...
	OK [linked cache]
Installing DESeq2 [1.30.1] ...
	OK [linked cache]
Installing DO.db [2.9] ...
	OK [linked cache]
Installing fgsea [1.16.0] ...
	OK [linked cache]
Installing GO.db [3.12.1] ...
	OK [linked cache]
Installing GOSemSim [2.16.1] ...
	OK [linked cache]
Installing qvalue [2.22.0] ...
	OK [linked cache]
Installing DOSE [3.16.0] ...
	OK [linked cache]
Installing RMySQL [0.10.22] ...
	OK [linked cache]
Installing RSpectra [0.16-0] ...
	OK [linked cache]
Installing RcppAnnoy [0.0.19] ...
	OK [linked cache]
Installing scattermore [0.7] ...
	OK [linked cache]
Installing SeuratObject [4.0.3] ...
	OK [linked cache]
Installing spatstat.sparse [2.0-0] ...
	OK [linked cache]
Installing spatstat.geom [2.3-0] ...
	OK [linked cache]
Installing spatstat.core [2.3-1] ...
	OK [linked cache]
Installing sitmo [2.0.2] ...
	OK [linked cache]
Installing dqrng [0.3.0] ...
	OK [linked cache]
Installing uwot [0.1.10] ...
	OK [linked cache]
Installing Seurat [4.0.5] ...
	OK [linked cache]
Installing ggfun [0.0.6] ...
	OK [linked cache]
Installing gridGraphics [0.5-1] ...
	OK [linked cache]
Installing yulab.utils [0.0.4] ...
	OK [linked cache]
Installing ggplotify [0.1.0] ...
	OK [linked cache]
Installing aplot [0.1.1] ...
	OK [linked cache]
Installing attempt [0.3.1] ...
	OK [linked cache]
Installing biologicSeqTools2 [0.0.0.9000] ...
	OK [linked cache]
Installing biologicToolsSC [0.0.0.9000] ...
	OK [linked cache]
Installing config [0.3.1] ...
	OK [linked cache]
Installing dockerfiler [0.1.4] ...
	OK [linked cache]
Installing here [1.0.1] ...
	OK [linked cache]
Installing golem [0.3.1] ...
	OK [linked cache]
Installing officer [0.4.1] ...
	OK [linked cache]
Installing systemfonts [1.0.3] ...
	OK [linked cache]
Installing gdtools [0.2.3] ...
	OK [linked cache]
Installing rvg [0.2.5] ...
	OK [linked cache]
Installing thinkr [0.15] ...
	OK [linked cache]
Installing biologicViewerSC [0.0.1.0] ...
	OK [linked cache]
Installing scatterpie [0.1.7] ...
	OK [linked cache]
Installing shadowtext [0.0.9] ...
	OK [linked cache]
Installing enrichplot [1.10.2] ...
	OK [linked cache]
Installing rvcheck [0.2.1] ...
	OK [linked cache]
Installing clusterProfiler [3.18.1] ...
	OK [linked cache]
Installing ggsankey [0.0.99999] ...
	OK [linked cache]
Installing tidytree [0.3.9] ...
	OK [linked cache]
Installing treeio [1.19.1] ...
	OK [linked cache]
Installing ggtree [3.5.0.900] ...
	OK [linked cache]
Installing jrc [0.5.0] ...
	OK [linked cache]
Installing sleepwalk [0.3.2] ...
	OK [linked cache]
The following package(s) have been updated:

	BiocManager [installed version 1.30.16 != loaded version 1.30.19]

Consider restarting the R session and loading the newly-installed packages.

Collecting package metadata (current_repodata.json): ...working... done
Solving environment: ...working... done

## Package Plan ##

  environment location: /camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/533_scRNAseq_subclustering_of_496_C7_C11_C12_subset3/scripts/scRNAseq/analyses/Main_Analysis/renv/python/r-reticulate

  added / updated specs:
    - numpy
    - python=3.6


The following NEW packages will be INSTALLED:

  _libgcc_mutex      conda-forge/linux-64::_libgcc_mutex-0.1-conda_forge
  _openmp_mutex      conda-forge/linux-64::_openmp_mutex-4.5-2_gnu
  ca-certificates    conda-forge/linux-64::ca-certificates-2022.12.7-ha878542_0
  certifi            conda-forge/linux-64::certifi-2016.9.26-py36_0
  ld_impl_linux-64   conda-forge/linux-64::ld_impl_linux-64-2.40-h41732ed_0
  libblas            conda-forge/linux-64::libblas-3.9.0-16_linux64_openblas
  libcblas           conda-forge/linux-64::libcblas-3.9.0-16_linux64_openblas
  libffi             conda-forge/linux-64::libffi-3.4.2-h7f98852_5
  libgcc-ng          conda-forge/linux-64::libgcc-ng-12.2.0-h65d4601_19
  libgfortran-ng     conda-forge/linux-64::libgfortran-ng-12.2.0-h69a702a_19
  libgfortran5       conda-forge/linux-64::libgfortran5-12.2.0-h337968e_19
  libgomp            conda-forge/linux-64::libgomp-12.2.0-h65d4601_19
  liblapack          conda-forge/linux-64::liblapack-3.9.0-16_linux64_openblas
  libnsl             conda-forge/linux-64::libnsl-2.0.0-h7f98852_0
  libopenblas        conda-forge/linux-64::libopenblas-0.3.21-pthreads_h78a6416_3
  libsqlite          conda-forge/linux-64::libsqlite-3.40.0-h753d276_0
  libstdcxx-ng       conda-forge/linux-64::libstdcxx-ng-12.2.0-h46fd767_19
  libzlib            conda-forge/linux-64::libzlib-1.2.13-h166bdaf_4
  ncurses            conda-forge/linux-64::ncurses-6.3-h27087fc_1
  numpy              conda-forge/linux-64::numpy-1.19.5-py36hfc0c790_2
  openssl            conda-forge/linux-64::openssl-1.1.1s-h0b41bf4_1
  pip                conda-forge/linux-64::pip-20.0.2-py36_1
  python             conda-forge/linux-64::python-3.6.15-hb7a2778_0_cpython
  python_abi         conda-forge/linux-64::python_abi-3.6-2_cp36m
  readline           conda-forge/linux-64::readline-8.1.2-h0f457ee_0
  setuptools         conda-forge/linux-64::setuptools-49.6.0-py36h5fab9bb_3
  sqlite             conda-forge/linux-64::sqlite-3.40.0-h4ff8645_0
  tk                 conda-forge/linux-64::tk-8.6.12-h27826a3_0
  wheel              conda-forge/linux-64::wheel-0.34.2-py36_0
  xz                 conda-forge/linux-64::xz-5.2.6-h166bdaf_0


Preparing transaction: ...working... done
Verifying transaction: ...working... done
Executing transaction: ...working... done


==> WARNING: A newer version of conda exists. <==
  current version: 4.8.3
  latest version: 22.11.1

Please update conda by running

    $ conda update -n base conda


#
# To activate this environment, use
#
#     $ conda activate /camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/533_scRNAseq_subclustering_of_496_C7_C11_C12_subset3/scripts/scRNAseq/analyses/Main_Analysis/renv/python/r-reticulate
#
# To deactivate an active environment, use
#
#     $ conda deactivate

Attaching SeuratObject
[1] "Sample in group Myc4W:Myc4W_A5, Myc4W_B5, Myc4W_A11, Myc4W_B11, Myc4W_C1, Myc4W_D1, Myc4W_A1, Myc4W_B1"
[1] "Myc4W: 2694"
[1] "Processed."
[1] "Sample in group Myc4W_2hrsOff:Myc4W_2hrsOff_C11, Myc4W_2hrsOff_D11, Myc4W_2hrsOff_E11, Myc4W_2hrsOff_C1, Myc4W_2hrsOff_D1, Myc4W_2hrsOff_E1"
[1] "Myc4W_2hrsOff: 908"
[1] "Processed."
[1] "Sample in group Myc4W_4hrsOff:Myc4W_4hrsOff_F11, Myc4W_4hrsOff_G11, Myc4W_4hrsOff_H11, Myc4W_4hrsOff_F1, Myc4W_4hrsOff_G1, Myc4W_4hrsOff_H1"
[1] "Myc4W_4hrsOff: 786"
[1] "Processed."
[1] "Sample in group Myc4W_8hrsOff:Myc4W_8hrsOff_C5, Myc4W_8hrsOff_D5, Myc4W_8hrsOff_A1, Myc4W_8hrsOff_B1"
[1] "Myc4W_8hrsOff: 6915"
[1] "Processed."
[1] "Sample in group Myc4W_12hrsOff:Myc4W_12hrsOff_E5, Myc4W_12hrsOff_F5, Myc4W_12hrsOff_G1, Myc4W_12hrsOff_H1"
[1] "Myc4W_12hrsOff: 3911"
[1] "Processed."
[1] "Sample in group Myc4W_24hrsOff:Myc4W_24hrsOff_G5, Myc4W_24hrsOff_H5, Myc4W_24hrsOff_E1, Myc4W_24hrsOff_F1"
[1] "Myc4W_24hrsOff: 1152"
[1] "Processed."
