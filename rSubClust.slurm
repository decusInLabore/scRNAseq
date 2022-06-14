# Bootstrapping renv 0.15.4-11 -----------------------------------------------
* Downloading renv 0.15.4-11 from GitHub ... OK
* Installing renv 0.15.4-11 ... Done!
* Successfully installed and loaded renv 0.15.4-11.
Retrieving 'https://cloud.r-project.org/src/contrib/BiocManager_1.30.18.tar.gz' ...
	OK [file is up to date]
Installing BiocManager [1.30.18] ...
	OK [built from source]
Moving BiocManager [1.30.18] into the cache ...
	OK [moved to cache in 1.8 milliseconds]
Warning messages:
1: This project is configured to use R version '4.0.3', but '3.6.0' is currently being used. 
2: This project is configured to use Bioconductor '3.12', which is not compatible with R '3.6.0'.
Use 'renv::init(bioconductor = TRUE)' to re-initialize this project with the latest Bioconductor release. 
Loading required namespace: remotes
Failed with error:  ‘there is no package called ‘remotes’’
Retrieving 'https://cloud.r-project.org/src/contrib/remotes_2.4.2.tar.gz' ...
	OK [file is up to date]
Installing remotes [2.4.2] ...
	OK [built from source]
Moving remotes [2.4.2] into the cache ...
	OK [moved to cache in 1.8 milliseconds]
Downloading GitHub repo rstudio/renv@HEAD
Running `R CMD build`...
* checking for file ‘/tmp/slurm_38898073/RtmpD8yI2m/remotes770759f9070f/rstudio-renv-b4403d9/DESCRIPTION’ ... OK
* preparing ‘renv’:
* checking DESCRIPTION meta-information ... OK
* checking for LF line-endings in source and make files and shell scripts
* checking for empty or unneeded directories
Removed empty directory ‘renv/tests/testthat/resources/dummy.Rmd’
* building ‘renv_0.15.5-17.tar.gz’
Installing package into ‘/camp/stp/babs/working/boeings/Projects/hillc/danielle.park/450A_sc_scPDAC_PRJCA001063full_Fibroblast_subclustering/workdir/scRNAseq/renv/library/R-3.6/x86_64-pc-linux-gnu’
(as ‘lib’ is unspecified)
Warning messages:
1: This project is configured to use R version '4.0.3', but '3.6.0' is currently being used. 
2: This project is configured to use Bioconductor '3.12', which is not compatible with R '3.6.0'.
Use 'renv::init(bioconductor = TRUE)' to re-initialize this project with the latest Bioconductor release. 
* installing *source* package ‘renv’ ...
** using staged installation
** R
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (renv)
The following package(s) will be updated:

# Bioconductor =======================
- AUCell                   [* -> 1.12.0]
- AnnotationDbi            [* -> 1.52.0]
- AnnotationHub            [* -> 2.22.1]
- Biobase                  [* -> 2.50.0]
- BiocFileCache            [* -> 1.14.0]
- BiocGenerics             [* -> 0.36.1]
- BiocParallel             [* -> 1.24.1]
- BiocVersion              [* -> 3.12.0]
- DESeq2                   [* -> 1.30.1]
- DO.db                    [* -> 2.9]
- DOSE                     [* -> 3.16.0]
- DelayedArray             [* -> 0.16.3]
- ExperimentHub            [* -> 1.16.1]
- GO.db                    [* -> 3.12.1]
- GOSemSim                 [* -> 2.16.1]
- GSEABase                 [* -> 1.52.1]
- GenomeInfoDb             [* -> 1.26.7]
- GenomeInfoDbData         [* -> 1.2.4]
- GenomicRanges            [* -> 1.42.0]
- IRanges                  [* -> 2.24.1]
- MatrixGenerics           [* -> 1.2.1]
- S4Vectors                [* -> 0.28.1]
- SingleCellExperiment     [* -> 1.12.0]
- SummarizedExperiment     [* -> 1.20.0]
- XVector                  [* -> 0.30.0]
- annotate                 [* -> 1.68.0]
- clusterProfiler          [* -> 3.18.1]
- clustifyr                [* -> 1.2.0]
- enrichplot               [* -> 1.10.2]
- fgsea                    [* -> 1.16.0]
- genefilter               [* -> 1.72.1]
- geneplotter              [* -> 1.68.0]
- graph                    [* -> 1.68.0]
- interactiveDisplayBase   [* -> 1.28.0]
- qvalue                   [* -> 2.22.0]
- zlibbioc                 [* -> 1.36.0]

# CRAN ===============================
- BiocManager              [1.30.18 -> 1.30.16]
- KernSmooth               [2.23-15 -> 2.23-17]
- MASS                     [7.3-51.4 -> 7.3-54]
- Matrix                   [1.2-17 -> 1.3-4]
- boot                     [1.3-22 -> 1.3-25]
- cluster                  [2.0.8 -> 2.1.0]
- codetools                [0.2-16 -> 0.2-18]
- foreign                  [0.8-71 -> 0.8-79]
- lattice                  [0.20-38 -> 0.20-45]
- mgcv                     [1.8-28 -> 1.8-38]
- nlme                     [3.1-139 -> 3.1-153]
- nnet                     [7.3-12 -> 7.3-14]
- remotes                  [2.4.2 -> 2.4.1]
- survival                 [2.44-1.1 -> 3.1-12]
- BH                       [* -> 1.72.0-3]
- DBI                      [* -> 1.1.1]
- DT                       [* -> 0.13]
- FNN                      [* -> 1.1.3]
- GlobalOptions            [* -> 0.1.1]
- MatrixModels             [* -> 0.4-1]
- R.methodsS3              [* -> 1.8.1]
- R.oo                     [* -> 1.24.0]
- R.utils                  [* -> 2.11.0]
- R6                       [* -> 2.5.1]
- RANN                     [* -> 2.6.1]
- RColorBrewer             [* -> 1.1-2]
- RCurl                    [* -> 1.98-1.5]
- RMySQL                   [* -> 0.10.22]
- ROCR                     [* -> 1.0-7]
- RSQLite                  [* -> 2.2.0]
- RSpectra                 [* -> 0.16-0]
- Rcpp                     [* -> 1.0.8.3]
- RcppAnnoy                [* -> 0.0.19]
- RcppArmadillo            [* -> 0.9.870.2.0]
- RcppEigen                [* -> 0.3.3.7.0]
- RcppProgress             [* -> 0.4.2]
- Rtsne                    [* -> 0.15]
- Seurat                   [* -> 4.0.5]
- SeuratObject             [* -> 4.0.3]
- SparseM                  [* -> 1.78]
- XML                      [* -> 3.99-0.3]
- abind                    [* -> 1.4-5]
- ape                      [* -> 5.5]
- aplot                    [* -> 0.1.1]
- askpass                  [* -> 1.1]
- assertthat               [* -> 0.2.1]
- attempt                  [* -> 0.3.1]
- backports                [* -> 1.1.6]
- base64enc                [* -> 0.1-3]
- bit                      [* -> 1.1-15.2]
- bit64                    [* -> 0.9-7]
- bitops                   [* -> 1.0-6]
- blob                     [* -> 1.2.1]
- brew                     [* -> 1.0-6]
- brio                     [* -> 1.1.2]
- broom                    [* -> 0.5.6]
- bslib                    [* -> 0.3.1]
- caTools                  [* -> 1.18.0]
- cachem                   [* -> 1.0.6]
- callr                    [* -> 3.7.0]
- car                      [* -> 3.0-7]
- carData                  [* -> 3.0-3]
- cellranger               [* -> 1.1.0]
- circlize                 [* -> 0.4.8]
- cli                      [* -> 3.1.0]
- clipr                    [* -> 0.7.1]
- colorspace               [* -> 2.0-2]
- colourpicker             [* -> 1.0]
- commonmark               [* -> 1.7]
- config                   [* -> 0.3.1]
- corrplot                 [* -> 0.84]
- covr                     [* -> 3.5.0]
- cowplot                  [* -> 1.1.1]
- cpp11                    [* -> 0.4.1]
- crayon                   [* -> 1.4.2]
- credentials              [* -> 1.3.1]
- crosstalk                [* -> 1.1.0.1]
- curl                     [* -> 4.3.2]
- data.table               [* -> 1.12.8]
- dbplyr                   [* -> 1.4.3]
- deldir                   [* -> 1.0-6]
- desc                     [* -> 1.4.0]
- devtools                 [* -> 2.3.0]
- diffobj                  [* -> 0.3.5]
- digest                   [* -> 0.6.28]
- dockerfiler              [* -> 0.1.4]
- dotCall64                [* -> 1.0-1]
- downloader               [* -> 0.4]
- dplyr                    [* -> 1.0.7]
- dqrng                    [* -> 0.3.0]
- ellipsis                 [* -> 0.3.2]
- enrichR                  [* -> 3.0]
- entropy                  [* -> 1.2.1]
- evaluate                 [* -> 0.14]
- fansi                    [* -> 0.5.0]
- farver                   [* -> 2.1.0]
- fastmap                  [* -> 1.1.0]
- fastmatch                [* -> 1.1-3]
- fields                   [* -> 13.3]
- fitdistrplus             [* -> 1.0-14]
- fontawesome              [* -> 0.2.2]
- forcats                  [* -> 0.5.0]
- foreach                  [* -> 1.5.0]
- formatR                  [* -> 1.11]
- fs                       [* -> 1.5.0]
- furrr                    [* -> 0.2.3]
- futile.logger            [* -> 1.4.3]
- futile.options           [* -> 1.0.1]
- future                   [* -> 1.19.1]
- future.apply             [* -> 1.4.0]
- gam                      [* -> 1.16.1]
- gdata                    [* -> 2.18.0]
- gdtools                  [* -> 0.2.3]
- generics                 [* -> 0.1.1]
- gert                     [* -> 1.4.3]
- gganimate                [* -> 1.0.7]
- ggdendro                 [* -> 0.1.22]
- ggforce                  [* -> 0.3.3]
- ggfun                    [* -> 0.0.6]
- ggplot2                  [* -> 3.3.5]
- ggplotify                [* -> 0.1.0]
- ggpubr                   [* -> 0.3.0]
- ggraph                   [* -> 2.0.5]
- ggrepel                  [* -> 0.8.2]
- ggridges                 [* -> 0.5.2]
- ggsci                    [* -> 2.9]
- ggsignif                 [* -> 0.6.0]
- gh                       [* -> 1.3.0]
- git2r                    [* -> 0.26.1]
- gitcreds                 [* -> 0.1.1]
- globals                  [* -> 0.14.0]
- glue                     [* -> 1.5.0]
- goftest                  [* -> 1.2-2]
- golem                    [* -> 0.3.1]
- gplots                   [* -> 3.0.3]
- graphlayouts             [* -> 0.7.1]
- gridExtra                [* -> 2.3]
- gridGraphics             [* -> 0.5-1]
- gtable                   [* -> 0.3.0]
- gtools                   [* -> 3.8.2]
- haven                    [* -> 2.2.0]
- hdf5r                    [* -> 1.3.2]
- here                     [* -> 1.0.1]
- hexbin                   [* -> 1.28.1]
- highr                    [* -> 0.9]
- hms                      [* -> 1.1.1]
- htmltools                [* -> 0.5.2]
- htmlwidgets              [* -> 1.5.1]
- httpuv                   [* -> 1.6.3]
- httr                     [* -> 1.4.2]
- ica                      [* -> 1.0-2]
- igraph                   [* -> 1.2.5]
- ini                      [* -> 0.3.1]
- irlba                    [* -> 2.3.3]
- isoband                  [* -> 0.2.5]
- iterators                [* -> 1.0.12]
- jquerylib                [* -> 0.1.4]
- jrc                      [* -> 0.5.0]
- jsonlite                 [* -> 1.7.2]
- kernlab                  [* -> 0.9-29]
- knitr                    [* -> 1.36]
- labeling                 [* -> 0.4.2]
- lambda.r                 [* -> 1.2.4]
- later                    [* -> 1.3.0]
- lazyeval                 [* -> 0.2.2]
- leiden                   [* -> 0.3.3]
- lifecycle                [* -> 1.0.1]
- listenv                  [* -> 0.8.0]
- lme4                     [* -> 1.1-23]
- lmtest                   [* -> 0.9-37]
- locfit                   [* -> 1.5-9.4]
- lsei                     [* -> 1.2-0]
- lubridate                [* -> 1.7.8]
- magrittr                 [* -> 2.0.1]
- maps                     [* -> 3.4.0]
- maptools                 [* -> 0.9-9]
- matrixStats              [* -> 0.61.0]
- memoise                  [* -> 1.1.0]
- mime                     [* -> 0.12]
- miniUI                   [* -> 0.1.1.1]
- minqa                    [* -> 1.2.4]
- mixtools                 [* -> 1.2.0]
- modelr                   [* -> 0.1.6]
- munsell                  [* -> 0.5.0]
- nloptr                   [* -> 1.2.2.1]
- npsurv                   [* -> 0.4-0]
- officer                  [* -> 0.4.1]
- openssl                  [* -> 1.4.5]
- openxlsx                 [* -> 4.1.4]
- parsedate                [* -> 1.2.1]
- patchwork                [* -> 1.0.0]
- pbapply                  [* -> 1.4-2]
- pbkrtest                 [* -> 0.4-8.6]
- pheatmap                 [* -> 1.0.12]
- pillar                   [* -> 1.6.4]
- pkgbuild                 [* -> 1.2.0]
- pkgconfig                [* -> 2.0.3]
- pkgload                  [* -> 1.2.3]
- plogr                    [* -> 0.2.0]
- plotly                   [* -> 4.9.2.1]
- plyr                     [* -> 1.8.6]
- png                      [* -> 0.1-7]
- polyclip                 [* -> 1.10-0]
- polynom                  [* -> 1.4-0]
- praise                   [* -> 1.0.0]
- prettyunits              [* -> 1.1.1]
- processx                 [* -> 3.5.2]
- progress                 [* -> 1.2.2]
- promises                 [* -> 1.2.0.1]
- proxy                    [* -> 0.4-24]
- ps                       [* -> 1.6.0]
- purrr                    [* -> 0.3.4]
- quantreg                 [* -> 5.55]
- rappdirs                 [* -> 0.3.3]
- rcmdcheck                [* -> 1.4.0]
- readr                    [* -> 1.3.1]
- readxl                   [* -> 1.3.1]
- rematch                  [* -> 1.0.1]
- rematch2                 [* -> 2.1.2]
- reprex                   [* -> 0.3.0]
- reshape2                 [* -> 1.4.4]
- reticulate               [* -> 1.15]
- rhub                     [* -> 1.1.1]
- rio                      [* -> 0.5.16]
- rjson                    [* -> 0.2.20]
- rlang                    [* -> 0.4.12]
- rmarkdown                [* -> 2.13]
- roxygen2                 [* -> 7.1.2]
- rprojroot                [* -> 2.0.2]
- rsample                  [* -> 0.1.1]
- rstatix                  [* -> 0.5.0]
- rstudioapi               [* -> 0.13]
- rversions                [* -> 2.0.1]
- rvest                    [* -> 0.3.5]
- rvg                      [* -> 0.2.5]
- sass                     [* -> 0.4.0]
- scales                   [* -> 1.1.1]
- scattermore              [* -> 0.7]
- scatterpie               [* -> 0.1.7]
- sctransform              [* -> 0.3.2]
- segmented                [* -> 1.1-0]
- selectr                  [* -> 0.4-2]
- sessioninfo              [* -> 1.2.1]
- shadowtext               [* -> 0.0.9]
- shape                    [* -> 1.4.4]
- shiny                    [* -> 1.7.1]
- shinyjs                  [* -> 1.1]
- sitmo                    [* -> 2.0.2]
- sleepwalk                [* -> 0.3.2]
- slider                   [* -> 0.2.2]
- snow                     [* -> 0.4-4]
- sourcetools              [* -> 0.1.7]
- sp                       [* -> 1.4-1]
- spam                     [* -> 2.7-0]
- spatstat.core            [* -> 2.3-1]
- spatstat.data            [* -> 2.1-0]
- spatstat.geom            [* -> 2.3-0]
- spatstat.sparse          [* -> 2.0-0]
- spatstat.utils           [* -> 2.2-0]
- statmod                  [* -> 1.4.34]
- stringi                  [* -> 1.7.5]
- stringr                  [* -> 1.4.0]
- sys                      [* -> 3.4]
- systemfonts              [* -> 1.0.3]
- tensor                   [* -> 1.5]
- testthat                 [* -> 3.1.0]
- thinkr                   [* -> 0.15]
- tibble                   [* -> 3.1.6]
- tidygraph                [* -> 1.2.0]
- tidyr                    [* -> 1.1.4]
- tidyselect               [* -> 1.1.1]
- tidytree                 [* -> 0.3.9]
- tidyverse                [* -> 1.3.0]
- tinytex                  [* -> 0.34]
- tweenr                   [* -> 1.0.2]
- usethis                  [* -> 2.1.3]
- utf8                     [* -> 1.2.2]
- uuid                     [* -> 1.0-3]
- uwot                     [* -> 0.1.10]
- vctrs                    [* -> 0.3.8]
- viridis                  [* -> 0.5.1]
- viridisLite              [* -> 0.4.0]
- waldo                    [* -> 0.3.1]
- warp                     [* -> 0.2.0]
- whisker                  [* -> 0.4]
- whoami                   [* -> 1.3.0]
- withr                    [* -> 2.4.2]
- xfun                     [* -> 0.30]
- xml2                     [* -> 1.3.2]
- xopen                    [* -> 1.0.0]
- xtable                   [* -> 1.8-4]
- yaml                     [* -> 2.2.1]
- yulab.utils              [* -> 0.0.4]
- zip                      [* -> 2.2.0]
- zoo                      [* -> 1.8-7]

# GitHub =============================
- renv                     [ver: 0.15.5-17 -> 0.15.4-11; ref: HEAD -> main; sha: b4403d91 -> 203d222e]
- DoubletFinder            [* -> chris-mcginnis-ucsf/DoubletFinder]
- biologicSeqTools         [* -> decusInLabore/biologicSeqTools@main]
- biologicSeqTools2        [* -> decusInLabore/biologicSeqTools2@main]
- biologicToolsSC          [* -> decusInLabore/biologicToolsSC@main]
- biologicViewerSC         [* -> decusinlabore/biologicViewerSC@main]
- ggtree                   [* -> YuLab-SMU/ggtree]
- rvcheck                  [* -> GuangchuangYu/rvcheck]
- treeio                   [* -> GuangchuangYu/treeio]

Retrieving 'https://api.github.com/repos/rstudio/renv/tarball/203d222e8d540f0ab29f251986fd47a3363907d7' ...
	OK [downloaded 465 Kb in 0.5 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/boot/boot_1.3-25.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/cluster/cluster_2.1.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/foreign/foreign_0.8-79.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/KernSmooth/KernSmooth_2.23-17.tar.gz' ...
	OK [downloaded 24.8 Kb in 0.2 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/nnet/nnet_7.3-14.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/survival/survival_3.1-12.tar.gz' ...
	OK [downloaded 7.1 Mb in 1.1 secs]
* Querying repositories for available source packages ... Done!
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/AUCell_1.12.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/data.table/data.table_1.12.8.tar.gz' ...
	OK [downloaded 4.7 Mb in 0.9 secs]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/GSEABase_1.52.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/AnnotationDbi_1.52.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/RSQLite/RSQLite_2.2.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/bit64/bit64_0.9-7.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/bit/bit_1.1-15.2.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/glue/glue_1.5.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/memoise/memoise_1.1.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Rcpp_1.0.8.3.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/BH/BH_1.72.0-3.tar.gz' ...
	OK [downloaded 11.3 Mb in 0.3 secs]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/S4Vectors_0.28.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/BiocGenerics_0.36.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/Biobase_2.50.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/IRanges_2.24.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/annotate_1.68.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/graph_1.68.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/segmented/segmented_1.1-0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/SummarizedExperiment_1.20.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/GenomeInfoDb_1.26.7.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/data/annotation/src/contrib/GenomeInfoDbData_1.2.4.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/DelayedArray_0.16.3.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/MatrixGenerics_1.2.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/GenomicRanges_1.42.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/XVector_0.30.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/zlibbioc_1.36.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/AnnotationHub_2.22.1.tar.gz' ...
	OK [downloaded 1.2 Mb in 0.1 secs]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/BiocVersion_3.12.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/interactiveDisplayBase_1.28.0.tar.gz' ...
	OK [downloaded 11.4 Kb in 0.1 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/DT/DT_0.13.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/htmlwidgets/htmlwidgets_1.5.1.tar.gz' ...
	OK [downloaded 846.9 Kb in 0.2 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/crosstalk/crosstalk_1.1.0.1.tar.gz' ...
	OK [downloaded 657.3 Kb in 0.6 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/generics/generics_0.1.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/BiocFileCache_1.14.0.tar.gz' ...
	OK [downloaded 282.5 Kb in 0.2 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/dbplyr/dbplyr_1.4.3.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/BiocParallel_1.24.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/snow_0.4-4.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/DESeq2_1.30.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/genefilter_1.72.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/geneplotter_1.68.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.9.870.2.0.tar.gz' ...
	OK [downloaded 1.6 Mb in 0.7 secs]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/DOSE_3.16.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/fgsea_1.16.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/GOSemSim_2.16.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/data/annotation/src/contrib/GO.db_3.12.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/qvalue_2.22.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://api.github.com/repos/chris-mcginnis-ucsf/DoubletFinder/tarball/67fb8b5808eb16167ead5f9b439677cc24837554' ...
	OK [downloaded 1.4 Mb in 0.4 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/fields_13.3.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/ROCR/ROCR_1.0-7.tar.gz' ...
	OK [downloaded 105.3 Kb in 0.4 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/gplots/gplots_3.0.3.tar.gz' ...
	OK [downloaded 504.1 Kb in 0.6 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/gdata/gdata_2.18.0.tar.gz' ...
	OK [downloaded 1 Mb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/caTools/caTools_1.18.0.tar.gz' ...
	OK [downloaded 63.1 Kb in 0.3 secs]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/ExperimentHub_1.16.1.tar.gz' ...
	OK [downloaded 467.8 Kb in 0.1 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/GlobalOptions/GlobalOptions_0.1.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/MatrixModels/MatrixModels_0.4-1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/RcppEigen/RcppEigen_0.3.3.7.0.tar.gz' ...
	OK [downloaded 1.6 Mb in 0.7 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/Seurat/Seurat_4.0.5.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/fitdistrplus/fitdistrplus_1.0-14.tar.gz' ...
	OK [downloaded 1.5 Mb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/npsurv/npsurv_0.4-0.tar.gz' ...
	OK [downloaded 28.6 Kb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/lsei/lsei_1.2-0.tar.gz' ...
	OK [downloaded 33.6 Kb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/future.apply/future.apply_1.4.0.tar.gz' ...
	OK [downloaded 56 Kb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/ggridges/ggridges_0.5.2.tar.gz' ...
	OK [downloaded 2.2 Mb in 0.6 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/igraph/igraph_1.2.5.tar.gz' ...
	OK [downloaded 2.6 Mb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/leiden/leiden_0.3.3.tar.gz' ...
	OK [downloaded 497.9 Kb in 0.2 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/reticulate/reticulate_1.15.tar.gz' ...
	OK [downloaded 1.4 Mb in 0.2 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/lmtest/lmtest_0.9-37.tar.gz' ...
	OK [downloaded 180.1 Kb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/zoo/zoo_1.8-7.tar.gz' ...
	OK [downloaded 841.1 Kb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/patchwork/patchwork_1.0.0.tar.gz' ...
	OK [downloaded 3 Mb in 0.8 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/pbapply/pbapply_1.4-2.tar.gz' ...
	OK [downloaded 14.7 Kb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/plotly/plotly_4.9.2.1.tar.gz' ...
	OK [downloaded 3.5 Mb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/cpp11/cpp11_0.4.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/hexbin/hexbin_1.28.1.tar.gz' ...
	OK [downloaded 643.7 Kb in 0.2 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/scattermore/scattermore_0.7.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_4.0.3.tar.gz' ...
	OK [downloaded 147.8 Kb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/spatstat.core/spatstat.core_2.3-1.tar.gz' ...
	OK [downloaded 1.5 Mb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/spatstat.utils/spatstat.utils_2.2-0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/spatstat.sparse/spatstat.sparse_2.0-0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/tensor_1.5.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/goftest/goftest_1.2-2.tar.gz' ...
	OK [downloaded 18.1 Kb in 0.2 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/spatstat.data/spatstat.data_2.1-0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/spatstat.geom/spatstat.geom_2.3-0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/deldir_1.0-6.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/SingleCellExperiment_1.12.0.tar.gz' ...
	OK [downloaded 622.3 Kb in 0.6 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/SparseM/SparseM_1.78.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/ggfun_0.0.6.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/backports/backports_1.1.6.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/openxlsx/openxlsx_4.1.4.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://api.github.com/repos/decusInLabore/biologicSeqTools2/tarball/d7f9a8c8e68f0fd2b97ca0b82cfa63fd01b3a430' ...
	OK [downloaded 107.4 Kb in 0.3 secs]
Retrieving 'https://api.github.com/repos/decusInLabore/biologicToolsSC/tarball/1d65b32df7d535512b36dd7d61d0e00f9c694899' ...
	OK [downloaded 94 Kb in 0.3 secs]
Retrieving 'https://api.github.com/repos/decusinlabore/biologicViewerSC/tarball/d3e057a7098ed1607b7f252b798cc1293b4ebf32' ...
	OK [downloaded 208.3 Kb in 0.3 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/usethis/usethis_2.1.3.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/gert/gert_1.4.3.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/xfun/xfun_0.30.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/devtools/devtools_2.3.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/covr/covr_3.5.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/git2r/git2r_0.26.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/sessioninfo/sessioninfo_1.2.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/rversions/rversions_2.0.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/lubridate/lubridate_1.7.8.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/officer/officer_0.4.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/uuid/uuid_1.0-3.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/readr/readr_1.3.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/colourpicker/colourpicker_1.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/shinyjs/shinyjs_1.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/broom/broom_0.5.6.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/car/car_3.0-7.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.4-8.6.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/lme4/lme4_1.1-23.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.2.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/statmod/statmod_1.4.34.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/quantreg/quantreg_5.55.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/maptools/maptools_0.9-9.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/sp/sp_1.4-1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/rio/rio_0.5.16.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/haven/haven_2.2.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/forcats/forcats_0.5.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/carData/carData_3.0-3.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/circlize/circlize_0.4.8.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/shape/shape_1.4.4.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/clusterProfiler_3.18.1.tar.gz' ...
	OK [downloaded 2.5 Mb in 0.7 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/downloader_0.4.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/enrichplot_1.10.2.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/scatterpie_0.1.7.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/shadowtext/shadowtext_0.0.9.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://api.github.com/repos/GuangchuangYu/rvcheck/tarball/1a504baf6e689ec890a45d48b9180ddb06b4990d' ...
	OK [downloaded 6 Kb in 0.3 secs]
Retrieving 'https://bioconductor.org/packages/3.12/bioc/src/contrib/clustifyr_1.2.0.tar.gz' ...
	OK [downloaded 3 Mb in 0.7 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/entropy/entropy_1.2.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/proxy/proxy_0.4-24.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/foreach/foreach_1.5.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/iterators/iterators_1.0.12.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/gam/gam_1.16.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/ggpubr/ggpubr_0.3.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/ggsignif/ggsignif_0.6.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/rstatix/rstatix_0.5.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://api.github.com/repos/YuLab-SMU/ggtree/tarball/84d26a5835da42a6d21fb1946ecb37fb0388f276' ...
	OK [downloaded 318.1 Kb in 0.4 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/tidytree_0.3.9.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://api.github.com/repos/GuangchuangYu/treeio/tarball/9d4e4bff464684240e74044a1206058692bf825c' ...
	OK [downloaded 627.9 Kb in 0.5 secs]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/hdf5r/hdf5r_1.3.2.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/modelr/modelr_0.1.6.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/reprex/reprex_0.3.0.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/rmarkdown/rmarkdown_2.13.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/rsample_0.1.1.tar.gz' ...
	OK [file is up to date]
Retrieving 'https://cloud.r-project.org/src/contrib/Archive/rvest/rvest_0.3.5.tar.gz' ...
	OK [file is up to date]
Installing BiocManager [1.30.16] ...
	OK [linked cache]
Installing remotes [2.4.1] ...
	OK [linked cache]
Installing renv [0.15.4-11] ...
	OK [built from source]
Moving renv [0.15.4-11] into the cache ...
	OK [moved to cache in 2.3 milliseconds]
Installing boot [1.3-25] ...
	OK [built from source]
Moving boot [1.3-25] into the cache ...
	OK [moved to cache in 1.6 milliseconds]
Installing cluster [2.1.0] ...
	OK [built from source]
Moving cluster [2.1.0] into the cache ...
	OK [moved to cache in 1.7 milliseconds]
Installing codetools [0.2-18] ...
	OK [linked cache]
Installing foreign [0.8-79] ...
	FAILED
Error installing package 'foreign':
===================================

ERROR: this R is version 3.6.0, package 'foreign' requires R >= 4.0.0
Error: install of package 'foreign' failed [error code 1]
Traceback (most recent calls last):
12: renv::restore(prompt = FALSE)
11: renv_restore_run_actions(project, diff, current, lockfile, rebuild)
10: renv_install_impl(records)
 9: renv_install_staged(records)
 8: renv_install_default(records)
 7: handler(package, renv_install_package(record))
 6: renv_install_package(record)
 5: withCallingHandlers(renv_install_package_impl(record), error = function(e) {
        vwritef("\tFAILED")
        writef(e$output)
    })
 4: renv_install_package_impl(record)
 3: r_cmd_install(package, path)
 2: r_exec_error(package, output, "install", status)
 1: stop(error)
Execution halted
