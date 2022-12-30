
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



processing file: sc_PartC_Analysis.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |....                                                                  |   5%
   inline R code fragments

  |                                                                              |.......                                                               |  10%
label: setup (with options) 
List of 1
 $ include: logi FALSE

  |                                                                              |..........                                                            |  15%
  ordinary text without R code

  |                                                                              |..............                                                        |  20%
label: hpc_notes (with options) 
List of 1
 $ include: logi FALSE

  |                                                                              |..................                                                    |  25%
  ordinary text without R code

  |                                                                              |.....................                                                 |  30%
label: populate_meta_data_database (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi FALSE
 $ results: logi FALSE

Attaching SeuratObject
Registered S3 method overwritten by 'cli':
  method     from         
  print.boxx spatstat.geom
── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──
✔ ggplot2 3.3.5     ✔ purrr   0.3.4
✔ tibble  3.1.6     ✔ dplyr   1.0.7
✔ tidyr   1.1.4     ✔ stringr 1.4.0
✔ readr   1.3.1     ✔ forcats 0.5.0
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
Loading required namespace: remotes
  |                                                                              |........................                                              |  35%
  ordinary text without R code

  |                                                                              |............................                                          |  40%
label: old_res (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi FALSE
 $ results: logi FALSE

Loading required package: DBI
  |                                                                              |................................                                      |  45%
  ordinary text without R code

  |                                                                              |...................................                                   |  50%
label: create_app_upload_tables (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi FALSE
 $ results: logi FALSE


Attaching package: 'biologicViewerSC'

The following object is masked from 'package:biologicToolsSC':

    addDf2seuratMetaData

The following objects are masked from 'package:biologicSeqTools2':

    assignDbUsersAndPrivileges, inferDBcategories, killDbConnections,
    upload.datatable.to.database, uploadDbTableInfile


Attaching package: 'scales'

The following object is masked from 'package:purrr':

    discard

The following object is masked from 'package:readr':

    col_factor

  |                                                                              |......................................                                |  55%
  ordinary text without R code

  |                                                                              |..........................................                            |  60%
label: unnamed-chunk-1 (with options) 
List of 3
 $ eval   : logi FALSE
 $ echo   : logi FALSE
 $ results: logi FALSE

  |                                                                              |..............................................                        |  65%
  ordinary text without R code

  |                                                                              |.................................................                     |  70%
label: populate_expr_database (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi FALSE
 $ results: logi FALSE


Attaching package: 'gplots'

The following object is masked from 'package:stats':

    lowess

  |                                                                              |....................................................                  |  75%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: create_ref (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi FALSE
 $ results: logi FALSE

  |                                                                              |............................................................          |  85%
  ordinary text without R code

  |                                                                              |...............................................................       |  90%
label: create_microwebsite (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi FALSE
 $ results: logi FALSE

  |                                                                              |..................................................................    |  95%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: documentation (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi TRUE
 $ results: logi TRUE


output file: sc_PartC_Analysis.knit.md

/camp/apps/eb/software/pandoc/2.2.3.2-foss-2016b/bin/pandoc +RTS -K512m -RTS sc_PartC_Analysis.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/528_scRNAseq_subclustering_of_496_C7_C12_subset1/html_local/sc_PartC_Analysis.html --lua-filter /camp/stp/babs/working/boeings/Projects/home_extra/R/x86_64-pc-linux-gnu-library/4.0/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/stp/babs/working/boeings/Projects/home_extra/R/x86_64-pc-linux-gnu-library/4.0/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 3 --variable toc_float=1 --variable toc_selectors=h1,h2,h3 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/stp/babs/working/boeings/Projects/home_extra/R/x86_64-pc-linux-gnu-library/4.0/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --variable theme=bootstrap --mathjax --variable 'mathjax-url=https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --include-in-header /tmp/slurm_47119208/RtmpL1QY3B/rmarkdown-strdbfb2c719d6a.html --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/528_scRNAseq_subclustering_of_496_C7_C12_subset1/html_local/sc_PartC_Analysis.html
Warning messages:
1: In png(..., res = dpi, units = "in") :
  unable to open connection to X11 display ''
2: replacing previous import 'RMySQL::isIdCurrent' by 'RSQLite::isIdCurrent' when loading 'biologicSeqTools2' 
3: replacing previous import 'RMySQL::isIdCurrent' by 'RSQLite::isIdCurrent' when loading 'biologicViewerSC' 
4: In heatmap.2(m.df1, trace = "none", dendrogram = "both", density.info = "none",  :
  Discrepancy: Colv is FALSE, while dendrogram is `both'. Omitting column dendogram.
5: In .local(conn, statement, ...) :
  Unsigned INTEGER in col 0 imported as numeric
6: In .local(conn, statement, ...) :
  unrecognized MySQL field type 7 in column 2 imported as character
7: In in_dir(input_dir(), expr) :
  You changed the working directory to /camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/528_scRNAseq_subclustering_of_496_C7_C12_subset1/workdir (probably via setwd()). It will be restored to /camp/stp/babs/working/boeings/Projects/evang/roderik.kortlever/528_scRNAseq_subclustering_of_496_C7_C12_subset1/scripts/scRNAseq/analyses/Main_Analysis. See the Note section in ?knitr::knit
