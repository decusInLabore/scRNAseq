

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
 $ echo   : symbol F
 $ results: symbol F

── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──
✔ ggplot2 3.3.3     ✔ purrr   0.3.4
✔ tibble  3.1.0     ✔ dplyr   1.0.5
✔ tidyr   1.1.3     ✔ stringr 1.4.0
✔ readr   1.4.0     ✔ forcats 0.5.1
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()

Attaching package: 'scales'

The following object is masked from 'package:purrr':

    discard

The following object is masked from 'package:readr':

    col_factor

Skipping install of 'biologicSeqTools' from a github remote, the SHA1 (3c916392) has not changed since last install.
  Use `force = TRUE` to force installation
Skipping install of 'biologicToolsSC' from a github remote, the SHA1 (c9ec492f) has not changed since last install.
  Use `force = TRUE` to force installation
  |                                                                              |........................                                              |  35%
  ordinary text without R code

  |                                                                              |............................                                          |  40%
label: old_res (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol F
 $ results: symbol F

Loading required package: DBI
  |                                                                              |................................                                      |  45%
  ordinary text without R code

  |                                                                              |...................................                                   |  50%
label: create_app_upload_tables (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol F
 $ results: symbol F

Skipping install of 'biologicSeqTools' from a github remote, the SHA1 (3c916392) has not changed since last install.
  Use `force = TRUE` to force installation
Skipping install of 'biologicViewerSC' from a github remote, the SHA1 (d0e93c7c) has not changed since last install.
  Use `force = TRUE` to force installation

Attaching package: 'biologicViewerSC'

The following object is masked from 'package:biologicToolsSC':

    addDf2seuratMetaData

  |                                                                              |......................................                                |  55%
  ordinary text without R code

  |                                                                              |..........................................                            |  60%
label: populate_expr_database_1 (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol F
 $ results: symbol F

  |                                                                              |..............................................                        |  65%
  ordinary text without R code

  |                                                                              |.................................................                     |  70%
label: populate_expr_database (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol F
 $ results: symbol F


Attaching package: 'gplots'

The following object is masked from 'package:stats':

    lowess

  |                                                                              |....................................................                  |  75%
  ordinary text without R code

  |                                                                              |........................................................              |  80%
label: create_ref (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol F
 $ results: symbol F

  |                                                                              |............................................................          |  85%
  ordinary text without R code

  |                                                                              |...............................................................       |  90%
label: create_microwebsite (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol F
 $ results: symbol F

  |                                                                              |..................................................................    |  95%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: documentation (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol T
 $ results: symbol T


output file: sc_PartC_Analysis.knit.md

/camp/apps/eb/dev/software/Anaconda2/2019.10/bin/pandoc +RTS -K512m -RTS sc_PartC_Analysis.utf8.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/446B_scRNAseq_subClustering_pbl359B_12clusters/workdir/sc_PartC_Analysis.html --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 3 --variable toc_float=1 --variable toc_selectors=h1,h2,h3 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --variable theme=bootstrap --include-in-header /tmp/RtmpXNixFV/rmarkdown-str4d3cd6b3a51.html --mathjax --variable 'mathjax-url:https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/446B_scRNAseq_subClustering_pbl359B_12clusters/workdir/sc_PartC_Analysis.html
Warning messages:
1: replacing previous import 'RMySQL::isIdCurrent' by 'RSQLite::isIdCurrent' when loading 'biologicSeqTools' 
2: replacing previous import 'RMySQL::isIdCurrent' by 'RSQLite::isIdCurrent' when loading 'biologicViewerSC' 
3: In heatmap.2(m.df1, trace = "none", dendrogram = "both", density.info = "none",  :
  Discrepancy: Colv is FALSE, while dendrogram is `both'. Omitting column dendogram.
4: In in_dir(input_dir(), evaluate(code, envir = env, new_device = FALSE,  :
  You changed the working directory to /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/446B_scRNAseq_subClustering_pbl359B_12clusters/workdir (probably via setwd()). It will be restored to /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/446B_scRNAseq_subClustering_pbl359B_12clusters/workdir/scRNAseq. See the Note section in ?knitr::knit
