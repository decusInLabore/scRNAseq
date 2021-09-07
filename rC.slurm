

processing file: sc_PartC_Analysis.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |....                                                                  |   6%
   inline R code fragments

  |                                                                              |.........                                                             |  12%
label: setup (with options) 
List of 1
 $ include: logi FALSE

  |                                                                              |.............                                                         |  19%
  ordinary text without R code

  |                                                                              |..................                                                    |  25%
label: hpc_notes (with options) 
List of 1
 $ include: logi FALSE

  |                                                                              |......................                                                |  31%
  ordinary text without R code

  |                                                                              |..........................                                            |  38%
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

Loading required package: S4Vectors
Loading required package: stats4
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: 'BiocGenerics'

The following objects are masked from 'package:parallel':

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from 'package:dplyr':

    combine, intersect, setdiff, union

The following objects are masked from 'package:stats':

    IQR, mad, sd, var, xtabs

The following objects are masked from 'package:base':

    anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    union, unique, unsplit, which, which.max, which.min


Attaching package: 'S4Vectors'

The following objects are masked from 'package:dplyr':

    first, rename

The following object is masked from 'package:tidyr':

    expand

The following object is masked from 'package:base':

    expand.grid

Loading required package: IRanges

Attaching package: 'IRanges'

The following objects are masked from 'package:dplyr':

    collapse, desc, slice

The following object is masked from 'package:purrr':

    reduce

Loading required package: GenomicRanges
Loading required package: GenomeInfoDb
Loading required package: SummarizedExperiment
Loading required package: Biobase
Welcome to Bioconductor

    Vignettes contain introductory material; view with
    'browseVignettes()'. To cite Bioconductor, see
    'citation("Biobase")', and for packages 'citation("pkgname")'.

Loading required package: DelayedArray
Loading required package: matrixStats

Attaching package: 'matrixStats'

The following objects are masked from 'package:Biobase':

    anyMissing, rowMedians

The following object is masked from 'package:dplyr':

    count

Loading required package: BiocParallel

Attaching package: 'DelayedArray'

The following objects are masked from 'package:matrixStats':

    colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

The following object is masked from 'package:purrr':

    simplify

The following objects are masked from 'package:base':

    aperm, apply, rowsum


Attaching package: 'SummarizedExperiment'

The following object is masked from 'package:Seurat':

    Assays

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

Loading required package: DBI
  |                                                                              |...............................                                       |  44%
  ordinary text without R code

  |                                                                              |...................................                                   |  50%
label: populate_expr_database_1 (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol F
 $ results: symbol F

  |                                                                              |.......................................                               |  56%
  ordinary text without R code

  |                                                                              |............................................                          |  62%
label: populate_expr_database (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol F
 $ results: symbol F


Attaching package: 'gplots'

The following object is masked from 'package:IRanges':

    space

The following object is masked from 'package:S4Vectors':

    space

The following object is masked from 'package:stats':

    lowess

  |                                                                              |................................................                      |  69%
  ordinary text without R code

  |                                                                              |....................................................                  |  75%
label: create_ref (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol F
 $ results: symbol F

  |                                                                              |.........................................................             |  81%
  ordinary text without R code

  |                                                                              |.............................................................         |  88%
label: create_microwebsite (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol F
 $ results: symbol F

  |                                                                              |..................................................................    |  94%
  ordinary text without R code

  |                                                                              |......................................................................| 100%
label: documentation (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : symbol T
 $ results: symbol T


output file: sc_PartC_Analysis.knit.md

/camp/apps/eb/dev/software/Anaconda2/2019.10/bin/pandoc +RTS -K512m -RTS sc_PartC_Analysis.utf8.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449D_subclustering_449C_SC20193_C1_2_3_4_5/workdir/sc_PartC_Analysis.html --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmarkdown/lua/latex-div.lua --self-contained --variable bs3=TRUE --standalone --section-divs --table-of-contents --toc-depth 3 --variable toc_float=1 --variable toc_selectors=h1,h2,h3 --variable toc_collapsed=1 --variable toc_smooth_scroll=1 --variable toc_print=1 --template /camp/apps/misc/stp/babs/manual/software/r/r-site-library/bioc-3.9/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --variable theme=bootstrap --include-in-header /tmp/RtmpMCBl6A/rmarkdown-strc30866b6431f.html --mathjax --variable 'mathjax-url:https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --variable code_folding=hide --variable code_menu=1 

Output created: /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449D_subclustering_449C_SC20193_C1_2_3_4_5/workdir/sc_PartC_Analysis.html
Warning messages:
1: replacing previous import 'RMySQL::isIdCurrent' by 'RSQLite::isIdCurrent' when loading 'biologicSeqTools' 
2: In heatmap.2(m.df1, trace = "none", dendrogram = "both", density.info = "none",  :
  Discrepancy: Colv is FALSE, while dendrogram is `both'. Omitting column dendogram.
3: In in_dir(input_dir(), evaluate(code, envir = env, new_device = FALSE,  :
  You changed the working directory to /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449D_subclustering_449C_SC20193_C1_2_3_4_5/workdir (probably via setwd()). It will be restored to /camp/stp/babs/working/boeings/Projects/bonfantip/roberta.ragazzini/449D_subclustering_449C_SC20193_C1_2_3_4_5/workdir/scRNAseq. See the Note section in ?knitr::knit
