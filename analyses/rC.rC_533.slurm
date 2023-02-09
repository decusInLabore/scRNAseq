
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
Quitting from lines 42-150 (sc_PartC_Analysis.Rmd) 
Error in eval(ei, envir) : 
  Check if a biologic object has been initiated and is stored in [projectDir]/data/biologic_active_object/. 
Calls: <Anonymous> ... eval -> eval -> source -> withVisible -> eval -> eval
In addition: Warning messages:
1: In png(..., res = dpi, units = "in") :
  unable to open connection to X11 display ''
2: replacing previous import 'RMySQL::isIdCurrent' by 'RSQLite::isIdCurrent' when loading 'biologicSeqTools2' 

Execution halted
