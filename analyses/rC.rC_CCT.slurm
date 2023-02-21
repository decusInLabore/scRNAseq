
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

Collecting package metadata (current_repodata.json): ...working... done
Solving environment: ...working... done

## Package Plan ##

  environment location: /camp/stp/babs/working/boeings/Projects/babs/deborah/cell_cycle_regression_test/scripts/scRNAseq/analyses/Main_Analysis/renv/python/r-reticulate

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
  openssl            conda-forge/linux-64::openssl-1.1.1t-h0b41bf4_0
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
#     $ conda activate /camp/stp/babs/working/boeings/Projects/babs/deborah/cell_cycle_regression_test/scripts/scRNAseq/analyses/Main_Analysis/renv/python/r-reticulate
#
# To deactivate an active environment, use
#
#     $ conda deactivate

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
Error in readChar(con, 5L, useBytes = TRUE) : cannot open the connection
Calls: <Anonymous> ... eval_with_user_handlers -> eval -> eval -> load -> readChar
In addition: Warning messages:
1: In png(..., res = dpi, units = "in") :
  unable to open connection to X11 display ''
2: replacing previous import 'RMySQL::isIdCurrent' by 'RSQLite::isIdCurrent' when loading 'biologicSeqTools2' 
3: In readChar(con, 5L, useBytes = TRUE) :
  cannot open file '/camp/stp/babs/working/boeings/Projects/babs/deborah/cell_cycle_regression_test/workdir/': it is a directory

Execution halted
