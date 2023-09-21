
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



processing file: DGE_Marker_Analysis.Rmd
  |                                                                              |                                                                      |   0%  |                                                                              |.                                                                     |   2%
  ordinary text without R code

  |                                                                              |...                                                                   |   4%
label: setup (with options) 
List of 1
 $ include: logi FALSE

  |                                                                              |....                                                                  |   6%
  ordinary text without R code

  |                                                                              |......                                                                |   8%
label: unnamed-chunk-1 (with options) 
List of 1
 $ include: logi FALSE

  |                                                                              |.......                                                               |  10%
  ordinary text without R code

  |                                                                              |.........                                                             |  12%
label: unnamed-chunk-2 (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi FALSE
 $ results: logi FALSE

Loading required namespace: remotes
  |                                                                              |..........                                                            |  14%
  ordinary text without R code

  |                                                                              |...........                                                           |  16%
label: unnamed-chunk-3 (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi TRUE
 $ results: logi FALSE

Collecting package metadata (current_repodata.json): ...working... done
Solving environment: ...working... done

## Package Plan ##

  environment location: /nemo/stp/babs/working/boeings/Projects/boeings/stefan.boeing/dev_scRNAseq_PBMCstim/scripts/scRNAseq/analyses/DGE/renv/python/r-reticulate

  added / updated specs:
    - numpy
    - python=3.6


The following packages will be downloaded:

    package                    |            build
    ---------------------------|-----------------
    libsqlite-3.43.0           |       h2797004_0         821 KB  conda-forge
    sqlite-3.43.0              |       h2c6b66d_0         812 KB  conda-forge
    ------------------------------------------------------------
                                           Total:         1.6 MB

The following NEW packages will be INSTALLED:

  _libgcc_mutex      conda-forge/linux-64::_libgcc_mutex-0.1-conda_forge
  _openmp_mutex      conda-forge/linux-64::_openmp_mutex-4.5-2_gnu
  ca-certificates    conda-forge/linux-64::ca-certificates-2023.7.22-hbcca054_0
  certifi            conda-forge/linux-64::certifi-2016.9.26-py36_0
  ld_impl_linux-64   conda-forge/linux-64::ld_impl_linux-64-2.40-h41732ed_0
  libblas            conda-forge/linux-64::libblas-3.9.0-17_linux64_openblas
  libcblas           conda-forge/linux-64::libcblas-3.9.0-17_linux64_openblas
  libffi             conda-forge/linux-64::libffi-3.4.2-h7f98852_5
  libgcc-ng          conda-forge/linux-64::libgcc-ng-13.1.0-he5830b7_0
  libgfortran-ng     conda-forge/linux-64::libgfortran-ng-13.1.0-h69a702a_0
  libgfortran5       conda-forge/linux-64::libgfortran5-13.1.0-h15d22d2_0
  libgomp            conda-forge/linux-64::libgomp-13.1.0-he5830b7_0
  liblapack          conda-forge/linux-64::liblapack-3.9.0-17_linux64_openblas
  libnsl             conda-forge/linux-64::libnsl-2.0.0-h7f98852_0
  libopenblas        conda-forge/linux-64::libopenblas-0.3.23-pthreads_h80387f5_0
  libsqlite          conda-forge/linux-64::libsqlite-3.43.0-h2797004_0
  libstdcxx-ng       conda-forge/linux-64::libstdcxx-ng-13.1.0-hfd8a6a1_0
  libzlib            conda-forge/linux-64::libzlib-1.2.13-hd590300_5
  ncurses            conda-forge/linux-64::ncurses-6.4-hcb278e6_0
  numpy              conda-forge/linux-64::numpy-1.19.5-py36hfc0c790_2
  openssl            conda-forge/linux-64::openssl-1.1.1v-hd590300_0
  pip                conda-forge/linux-64::pip-20.0.2-py36_1
  python             conda-forge/linux-64::python-3.6.15-hb7a2778_0_cpython
  python_abi         conda-forge/linux-64::python_abi-3.6-2_cp36m
  readline           conda-forge/linux-64::readline-8.2-h8228510_1
  setuptools         conda-forge/linux-64::setuptools-49.6.0-py36h5fab9bb_3
  sqlite             conda-forge/linux-64::sqlite-3.43.0-h2c6b66d_0
  tk                 conda-forge/linux-64::tk-8.6.12-h27826a3_0
  wheel              conda-forge/noarch::wheel-0.36.2-pyhd3deb0d_0
  xz                 conda-forge/linux-64::xz-5.2.6-h166bdaf_0



Downloading and Extracting Packages
sqlite-3.43.0        | 812 KB    |            |   0% sqlite-3.43.0        | 812 KB    | 1          |   2% sqlite-3.43.0        | 812 KB    | ########## | 100% 
libsqlite-3.43.0     | 821 KB    |            |   0% libsqlite-3.43.0     | 821 KB    | ########## | 100% 
Preparing transaction: ...working... done
Verifying transaction: ...working... done
Executing transaction: ...working... done


==> WARNING: A newer version of conda exists. <==
  current version: 4.8.3
  latest version: 23.7.3

Please update conda by running

    $ conda update -n base conda


#
# To activate this environment, use
#
#     $ conda activate /nemo/stp/babs/working/boeings/Projects/boeings/stefan.boeing/dev_scRNAseq_PBMCstim/scripts/scRNAseq/analyses/DGE/renv/python/r-reticulate
#
# To deactivate an active environment, use
#
#     $ conda deactivate

Attaching SeuratObject

Attaching package: 'dplyr'

The following objects are masked from 'package:stats':

    filter, lag

The following objects are masked from 'package:base':

    intersect, setdiff, setequal, union

  |                                                                              |.............                                                         |  18%
  ordinary text without R code

  |                                                                              |..............                                                        |  20%
label: unnamed-chunk-4 (with options) 
List of 4
 $ echo   : logi FALSE
 $ eval   : logi TRUE
 $ warning: logi FALSE
 $ results: chr "asis"

  |                                                                              |................                                                      |  22%
  ordinary text without R code

  |                                                                              |.................                                                     |  24%
label: unnamed-chunk-5 (with options) 
List of 3
 $ eval   : logi TRUE
 $ echo   : logi TRUE
 $ results: logi FALSE

The following grouping variables have 1 value and will be ignored: ident
All grouping variables have 1 value only. Computing across all cells.
Warning: The following arguments are not used: row.names
Centering and scaling data matrix
  |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  69%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  81%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%
Quitting from lines 127-736 (DGE_Marker_Analysis.Rmd) 
Error in loadNamespace(name) : there is no package called 'ggrastr'
Calls: <Anonymous> ... loadNamespace -> withRestarts -> withOneRestart -> doWithOneRestart
In addition: Warning messages:
1: In png(..., res = dpi, units = "in") :
  unable to open connection to X11 display ''
2: replacing previous import 'RMySQL::isIdCurrent' by 'RSQLite::isIdCurrent' when loading 'biologicSeqTools2' 
3: In eng_r(options) :
  Failed to tidy R code in chunk 'unnamed-chunk-5'. Reason:
Error in base::parse(text = code, keep.source = FALSE) : 
  <text>:186:1: unexpected symbol
185: 
186: sort_by
     ^


Execution halted
