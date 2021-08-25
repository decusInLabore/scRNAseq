ml purge
ml Anaconda3/2019.07

module load Singularity/3.6.4

## Reset chachedir 
## https://github.com/COMBINE-lab/alevin-fry#a-quick-start-run-through-on-sample-data

export SINGULARITY_CACHEDIR=/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/FASTQ_files/alevin/cache
singularity pull docker://combinelab/usefulaf:latest


conda install -c bioconda alevin-fry 




###########################################
##

## Organise files 
cd $AF_SAMPLE_DIR
$ mkdir -p human_CR_3.0/fasta
$ mkdir -p human_CR_3.0/genes
$ wget -v -O human_CR_3.0/fasta/genome.fa -L https://umd.box.com/shared/static/3kuh1lc03bxg1d3hi1jfloez7zoutfjc
$ wget -v -O human_CR_3.0/genes/genes.gtf -L https://umd.box.com/shared/static/tvyg43710ufuuvp8mnuoanowm6xmkbjk
$ mkdir -p data/pbmc_1k_v3_fastqs
$ wget -v -O data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz -L https://umd.box.com/shared/static/bmhtt9db8ojhmbkb6d98mt7fnsdhsymm
$ wget -v -O data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz -L https://umd.box.com/shared/static/h8ymvs2njqiygfsu50jla2uce6p6etke
$ wget -v -O data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz -L https://umd.box.com/shared/static/hi8mkx1yltmhnl9kn22n96xtic2wqm5i
$ wget -v -O data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz -L https://umd.box.com/shared/static/4sn4pio63kk7pku52eo3xg9ztf5tq1ul


## build splicing index 
$ singularity exec --cleanenv \
--bind $AF_SAMPLE_DIR:/workdir \
--pwd /usefulaf/bash usefulaf.sif \
./simpleaf index \
-f /workdir/human_CR_3.0/fasta/genome.fa \
-g /workdir/human_CR_3.0/genes/genes.gtf \
-l 91 -t 16 -o /workdir/human_CR_3.0_splici

## Quantify sample ##
$ singularity exec --cleanenv \
--bind $AF_SAMPLE_DIR:/workdir \
--pwd /usefulaf/bash usefulaf.sif \
./simpleaf quant \
-1 /workdir/data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,/workdir/data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz \
-2 /workdir/data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,/workdir/data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz \
-i /workdir/human_CR_3.0_splici/index \
-o /workdir/quants/pbmc1k_v3 \
-f u -c v3 -r cr-like \
-m /workdir/human_CR_3.0_splici/ref/transcriptome_splici_fl86_t2g_3col.tsv \
-t 16

##
###########################################

###########################################
## Try on A4 sample 

## Pull docker container 
module load Singularity/3.6.4



## Reset chachedir 
## https://github.com/COMBINE-lab/alevin-fry#a-quick-start-run-through-on-sample-data

export SINGULARITY_CACHEDIR=/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/FASTQ_files/alevin/cache
export AF_SAMPLE_DIR=$PWD/alevin

singularity pull docker://combinelab/usefulaf:latest



## path to files

## GTF file
"/camp/svc/reference/Genomics/10x/10x_transcriptomes/GRCh38-3.0.0_GFP-tdT/genes/genes.gtf"

## FA file 
"/camp/svc/reference/Genomics/10x/10x_transcriptomes/GRCh38-3.0.0_GFP-tdT/fasta/genome.fa"

## A4 FASTQ-dir
/workdir/FASTQ_files/CAM441A4/CAM441A4*

  
## R1s
/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L001_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L002_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L003_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L004_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L005_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L006_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L007_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L008_R1_001.fastq.gz

## R2s
/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L001_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L002_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L003_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L004_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L005_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L006_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L007_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L008_R2_001.fastq.gz

## build splicing index 

## Build splice index in 
singularity exec --cleanenv \
--bind $AF_SAMPLE_DIR:/workdir \
--pwd /usefulaf/bash usefulaf.sif \
./simpleaf index \
-f /workdir/human_10X_tdt/fasta/genome.fa \
-g /workdir/human_10X_tdt/genes/genes.gtf \
-l 91 -t 16 -o /workdir/human10X_tdt_splici


## Quantify sample ##
## FASTQ files need to be copied into a directory downstream of /workdir

singularity exec --cleanenv \
--bind $AF_SAMPLE_DIR:/workdir \
--pwd /usefulaf/bash usefulaf.sif \
./simpleaf quant \
-1 /workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L001_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L002_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L003_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L004_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L005_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L006_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L007_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L008_R1_001.fastq.gz \
-2 /workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L001_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L002_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L003_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L004_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L005_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L006_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L007_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A4/CAM441A4_S4_L008_R2_001.fastq.gz \
-i /workdir/human10X_tdt_splici/index \
-o /workdir/quants/CAM441A4 \
-f u -c v3 -r cr-like \
-m /workdir/human10X_tdt_splici/ref/transcriptome_splici_fl86_t2g_3col.tsv \
-t 16

## A4 sample complete ##

#################################################
## Calculate A5 sample 
mkdir FASTQ_files/CAM441A5

cp /camp/stp/babs/inputs/sequencing/data/bonfantip/sara.campinoti/SC19153/primary_data/210315_K00102_0567_BHKTFKBBXY/fastq/CAM441A5* .

R1s
/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L001_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L005_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L002_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L006_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L003_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L007_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L004_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L008_R1_001.fastq.gz

R2s
/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L001_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L005_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L002_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L006_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L003_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L007_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L004_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L008_R2_001.fastq.gz

## Copy fastq files into local FASTQ-folder

singularity exec --cleanenv \
--bind $AF_SAMPLE_DIR:/workdir \
--pwd /usefulaf/bash usefulaf.sif \
./simpleaf quant \
-1 /workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L001_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L005_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L002_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L006_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L003_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L007_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L004_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L008_R1_001.fastq.gz \
-2 /workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L001_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L005_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L002_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L006_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L003_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L007_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L004_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A5/CAM441A5_S5_L008_R2_001.fastq.gz \
-i /workdir/human10X_tdt_splici/index \
-o /workdir/quants/CAM441A5 \
-f u -c v3 -r cr-like \
-m /workdir/human10X_tdt_splici/ref/transcriptome_splici_fl86_t2g_3col.tsv \
-t 16

##
#################################################

#################################################
## Calculate A6 sample 
mkdir FASTQ_files/CAM441A6

cp /camp/stp/babs/inputs/sequencing/data/bonfantip/sara.campinoti/SC19153/primary_data/210315_K00102_0567_BHKTFKBBXY/fastq/CAM441A6* .
#cp /camp/stp/babs/inputs/sequencing/data/bonfantip/sara.campinoti/SC19153/primary_data/190820_K00102_0382_BH7MGCBBXY/fastq/SC19153/CAM441A6* .
#cp /camp/stp/babs/inputs/sequencing/data/bonfantip/sara.campinoti/SC19153/primary_data/191025_K00102_0414_BHFFWMBBXY/fastq/CAM441A6* .
# cp /camp/stp/babs/inputs/sequencing/data/bonfantip/sara.campinoti/SC19153/primary_data/190927_K00102_0401_BHF357BBXY/fastq/CAM441A6* .


R1s
/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L001_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L004_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L007_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L002_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L005_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L008_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L003_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L006_R1_001.fastq.gz

R2s
/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L001_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L004_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L007_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L002_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L005_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L008_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L003_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L006_R2_001.fastq.gz

## Copy fastq files into local FASTQ-folder

singularity exec --cleanenv \
--bind $AF_SAMPLE_DIR:/workdir \
--pwd /usefulaf/bash usefulaf.sif \
./simpleaf quant \
-1 /workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L001_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L004_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L007_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L002_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L005_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L008_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L003_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L006_R1_001.fastq.gz \
-2 /workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L001_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L004_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L007_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L002_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L005_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L008_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L003_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A6/CAM441A6_S6_L006_R2_001.fastq.gz \
-i /workdir/human10X_tdt_splici/index \
-o /workdir/quants/CAM441A6 \
-f u -c v3 -r cr-like \
-m /workdir/human10X_tdt_splici/ref/transcriptome_splici_fl86_t2g_3col.tsv \
-t 16

##
#################################################

#################################################
## Calculate A7 sample 
mkdir FASTQ_files/CAM441A7

cp /camp/stp/babs/inputs/sequencing/data/bonfantip/sara.campinoti/SC19153/primary_data/210315_K00102_0567_BHKTFKBBXY/fastq/CAM441A7* .

R1s
/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L001_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L004_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L007_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L002_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L005_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L008_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L003_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L006_R1_001.fastq.gz

R2s
/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L001_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L004_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L007_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L002_R2_001.fastq.gz  /workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L005_R2_001.fastq.gz  /workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L008_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L003_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L006_R2_001.fastq.gz
## Copy fastq files into local FASTQ-folder

singularity exec --cleanenv \
--bind $AF_SAMPLE_DIR:/workdir \
--pwd /usefulaf/bash usefulaf.sif \
./simpleaf quant \
-1 /workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L001_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L004_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L007_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L002_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L005_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L008_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L003_R1_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L006_R1_001.fastq.gz \
-2 /workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L001_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L004_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L007_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L002_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L005_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L008_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L003_R2_001.fastq.gz,/workdir/FASTQ_files/CAM441A7/CAM441A7_S7_L006_R2_001.fastq.gz \
-i /workdir/human10X_tdt_splici/index \
-o /workdir/quants/CAM441A7 \
-f u -c v3 -r cr-like \
-m /workdir/human10X_tdt_splici/ref/transcriptome_splici_fl86_t2g_3col.tsv \
-t 16

##
#################################################


###############################################################################
## Create combined spliced and unspliced matrix in sce objec                 ##

install.packages("devtools")
devtools::install_github("decusInLabore/biologicSeqTools")
library(biologicSeqTools)


## Load reference Seurat object ##
FNs <- "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/workdir/pbl359B.Seurat.Robj"
load(FNs)

dfMeta <- OsC@meta.data
dfMeta <- dfMeta[,c("cellID", "sampleID", "sampleName")]
dfMeta[["extension"]] <- sapply(dfMeta$cellID, function(x) unlist(strsplit(x, "_"))[2])

FNo <- "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/workdir/pbl359B.bioLOGIC.Robj"
load(FNo)

dfAnno <- Obio@dfGeneAnnotation
dfAnno <- unique(dfAnno[,c("ENSG", "hgnc_symbol")])

# Epi6 - CAM441A4
# Epi7 - CAM441A5
# Epi10 - CAM441A6
# Epi14 - CAM441A7

fileList <- list(
  "Epip6" = "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/af_test_workdir/quants/CAM441A4/quant",
  "Epip7" = "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/af_test_workdir/quants/CAM441A5/quant",
  "Epip10" = "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/af_test_workdir/quants/CAM441A6/quant",
  "Epip14" = "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/af_test_workdir/quants/CAM441A7/quant"
)

splicedList <- list()
unsplicedList <- list()

for (i in 1:length(fileList)){
  
    frydir <- fileList[[i]]
  
    suppressPackageStartupMessages({
      library(rjson)
      library(Matrix)
      library(SingleCellExperiment)
    })
  
    qfile <- file.path(frydir, "quant.json")
    if (!file.exists(qfile)) {
      qfile <- file.path(frydir, "meta_info.json")
    }
    
    meta_info <- fromJSON(file = qfile)
    ng <- meta_info$num_genes
    usa_mode <- meta_info$usa_mode
    
    # read in count matrix
    af_raw <- readMM(file = file.path(frydir, "alevin", "quants_mat.mtx"))
    # if usa mode, each gene gets 3 rows, so the actual number of genes is ng/3
    if (usa_mode) {
      if (ng %% 3 != 0) {
        stop("The number of quantified targets is not a multiple of 3")
      }
      ng <- as.integer(ng/3)
    }
    
    # read in gene name file and cell barcode file
    afg <- read.csv(file.path(frydir, "alevin", "quants_mat_cols.txt"), 
                    strip.white = TRUE, header = FALSE, nrows = ng #, 
                    #col.names = c("gene_ids"), 
                    #row.names = 1
    )[1:ng,1]
    afc <- read.csv(file.path(frydir, "alevin", "quants_mat_rows.txt"), 
                    strip.white = TRUE, header = FALSE #,
                    #col.names = c("barcodes"), row.names = 1
    )
    
    ## Add id extension based on integrated experiment ##
    extension <- as.vector(dfMeta[dfMeta$sampleName == names(fileList)[i], "extension"])
    
    afc <- as.vector(paste0(afc[,1], "_", extension))
    
    
    
    rd <- list("S" = seq(1, ng), "U" =  seq(ng + 1, 2 * ng),
               "A" =  seq(2 * ng + 1, 3 * ng))
    
    ## spliced counts matrix ##
    which_counts <- c("S","A")
    splicedMat <- af_raw[, rd[[which_counts[1]]], drop = FALSE]
    for (wc in which_counts[-1]) {
      splicedMat <- splicedMat + af_raw[, rd[[wc]], drop = FALSE]
    }
    
    ############################
    ## Rename afg
    
    h <- dfAnno$ENSG[dfAnno$ENSG %in% afg]
    dfTempAnno <- dfAnno[dfAnno$ENSG %in% h, ]
    
    dfConversion <- data.frame(ENSG = afg, hgnc_symbol = afg)
    row.names(dfConversion) <- dfConversion$ENSG 
    
    
    dfConversion$hgnc_symbol <- as.character(dfConversion$hgnc_symbol)
    dfConversion[dfTempAnno$ENSG, "hgnc_symbol"] <- dfTempAnno$hgnc_symbol
    dfConversion <- dfConversion[afg, ]
    
    afgGN <- make.unique(dfConversion$hgnc_symbol)
    # sum(duplicated(afgGN))
    
    
    
    ##
    ############################
    
    ## unspliced matrix ##
    which_counts <- c("U")
    unsplicedMat <- af_raw[, rd[[which_counts[1]]], drop = FALSE]
    for (wc in which_counts[-1]) {
      unsplicedMat <- unsplicedMat + af_raw[, rd[[wc]], drop = FALSE]
    }
    
    
    # create SingleCellExperiment object
    spliced <- t(splicedMat)
    row.names(spliced) <- afgGN
    colnames(spliced) <- as.vector(afc)
    unspliced = t(unsplicedMat)
    row.names(unspliced) <- afgGN
    colnames(unspliced) <- as.vector(afc)
    
    
    
    #############################################################
    ## Translate gene IDs into gene names and make unique
    splicedList[[names(fileList)[i]]] <- spliced
    
    unsplicedList[[names(fileList)[i]]] <- unspliced
    ##
    #############################################################
}


###############################################################################
## Create one spliced and one unspliced matrix                               ##

geneOrder <- row.names(splicedList[[1]])

unityMatSpliced <- list()
unityMatUnSpliced <- list()
  
for (i in 1:length(splicedList)){
    splicedList[[i]] <- splicedList[[i]][order(geneOrder), ]
    unsplicedList[[i]] <- unsplicedList[[i]][order(geneOrder), ]
    
    if (i ==1){
      unityMatSpliced <- splicedList[[i]]
    } else {
      unityMatSpliced <- cbind(unityMatSpliced, splicedList[[i]])
    }
    
    if (i ==1){
      unityMatUnSpliced <- unsplicedList[[i]]
    } else {
      unityMatUnSpliced <- cbind(unityMatUnSpliced, unsplicedList[[i]])
    }
}




## Done                                                                      ##
###############################################################################


###############################################################################
## Cells in experiment                                                       ##

cellsInExperiment <- OsC@meta.data$cellID

unityMatSpliced <- unityMatSpliced[,cellsInExperiment]
#unityMatSpliced <- unityMatSpliced[rowSums(unityMatSpliced) > 0,]

unityMatUnSpliced <- unityMatUnSpliced[,cellsInExperiment]
#unityMatUnSpliced <- unityMatUnSpliced[row.names(unityMatSpliced),]

## Create single-cell object ##
# sce <- SingleCellExperiment(list(spliced = unityMatSpliced, unspliced = unityMatUnSpliced),
#                             colData =DataFrame(cells=colnames(unityMatSpliced))   ,
#                             rowData = DataFrame(gene=row.names(unityMatSpliced))
# )
# 
# 
# sce <- sce[sce$sizefactor > 0] 

###############################################################################
## python https://combine-lab.github.io/alevin-fry-tutorials/2021/alevin-fry-velocity/

conda deactivate
conda activate scvelo-0.2.2

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scvelo as scv
import matplotlib
import scipy
import json
import os


fileList <- list(
  "Epip6" = "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/af_test_workdir/quants/CAM441A4/quant",
  "Epip7" = "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/af_test_workdir/quants/CAM441A5/quant",
  "Epip10" = "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/af_test_workdir/quants/CAM441A6/quant",
  "Epip14" = "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/af_test_workdir/quants/CAM441A7/quant"
)

frydir = "/camp/stp/babs/working/boeings/Projects/bonfantip/sara.campinoti/359B_PBL_SC_subset_P6_P7_P10_P14_scRNAseq_lung_thymus_SC19153/af_test_workdir/quants/CAM441A4/quant"
e2n_path = "data/geneid_to_name.txt"
fpath = os.path.sep.join([frydir, "quant.json"])
if !os.path.exists(fpath):
  fpath = os.path.sep.join([frydir, "meta_info.json"])

meta_info = json.load(open(fpath))
ng = meta_info['num_genes']
usa_mode = meta_info['usa_mode']

if usa_mode:
  print("processing input in USA mode, will return A+S as the spliced count, and U as the unspliced count")
else:
  print("please follow previous steps to generate the ount matrix in the USA mode")
assert(False)

af_raw = sc.read_mtx(os.path.sep.join([frydir, "alevin", "quants_mat.mtx"]))
ng = int(ng/3)
e2n = dict([ l.rstrip().split() for l in open(e2n_path).readlines()])
var_names = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"])).readlines()][:ng]
var_names = [e2n[e] for e in var_names]

obs_names = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_rows.txt"])).readlines() ]

x = af_raw.X
spliced = x[:,range(0,ng)] + x[:,range(2*ng,3*ng)]
unspliced = x[:,range(ng, 2*ng)]





###############################################################################
## Save as matrix to read into python                                        ##

scvelo_combined.obs <- as.data.frame(cbind(OsC@meta.data$sampleName, OsC@meta.data$seurat_clusters)) #where combined is the seurat object created from combined.emat
scvelo_combined.obsm <- as.data.frame(Seurat::Embeddings(OsC, reduction = "umap"))

writeMM(t(unityMatSpliced), file = "combined.spliced.mtx")
writeMM(t(unityMatUnSpliced), file = "combined.unspliced.mtx")
write.csv(scvelo_combined.obs, file = "combined.obs.csv", row.names = FALSE)
write.csv(scvelo_combined.obsm, file = "combined.obsm.csv", row.names = FALSE)

## Done                                                                      ##
###############################################################################


combined.emat <- unityMatSpliced #where combined.cm is the .rds output of dropEst
combined.nmat <- unityMatUnSpliced
#combined.smat <- combined.cm$spanning

scvelo_combined.obs <- as.data.frame(cbind(OsC@meta.data$sampleName, OsC@meta.data$seurat_clusters)) #where combined is the seurat object created from combined.emat
colnames(scvelo_combined.obs) <- c("orig.ident", "clusters")
scvelo_combined.obsm <- as.data.frame(Embeddings(combined, reduction = "umap"))

writeMM(t(scvelo_combined.emat), file = "~/combined.emat.mtx")
writeMM(t(scvelo_combined.nmat), file = "~/combined.nmat.mtx")
writeMM(t(scvelo_combined.smat), file = "~/combined.smat.mtx")
write.csv(scvelo_combined.obs, file = "~/combined.obs.csv", row.names = FALSE)
write.csv(scvelo_combined.obsm, file = "~/combined.obsm.csv", row.names = FALSE)


#########################################
## Now on in python

conda deactivate
conda activate scvelo-0.2.2

python

import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib




#import anndata
import scvelo as scv
#import matplotlib
#import scipy
import json
import os



And here is my code in python

sdata = ad.read_mtx("combined.spliced.mtx")
sdata.layers['spliced'] = sdata.X
ndata = ad.read_mtx("combined.unspliced.mtx")
sdata.layers['unspliced'] = ndata.X
#mdata = ad.read_mtx("~/combined.smat.mtx")
#sdata.layers['ambiguous'] = mdata.X
pd_obs = pd.read_csv("combined.obs.csv")
pd_obsm = pd.read_csv("combined.obsm.csv")
sdata.obs = pd_obs
sdata.obsm["X_umap"] = pd_obsm.values
sdata
AnnData object with n_obs × n_vars = 17829 × 14152 
obs: 'orig.ident', 'clusters'
obsm: 'X_umap'
layers: 'spliced', 'unspliced', 'ambiguous'

## Following https://scvelo.readthedocs.io/getting_started/
import scvelo as scv
#scv.set_figure_params()
adata = sdata

## from https://combine-lab.github.io/alevin-fry-tutorials/2021/alevin-fry-velocity/
adata.var_names_make_unique()

# get the proportion of spliced and unspliced count
scv.utils.show_proportions(adata)

# filter cells and genes, then normalize expression values
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000,enforce=True)

sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.tsne(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
#scv.pp.moments(adata)

## ref
adata2 = scv.datasets.pancreas()
sc.tl.pca(adata)
matplotlib.use('AGG')
scv.settings.set_figure_params('scvelo')

scv.utils.show_proportions(adata)

velo.out <- velociraptor::scvelo(sce, subset.row=top.hvgs, assay.X="spliced")


################################
## save as hd5
library(zellkonverter)

out_path <- tempfile(pattern = ".h5ad")
out_path <- "sce.h5ad"
zellkonverter::writeH5AD(sce, file = out_path)

library(scRNAseq)

sce_zeisel <- ZeiselBrainData()
out_path <- tempfile(pattern = ".h5ad")
writeH5AD(sce_zeisel, file = out_path)


##
################################

##############################
## save as anndata

library(reticulate)
anndata <- import("anndata", convert = FALSE)
adata <- anndata$AnnData(
  X = t(counts(sce_object)),
  obs = data.frame(colData(sce_object)),
  obsm  = list(
    "X_emb1" = as.matrix(reducedDim(sce_object, "emb1")),
    "X_emb2" = as.matrix(reducedDim(sce_object, "emb2"))
  )
)
anndata$AnnData$write(adata, 'filename.h5ad')

## done
#################################


## Preprocessing 
# https://github.com/edward130603/BayesSpace/issues/28
# Computing the size factors.


## Make 


txis <- sce
sizeFactors(sce) <- 2^rnorm(ncol(sce))
assays(txis) <- list(
  counts = mCounts,
  spliced = mSpliced,
  unspliced = mUnspliced
)

txis <- scater::logNormCounts(txis)
txis <- scater::runPCA(txis)
txis <- scater::runTSNE(txis, dimred = "PCA")

b <- data.matrix(OsC@reductions$umap@cell.embeddings)

reducedDims(txis) <- list(
    umap=data.matrix(OsC@reductions$umap@cell.embeddings),
    pca=data.matrix(OsC@reductions$pca@cell.embeddings),
    tsne=data.matrix(OsC@reductions$tsne@cell.embeddings)
)

# https://satijalab.org/loomr/loomr_tutorial

devtools::install_github(repo = "hhoeflin/hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

# Load loomR
library(loomR)

## write hd5 file 
zellkonverter::writeH5AD(sce, file = "test.loom")

scvelo(
  txis,
  assay.X = "counts",
  assay.spliced = "spliced",
  assay.unspliced = "unspliced"
)

txis <- txis[, colSums(counts(txis)) > 0]
velo.out <- scvelo(txis, subset.row=top.hvgs, assay.X="spliced")
velo.out2 <- scvelo(txis, assay.X=1, subset.row=top.hvgs, use.dimred="umap") 





##
###############################################################################

#############
## Next try 
# https://jef.works/blog/2020/08/25/using-scvelo-in-R-using-reticulate/
use_condaenv("scvelo-0.2.2", required = TRUE)
scv <- import("scvelo")
scv$logging$print_version()

adata <- zellkonverter::SCE2AnnData(sce, X_name = NULL, skip_assays = FALSE)
###############################################################################
## Run velocityraptor 

library(velociraptor)

sizeFactors(txis) <- scran::computeSumFactors(txis)
sizeFactors(txis) <- 2^rnorm(ncol(txis))
velo.out <- scvelo(txis, assay.X="spliced")
velo.out


txis <- scuttle::logNormCounts(txis, assay.type=1)

library(scran)
dec <- modelGeneVar(txis)
top.hvgs <- getTopHVGs(dec, n=2000)

library(velociraptor)
velo.out <- velociraptor::scvelo(sce, subset.row=top.hvgs, assay.X="spliced")
velo.out

## done 
##############33


#############################################################
## Read results and create spliced, unspliced loom files 


fry_to_spliced_unspliced_loom <- function(frydir, verbose = FALSE) {
  suppressPackageStartupMessages({
    library(rjson)
    library(Matrix)
    library(SingleCellExperiment)
  })
  
  # read in metadata
  qfile <- file.path(frydir, "quant.json")
  if (!file.exists(qfile)) {
    qfile <- file.path(frydir, "meta_info.json")
  }
  
  meta_info <- fromJSON(file = qfile)
  ng <- meta_info$num_genes
  usa_mode <- meta_info$usa_mode
  
  # if (usa_mode) {
  #   if (length(which_counts) == 0) {
  #     stop("Please at least provide one status in 'U' 'S' 'A' ")
  #   }
  #   if (verbose) {
  #     message("processing input in USA mode, will return ", paste(which_counts, collapse = '+'))
  #   }
  # } else if (verbose) {
  #   message("processing input in standard mode, will return spliced count")
  # }
  
  # read in count matrix
  af_raw <- readMM(file = file.path(frydir, "alevin", "quants_mat.mtx"))
  # if usa mode, each gene gets 3 rows, so the actual number of genes is ng/3
  if (usa_mode) {
    if (ng %% 3 != 0) {
      stop("The number of quantified targets is not a multiple of 3")
    }
    ng <- as.integer(ng/3)
  }
  
  # read in gene name file and cell barcode file
  afg <- read.csv(file.path(frydir, "alevin", "quants_mat_cols.txt"), 
                  strip.white = TRUE, header = FALSE, nrows = ng #, 
                  #col.names = c("gene_ids"), 
                  #row.names = 1
  )[1:ng,1]
  afc <- read.csv(file.path(frydir, "alevin", "quants_mat_rows.txt"), 
                  strip.white = TRUE, header = FALSE #,
                  #col.names = c("barcodes"), row.names = 1
  )
  
  
  rd <- list("S" = seq(1, ng), "U" =  seq(ng + 1, 2 * ng),
             "A" =  seq(2 * ng + 1, 3 * ng))
  
  ## spliced counts matrix ##
  which_counts <- c("S","A")
  splicedMat <- af_raw[, rd[[which_counts[1]]], drop = FALSE]
  for (wc in which_counts[-1]) {
    splicedMat <- splicedMat + af_raw[, rd[[wc]], drop = FALSE]
  }
  
  ## unspliced matrix ##
  which_counts <- c("U")
  unsplicedMat <- af_raw[, rd[[which_counts[1]]], drop = FALSE]
  for (wc in which_counts[-1]) {
    unsplicedMat <- unsplicedMat + af_raw[, rd[[wc]], drop = FALSE]
  }
  
  
  # create SingleCellExperiment object
  spliced <- t(splicedMat)
  row.names(spliced) <- afg
  colnames(spliced) <- as.vector(afc[,1])
  unspliced = t(unsplicedMat)
  row.names(spliced) <- afg
  colnames(spliced) <- as.vector(afc[,1])
  
  sce <- SingleCellExperiment(list(spliced = spliced, unspliced = unspliced),
                              colData =DataFrame(label=as.vector(afc[,1]))   ,
                              rowData = DataFrame(length=afg)
  )
  sce
}

dir.create("quants/CAM441A4/loomFile")

frydir <- "quants/CAM441A4/quant"
sce <- fry_to_spliced_unspliced_loom(frydir, verbose = FALSE)

sceasy::convertFormat(sce_object = sce, from="sce", to="loom",
                      outFile='quants/CAM441A4/loomFile/CAM441A4.loom')


## done 
#############################################################

#############################################################
## Add cell IDs based on Seurat object 

##
###########################################################

#############################################################
## Merge all sce objects into one

##
##############################################################

##################################################
## Run scvelo in velocitiraptor
# https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/

library(velociraptor)
velo.out <- scvelo(sce, subset.row=top.hvgs, assay.X="spliced")
velo.out <- scvelo(sce, assay.X="spliced")
velo.out


## 
###################################################

###
#Process splice matrix in R
conda deactivate

## List environments
conda info --envs

## Get spliced/unspliced matrices into R
conda activate R-4.0.2-BABS
# https://github.com/COMBINE-lab/alevin-fry#a-quick-start-run-through-on-sample-data

###############################################################################
## Function source                                                           ##

load_fry <- function(frydir, which_counts = c('S', 'A'), verbose = FALSE) {
  suppressPackageStartupMessages({
    library(rjson)
    library(Matrix)
    library(SingleCellExperiment)
  })
  
  # read in metadata
  qfile <- file.path(frydir, "quant.json")
  if (!file.exists(qfile)) {
    qfile <- file.path(frydir, "meta_info.json")
  }
  
  meta_info <- fromJSON(file = qfile)
  ng <- meta_info$num_genes
  usa_mode <- meta_info$usa_mode
  
  if (usa_mode) {
    if (length(which_counts) == 0) {
      stop("Please at least provide one status in 'U' 'S' 'A' ")
    }
    if (verbose) {
      message("processing input in USA mode, will return ", paste(which_counts, collapse = '+'))
    }
  } else if (verbose) {
    message("processing input in standard mode, will return spliced count")
  }
  
  # read in count matrix
  af_raw <- readMM(file = file.path(frydir, "alevin", "quants_mat.mtx"))
  # if usa mode, each gene gets 3 rows, so the actual number of genes is ng/3
  if (usa_mode) {
    if (ng %% 3 != 0) {
      stop("The number of quantified targets is not a multiple of 3")
    }
    ng <- as.integer(ng/3)
  }
  
  # read in gene name file and cell barcode file
  afg <- read.csv(file.path(frydir, "alevin", "quants_mat_cols.txt"), 
                  strip.white = TRUE, header = FALSE, nrows = ng #, 
                  #col.names = c("gene_ids"), 
                  #row.names = 1
                  )
  afc <- read.csv(file.path(frydir, "alevin", "quants_mat_rows.txt"), 
                  strip.white = TRUE, header = FALSE #,
                  #col.names = c("barcodes"), row.names = 1
                  )
  
  # if in usa_mode, sum up counts in different status according to which_counts
  if (usa_mode) {
    rd <- list("S" = seq(1, ng), "U" =  seq(ng + 1, 2 * ng),
               "A" =  seq(2 * ng + 1, 3 * ng))
    o <- af_raw[, rd[[which_counts[1]]], drop = FALSE]
    for (wc in which_counts[-1]) {
      o <- o + af_raw[, rd[[wc]], drop = FALSE]
    }
  } else {
    o <- af_raw
  }
  
  # create SingleCellExperiment object
  sce <- SingleCellExperiment(list(counts = t(o)),
                              colData = afc,
                              rowData = afg
  )
  sce
}

##                                                                           ##
###############################################################################

m <- load_fry("$AF_SAMPLE_DIR/quants/pbmc1k_v3/quant", which_counts=c('S', 'A'))
mSpliced <- load_fry("quants/CAM441A4/quant", which_counts=c('S', 'A'))
mUnspliced <- load_fry("quants/CAM441A4/quant", which_counts=c('U'))




###################################################
## Use scVelo/python
conda activate scvelo-0.2.2
python
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scvelo as scv
import matplotlib
import scipy
import json
import os

frydir = "quants/CAM441A4/quant"
e2n_path = "data/geneid_to_name.txt"
fpath = os.path.sep.join([frydir, "quant.json"])
if !os.path.exists(fpath):
  fpath = os.path.sep.join([frydir, "meta_info.json"])

meta_info = json.load(open(fpath))
ng = meta_info['num_genes']
usa_mode = meta_info['usa_mode']

if usa_mode:
  print("processing input in USA mode, will return A+S as the spliced count, and U as the unspliced count")
else:
  print("please follow previous steps to generate the ount matrix in the USA mode")
assert(False)

af_raw = sc.read_mtx(os.path.sep.join([frydir, "alevin", "quants_mat.mtx"]))
ng = int(ng/3)
e2n = dict([ l.rstrip().split() for l in open(e2n_path).readlines()])
var_names = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_cols.txt"])).readlines()][:ng]
var_names = [e2n[e] for e in var_names]

obs_names = [ l.rstrip() for l in open(os.path.sep.join([frydir, "alevin", "quants_mat_rows.txt"])).readlines() ]

x = af_raw.X
spliced = x[:,range(0,ng)] + x[:,range(2*ng,3*ng)]
unspliced = x[:,range(ng, 2*ng)]

## Running scVelo based on a Seurat object 
# https://github.com/basilkhuder/Seurat-to-RNA-Velocity#generating-loom-files

## Run scVelo
# create AnnData using spliced and unspliced count matrix
adata = anndata.AnnData(X = spliced,
                        layers = dict(spliced = spliced,
                                      unspliced = unspliced))
adata.obs_names = obs_names
adata.var_names = var_names
adata.var_names_make_unique()


# create AnnData using spliced and unspliced count matrix
# adata = anndata.AnnData(X = spliced,
#                         layers = dict(spliced = spliced,
#                                       unspliced = unspliced))
# adata.obs_names = obs_names
# adata.var_names = var_names
# adata.var_names_make_unique()

# get embeddings
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.tsne(adata)
#sc.tl.umap(adata, n_components = 2)

# housekeeping
matplotlib.use('AGG')
#scv.settings.set_figure_params('scvelo')

# get the proportion of spliced and unspliced count
scv.utils.show_proportions(adata)

# filter cells and genes, then normalize expression values
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000,enforce=True)

# scVelo pipeline
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata, n_jobs = 11)
scv.tl.velocity(adata, mode = 'dynamical')
scv.tl.velocity_graph(adata)
adata.write('test_full_dim_scvelo.h5ad', compression='gzip')
scv.pl.velocity_embedding_stream(adata, basis='X_umap', save="pancreas_full_dim.png")



###############################################
## Export adata object as Seurat object 
# source https://theislab.github.io/scanpy-in-R/
exprs <- t(py$adata$X)
colnames(exprs) <- py$adata$obs_names$to_list()
rownames(exprs) <- py$adata$var_names$to_list()
# Create the Seurat object
seurat <- CreateSeuratObject(exprs)
# Set the expression assay
seurat <- SetAssayData(seurat, "data", exprs)
# Add observation metadata
seurat <- AddMetaData(seurat, py$adata$obs)
# Add fetaure metadata
seurat[["RNA"]][["n_cells"]] <- py$adata$var["n_cells"]
# Add embedding
embedding <- py$adata$obsm["X_umap"]
rownames(embedding) <- py$adata$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")
seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")

## Done 
###########################################

########################################
## Create singel cell experiment 

sce <- SingleCellExperiment(
  assays      = list(logcounts = t(py$adata$X)),
  colData     = py$adata$obs,
  rowData     = py$adata$var,
  reducedDims = list(umap = py$adata$obsm["X_umap"])
)
sce

## Done
########################################


# get embeddings
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.tsne(adata)
sc.tl.umap(adata, n_components = 2)

# housekeeping
matplotlib.use('AGG')
scv.settings.set_figure_params('scvelo')

# get the proportion of spliced and unspliced count
scv.utils.show_proportions(adata)

# filter cells and genes, then normalize expression values
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000,enforce=True)

# scVelo pipeline
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata, n_jobs = 11)
scv.tl.velocity(adata, mode = 'dynamical')
scv.tl.velocity_graph(adata)
adata.write('pancreas_full_dim_scvelo.h5ad', compression='gzip')
scv.pl.velocity_embedding_stream(adata, basis='X_umap', save="pancreas_full_dim.png")




## Continue in python
module load Python/3.9.5-GCCcore-10.3.0
pip install scanpy


$ singularity exec --cleanenv \
--bind $AF_SAMPLE_DIR:/workdir \
--pwd /usefulaf/bash usefulaf.sif \
./simpleaf quant \
-1 /workdir/data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,/workdir/data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz \
-2 /workdir/data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,/workdir/data/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz \
-i /workdir/human_CR_3.0_splici/index \
-o /workdir/quants/pbmc1k_v3 \
-f u -c v3 -r cr-like \
-m /workdir/human_CR_3.0_splici/ref/transcriptome_splici_fl86_t2g_3col.tsv \
-t 16

  
##
###########################################


conda env list

conda activate R-4.0.2-BABS

###################################
## R

make_splici_txome <- function(gtf_path,
                              genome_path,
                              read_length,
                              flank_trim_length = 5,
                              output_dir,
                              extra_spliced=NULL,
                              extra_unspliced=NULL,
                              dedup_seqs=FALSE) {
  # if you get some error from .get_cds_IDX, please try to rerun the code again
  # read length is the scRNA-seq read length
  # flank trim length is used to avoid marginal case when dealing with junctional reads
  # assumes the following packages have been imported
  # eisaR, Biostrings, BSgenome, dplyr, stringr
  
  ########################################################################################################
  # Preprocessing
  ########################################################################################################
  
  suppressPackageStartupMessages({
    library(eisaR)
    library(Biostrings)
    library(BSgenome)
    library(stringr)
    library(GenomicFeatures)
  })
  
  if (!dir.exists(output_dir)) {
    dir.create(file.path(output_dir),recursive = TRUE, showWarnings = FALSE)
  }
  # make sure flank_length makes sense
  flank_length = read_length - flank_trim_length
  if (flank_length < 0 ){
    stop("flank trim length is larger than read length!")
  }
  # make sure gtf file exists
  if (!file.exists(gtf_path)) {
    stop("The following file does not exist: \n", gtf_path)
  }
  
  # make sure fasta file exists
  if (!file.exists(genome_path)) {
    stop("The following file does not exist: \n", genome_path)
  }
  
  # output file names
  file_name_prefix = paste0("transcriptome_splici_fl", flank_length)
  out_fa <- file.path(output_dir, paste0(file_name_prefix, ".fa"))
  out_t2g <- file.path(output_dir, paste0(file_name_prefix, "_t2g.tsv"))
  out_t2g3col <- file.path(output_dir, paste0(file_name_prefix, "_t2g_3col.tsv"))
  
  
  
  #########################################################################################################
  # Process gtf to get spliced and introns
  #########################################################################################################
  message("============processing gtf to get spliced and introns============")
  # fl is the flank length, here we set it to
  # the read length - 5
  grl <- suppressWarnings(getFeatureRanges(
    gtf = file.path(gtf_path),
    featureType = c("spliced", "intron"),
    intronType = "separate",
    flankLength = flank_length,
    joinOverlappingIntrons = TRUE,
    verbose = TRUE
  ))
  
  #########################################################################################################
  # Get spliced related stuffs
  #########################################################################################################
  
  # spliced ranges has no dash in it
  spliced_grl = grl[str_detect(names(grl), "-", negate = TRUE)]
  
  #########################################################################################################
  # Get reduced introns
  #########################################################################################################
  
  # identify all introns and convert to GRanges
  intron_gr = unlist(grl[str_detect(names(grl), "-")])
  # group introns by gene, then collapse ovelaping ranges!
  intron_grl = reduce(split(intron_gr, intron_gr$gene_id))
  
  # clean txp names and gene names
  intron_gr <- BiocGenerics::unlist(intron_grl)
  intron_gr$exon_rank <- 1L
  intron_gr$transcript_id <- word(names(intron_gr), 1, sep = '-')
  intron_gr$gene_id <- intron_gr$transcript_id
  intron_gr$type <- "exon"
  intron_gr$transcript_id <- make.unique(paste0(intron_gr$transcript_id, "-I"), sep = '')
  intron_gr$gene_id <- paste0(intron_gr$gene_id, "-I")
  intron_gr$exon_id <- intron_gr$transcript_id
  names(intron_gr) <- NULL
  mcols(intron_gr) <-
    S4Vectors::mcols(intron_gr)[, c("exon_id", "exon_rank",
                                    "transcript_id", "gene_id", "type")]
  # remake intron GRangesList
  intron_grl <- BiocGenerics::relist(intron_gr, lapply(
    structure(seq_along(intron_gr),
              names = intron_gr$transcript_id), function(i) i))
  
  
  #########################################################################################################
  # extract sequences from genome
  #########################################################################################################
  
  message("============extracting spliced and intron sequences from genome============")
  
  # load the genome sequence
  x <- Biostrings::readDNAStringSet(file.path(genome_path))
  # get the first word as the name
  names(x) <- word(names(x), 1)
  
  
  grl = c(spliced_grl, intron_grl)
  
  # make sure introns don't out of boundary
  seqlevels(grl) <- seqlevels(x)
  seqlengths(grl) <- suppressWarnings(seqlengths(x))
  grl <- trim(grl)
  
  seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = x,
    transcripts = grl
  )
  
  # If having duplicated sequences, only keep one
  if(dedup_seqs) {
    seqs = unique(seqs)
    grl = grl[names(seqs)]
  }
  
  
  # save some space
  rm(x)
  #########################################################################################################
  # process final outputs
  #########################################################################################################
  message("Writing outputs...")
  
  df <- getTx2Gene(grl)
  write.table(df, out_t2g, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  df <- df %>%
    dplyr::mutate(gene_id = word(gene_id, 1, sep = '-'),
                  status = ifelse(str_detect(transcript_id, '-'), 'U', 'S'))
  
  writeXStringSet(seqs, out_fa, format = "fasta")
  write.table(df, file.path(output_dir, paste0(file_name_prefix, "_t2g_3col.tsv")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  # optional: adding extra spliced and unspliced sequences from an fasta file
  if (!is.null(extra_spliced)) {
    if (!file.exists(extra_spliced)) {
      warning("provided extra_sequences file does not exist, will ignore it")
    } else {
      fa = file(extra_spliced, open="r")
      lns = readLines(fa)
      close(fa)
      for (ln in lns) {
        if (startsWith(ln, ">")) {
          # it is a header, write to t2g files and fasta file
          txp_name = gsub(">", "", ln)
          write.table(matrix(c(txp_name, txp_name), nrow = 1), file = out_t2g, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(matrix(c(txp_name, txp_name, "S"), nrow = 1), file = out_t2g3col, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        } else {
          # if not a header, just write to fasta file
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        }
      }
    }
  }
  
  if (!is.null(extra_unspliced)) {
    if (!file.exists(extra_unspliced)) {
      warning("provided extra_sequences file does not exist, will ignore it")
    } else {
      fa = file(extra_unspliced, open="r")
      lns = readLines(fa)
      close(fa)
      for (ln in lns) {
        if (startsWith(ln, ">")) {
          # it is a header, write to t2g file and fasta file
          txp_name = gsub(">", "", ln)
          write.table(matrix(c(txp_name, paste0(txp_name, "-U")), nrow = 1), file = out_t2g, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(matrix(c(txp_name, txp_name, "U"), nrow = 1), file = out_t2g3col, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        } else {
          # if not a header, just write to fasta file
          write.table(ln, file = out_fa, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
        }
      }
    }
  }
  
  
  
  
  message("Done.")
}

##################
##
##


suppressPackageStartupMessages({
  library(eisaR)
  library(Biostrings)
  library(BSgenome)
  library(stringr)
  library(GenomicFeatures)
})

ref_path <- file.path("/camp/svc/reference/Genomics/10x/10x_transcriptomes/GRCh38-3.0.0_GFP-tdT")
gtf_path = file.path(ref_path, "genes/genes.gtf")
genome_path = file.path(ref_path, "fasta/genome.fa")
read_length = 151
flank_trim_length = 5

output_dir = paste0("transcriptome_splici_fl", read_length - flank_trim_length)

make_splici_txome(gtf_path=gtf_path, 
                  genome_path=genome_path, 
                  read_length=read_length, 
                  flank_trim_length=flank_trim_length, 
                  output_dir=output_dir)

gffread refdata-cellranger-mm10-2.1.0/genes/genes.gtf -o refdata-cellranger-mm10-2.1.0/genes/genes.gff
$ grep "gene_name" refdata-cellranger-mm10-2.1.0/genes/genes.gff | cut -f9 | cut -d';' -f2,3 | sed 's/=/ /g' | sed 's/;/ /g' | cut -d' ' -f2,4 | sort | uniq > geneid_to_name.txt


## old below 
suppressPackageStartupMessages({
  library(eisaR)
  library(Biostrings)
  library(BSgenome)
  library(stringr)
  library(GenomicFeatures)
})

gtf_path = file.path( "/camp/svc/reference/Genomics/10x/10x_transcriptomes/GRCh38-3.0.0_GFP-tdT/genes/genes.gtf")
genome_path = file.path( "/camp/svc/reference/Genomics/10x/10x_transcriptomes/GRCh38-3.0.0_GFP-tdT/fasta/genome.fa")
read_length = 91
flank_trim_length = 5
output_dir = paste0("transcriptome_splici_fl", read_length - flank_trim_length)

make_splici_txome(gtf_path=gtf_path, 
                  genome_path=genome_path, 
                  read_length=read_length, 
                  flank_trim_length=flank_trim_length, 
                  output_dir=output_dir)

##
#####################



wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

ml purge
ml Anaconda3/2019.07
conda activate R-4.0.2-BABS

library(Biostrings)
library(BSgenome)
library(eisaR)
library(GenomicFeatures)
library(SummarizedExperiment)
library(tximeta)

library(rjson)
library(reticulate)
library(SingleCellExperiment)
library(scater)

##
gtf <- "/camp/svc/reference/Genomics/10x/10x_transcriptomes/GRCh38-3.0.0_GFP-tdT/genes/genes.gtf"

faFile <- "/camp/svc/reference/Genomics/10x/10x_transcriptomes/GRCh38-3.0.0_GFP-tdT/fasta/genome.fa"
#faFile <- "GRCh38.primary_assembly.genome.fa"
#cat GRCh38.primary_assembly.genome.fa | sed -r 's/(N[CTW]_[0-9]*)/chr\1/g' > genome.mod.modified.fa
#faFile <- "genome.mod.modified.fa"

cat /camp/svc/reference/Genomics/10x/10x_transcriptomes/GRCh38-3.0.0_GFP-tdT/genes/genes.gtf | sed -r 's/(N[CTW]_[0-9]*)/chr\1/g' > genome.modified.gtf

gtf <- "genome.modified.gtf"

grl <- eisaR::getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "intron"), 
  intronType = "separate", 
  flankLength = 90L, 
  joinOverlappingIntrons = FALSE, 
  verbose = TRUE
)



genome <- Biostrings::readDNAStringSet(
  faFile
)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl
)
Biostrings::writeXStringSet(
  seqs, filepath = "gencode.vM24.annotation.expanded.fa"
)