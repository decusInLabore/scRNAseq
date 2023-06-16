wait_on_lsf() { ## wait on jobs{
sleep 300
n=`squeue --name=$project  | wc -l`
while [ $n -ne 1 ]
do
n=`squeue --name=$project  | wc -l`
((z=$n-1))
#number of running
echo "$project jobs running: $z"
#number of pending
sleep 300
done
}

project_id=511
project=rA_$project_id
cd ./analyses/QC/
#sbatch --time=24:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r QC.Rmd" --job-name=$project -p hmem --mem=600G -o ../../../../workdir/rA.$project.slurm >> ../commands.txt
sbatch --time=24:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r QC.Rmd" --job-name=$project --mem=200G -o ../../../../workdir/rA.$project.slurm >> ../commands.txt
cd -

wait_on_lsf
project=rB_$project_id
cd ./analyses/Main_Analysis/
#sbatch --time=72:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r Main_Analysis.Rmd" --job-name=$project -p hmem --mem=1500G -o ../../../../workdir/rB.$project.slurm >> ../commands.txt
sbatch --time=72:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r Main_Analysis.Rmd" --job-name=$project --mem=200G -o ../../../../workdir/rB.$project.slurm >> ../commands.txt
cd -

wait_on_lsf
cd ./analyses/Cluster_Finder/
project=rCF_$project_id
sbatch --time=24:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r ClusterFinder.Rmd" --job-name=$project --mem=200G -o ../../../../workdir/rCF.$project.tslurm >> ../commands.txt
cd -

wait_on_lsf
project=rCTA_$project_id
cd ./analyses/Celltype_Assignment/
sbatch --time=24:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r Celltype_Assignment.Rmd" --job-name=$project --mem=200G -o ../../../../workdir/rCTA.$project.tslurm >> ../commands.txt
cd -

wait_on_lsf
project=rC_$project_id
cd ./analyses/Main_Analysis/
# sbatch --time=24:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r sc_PartC_Analysis.Rmd " --job-name=$project --mem=200G -o ../rC.$project.slurm >> commands.txt
sbatch --time=24:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r sc_PartC_Analysis.Rmd " --job-name=$project --mem=200G -o ../../../../workdir/rC.$project.slurm >> commands.txt
cd -

project=rIntTest_$project_id
cd ./analyses/Integration_Benchmarking
sbatch --time=24:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r Integration_Benchmarking.Rmd " --job-name=$project --mem=200G -o ../rC.$project.slurm >> commands.txt
cd -

wait_on_lsf
project=rDGE_$project_id
cd ./analyses/DGE
sbatch --time=12:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r DGE_Analysis.Rmd " --job-name=$project --mem=200G -o ../rDGE.$project.slurm >> commands.txt
cd -

wait_on_lsf
project=rT_$project_id
sbatch --time=24:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript runRmd.r sc_Trajectory_Analysis_Module.Rmd " --job-name=$project -p hmem --mem=300G -o ../rT.$project.slurm >> commands.txt


sbatch --time=18:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;Rscript runCM.r" --job-name=$project -p hmem --mem=300G -o rCN.slurm >> commands.txt

project=rB446
#sbatch --time=18:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;Rscript runB.r" --job-name=$project -p hmem --mem=300G -o rB.slurm >> commands.txt
sbatch --time=12:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;Rscript runRmd.r sc_PartB_Analysis.Rmd" --job-name=$project --mem=200G -o $project.slurm >> commands.txt


project=r359sub7
sbatch --time=18:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;Rscript runDGE.r" --job-name=$project -p hmem --mem=300G -o rDGE.slurm >> commands.txt

#project=rC409
#sbatch --time=12:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;Rscript runC.r" --job-name=$project  --mem=100G -o rC.slurm >> commands.txt

sbatch --time=12:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;Rscript runC.r" --job-name=$project -p hmem --mem=300G -o rC.slurm >> commands.txt

project=rT446
sbatch --time=12:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;Rscript runT.r" --job-name=$project -p hmem --mem=300G -o rA.slurm >> commands.txt
sbatch --time=12:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/3.6.0-foss-2016b-BABS;Rscript runRmd.r sc_Trajectory_Analysis.Rmd" --job-name=$project -p hmem --mem=300G -o rT.slurm >> commands.txt

## Run subsetting script
project=subSet
sbatch --time=04:00:00 --wrap "module purge;source /camp/stp/babs/working/software/modulepath_new_software_tree_2018-08-13;module load pandoc/2.2.3.2-foss-2016b;ml R/4.0.3-foss-2020a;Rscript subsetting.dataset.for.subclustering.R" --job-name=$project --mem=200G -o rSub.slurm >> commands.txt
