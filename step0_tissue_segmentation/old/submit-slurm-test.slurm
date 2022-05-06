#!/bin/bash
#! SLURM job script for alma

#!##########################HOW TO RUN THIS SCRIPT #######################################################
#!#### To run on the first slide : sbatch --array=0 submit-slurm-test.slurm                         ######
#!#### To run on the 3rd slide : sbatch --array=2 submit-slurm-test.slurm                           ######
#!#### To run on slides 0-200 in parallel: sbatch --array=0-200 submit-slurm-test.slurm             ######
#!#### To run on slides 20-30 in parallel: sbatch --array=20-30 submit-slurm-test.slurm             ######
#!##########################HOR TO RUN THIS SCRIPT #######################################################

#!################################################################
#!#### Modify the options in this section as appropriate    ######
#!#### Change path to log files to your desired location	######
#!#### rds-cache path: /data/rds-cache/DMP/DUDMP/COPAINGE/  ######
#!################################################################

#SBATCH --job-name=til
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=7:00:00
#SBATCH --mail-user=yeman.hagos@icr.ac.uk
#SBATCH --no-requeue
#SBATCH --output=/data/scratch/DMP/DUDMP/COPAINGE/yhagos/Logs/%x_slurm.%A.%a.out
#SBATCH --error=/data/scratch/DMP/DUDMP/COPAINGE/yhagos/Logs/%x_slurm.%A.%a.err

# environment variables for singularity [NO NEED TO CHANGE]
#export SINGULARITY_CACHEDIR=/data/scratch/DMP/DUDMP/COPAINGE/.singularity/cache
#export SINGULARITY_TMPDIR=/data/scratch/DMP/DUDMP/COPAINGE/.singularity/tmp

#!###### you need to change this part and output/error Log paths (ABOVE) ONLY [No change in main.py] #####################
# code can be anywhere accessible to the cluster. This .slurm file is not necessary to be in same directory with code,
# as as long as you provide the path to the code. So, you can copy this .slurm file any where (your home dir, scratch or
# rds-cache) without the src code and refer the path for the code to run in "code_dir" bellow.
#Log path should be absolute path in alma, not in docker (make sure the the path you provide exists)
code_dir=/data/scratch/DMP/DUDMP/COPAINGE/yhagos/Vectra_Validation/mIF_pipeline/tissue_segmentation
data=/data/scratch/DMP/DUDMP/COPAINGE/yhagos/Vectra_Validation/data/Macrophages
output_dir=/data/scratch/DMP/DUDMP/COPAINGE/yhagos/Vectra_Validation/data/results/Macrophages/tissue_segmentation
mkdir -p $output_dir
#file_name_pattern='*'
#!#########################################################################################################################

# get a file to run based on file index from slurm array task id
#slides=$(find ${data} -maxdepth 1 -type d -name ${file_name_pattern})
#slides=$(find ${data} -maxdepth 1 -type f -name ${file_name_pattern})
slides=$(find ${data} -maxdepth 1 -type d -links 2)
echo $data
#echo $file_name_pattern
echo $slides
array=($slides)
full_path=${array[${SLURM_ARRAY_TASK_ID}]}
file_name=$(basename "${full_path}")

# code, input data and output data arguments should be absolute path in docker, not in alma (SEE --BIND below, how
# the data, output and code directories in alma are mounted in  /input, /output and /code directories in docker, respectively)
result_dir=/output
data_dir=/input/${file_name}

# display some details (can help for debugging)
echo "number of slides: ${#array[@]}"
echo "evaluating slide index:" ${SLURM_ARRAY_TASK_ID}
echo "evaluating slide name: ${file_name}"
echo "all slides: ${slides}"
echo "SLURM_ARRAY_TASK_ID:" ${SLURM_ARRAY_TASK_ID}
echo "current working directory:" ${PWD}
echo "created by:$USER"
echo "created on: $(date)"
echo "code_dir:$code_dir"
echo "data:$data"

echo "########################### paths in docker ####################"
echo "output_dir:$result_dir"
echo "slide evaluated:$data_dir"
echo "########################### paths in docker ####################"

singularity exec \
	--bind "$code_dir:/code" \
	--bind "$data:/input" \
	--bind "$output_dir:/output" \
	/data/scratch/DMP/DUDMP/COPAINGE/yhagos/Glioblastoma/mIF_pipeline/docker/mif.sif \
	/usr/local/bin/python /code/main.py -d $data_dir -o $result_dir

echo "DONE!!"
