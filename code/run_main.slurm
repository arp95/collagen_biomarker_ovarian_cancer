#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=03-15:00:00
#SBATCH --mail-user=$USER@case.edu
#SBATCH --job-name="features9"
#SBATCH --output="log.pipeline.features9"


#=========loading libraries==========#
module swap intel gcc
module load python/3.7.0
module load matlab
#pip3.6 install --upgrade setuptools pip --user
#pip3.6 install open slide-python --user
#pip3.6 install opencv-python --user


cd /mnt/rstor/CSE_BME_AXM788/home/axa1399/collagen_biomarker_ovarian_cancer/code/
time matlab -nodisplay -r main_2