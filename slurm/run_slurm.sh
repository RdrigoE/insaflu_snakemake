#SBATCH --output /path/to/insaflu_snakemake/profile-%j.out
#SBATCH --error /path/to/insaflu_snakemake/profile-%j.err
#SBATCH -p {partition}
#SBATCH -q {queue}

# activate conda in general
source path/to/.bashrc # to conda init setting

# activate a specific conda environment, if you so choose
conda activate /path/to/mambaforge/envs/insaflu

# go to a particular directory
cd /path/to/insaflu_snakemake

# make things fail on errors
set -o nounset
set -o errexit
set -x

### run your commands here!

snakemake --profile slurm
