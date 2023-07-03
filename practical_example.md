# INSAFLU Snakemake Software Documentation

## Table of Contents

- [Installation Process](#installation-process)
- [Activate Environment](#activate-environment)
- [Explore Folder Structure](#explore-folder-structure)
- [Execution](#execution)
- [SLURM Execution](#execution-in-slurm-cluster)
- [Exploring the Results](#exploring-the-results)

## Installation Process

Download and install Mambaforge:

```shell
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```
## Activate Environment

Activate the Conda environment and create a new environment for INSAFLU:

```shell
conda activate base
mamba env create --name insaflu --file config/insaflu.yaml
conda activate insaflu
```
## Explore Folder Structure

Go to the user/reads directory and verify the presence of the files:
``` 
        demo_SARSCoV2_001_1P.fastq.gz
        demo_SARSCoV2_002_1P.fastq.gz
        demo_SARSCoV2_100.fastq.gz
```
Verify the configuration files for the workflow in user/samples.csv. Note that the fastq1 and fastq2 entries should not have the .fastq.gz extension. The sample names must be unique.

Check the project configuration in user/project.yaml. Here, you can change the project name, choose the type of run, change references, and modify SNP calling and consensus assembly software.

Examine user/parameters.yaml for specific parameters in various steps of the pipeline. Only change them if you know what you're doing.

Visit config/threads.yaml to change the number of threads used in each rule.

Go to config/max_memory.yaml to define the necessary memory for each rule.

Check constants.yaml. This file ties everything together by routing all the files described above to the workflow.

## Execution

To execute the workflow locally, run the following command:

```shell
snakemake -c 8 --use-conda
```

## Execution in SLURM cluster

Change the execution to the SLURM execution by using two files:

- Edit the settings for SLURM execution in slurm/config.yaml.
- Edit your paths in slurm/run_slurm.sh and run with:

```shell
sbatch slurm/run_slurm.sh
```

## Exploring the Results

The results are in the results directory.
```shell
samples
align_samples
projects
```


results/samples contains directories for each sample with FastQC results and trimmed reads.

results/align_samples contains the results of SNP variant calling and consensus assembly.

results/projects contains project-related files and align_sample files. The main files are in the main_result folder, which includes alignments, variants, and coverage files.


