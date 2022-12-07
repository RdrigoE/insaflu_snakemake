# insaflu_snakemake

## Create environment 

### install mamba forge:
    
   `$ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh`
    
   `$ bash Mambaforge-Linux-x86_64.sh`

### activate conda:
    
   `$ conda activate base`

### create environment:
    
   `$ mamba env create --name insaflu --file config/insaflu.yaml`

## Activate Environment 
   `$ conda activate insaflu`

## Deactivate Enviroment
   `$ conda deactivate insaflu`

## Create new folder named 'user_data'
   `$ mkdir data`
   
   `$ cp path/to/your/files/your_files path/to/data`
## Create a csv in the folder config_user called sample_info.csv
 This csv has to have the following format: 
   sample_name,fastq1,fastq2,tech
   your_sample_name,file_with_fastq1, file_with_fastq2,illumina/ont

 All samples much be in fastq.gz format (do not append .fastq.gz in the fastq1 and fastq2 fields)
## Go to user_metadata/sample_metadata.yaml and fill the following fields
 - only_samples
 - project_name
 - fasta_reference
 - gb_reference

## Change the parameters in the run in user_metadata/parameters.yaml
   - Here is possible to change the parameters in some softwares, change at your own responsability. A copy of the default parameters are in config/default_parameters.yaml
   - Each project will recive a copy of this file inside projects/project_name_folder/

## You must keep the fasta identifier equal to the GenBank Locus
">ON563414" == "LOCUS       ON563414"

## Run your project or analyse your samples with the following command:
   `$ snakemake -c {threads} --use-conda`

