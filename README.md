# INSaFLU snakemake - UNDER DEVELOPMENT

## Create Environment 

### Install Mamba Forge:
    
   ```$ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh```
    
   ```$ bash Mambaforge-Linux-x86_64.sh```

### Activate Conda:
    
   ```$ conda activate base```

### Create Environment:
    
   ```$ mamba env create --name insaflu --file config/insaflu.yaml```

## Activate Environment 
   ```$ conda activate insaflu```

## Deactivate Enviroment
   ```$ conda deactivate insaflu```

## Change directory to user, create reads directory and move your reads to this directory
   ```
   $ cd user
   $ mkdir reads
   $ cp path/to/your/files/your_files.fastq.gz path/user/reads/
   $ cd ..
   ```
## Edit the csv in the folder user called sample_metadata.csv
   - This csv has to have the following format: 
      sample_name,fastq1,fastq2,tech
      your_sample_name,file_with_fastq1, file_with_fastq2,illumina or ont

   - All samples much be in fastq.gz format (do not append .fastq.gz in the fastq1 and fastq2 fields)
## Go to user/project_info.yaml and fill the following fields
 - only_samples
 - project_name
 - fasta_reference
 - gb_reference

## Change the parameters in the run in user/parameters.yaml
   - Here is possible to change the parameters in some softwares, change at your own responsability. A copy of the default parameters are in config/default_parameters.yaml
   - Each project will recive a copy of this file inside projects/project_name_folder/
## The references are in user/references folder
To add new ones just copy your fasta and genbank file.
### You must keep the fasta identifier equal to the GenBank Locus
   ">ON563414" == "LOCUS       ON563414"

## Run your project or analyse your samples with the following command:
   `$ cd workflow`
   
   `$ snakemake -c {threads} --use-conda`
   If it doesn't work there is an alternative:
   `$ snakemake -s workflow/Snakefile --use-conda --cores {threads}`
