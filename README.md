# insaflu_snakemake on linux

## Create environment 
install mamba forge:
    $ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
    $ bash Mambaforge-Linux-x86_64.sh
activate conda:
    $ conda activate base
create environment:
    $ mamba env create --name insaflu --file config/insaflu.yaml
## Activate Environment 
    $ conda activate insaflu

## Deactivate Enviroment
    $ conda deactivate insaflu