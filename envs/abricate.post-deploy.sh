#!/bin/bash 
set -x
mkdir $CONDA_PREFIX/db/insaflu/
cp ./db/db_influenza_typing_v8.fasta $CONDA_PREFIX/db/insaflu/sequences
conda activate $CONDA_PREFIX
abricate --setupdb