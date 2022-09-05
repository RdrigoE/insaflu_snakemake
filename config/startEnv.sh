#!/bin/bash

conda activate read-quality-analysis
snakemake --dag | dot -Tpng > dag.png