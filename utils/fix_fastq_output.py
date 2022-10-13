
import sys
import yaml 
import os
sample = sys.argv[1]


with open('config/config_run.yaml') as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

fastq1 = config_user['samples'][sample]['fastq1']
fastq2 = config_user['samples'][sample]['fastq2']

if fastq1 and fastq2:
    run = f'mv ./samples/{sample}/raw_fastqc/{fastq1}_fastqc.html ./samples/{sample}/raw_fastqc/{sample}_1_fastqc.html'
    os.system(f'{run}')
    run = f'mv ./samples/{sample}/raw_fastqc/{fastq2}_fastqc.html ./samples/{sample}/raw_fastqc/{sample}_2_fastqc.html'
    os.system(f'{run}')
    run = f'mv ./samples/{sample}/raw_fastqc/{fastq1}_fastqc.zip ./samples/{sample}/raw_fastqc/{sample}_1_fastqc.zip'
    os.system(f'{run}')
    run = f'mv ./samples/{sample}/raw_fastqc/{fastq2}_fastqc.zip ./samples/{sample}/raw_fastqc/{sample}_2_fastqc.zip'
    os.system(f'{run}')
elif fastq1:
    run = f'mv ./samples/{sample}/raw_fastqc/{fastq1}_fastqc.html ./samples/{sample}/raw_fastqc/{sample}_fastqc.html'
    os.system(f'{run}')
    run = f'mv ./samples/{sample}/raw_fastqc/{fastq1}_fastqc.zip ./samples/{sample}/raw_fastqc/{sample}_fastqc.zip'
    os.system(f'{run}')