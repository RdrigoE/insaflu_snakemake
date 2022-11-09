"""Handle errors from quality analysis step with fastqc"""
import sys
import os
import yaml

SAMPLE = sys.argv[1]


with open("config/config_run.yaml", encoding="UTF8") as file:
    config_user = yaml.load(file, Loader=yaml.FullLoader)

fastq1 = config_user["samples"][SAMPLE]["fastq1"]
fastq2 = config_user["samples"][SAMPLE]["fastq2"]

if fastq1 and fastq2:
    run = f"mv ./samples/{SAMPLE}/raw_fastqc/{fastq1}_fastqc.html ./samples/{SAMPLE}/raw_fastqc/{SAMPLE}_1_fastqc.html"
    os.system(f"{run}")
    run = f"mv ./samples/{SAMPLE}/raw_fastqc/{fastq2}_fastqc.html ./samples/{SAMPLE}/raw_fastqc/{SAMPLE}_2_fastqc.html"
    os.system(f"{run}")
    run = f"mv ./samples/{SAMPLE}/raw_fastqc/{fastq1}_fastqc.zip ./samples/{SAMPLE}/raw_fastqc/{SAMPLE}_1_fastqc.zip"
    os.system(f"{run}")
    run = f"mv ./samples/{SAMPLE}/raw_fastqc/{fastq2}_fastqc.zip ./samples/{SAMPLE}/raw_fastqc/{SAMPLE}_2_fastqc.zip"
    os.system(f"{run}")
elif fastq1:
    run = f"mv ./samples/{SAMPLE}/raw_fastqc/{fastq1}_fastqc.html ./samples/{SAMPLE}/raw_fastqc/{SAMPLE}_fastqc.html"
    os.system(f"{run}")
    run = f"mv ./samples/{SAMPLE}/raw_fastqc/{fastq1}_fastqc.zip ./samples/{SAMPLE}/raw_fastqc/{SAMPLE}_fastqc.zip"
    os.system(f"{run}")
