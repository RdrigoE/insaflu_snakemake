"""Handle errors from quality analysis step with fastqc"""
import sys
import os
from yaml_io import read_yaml

if __name__ == "__main__":
    dic_directory: dict[str, str] = read_yaml("../config/constants.yaml")

    SAMPLE = sys.argv[1]

    config_user = read_yaml(f"{dic_directory['config_file']}")

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
