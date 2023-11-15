"""Import sample data from user"""
from typing import Literal, Optional

import pandas as pd
from scripts.yaml_io import read_yaml

dic_directory = read_yaml("config/constants.yaml")


class Data:
    """Object with sample data information"""

    def __init__(
        self, file: str, consensus_illumina: Literal["iVar"] | Literal["snippy"]
    ) -> None:
        self.user_df = pd.read_csv(file)
        self.consensus_illumina = consensus_illumina

    def get_sample_names(self) -> list[str]:
        return list(self.user_df["sample_name"])

    def get_sample_1(self) -> list[str]:
        return list(self.user_df["fastq1"])

    def get_sample_2(self) -> list[Optional[str]]:
        return list(self.user_df["fastq2"])

    def get_sample_type(self) -> dict[str, str]:
        names = self.get_sample_names()
        tech_type = self.user_df["tech"]
        type_dic: dict[str, str] = {}
        for index, _ in enumerate(names):
            type_dic[names[index]] = (
                "medaka" if tech_type[index] == "ont" else self.consensus_illumina
            )
        return type_dic

    def get_dic(self) -> dict[str, dict[str, Optional[str]]]:
        dic: dict = {}
        names = self.get_sample_names()
        for name in names:
            dic[name] = {}

        for idx, sample_name_1 in enumerate(self.get_sample_1()):
            dic[names[idx]]["fastq1"] = (
                sample_name_1 if isinstance(sample_name_1, str) else None
            )

        for idx, sample_name_2 in enumerate(self.get_sample_2()):
            dic[names[idx]]["fastq2"] = (
                sample_name_2 if isinstance(sample_name_2, str) else None
            )

        for idx, tech in enumerate(list(self.user_df["tech"])):
            dic[names[idx]]["tech"] = tech if isinstance(tech, str) else None

        return dic

    def get_options(self) -> list[Optional[dict[str, dict[str, Optional[str]]]]]:
        project_dic = self.get_dic()
        final_dictionary: Optional[dict[str, dict[str, Optional[str]]]] = {}
        single_illumina: Optional[dict[str, dict[str, Optional[str]]]] = {}
        paired_illumina: Optional[dict[str, dict[str, Optional[str]]]] = {}
        ont_samples = {}

        for names, _ in project_dic.items():
            if (
                project_dic[names]["fastq1"]
                and project_dic[names]["fastq2"]
                and project_dic[names]["tech"] == "illumina"
            ):
                paired_illumina[names] = project_dic[names]
                final_dictionary[names] = project_dic[names]
            elif (
                project_dic[names]["fastq1"]
                and not project_dic[names]["fastq2"]
                and project_dic[names]["tech"] == "illumina"
            ):
                single_illumina[names] = project_dic[names]
                final_dictionary[names] = project_dic[names]
            elif (
                project_dic[names]["fastq1"]
                and not project_dic[names]["fastq2"]
                and project_dic[names]["tech"] == "ont"
            ):
                ont_samples[names] = project_dic[names]
                final_dictionary[names] = project_dic[names]
        return [
            paired_illumina,
            single_illumina,
            ont_samples,
            final_dictionary,
        ]

    def get_assembler(self) -> str:
        return self.consensus_illumina


def get_data_in_align_form(
    assembler: str,
    single_illumina: dict[str, Optional[str]],
    paired_illumina: dict[str, Optional[str]],
) -> list[dict[str, Optional[str]]]:
    paired_illumina_snippy: dict[str, Optional[str]] = {}
    single_illumina_snippy: dict[str, Optional[str]] = {}
    paired_illumina_ivar: dict[str, Optional[str]] = {}
    single_illumina_ivar: dict[str, Optional[str]] = {}

    if assembler == "snippy":
        paired_illumina_snippy = paired_illumina
        single_illumina_snippy = single_illumina
        paired_illumina_ivar = {}
        single_illumina_ivar = {}
    elif assembler == "iVar":
        paired_illumina_ivar = paired_illumina
        single_illumina_ivar = single_illumina
        paired_illumina_snippy = {}
        single_illumina_snippy = {}
    return [
        paired_illumina_snippy,
        single_illumina_snippy,
        paired_illumina_ivar,
        single_illumina_ivar,
    ]
