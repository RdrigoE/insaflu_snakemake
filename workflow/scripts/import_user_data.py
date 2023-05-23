"""Import sample data from user"""
from typing import Optional
import pandas
from scripts.yaml_io import read_yaml

dic_directory = read_yaml("config/constants.yaml")


class Data:
    """Object with sample data information"""

    def __init__(self, file: str, consensus_illumina) -> None:
        self.user_df: pandas.DataFrame = pandas.read_csv(file)
        self.consensus_illumina = consensus_illumina

    def get_sample_names(self) -> list[str]:
        """
        The get_sample_names function returns a list of sample names from the
        user_df dataframe.


        :param self: Access the attributes and methods of the class in python
        :return: A list of the sample names in the user_df
        :doc-author: Trelent
        """
        return list(self.user_df["sample_name"])

    def get_sample_1(self) -> list[str]:
        """
        The get_sample_1 function returns a list of the fastq files associated with
        the sample.  This is accomplished by returning the values in the &quot;fastq&quot; column
        of a pandas dataframe.

        :param self: Access the attributes and methods of the class in python
        :return: A list of the fastq file names from the user_df dataframe
        :doc-author: Trelent
        """
        return list(self.user_df["fastq1"])

    def get_sample_2(self) -> list[Optional[str]]:
        """
        The get_sample_2 function returns a list of the second fastq file for each sample.
        This is useful when you want to run a program that requires paired-end reads.

        :param self: Access the attributes and methods of the class in python
        :return: A list of the values in the column &quot;fastq2&quot; from the user_df
        :doc-author: Trelent
        """
        return list(self.user_df["fastq2"])

    def get_sample_type(self) -> dict[str, str]:
        """
        The get_sample_type function returns a dictionary of sample names and their corresponding
        sequencing technology. The function is called by the get_sample_metadata function to create
        a dictionary of metadata for each sample.

        :param self: Reference the class instance
        :return: A dictionary with sample names as keys and the corresponding assembler type as values
        :doc-author: Trelent
        """
        names = self.get_sample_names()
        tech_type: pandas.Series = self.user_df["tech"]
        type_dic: dict[str, str] = {}
        for index, _ in enumerate(names):
            type_dic[names[index]] = (
                "medaka" if tech_type[index] == "ont" else self.consensus_illumina
            )
        return type_dic

    def get_dic(self) -> dict[str, dict[str, Optional[str]]]:
        """
        The get_dic function creates a dictionary of dictionaries. The outer dictionary is keyed
        by sample names, and the inner
        dictionary contains keys for each column in the user_df dataframe. For example:
        {'sample_name': {'fastq2': 'path/to/file', 'fastq2': 'path/to/file', ...}, ...}

        :param self: Refer to the object that is calling the function
        :return: A dictionary with the sample names as keys and a subdictionary of fastq files,
        tech, and species
        :doc-author: Trelent
        """
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
        """
        The get_options function returns a list of dictionaries.
        The first dictionary is the paired illumina samples, the second is th
        single illumina samples, and the third is ont_samples. The fourth
        dictionary contains all of the information for each sample in a
        project and will be used to create yaml files for Snakemake.
        Lastly, it :param self: Access the attributes and methods of
        the class in python
        :return: A list of dictionaries
        :doc-author: Trelent
        """
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
        """
        The get_assembler function returns the assembler used in the pipeline.


        :param self: Access the attributes and methods of the class in python
        :return: The assembler that is being used in the current instance of the class
        :doc-author: Trelent
        """
        return self.consensus_illumina


def get_data_in_align_form(
    assembler: str,
    single_illumina: dict[str, Optional[str]],
    paired_illumina: dict[str, Optional[str]],
) -> list[dict[str, Optional[str]]]:
    """
    The get_data_in_align_form function takes in the assembler, single_illumina
    and paired_illumina dictionaries. It then checks to see if the assembler is
    snippy or iVar. If it is snippy, it returns the paired illumina dictionary 
    as well as both single illumina dictionaries.
    If it is iVar, then it returns both paired illumina dictionaries and both 
    single illumina dictionaries.

    :param assembler:str: Determine which assembler to use
    :param single_illumina:dict[str: Store the name of the file that contains the single-end
    illumina reads
    :param Optional[str]]: Indicate that the parameter is optional
    :param paired_illumina:dict[str: Store the paired illumina reads
    :param Optional[str]]: Indicate that the value can be none
    :param : Determine which assembler to use
    :return: A list of four dictionaries
    :doc-author: Trelent
    """
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
