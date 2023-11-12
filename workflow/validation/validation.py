# only_samples: bool
# project_name: str (no_spaces)
# fasta_reference: valid_fasta
# gb_reference: valid_genbank
# abricate: bool
# illumina_consensus: snippy | iVar
# primers_fasta: valid fasta and check if there is the other file
# get_minor_variants: bool

from pathlib import Path
from Bio import SeqIO


class Validator:
    @staticmethod
    def string(value: str) -> bool:
        """
        Validates a string checking if it has not empty characters
        and checks if there are spaces
        """
        value = value.strip()
        if len(value) == 0 or value.find(" ") != -1:
            return False
        return True

    @staticmethod
    def includes(value: str, value_list: list[str]) -> bool:
        """
        Checks if a value is equal to a string entry in a list
        """
        for hit in value_list:
            if value == hit:
                return True
        return False

    @staticmethod
    def is_fasta(path: str) -> bool:
        """
        The is_fasta function checks if a file is in FASTA format.

        :param filename: str: Specify the file to be checked
        :return: True if the file is a fasta file, and false otherwise
        :doc-author: Trelent
        """
        with open(path, "r", encoding="utf-8") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)

    @staticmethod
    def is_genbank(path: str) -> bool:
        """
        The is_gbk function checks if a file is in GenBank format.

        :param filename: str: Specify the file to be checked
        :return: True if the file is a genbank file, and false otherwise
        :doc-author: Trelent
        """
        with open(path, "r", encoding="utf-8") as handle:
            gbk = SeqIO.parse(handle, "genbank")
            return any(gbk)

    @staticmethod
    def same_identifiers(fasta_path: str, gbk_path: str) -> bool:
        """
        The same_identifiers function checks whether the IDs in a FASTA file match those
        in a GenBank file. It does this by opening both files and parsing them with
        the SeqIO module from Biopython. The function returns True if the lists of IDs
        are identical, and False otherwise.

        :param fasta_file: str: Specify the path to the fasta file
        :param gbk_file: str: Specify the path to the genbank file
        :return: True if the order of the ids in both files are equal
        """
        with open(fasta_path, "r", encoding="utf-8") as handle:
            fasta = SeqIO.parse(handle, "fasta")  # type: ignore
            fasta_ids: list[str] = [record.id for record in fasta]

        with open(gbk_path, "r", encoding="utf-8") as handle:
            gbk = SeqIO.parse(handle, "genbank")
            gbk_locus: list[str] = [record.name for record in gbk]

        return sorted(fasta_ids) == sorted(gbk_locus)

    @staticmethod
    def check_file_in_the_directory(directory: str, file_name: str) -> bool:
        path = Path(directory)
        return not path.is_relative_to(Path(directory + file_name))
