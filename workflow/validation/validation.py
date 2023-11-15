import os

from pathlib import Path
from Bio import SeqIO


class Validator:
    @staticmethod
    def string(value: str) -> bool:
        """
        Validates a string checking if it has not empty characters
        and checks if there are spaces
        """
        if not isinstance(value, str):
            return False
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
    def is_file(path: str) -> bool:
        return os.path.exists(path)

    @staticmethod
    def is_fasta(path: str) -> bool:
        if not Validator.is_file(path):
            return False
        with open(path, "r", encoding="utf-8") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)

    @staticmethod
    def is_genbank(path: str) -> bool:
        if not Validator.is_file(path):
            return False
        with open(path, "r", encoding="utf-8") as handle:
            gbk = SeqIO.parse(handle, "genbank")
            return any(gbk)

    @staticmethod
    def same_identifiers(fasta_path: str, gbk_path: str) -> bool:
        if not Validator.is_file(fasta_path) and not Validator.is_file(gbk_path):
            return False
        with open(fasta_path, "r", encoding="utf-8") as handle:
            fasta = SeqIO.parse(handle, "fasta")  # type: ignore
            fasta_ids: list[str] = [record.id for record in fasta]

        with open(gbk_path, "r", encoding="utf-8") as handle:
            gbk = SeqIO.parse(handle, "genbank")
            gbk_locus: list[str] = [record.name for record in gbk]

        return sorted(fasta_ids) == sorted(gbk_locus)

    @staticmethod
    def check_file_in_the_directory(directory: str, file_name: str) -> bool:
        if not Validator.is_file(directory + file_name):
            return False

        path = Path(directory)

        return not path.is_relative_to(Path(directory + file_name))
