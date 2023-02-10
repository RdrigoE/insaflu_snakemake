"""Yaml interactions"""
import yaml


def read_yaml(yaml_file_path: str) -> dict:
    """
    The read_yaml function reads a yaml file and returns the contents as a dictionary.
    The read_yaml function accepts one argument, which is the path to the yaml file.

    :param yaml_file_path: Specify the path to the yaml file that is going to be read
    :return: A dictionary
    :doc-author: Trelent
    """
    with open(yaml_file_path, encoding="utf-8") as yaml_file:
        return yaml.load(yaml_file, Loader=yaml.FullLoader)


def write_yaml(yaml_file_path: str, dump_dict: dict) -> None:
    """
    The write_yaml function writes a dictionary to a yaml file.

    :param yaml_file_path:str: Specify the path to the yaml file that will be written
    :param dump_dict:dict: Specify the dictionary that is to be written to a yaml file
    :return: None
    :doc-author: Trelent
    """
    with open(yaml_file_path, "w", encoding="utf-8") as file:
        yaml.dump(dump_dict, file)
