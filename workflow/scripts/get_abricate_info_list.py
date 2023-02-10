"""Get info from abricate"""
import csv
import re
from collections import Counter
import sys


def get_abricate_info(path):
    """
    The get_abricate_info function takes in a path to an abricate file and returns
    the most common species and genus found in that file. The function first opens
    the csv file, reads it into a list of lists, then iterates through each entry
    to find all entries with "genus" or "species" as their header. If they have
    "genus", then they are added to a list of genera. If they have "species", then
    they are added to a list of species. After this is done for every
    entry, we use Counter from collections (imported at top)

    :param path:str: Specify the path to the abricate file
    :return: A list of two strings: the most common species and genus found in a given abricate file
    """

    with open(path, encoding="utf-8") as csv_handler:
        reader = csv.reader(csv_handler)
        info_list = list(reader)

    genus_list = []
    species_list = []

    for entry_list in info_list:
        entry = entry_list[0]
        if "genus" in entry:
            genus = re.findall("(?<=genus~~~)(.*?)(?=~~~)", entry)[0]
            genus_list.append(genus)

        elif "species" in entry:
            species = re.findall("(?<=species~~~)(.*?)(?=~~~)", entry)[0]
            species_list.append(species)
    if genus_list and species_list:
        final_genus = set(genus_list)
        final_species = set(species_list)
        return [final_species, final_genus]
    elif genus_list and not species_list:
        final_genus = set(genus_list)
        return ["", final_genus]
    elif not genus_list and species_list:
        final_species = set(species_list)
        return [final_species, ""]
    else:
        return ["", ""]


def write_to_file(output_path, info_list) -> None:
    """
    The write_to_file function writes the species and genus names to a file.
    The function takes two arguments: output_path, which is the path of the file we want to write
    to;
    and info_list, which is a list containing both species and genus names. The function returns
    None.

    :param output_path:str: Specify the path to the output file
    :param info_list:list[str]: Store the information that will be written to the file
    :return: None
    """
    species = info_list[0]
    genus = info_list[1]

    with open(output_path, "w", encoding="utf-8") as output_handler:
        output_handler.write(f"Species: {species}\n")
        output_handler.write(f"Genus: {genus}")


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    write_to_file(output_file, get_abricate_info(input_file))
