"""Create coverage report"""
import re
import sys
import csv
from yaml_io import read_yaml
from extract_gb_info import get_locus, get_locus_len

dic_directory = read_yaml("../config/constants.yaml")


def merge_coverage(file_list, output_file):
    """
    The merge_coverage function takes a list of coverage files and merges them into one file.
    The first argument is the list of files, the second argument is the name for the output file.

    :param file_list: Store the list of files that will be merged
    :return: A list of lists with the coverage information for each sample
    :doc-author: Trelent
    """

    with open(file_list[0], "r", encoding="utf8") as file_handler:
        length = file_handler.readlines()[3].split()[1]

    species = get_locus(REFERENCE_GB)
    length = get_locus_len(REFERENCE_GB)
    super_header = ["Name", species, "", ""]
    super_header_1 = ["Length", length, "", ""]

    software_parameters = read_yaml(dic_directory["software_parameters"])

    coverage_value = software_parameters["min_coverage_consensus"]
    header = [
        "SAMPLES",
        "Mean depth of coverage",
        f"{'%'} of size covered by at least 1-fold",
        f"% of size covered by at least {coverage_value}-fold",
    ]
    # colocar os segmentos para fazer sentido com a flu
    final_output = [super_header, super_header_1, header]

    for file in file_list:
        with open(file, "r", encoding="utf8") as f:
            sample_name = re.findall("(?<=/)(.*?)(?=_coverage.csv)", file)[0].split(
                "/"
            )[0]
            info = f.readlines()[6].split()
            info[0] = sample_name
            final_output.append(info)

    with open(output_file, mode="w", encoding="utf8") as f:
        f_writer = csv.writer(
            f, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
        )
        f_writer.writerows(final_output)
        f.close()


def get_coverage(filename: str, n_locus: int):
    """
    The get_coverage function takes a filename as an argument and returns the coverage of each
    locus in that file.
    The function also takes n_locus as an argument, which is the number of loci in the file.

    :param filename: Specify the file to be read
    :param n_locus: Specify the number of loci in the file
    :return: A list of the coverage for each locus in a given sample
    :doc-author: Trelent
    """
    with open(filename, newline="", encoding="UTF8") as csvfile:
        csv_reader = csv.reader(csvfile, delimiter="\t")
        coverage_list = []
        sample_name = re.findall("(?<=/)(.*?)(?=_coverage.csv)", filename)[0].split(
            "/"
        )[0]
        for i in csv_reader:
            coverage_list.append(i)
        coverage_list = coverage_list[-1][-n_locus:]
        coverage_list.insert(0, sample_name)
        csvfile.close()
    return coverage_list


def create_machine_readable_coverage_file(coverage_files, output, n_locus):
    """
    The create_machine_readable_coverage_file function takes a list of coverage files and outputs
    a machine readable coverage file. The output is the name of the file you want to create.
    The last argument is n_locus which specifies how many loci are in each coverage file.

    :param coverage_files: Specify the files that contain the coverage data
    :param output: Specify the name of the output file
    :param n_locus: Indicate the number of locus in the genome
    :return: A list of lists, where each sublist contains the coverage each locus
    :doc-author: Trelent
    """
    coverage_list = coverage_files

    n_locus = int(n_locus)
    coverage_translate = []
    for filename in coverage_list:
        coverage_translate.append(get_coverage(filename, n_locus))

    with open(output, mode="w", encoding="UTF8") as output_file:
        f_writer = csv.writer(output_file, delimiter=",")
        f_writer.writerows(coverage_translate)
        output_file.close()


if __name__ == "__main__":
    COVERAGE_FILES = sys.argv[1].split()
    OUTPUT_FILE_1 = sys.argv[2]
    OUTPUT_FILE_2 = sys.argv[3]
    REFERENCE_GB = sys.argv[4]
    NUMBER_OF_LOCUS = len(get_locus(REFERENCE_GB))
    merge_coverage(COVERAGE_FILES, OUTPUT_FILE_1)
    create_machine_readable_coverage_file(
        COVERAGE_FILES, OUTPUT_FILE_2, NUMBER_OF_LOCUS
    )
