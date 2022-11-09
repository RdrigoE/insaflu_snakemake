"""Create coverage report"""
import re
import sys
import csv
from get_locus import get_locus


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
        sample_name = re.findall("(?<=main_result/)(.*?)(?=_coverage.csv)", filename)[0]
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
    coverage_list = coverage_files.split()

    n_locus = int(n_locus)
    coverage_translate = []
    for filename in coverage_list:
        print(get_coverage(filename, n_locus))
        coverage_translate.append(get_coverage(filename, n_locus))

    with open(output, mode="w", encoding="UTF8") as output_file:
        f_writer = csv.writer(output_file, delimiter=",")
        f_writer.writerows(coverage_translate)
        output_file.close()


if __name__ == "__main__":
    COVERAGE_FILES = sys.argv[1]
    OUTPUT_FILE = sys.argv[2]
    REFERENCE_GB = sys.argv[3]
    NUMBER_OF_LOCUS = len(get_locus(REFERENCE_GB))
    create_machine_readable_coverage_file(COVERAGE_FILES, OUTPUT_FILE, NUMBER_OF_LOCUS)
