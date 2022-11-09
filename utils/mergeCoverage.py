import re
import sys
import csv
import yaml


def merge_coverage(file_list):
    """
    The merge_coverage function takes a list of coverage files and merges them into one file.
    The first argument is the list of files, the second argument is the name for the output file.

    :param file_list: Store the list of files that will be merged
    :return: A list of lists with the coverage information for each sample
    :doc-author: Trelent
    """

    with open(file_list[0], "r") as f:
        length = f.readlines()[3].split()[1]

    species = "SARS_CoV_2"
    super_header = ["Name", species, "", ""]
    super_header_1 = ["Length", length, "", ""]
    with open("config_user/parameters.yaml") as file:
        software_parameters = yaml.load(file, Loader=yaml.FullLoader)
    coverage_value = software_parameters["min_coverage_consensus"]
    header = [
        "SAMPLES",
        "Mean depth of coverage",
        f"% of size covered by at least 1-fold",
        f"% of size covered by at least {coverage_value}-fold",
    ]
    # colocar os segmentos para fazer sentido com a flu
    final_output = [super_header, super_header_1, header]

    for file in file_list:
        with open(file, "r") as f:
            x = re.findall("(?<=/main_result/)(.*?)(?=_coverage.csv)", file)
            print("This is the re result: ", x)
            info = f.readlines()[6].split()
            info[0] = x[0]
            final_output.append(info)

    with open(sys.argv[2], mode="w") as f:
        f_writer = csv.writer(
            f, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
        )
        f_writer.writerows(final_output)
        f.close()


if __name__ == "__main__":
    file_list = sys.argv[1].split()
    merge_coverage(file_list)
