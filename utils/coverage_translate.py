import re
import sys
import csv
from get_locus import get_locus

def get_coverage(filename,n_locus):
    """
    The get_coverage function takes a filename as an argument and returns the coverage of each locus in that file.
    The function also takes n_locus as an argument, which is the number of loci in the file.
    
    :param filename: Specify the file to be read
    :param n_locus: Specify the number of loci in the file
    :return: A list of the coverage for each locus in a given sample
    :doc-author: Trelent
    """
    with open(filename, newline='') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter='\t')
        coverage_list = []
        x = re.findall("(?<=main_result/)(.*?)(?=_coverage.csv)",filename)
        for i in csv_reader:
            coverage_list.append(i)
        coverage_list = coverage_list[-1][-n_locus:]
        coverage_list.insert(0,x[0])
        csvfile.close()
    return coverage_list

def create_machine_readable_coverage_file(coverage_files,output,n_locus):
    """
    The create_machine_readable_coverage_file function takes a list of coverage files and outputs a machine readable
    coverage file. The output is the name of the file you want to create. The last argument is n_locus which specifies 
    how many loci are in each coverage file.
    
    :param coverage_files: Specify the files that contain the coverage data
    :param output: Specify the name of the output file
    :param n_locus: Indicate the number of locus in the genome
    :return: A list of lists, where each sublist contains the coverage each locus
    :doc-author: Trelent
    """
    list = coverage_files.split()

    n_locus = int(n_locus)
    coverage_translate = []
    for filename in list:
        print(get_coverage(filename,n_locus))
        coverage_translate.append(get_coverage(filename,n_locus))


    with open(output, mode='w') as f:
        f_writer = csv.writer(f, delimiter=',')
        f_writer.writerows(coverage_translate)
        f.close()

if __name__ == '__main__':
    coverage_files = sys.argv[1]
    output = sys.argv[2]
    reference_gb = sys.argv[3]
    n_locus = len(get_locus(reference_gb))
    create_machine_readable_coverage_file(coverage_files,output,n_locus )
