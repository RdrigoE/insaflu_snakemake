import re
import sys
import csv

def get_coverage(filename,n_locus):
    with open(filename, newline='') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter='\t')
        coverage_list = []
        x = re.findall("(?<=coverage/)(.*?)(?=_coverage.csv)",filename)
        for i in csv_reader:
            coverage_list.append(i)
        coverage_list = coverage_list[-1][-n_locus:]
        coverage_list.insert(0,x[0])
        csvfile.close()
    return coverage_list


list = sys.argv[1].split()

n_locus = int(sys.argv[3])
coverage_translate = []
for filename in list:
    print(get_coverage(filename,n_locus))
    coverage_translate.append(get_coverage(filename,n_locus))


with open(sys.argv[2], mode='w') as f:
    f_writer = csv.writer(f, delimiter=',')
    f_writer.writerows(coverage_translate)
    f.close()