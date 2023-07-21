#!/usr/bin/env python

'''
creator: dsobral
'''


from Bio import SeqIO
from optparse import OptionParser
import sys
import random
import os
import gzip
import glob


def get_locus_name_len(genbank_file):
    locus_name = []
    locus_len = []
    with open(genbank_file, encoding="utf-8") as handle_gb:
        for record in SeqIO.parse(handle_gb, "genbank"):
            locus_name.append(record.name)
            locus_len.append(str(len(record.seq)))
    return locus_name, locus_len


class Util(object):
    '''
    classdocs
    '''
    FORMAT_FASTA = "fasta"
    FORMAT_FASTQ = "fastq"
    EXTENSION_ZIP = ".gz"
    TEMP_DIRECTORY = "/tmp"
    COVERAGE_TEMP_DIRECTORY = "getCoverage"
    COUNT_DNA_TEMP_DIRECTORY = "count_dna_temp_directory"

    def __init__(self):
        '''
        Constructor
        '''
        pass

    def get_temp_file(self, file_name, sz_type):
        main_path = os.path.join(
            self.TEMP_DIRECTORY, self.COUNT_DNA_TEMP_DIRECTORY)
        if (not os.path.exists(main_path)):
            os.makedirs(main_path)
        while 1:
            return_file = os.path.join(main_path, "count_dna_" + file_name + "_" + str(
                random.randrange(100000, 999999, 10)) + "_file" + sz_type)
            if (not os.path.exists(return_file)):
                return return_file

    def is_integer(self, n_value):
        try:
            int(n_value)
            return True
        except ValueError:
            return False

    def is_gzip(self, file_name): return True if (
        file_name.rfind(".gz") == len(file_name) - 3) else False

    def __get_temp_file__(self, index_file_to_process, sz_type):
        main_path = os.path.join(
            self.TEMP_DIRECTORY, self.COVERAGE_TEMP_DIRECTORY)
        if (not os.path.exists(main_path)):
            os.makedirs(main_path)
        while 1:
            return_file = os.path.join(main_path, "seq_dna_" + str(index_file_to_process) + "_" + str(
                random.randrange(10000, 99999, 10)) + "_file." + sz_type)
            if (not os.path.exists(return_file)):
                return return_file


class CoverageElement(object):
    """
    Only have the number of reads and average
    """

    def __init__(self, element):
        self.element = element
        self.dt_data = {}

    def add_coverage(self, type_coverage, coverage):
        self.dt_data[type_coverage] = coverage

    def get_coverage(self, type_coverage):
        return self.dt_data.get(type_coverage, None)

    def __str__(self):
        return "Element: {}  {}".format(self.element, self.dt_data)


class Coverage(object):
    """
    Only have the number of reads and average
    """
    COVERAGE_ALL = "CoverageAll"
    COVERAGE_MORE_DEFINED_BY_USER = "CoverageMoreDefinedByUser"

    def __init__(self, limit_defined_by_user=10):
        self.limit_defined_by_user = limit_defined_by_user
        self.dt_data = {}
        self.ratio_value_defined_by_user = 0

    def get_dict_data(self): return self.dt_data

    def add_coverage(self, element, type_coverage, coverage):
        if (element in self.dt_data):
            self.dt_data[element].add_coverage(type_coverage, coverage)
        else:
            self.dt_data[element] = CoverageElement(element)
            self.dt_data[element].add_coverage(type_coverage, coverage)

    def get_coverage(self, element, type_coverage):
        if (element in self.dt_data):
            return self.dt_data[element].get_coverage(type_coverage)
        raise Exception("Error: there's no key like this: " + element)

    def __str__(self):
        sz_return = "snps\tLOCUS\t"
        for key in self.dt_data:
            sz_return += "{}\t".format(
                self.get_coverage(
                    key, Coverage.COVERAGE_ALL))
            print(sz_return)

        for key in self.dt_data:
            sz_return += "{}\t".format(
                self.get_coverage(
                    key, Coverage.COVERAGE_MORE_DEFINED_BY_USER))

        return sz_return


class ParseFile(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.data_file = None
        self.reference_dict = {}
        self.vect_reference = []
        self.util = Util()

    def is_gzip(self, file_name): return True if (
        file_name.rfind(".gz") == len(file_name) - 3) else False

    def parse_file(self, file_name):
        """
        """
        self.data_file = DataFile(file_name)
        if (self.util.is_gzip(file_name)):
            handle = gzip.open(file_name, mode='rt')
        else:
            handle = open(file_name)
        for line in handle:
            sz_temp = line.strip().lower()
            if (len(sz_temp) == 0 or sz_temp[0] == '#'):
                continue
            self.data_file.add_data(line)
        handle.close()
        return self.data_file

    def read_reference_fasta(self, reference_file):
        """
        test if the reference_file and ge the handle
        """
        if (not os.path.exists(reference_file)):
            raise Exception(
                "Can't locate the reference file: '" + reference_file + "'")

        # set temp file name
        temp_file_name = reference_file

        # create temp file
        b_temp_file = False
        if self.util.is_gzip(reference_file):
            b_temp_file = True
            temp_file_name = self.util.get_temp_file(
                "reference_file_", ".fasta")
            cmd = "gzip -cd " + reference_file + " > " + temp_file_name
            os.system(cmd)

        for rec in SeqIO.parse(temp_file_name, 'fasta'):
            self.reference_dict[rec.id] = len(str(rec.seq))
            self.vect_reference.append(rec.id)

        # remove temp file if necessary
        if b_temp_file and temp_file_name:
            os.remove(temp_file_name)


class DataFile(object):
    '''
    classdocs
    '''
    util = Util()

    def __init__(self, file_name):
        '''
        Constructor
        '''
        self.file_name = file_name
        self.vect_chromosomes = []
        self.dict_data = {}
        self.dict_data_coverage = {}
        self.previous_position = -1

    def get_vect_chromosomes(self): return self.vect_chromosomes
    def get_dict_data(self): return self.dict_data

    def add_data(self, line):
        if (len(line) == 0 or line[0] == '#'):
            return
        vect_data = line.split()
        if (len(vect_data) != 3):
            raise Exception("File: " + self.file_name +
                            "\nThis line must have three values '" + line + "'")
        if (not self.util.is_integer(vect_data[1])):
            raise Exception("File: " + self.file_name + "\nLine: '" +
                            line + "'\nThe locus need to be integer")
        if (not self.util.is_integer(vect_data[2])):
            raise Exception("File: " + self.file_name + "\nLine: '" +
                            line + "'\nThe coverage need to be integer")
        if (vect_data[0] in self.dict_data):
            if (int(vect_data[1]) <= (self.previous_position)):
                raise Exception("File: " + self.file_name + "\nLine: '" + line +
                                "'\nThe locus need to be greater than the predecessor in the file")
            self.dict_data[vect_data[0]].append([vect_data[1], vect_data[2]])
            self.previous_position = int(vect_data[1])
        else:
            self.vect_chromosomes.append(vect_data[0])
            self.dict_data[vect_data[0]] = [[vect_data[1], vect_data[2]]]
            self.previous_position = int(vect_data[1])

    def get_coverage(self, sz_chromosome, length_chromosome):
        if (sz_chromosome not in self.dict_data):
            return 0
        if (sz_chromosome in self.dict_data_coverage):
            return self.dict_data_coverage[sz_chromosome]
        if (length_chromosome == 0):
            return 0
        # medaka sometimes creates bigger references than the original, difference 2 or 3 number of bases
        if (len(self.dict_data[sz_chromosome]) > (length_chromosome * 1.10) or
                len(self.dict_data[sz_chromosome]) < (length_chromosome - (length_chromosome * 0.10))):
            raise Exception("Chromosome '%s' has different sizes. Coverage: %d; Reference: %d" % (
                sz_chromosome, len(self.dict_data[sz_chromosome]), length_chromosome))
        sum_total = 0
        for data_ in self.dict_data[sz_chromosome]:
            sum_total += int(data_[1])
        self.dict_data_coverage[sz_chromosome] = sum_total / \
            float(length_chromosome)
        return self.dict_data_coverage[sz_chromosome]

    def get_ratio_more_than(self, sz_chromosome, length_chromosome, value):
        if (sz_chromosome not in self.dict_data):
            return 0
        if (length_chromosome == 0):
            return 0
        # medaka sometimes creates bigger references than the original, difference 2 or 3 number of bases
        if (len(self.dict_data[sz_chromosome]) > (length_chromosome * 1.10) or
                len(self.dict_data[sz_chromosome]) < (length_chromosome - (length_chromosome * 0.10))):
            raise Exception("Chromosome '%s' has different sizes. Coverage: %d; Reference: %d" % (
                sz_chromosome, len(self.dict_data[sz_chromosome]), length_chromosome))
        sum_total = 0
        for data_ in self.dict_data[sz_chromosome]:
            sum_total += (1 if (int(data_[1]) > value) else 0)
        return sum_total / float(length_chromosome)


class GetCoverage(object):
    """
    get coverage from deep.gz file
    need deep.gz file and reference
    """
    util = Util()

    def __init__(self):
        self.vect_files_processed = []
        self.reference_dict = {}
        self.vect_reference = []

    def get_dict_reference(self): return self.reference_dict
    def get_vect_reference(self): return self.vect_reference
    def get_vect_files_processed(self): return self.vect_files_processed

    def process_input_files(self, input_path):
        """
                test if input_file is directory or file and read all files in the directory
        """

        print("Collecting files in: " + input_path)
        if (os.path.isfile(input_path)):
            if (os.path.exists(input_path)):
                self.vect_files_processed.append(input_path)
        else:
            self.vect_files_processed = glob.glob(input_path)

    def get_dict_with_coverage(self, deep_file):
        """
        get a dictonary of elements with coverage
        """
        self.reference_dict = {}
        self.vect_reference = []

        parse_file = ParseFile()
        data_file = parse_file.parse_file(deep_file)
        dt_out = {}
        for key in data_file.get_dict_data().keys():
            dt_out[key] = [int(value_[1])
                           for value_ in data_file.get_dict_data()[key]]
        return dt_out

    # def get_coverage(self, deep_file, reference, limit_defined_by_user = None, limit_defined_to_project = None):
    def get_coverage(self, deep_file, reference, limit_defined_by_user=10):
        """
        get an instance of coverage
        """
        self.reference_dict = {}
        self.vect_reference = []

        parse_file = ParseFile()
        data_file = parse_file.parse_file(deep_file)
        self.read_reference_fasta(reference)

        coverage = Coverage()
        for chromosome in self.vect_reference:
            if (chromosome not in self.reference_dict):
                raise Exception("Can't locate the chromosome '" +
                                chromosome + "' in reference file")
            coverage.add_coverage(chromosome, Coverage.COVERAGE_ALL, "%.1f" % (
                data_file.get_coverage(chromosome, self.reference_dict[chromosome])))
            coverage.add_coverage(chromosome, Coverage.COVERAGE_MORE_DEFINED_BY_USER,
                                  "%.1f" % (data_file.get_ratio_more_than(chromosome, self.reference_dict[chromosome], limit_defined_by_user - 1) * 100))
        return coverage

    def read_reference_fasta(self, reference_file):
        """
        test if the reference_file and ge the handle
        """
        if (not os.path.exists(reference_file)):
            raise Exception(
                "Can't locate the reference file: '" + reference_file + "'")

        # set temp file name
        temp_file_name = reference_file

        # create temp file
        b_temp_file = False
        if self.util.is_gzip(reference_file):
            b_temp_file = True
            temp_file_name = self.util.get_temp_file(
                "reference_file_", ".fasta")
            cmd = "gzip -cd " + reference_file + " > " + temp_file_name
            os.system(cmd)

        for rec in SeqIO.parse(temp_file_name, 'fasta'):
            self.reference_dict[rec.id] = len(str(rec.seq))
            self.vect_reference.append(rec.id)

        ###
        if b_temp_file and temp_file_name:
            os.remove(temp_file_name)


if __name__ == '__main__':

    """
    V0.6 release 21/11/2017
            add - 	ratio as a parameter
    V0.5 release 30/09/2017
            Fix - 	length chromosome
    V0.4 release 30/09/2017
            Fix - 	when coverage doesn't have at all coverage doesn't appear in the results
    V0.3 release 30/09/2017
            Fix - 	get coverage 0 for chromosomes that are missing in coverage file
    V0.2 release 30/09/2017
            Add - 	need the reference file
                            check if all chromosomes that are in the coverage are the ones that are in the reference.
                            get the size of the chromosomes from the reference
                            add chromosome length to the report
    V0.1 release 12/09/2017
            Add - coverage average base on the files <chromosome> <position> <deep coverage>:
    """

    b_debug = False
    if (b_debug):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        input_file = os.path.join(dir_path, "test/files/*.depth")
        output_file = "tmp_out_get_coverage.csv"
        reference = os.path.join(dir_path, "test/files/ref/ref_H3.fasta")
        threshold = 10
        get_coverage = GetCoverage()
        get_coverage.process_input_files(input_file)
        handle = open(output_file, "w")
        for file_to_process in get_coverage.get_vect_files_processed():
            coverage = get_coverage.get_coverage(deep_file=file_to_process,
                                                 reference=reference, limit_defined_by_user=threshold)
            print(coverage)
            handle.write("{},{}\n".format(

                str(os.path.basename(file_to_process)).split(".")[0], coverage))
            # print("{},{}\n".format(str(os.path.basename(file_to_process)).split(".")[0],coverage))
        handle.close()
        sys.exit(1)
    else:
        parser = OptionParser(
            usage="%prog [-h] -i -r -o [-t]", version="%prog 0.6", add_help_option=False)
        parser.add_option("-i", "--input", type="string", dest="input",
                          help="Input file or path with coverage files. Can be zipped.", metavar="IN_FILE")
        parser.add_option("-r", "--reference", type="string", dest="reference",
                          help="Reference file to get the length of the chromosomes to check his name.", metavar="REF_FILE")
        parser.add_option("-g", "--genbank", type="string", dest="genbank",
                          help="Genbank file")

        parser.add_option("-o", "--output", type="string",
                          dest="output", help="Output file name", metavar="OUT_FILE")
        parser.add_option("-t", "--threshold", type="int", dest="threshold",
                          help="Minimum coverage per position (default 10)", metavar="OUT_FILE")
        parser.add_option('-h', '--help', dest='help',
                          action='store_true', help='show this help message and exit')

        (options, args) = parser.parse_args()

        if (options.help):
            parser.print_help()
            print("")
            print("\tCreate an output file with several averages about the coverage.")
            print("\tOnly runs in linux or mac.")
            print(
                "\texample: python getCoverage -i '/usr/local/zpto/*.gz' -r reference.fasta -o resultsOut.csv")
            print("\texample: python getCoverage -i '/usr/local/zpto/*.gz' -r reference.fasta -o resultsOut.csv -t 30")
            print("")
            print(
                "\tThe input coverage files must be in this format '<chromosome> <position> <deep coverage>'")
            sys.exit(0)

        if (len(args) != 0):
            parser.error(
                "incorrect number of arguments, no of arguments: " + str(len(args)))

        if not options.input:   #
            parser.error('Name of the input files/path is not specified')

        if not options.genbank:   #
            parser.error('Name of the input files/path is not specified')

        if not options.reference:   #
            parser.error('Name of the reference file not specified')

        if not options.output:   #
            parser.error('Output file is not specified')

        threshold = 10
        if (options.threshold):
            threshold = options.threshold
        locus_name, locus_len = get_locus_name_len(options.genbank)
        get_coverage = GetCoverage()
        get_coverage.process_input_files(options.input)
        handle = open(options.output, "w")
        # print(get_coverage.get_vect_files_processed())
        for file_to_process in get_coverage.get_vect_files_processed():
            coverage = get_coverage.get_coverage(deep_file=file_to_process,
                                                 reference=options.reference,
                                                 limit_defined_by_user=threshold)
            print(file_to_process)
            print(str(os.path.basename(file_to_process)).split(".")[0])
            handle.write("\n")
            handle.write("Chromosome\n")
            print(locus_name, locus_len)
            name = "Name\t" + "\t".join(locus_name)
            handle.write(name.strip() + "\n")
            lenght = "Lenght\t" + "\t".join(locus_len)
            handle.write(lenght.strip() + "\n")
            handle.write("\n")
            handle.write("Coverage\tRatio>0\tRatio>9\n")
            handle.write("{}\n".format(
                str(coverage).strip()))

            # print("{},{}".format(str(os.path.basename(file_to_process)).split(".")[0],coverage))
            handle.close()
