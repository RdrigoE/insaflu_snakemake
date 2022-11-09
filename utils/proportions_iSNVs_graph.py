import re
import csv
import sys


def get_info_dic(info_list):
    """
    The get_info_dic function takes a list of strings and returns a dictionary with the first item in each string as the key and the second as its value.
    For example, if info_list = ['ID=gene:Solyc00g005000.2.3', 'Name=Potri.001G000100'] then get_info_dic(info_list) will return {'ID':'gene:Solyc00g005000.2.3', 'Name':'Potri001G000100'}.


    :param info_list: Create a dictionary of the information contained in the info column
    :return: A dictionary of the information in each line
    :doc-author: Trelent
    """
    info_dic = {}
    for item in info_list:
        pieces = item.split("=")
        info_dic[pieces[0]] = pieces[1]
    return info_dic


def get_row_dict(row):
    """
    The get_row_dict function takes a row from the vcf file and returns a dictionary with the following keys:
        CHROM, POS, ID, REF, ALT.


    :param row: Pass the row of data from the vcf file to the function
    :return: A dictionary with the chrom, pos, id, ref and alt values from a row in the vcf file
    :doc-author: Trelent
    """
    row_dic = {
        "CHROM": row[0],
        "POS": row[1],
        "ID": row[2],
        "REF": row[3],
        "ALT": row[4],
        "QUAL": row[5],
    }
    return row_dic


def create_graph_csv(files_path, output_file):
    """
    The create_graph_csv function creates a csv file with the following information:
        - Sample name
        - Number of variants with frequency less than 50% (less_50)
        - Number of variants between 90 and 50% (between_90_50)

    :param files_path: Pass the path of the files that will be used to create a csv file
    :param output_file: Specify the name of the output file
    :return: A list of lists with the number of variants that are less than 50% frequency and between 90% and 50%
    :doc-author: Trelent
    """
    file_final = [["Sample", "Less 50", "Between 90 and 50"]]
    for file in files_path:
        count_less = 0
        count_more_50_less_90 = 0
        with open(file, "r") as invcf:
            for line in invcf:
                try:
                    if line.startswith("#"):
                        continue
                    new_entry = []
                    line = line.strip().split()
                    info = line[7].split(";")
                    info_dic = get_info_dic(info)
                    row_dic = get_row_dict(line)
                    if "ANN" not in info_dic:
                        info_dic["ANN"] = "|||||||||||||||||||||||||||||||||||||||"
                        info_dic["FTYPE"] = ""
                    else:
                        info_dic["FTYPE"] = "CDS"

                    if (
                        "snp" in info_dic["TYPE"]
                        or "del" in info_dic["TYPE"]
                        or "ins" in info_dic["TYPE"]
                    ):
                        freq = round(float(info_dic["AO"]) / float(info_dic["DP"]), 2)
                        if freq < 0.50:
                            count_less += 1
                        elif freq > 0.50 and freq < 0.90:
                            count_more_50_less_90 += 1
                except:
                    continue
        file_final.append(
            [
                re.findall("(?<=/)(.*?)(?=_snpeff.vcf)", file)[0],
                count_less,
                count_more_50_less_90,
            ]
        )  # type: ignore

    with open(output_file, mode="w") as f:
        f_writer = csv.writer(f, delimiter=",")
        f_writer.writerows(file_final)
        f.close()


# def create_graph(output_file, out_img):
#     """
#     The create_graph function creates a bar graph of the number of SNVs at frequency 1-50% (minor iSNVs) and 50-90% for each sample.
#     The function takes two arguments: output_file, out_img. The output file is the csv file that contains the data to be graphed, while out_img is
#     the name of the image that will be created by this function.

#     :param output_file: Specify the path to the file containing the data that will be used for creating a graph
#     :param out_img: Specify the name of the output image
#     :return: The graph that we want to create
#     :doc-author: Trelent
#     """
#     with open(output_file,'r') as file:
#         rows = file.readlines()
#         dic = {
#             'samples' : [],
#             'bellow_fifty' : [],
#             'more_fifty' : []
#         }
#         for row in rows[1:]:
#             row = row[:-1]
#             row = row.split(',')
#             dic['samples'].append(row[0])
#             dic['bellow_fifty'].append(int(row[1]))
#             dic['more_fifty'].append(int(row[2]))
#         print(dic)
#     df = pd.DataFrame({'1-50% (minor iSNVs)':dic['bellow_fifty'],
#                     '50-90%':dic['more_fifty']},index=dic['samples'])

#     x = df.plot.barh(stacked=True,color = ["#ffc4de","#81cdff"], )
#     x.set_xlabel('SNVs at frequency 1-50% (minor iSNVs) and 50-90%')
#     x.set_ylabel('Samples')

#     plt.savefig(out_img, format = 'pdf')

if __name__ == "__main__":
    files_path = sys.argv[1].split(" ")
    output_file = sys.argv[2]
    # out_img = sys.argv[3]

    create_graph_csv(files_path, output_file)
    # create_graph(output_file, out_img)
