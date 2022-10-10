import re    
import csv
import sys



def get_info_dic(info_list):
    """
    The get_info_dic function takes a list of strings and returns a dictionary with the first item in each string as the key and the second as its value.
    For example, if info_list = ['ID=ADAR0001', 'TYPE=snRNA'], then get_info_dic(info_list) will return {'ID': 'ADAR0001', 'TYPE': 'snRNA'}
    
    :param info_list: Create a dictionary of the information contained in the info column
    :return: A dictionary of the key/value pairs in the info field
    :doc-author: Trelent
    """
    info_dic = {}
    for item in info_list:
        pieces = item.split('=')
        info_dic[pieces[0]] = pieces[1]
    return info_dic

def get_row_dict(row):
    """
    The get_row_dict function takes a row from the VCF file and returns a dictionary with the following keys:
        'CHROM' : row[0],
        'POS' : row[1],
        'ID' : row[2],
        'REF' : row[3],
        'ALT' : row[4],
        'QUAL' : row[5]
    
    :param row: Pass the row of data from the vcf file to the get_row_dict function
    :return: A dictionary of the row elements
    :doc-author: Trelent
    """
    row_dic = {
        'CHROM' : row[0],
        'POS' : row[1],
        'ID' : row[2],
        'REF' : row[3],
        'ALT' : row[4],
        'QUAL' : row[5]
    }
    return row_dic



def create_minor_iSNVs_inc_indels(files_path, output_file):
    """
    The create_minor_iSNVs_inc_indels function creates a csv file containing all the minor iSNVs and indels in the 
        given directory. The function takes two arguments: files_path, which is a list of paths to vcf files, and output_file,
        which is the name of the csv file that will be created by this function. The create_minor_iSNVs_inc_indels function 
        returns nothing.
    
    :param files_path: Pass the path to the files that are used in this function
    :param output_file: Specify the name of the file that will be created
    :return: A csv file with all the variants that are not major and do not have a locus tag
    :doc-author: Trelent
    """
    vcf_final = [['ID',
    'CHROM',
    'POS',
    'TYPE',
    'REF',
    'ALT',
    'FREQ',
    'COVERAGE',
    'EVIDENCE',
    'FTYPE',
    'STRAND',
    'NT_POS',
    'AA_POS',
    'EFFECT',
    'NT CHANGE',
    'AA CHANGE',
    'AA CHANGE ALT',
    'LOCUS_TAG',
    'GENE',
    'PRODUCT',
    'VARIANTS IN INCOMPLETE LOCUS']]
    for file in files_path:
        with open(file, 'r') as invcf:
            for line in invcf:
                try:
                    if line.startswith('#'):
                        continue
                    new_entry = []
                    line = line.strip().split()
                    info = line[7].split(";")
                    info_dic = get_info_dic(info)
                    row_dic = get_row_dict(line)
                    if 'ANN' not in info_dic:
                        info_dic['ANN'] = '|||||||||||||||||||||||||||||||||||||||'
                        info_dic['FTYPE'] = ''
                    else:
                        info_dic['FTYPE'] = 'CDS'
                    if 'snp' in info_dic['TYPE']:
                        if round(float(info_dic['AO'].replace(',','.'))/float(info_dic['DP'].replace(',','.')),2) < 0.50:
                            ann_list = info_dic['ANN'].split('|')
                            new_entry.append(re.findall("(?<=/freebayes/)(.*?)(?=_snpeff.vcf)",file)[0]) #ID
                            new_entry.append(row_dic['CHROM']) #CHROM
                            new_entry.append(row_dic['POS']) #POS
                            new_entry.append(info_dic['TYPE']) #TYPE
                            new_entry.append(row_dic['REF']) #REF
                            new_entry.append(row_dic['ALT']) #ALT
                            new_entry.append(round(float(info_dic['AO'].replace(',','.'))/float(info_dic['DP'].replace(',','.')),2)) #FREQ
                            new_entry.append(0) #COVERAGE
                            new_entry.append(0) #EVIDENCE
                            new_entry.append(info_dic['FTYPE']) #FTYPE
                            new_entry.append(0) #STRAND
                            new_entry.append(ann_list[11]) #NT_POS
                            new_entry.append(ann_list[13]) #AA_POS
                            new_entry.append(ann_list[1]) #EFFECT
                            new_entry.append(ann_list[9]) #NT CHANGE
                            new_entry.append(ann_list[10]) #AA CHANGE
                            new_entry.append(0) #AA CHANGE ALT
                            new_entry.append(0) #LOCUS_TAG
                            new_entry.append(ann_list[3]) #GENE
                            new_entry.append(0) #PRODUCT
                            new_entry.append(0) #VARIANTS IN INCOMPLETE LOCUS
                            vcf_final.append(new_entry)
                            print(new_entry)
                except:
                    continue
                
                    


    with open(output_file, mode='w') as f:
        f_writer = csv.writer(f, delimiter=',')
        f_writer.writerows(vcf_final)
        f.close()

if __name__ == '__main__':
    files_path = sys.argv[1].split(" ")
    output_file = sys.argv[2]
    create_minor_iSNVs_inc_indels(files_path, output_file)