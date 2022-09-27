import re
import csv
import sys

def get_info_dic(info_list):
    """
    The get_info_dic function takes a list of strings and returns a dictionary
    with the first item in each string as the key and the second as its value.
    
    
    :param info_list: Store the information of each line in the vcf file
    :return: A dictionary where the keys are the first item in each list, and the values are
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
        'ALT' : row[4]
    
    :param row: Store the current row in the iteration
    :return: A dictionary with the keys and values for each row
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

def validated_minor_iSNVs(files_path):
    """
    The validated_minor_iSNVs function takes a list of files and creates a csv file with the following columns:
        ID, CHROM, POS, TYPE, REF, ALT, FREQ (frequency), COVERAGE (coverage), EVIDENCE (evidence level), FTYPE
        STRAND(strand bias) NT_POS(nucleotide position in codon) AA_POS(amino acid position in protein sequence) 
        EFFECT(effect type of variant on gene product), NT CHANGE (nucleotide change from reference genome to sample genome) 
        AA CHANGE (amino acid change from reference
    
    :param files_path: Pass the path to the folder containing all of the vcf files
    :return: A list of lists
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
                        if round(float(info_dic['AO'])/float(info_dic['DP']),2) < 0.50:
                            new_entry.append(re.findall("(?<=/main_result/snpeff/)(.*?)(?=_snpeff.vcf)",file)[0]) #ID
                            new_entry.append(row_dic['CHROM']) #CHROM
                            new_entry.append(row_dic['POS']) #POS
                            new_entry.append(info_dic['TYPE']) #TYPE
                            new_entry.append(row_dic['REF']) #REF
                            new_entry.append(row_dic['ALT']) #ALT
                            new_entry.append(round(float(info_dic['AO'])/float(info_dic['DP']),2)) #FREQ
                            new_entry.append(0) #COVERAGE
                            new_entry.append(0) #EVIDENCE
                            new_entry.append(info_dic['FTYPE']) #FTYPE
                            new_entry.append(0) #STRAND
                            new_entry.append(info_dic['ANN'].split('|')[11]) #NT_POS
                            new_entry.append(info_dic['ANN'].split('|')[13]) #AA_POS
                            new_entry.append(info_dic['ANN'].split('|')[1]) #EFFECT
                            new_entry.append(info_dic['ANN'].split('|')[9]) #NT CHANGE
                            new_entry.append(info_dic['ANN'].split('|')[10]) #AA CHANGE
                            new_entry.append(0) #AA CHANGE ALT
                            new_entry.append(0) #LOCUS_TAG
                            new_entry.append(info_dic['ANN'].split('|')[3]) #GENE
                            new_entry.append(0) #PRODUCT
                            new_entry.append(0) #VARIANTS IN INCOMPLETE LOCUS
                            vcf_final.append(new_entry)
                except:
                    continue

                    
    with open(output_file, mode='w') as f:
        f_writer = csv.writer(f, delimiter=',')
        f_writer.writerows(vcf_final)
        f.close()


if __name__ == '__main__':
    files_path = sys.argv[1].split(" ")
    output_file = sys.argv[2]
    validated_minor_iSNVs(files_path)