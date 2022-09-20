    
import csv
import sys

def get_info_dic(info_list):
    info_dic = {}
    for item in info_list:
        pieces = item.split('=')
        info_dic[pieces[0]] = pieces[1]
    return info_dic

def get_row_dict(row):
    row_dic = {
        'CHROM' : row[0],
        'POS' : row[1],
        'ID' : row[2],
        'REF' : row[3],
        'ALT' : row[4],
        'QUAL' : row[5]
    }
    return row_dic


files_path = sys.argv[1].split(" ")
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

                if 'snp' in info_dic['TYPE'] or 'del' in info_dic['TYPE'] or 'ins' in info_dic['TYPE']:
                    
                    new_entry.append(file[0]) #ID
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
                    print(new_entry)
                    vcf_final.append(new_entry)
            except:
                continue

                
output_file = sys.argv[2]
with open(output_file, mode='w') as f:
    f_writer = csv.writer(f, delimiter=',')
    f_writer.writerows(vcf_final)
    f.close()