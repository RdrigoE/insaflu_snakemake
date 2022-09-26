import re 
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
file_final = [['Sample','Less 50','Between 90 and 50']]



for file in files_path:
    count_less = 0
    count_more_50_less_90 = 0
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
                    freq = round(float(info_dic['AO'])/float(info_dic['DP']),2)
                    if freq < 0.50:
                        count_less+=1
                    elif freq > 0.50 and freq < 0.90:
                        count_more_50_less_90 += 1
            except:
                continue
    file_final.append([re.findall("(?<=/)(.*?)(?=_snpeff.vcf)",file)[0],count_less,count_more_50_less_90])
    

                
output_file = sys.argv[2]
with open(output_file, mode='w') as f:
    f_writer = csv.writer(f, delimiter=',')
    f_writer.writerows(file_final)
    f.close()

import pandas as pd
import matplotlib.pyplot as plt

with open(output_file,'r') as file:
    rows = file.readlines()
    dic = {
        'samples' : [],
        'bellow_fifty' : [],
        'more_fifty' : []
    }
    for row in rows[1:]:
        row = row[:-1]
        row = row.split(',')
        dic['samples'].append(row[0])
        dic['bellow_fifty'].append(int(row[1]))
        dic['more_fifty'].append(int(row[2]))
    print(dic)
df = pd.DataFrame({'1-50% (minor iSNVs)':dic['bellow_fifty'],
                   '50-90%':dic['more_fifty']},index=dic['samples'])

x = df.plot.barh(stacked=True,color = ["#ffc4de","#81cdff"], )
x.set_xlabel('SNVs at frequency 1-50% (minor iSNVs) and 50-90%')
x.set_ylabel('Samples')

out_img = sys.argv[3]
plt.savefig(out_img, format = 'pdf')