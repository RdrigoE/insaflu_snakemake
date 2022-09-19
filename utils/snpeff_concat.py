import csv
import sys

files_path = sys.argv[1].split(" ")
vcf_final = [['ID',
 "CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","snps"
]]
for file in files_path:
    with open(file, 'r') as invcf:
        for line in invcf:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            line.insert(0,file)
            vcf_final.append(line)

output_file = sys.argv[2]
with open(output_file, mode='w') as f:
    f_writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    f_writer.writerows(vcf_final)
    f.close()