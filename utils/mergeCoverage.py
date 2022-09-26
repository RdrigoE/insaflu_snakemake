import re
import sys
import csv



list = sys.argv[1].split()
with open(list[0],"r") as f:
    length = (f.readlines()[3].split()[1])

species = "SARS_CoV_2"
super_header = ["Name",species,"",""]
super_header_1 = ["Length",length,"",""]
header = ["SAMPLES", "Mean depth of coverage", f"% of size covered by at least 1-fold", f"% of size covered by at least {90}-fold"]
#colocar os segmentos para fazer sentido com a flu
final_output = [super_header,super_header_1, header]

for file in list:
    with open(file, "r") as f:
        x = re.findall("(?<=/main_result/)(.*?)(?=_coverage.csv)",file)
        print("This is the re result: ",x)
        info = f.readlines()[6].split()
        info[0] = x[0]
        final_output.append(info)

with open(sys.argv[2], mode='w') as f:
    f_writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    f_writer.writerows(final_output)
    f.close()