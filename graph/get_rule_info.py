"""get_rule_graph"""
import os
import sys
from time import sleep


script_dir = sys.argv[1]
rules_path = sys.argv[2]
rules_file_path = sys.argv[3]
crawl_folder = sys.argv[4]
time_path = sys.argv[5]
rulegraph = sys.argv[6]
rulegraph_output = sys.argv[7]
# get rules with get_rules_benchmark.py
os.system(
    f"python {script_dir}get_rules_benchmark.py {rules_path} {rules_file_path}"
)
# get rules crawler.py
os.system(f"sed -i '{r's/{[a-zA-Z_]*}/*/g'}' {rules_file_path}")
os.system(f"sed -i 's/benchmark/./g' {rules_file_path}")
os.system(
    f"python {script_dir}crawler.py {crawl_folder} {rules_file_path} {time_path}"
)
# generate rulegraph
# with open(rulegraph, "w") as f:
#     subprocess.Popen(["snakemake", "--rulegraph"], stdout=f, shell=True)
os.system(f"snakemake --rulegraph > {rulegraph}")
# Change the file to feed to color_python.py
sleep(2)
os.system(
    f"python3 {script_dir}color_python.py {time_path} {rulegraph} {rulegraph_output}"
)
