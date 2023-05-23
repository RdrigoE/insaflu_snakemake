# How to do dag by color?
python ./graph/get_rule_info.py graph/ workflow/rules rules.csv results/benchmark/ time.csv rule_graph rule_graph.png

# Getting the information from benchmarks

Each rule has a benchmark file that saves the following information:

- s	
- h:m:s
- max_rss
- max_vms
- max_uss
- max_pss
- io_in
- io_out
- mean_load
- cpu_time

The intersting part is the second that the rule take from being called to being executed.

## Get every benchmark file pattern

The frist step is to go to every rule file (workflow/rules/*.smk) and search for the terms "rule " or "checkpoint ". In this stage we will save the file where we are, the name of the rule and the name of the file where the benchamark will go.

Then we will save to a file and apply some changes.

First there is the need to replace "benchmark" with "." : $s/benchmark/./g
Second replace all entries with {[a-zA-z]} with "*" : %s/[{[a-zA-z]*}]*/*/g

Then it is ready to feed to the crawler.py script.

## Crawling benchmark folder
1. Load in the rules and their patterns.
2. Compare each file inside benchmark folder to the pattern of each rule.
3. Save the file name in a dictionary where the key is the rule name and the values are a list of files.
4. Go thru the latter dictionary and for each rule check each file and save the value of the property in interest.
5. Find the max value of the list and save in a file with rule name and the value.

## Converting values into a range from red to green

1. Convert values from numeric to percentage
2. From percentage to rgb
3. From rbg to hex

## Get rulegraph from snakemake
snakemake --rulegraph > rulegraph.txt

## Replace default colors with new colors

python3 color_python.py results/benchmark/time.csv rulegraph.txt new_rulegraph.txt

## Generate png

cat new_rulegraph.txt | dot -Tpng > rulegraph_new.png
