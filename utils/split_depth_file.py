import sys

file_path = sys.argv[1]


last_pos = len(file_path) - file_path[::-1].index("/")
filename = file_path[last_pos:]
new_file = []

with open(file_path, "r") as f:
    for line in f.readlines():
        if line.split()[0] == filename[-7]:
            new_file.append(line)
with open(file_path, "w") as f:
    f.writelines(new_file)   
