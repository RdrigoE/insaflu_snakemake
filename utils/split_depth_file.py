import sys

file_path = sys.argv[1]
last_pos = len(file_path) - file_path[::-1].index("/")
print(file_path[last_pos:])
filename = file_path[last_pos:]
new_file = []
with open(file_path, "r") as f:
    for line in f.readlines():
        new_file.append(f"{filename[:-6]}__{line}")
with open(file_path, "w") as f:
    f.writelines(new_file)   
