from pathlib import Path

for item in Path("./").iterdir():
    if "post-deploy" not in item.name:
        with open(item.absolute()) as handler:
            lines = handler.readlines()
            lines = [line for line in lines if item.name[:-5] in line]
            if lines:
                print(lines[0].strip())
