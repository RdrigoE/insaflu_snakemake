from pathlib import Path
import csv
from dataclasses import dataclass
import sys


@dataclass
class Rule:
    owner: str
    name: str
    benchmark: str


def get_rules(filename: str) -> list[Rule]:
    rule_list = []
    with open(filename) as handler:
        a, b, c = False, False, False
        new_object: list[str] = []
        for line in handler.readlines():
            if a is False:
                new_object.append(filename)
                a = True
            if line.startswith("rule "):
                if b is False:
                    new_object.append(line.strip()[5:-1])
                    b = True
            elif line.startswith("checkpoint "):
                if b is False:
                    new_object.append(line.strip()[10:-1].strip())
                    b = True
            if "benchmark/" in line:
                if c is False:
                    new_object.append(line.strip())
                    b = False
                    a = False
                    rule_list.append(
                        Rule(new_object[0], new_object[1], new_object[2])
                    )
                    new_object = []
    return rule_list


def get_quotes(string: str) -> str:
    cut_start = string.index("benchmark")
    first_sliced = string[cut_start:]
    cut_end = first_sliced[::-1].index("vst.")
    return string[cut_start : len(string) - cut_end]


if __name__ == "__main__":
    path = Path(sys.argv[1])
    rule_list: list[Rule] = []
    for dir_item in path.iterdir():
        fresh_rules = get_rules(str(dir_item.absolute()))
        for rule in fresh_rules:
            if rule not in rule_list:
                rule_list.append(rule)

    new_list = list(
        map(lambda x: [x.owner, x.name, get_quotes(x.benchmark)], rule_list)
    )

    with open(sys.argv[2], "w", encoding="utf-8") as handler:
        writer = csv.writer(handler)
        writer.writerows(new_list)
