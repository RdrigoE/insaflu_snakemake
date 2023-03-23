import csv
import glob
from pathlib import Path
import sys
from typing import NamedTuple, TypedDict


class rule(NamedTuple):
    file_name: str
    rule: str
    path: str


def read_rules(rule_path):
    with open(rule_path) as handler:
        lines = list(csv.reader(handler))

    rule_dict: dict[str, rule] = {}

    for line in lines:
        entry = rule(file_name=line[0], rule=line[1], path=line[2])
        rule_dict[entry.rule] = entry
    return rule_dict


def compare_with_rules(string, rules):
    for _, rule in rules.items():
        if sys.argv[1] + string[1:] in glob.glob(sys.argv[1] + rule.path[1:]):
            print(
                sys.argv[1] + string[2:], glob.glob(sys.argv[1] + rule.path[1:])
            )
            return rule


rule_dict = {}


def crawl_folder(folder_name):
    for item in Path(folder_name).iterdir():
        if item.is_file():
            full_path = str(item.absolute())
            matched_rule = compare_with_rules(path, rules)
            if matched_rule is not None:
                if not rule_dict.get(matched_rule.rule, False):
                    rule_dict[matched_rule.rule] = [path]
                else:
                    rule_dict[matched_rule.rule].append(path)
        else:
            crawl_folder(item.absolute())


class benchmark(NamedTuple):
    seconds: float
    hours: str
    max_rss: float
    max_vms: float
    max_uss: float
    max_pss: float
    io_in: float
    io_out: float
    mean_load: float
    cpu_time: float


def get_time(path) -> benchmark:
    with open(path) as handler:
        line = list(csv.reader(handler, delimiter="\t"))[1]
    return benchmark(
        seconds=float(line[0]),
        hours=line[1],
        max_rss=float(line[2]),
        max_vms=float(line[3]),
        max_uss=float(line[4]),
        max_pss=float(line[5]),
        io_in=float(line[6]),
        io_out=float(line[7]),
        mean_load=float(line[8]),
        cpu_time=float(line[9]),
    )


os.system("cat hello")
rules = read_rules(sys.argv[2])
crawl_folder(sys.argv[1])
rule_seconds = {}

for k, v in rule_dict.items():
    max_value = 0
    for path in v:
        time = get_time(sys.argv[1] + path[1:]).max_rss
        if time > max_value:
            max_value = time
    rule_seconds[k] = max_value
# :%s/[{[a-zA-z]*}]*/*/g

with open(sys.argv[3], "w", encoding="utf-8") as handler:
    write = csv.writer(handler)
    write.writerow(["rule", "seconds"])
    for k_rule, time in rule_seconds.items():
        write.writerow([k_rule, time])
