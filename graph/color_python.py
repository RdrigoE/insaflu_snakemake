"""Write colors to dag"""

import csv
import os
import sys


def rgb_to_hex(r: int, g: int, b: int) -> str:
    """convert rbg to hex"""
    return f"#{r:02x}{g:02x}{b:02x}"


def read_rules(csv_path):
    with open(csv_path) as handler:
        lines = list(csv.reader(handler))[1:]
    rules = {}
    for line in lines:
        rules[line[0]] = float(line[1])
    return rules


def read_graph(graph_file):
    with open(graph_file) as handler:
        return handler.readlines()


def get_rule_colors(rules):
    list_rules = sorted(rules.items(), key=lambda x: x[1], reverse=True)
    first_too_big = list_rules[0][1] > list_rules[1][1] * 2
    new_dict = {}
    if first_too_big:
        max_value = list_rules[1][1]
        new_dict[list_rules[0][0]] = rgb_to_hex(0, 0, 0)
        new_dict[list_rules[1][0]] = rgb_to_hex(255, 0, 0)
        for rule, value in list_rules[2:]:
            percentage = value * 100 // max_value
            rgb = (0, 0, 0)
            if percentage > 50:
                rgb = (255, 255 - int(((percentage - 50) * 2) * 255 / 100), 0)
            elif percentage < 50:
                rgb = (int(((percentage) * 2) * 255 / 100), 255, 0)
            elif percentage == 50:
                rgb = (255, 255, 0)
            # print(rule, percentage, rgb)
            new_dict[rule] = rgb_to_hex(rgb[0], rgb[1], rgb[2])

    else:
        max_value = list_rules[0][1]
        new_dict[list_rules[0][0]] = rgb_to_hex(255, 0, 0)
        for rule, value in list_rules[1:]:
            percentage = value * 100 // max_value
            rgb = (0, 0, 0)
            if percentage > 50:
                rgb = (255, int(((percentage - 50) * 2) * 255 / 100), 0)
            elif percentage < 50:
                rgb = (255, int(((percentage) * 2) * 255 / 100), 0)
            elif percentage == 50:
                rgb = (255, 255, 0)
            new_dict[rule] = rgb_to_hex(rgb[0], rgb[1], rgb[2])

    return new_dict


def main():
    rules = sys.argv[1]
    graph = sys.argv[2]
    output = sys.argv[3]
    rule_dict = read_rules(rules)
    graph_lines = read_graph(graph)
    rule_colors = get_rule_colors(rule_dict)
    new_graph_lines = []
    for line in graph_lines:
        curr_rule = ""
        found_rule = False
        for rule in rule_dict.keys():
            if '"' + rule + '"' in line:
                found_rule = True
                curr_rule = rule
                break
        if found_rule:
            color = line[line.index("color") :]
            end = color.index(",")
            try:
                if rule_colors[curr_rule] == "#000000":
                    line = f'{line[: line.index("color")]} color = "{rule_colors[curr_rule]}" fontcolor="#ffffff" {line[line.index("color") + end :]}'
                else:
                    line = (
                        line[: line.index("color")]
                        + f'color = "{rule_colors[curr_rule]}"'
                        + line[line.index("color") + end :]
                    )
            except KeyError:
                line = (
                    line[: line.index("color")]
                    + f'color = "{rgb_to_hex(0,0,255)}"'
                    + line[line.index("color") + end :]
                )
            found_rule = False
        elif "color" in line and "edge" not in line and "graph" not in line:
            color = line[line.index("color") :]
            end = color.index(",")
            line = (
                line[: line.index("color")]
                + f'color = "{rgb_to_hex(0,0,255)}"'
                + line[line.index("color") + end :]
            )
        new_graph_lines.append(line)
    with open(output, "w", encoding="utf-8") as handler:
        for line in new_graph_lines:
            handler.write(line)
    os.system(f"sed -i 's/rounded/'rounded,filled'/g' {output}")
    os.system(f"sed -i 's/=rounded,filled/=rounded/g' {output}")
    os.system(f"cat {output} | dot -Tpdf > {output[:-3]}pdf")
    os.system(f"cat {output} | dot -Tpng > {output[:-3]}png")


if __name__ == "__main__":
    main()
    # print(rgb_to_hex(255, 10, 1))
