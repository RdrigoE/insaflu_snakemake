import yaml

with open("config/config.yaml", 'r') as stream:
    try:
        config=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


def is_single_end():
    return config["single_end"]

def is_sample_table():
    if config["sample_table"]["url"] == "":
        return False
    return True

def is_unique_file():
    if is_sample_table() or config["fastq2"]["path"] != "":
        return True
    return False

print(is_unique_file())
