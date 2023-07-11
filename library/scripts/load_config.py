import yaml
from log_entry import log_entry


def load_config(log_filename):
    with open("library/config.txt", "r") as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            print("Error when loading config.txt.")
    log_entry("MATseq configuration loaded.", True, log_filename)
    return config
