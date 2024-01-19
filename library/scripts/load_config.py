import yaml
import os
from mylogger import get_logger

logger = get_logger(__name__)


def load_config():
    # Sets path to ~/MATseq
    path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    file_path = os.path.join(path, "library", "config.txt")

    try:
        with open(file_path, "r") as stream:
            try:
                config = yaml.safe_load(stream)
                logger.info("MATseq configuration loaded.")
                return config
            except yaml.YAMLError as e:
                logger.error(f"Error loading config.txt from {file_path}: {e}")
    except FileNotFoundError:
        logger.error(f"config.txt not found in {file_path}")
    except Exception as e:
        logger.error(f"An error occurred while loading config.txt: {e}")
    return None


if __name__ == "__main__":
    config = load_config()
    if config:
        print(config)
    else:
        print("Failed to load MATseq configuration.")
