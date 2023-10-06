import logging
import os
from datetime import datetime


def get_logger(name):
    # Sets path to ~/MATseq
    path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if not os.path.isdir("logs"):
        os.mkdir("logs")

    LOG_FILENAME = datetime.now().strftime("%Y%m%d_%H%M%S") + ".log"
    log_file_path = os.path.join(path, "logs", LOG_FILENAME)

    file_formatter = logging.Formatter(
        "%(asctime)s~%(levelname)s~%(message)s~module:%(module)s~function:%(module)s"
    )
    file_info_formatter = logging.Formatter("%(asctime)s~%(levelname)s~%(message)s")
    console_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(file_formatter)

    file_handler_info = logging.FileHandler(log_file_path)
    file_handler_info.setLevel(logging.INFO)
    file_handler_info.setFormatter(file_info_formatter)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(console_formatter)

    logger = logging.getLogger(name)
    logger.addHandler(file_handler)
    logger.addHandler(file_handler_info)
    logger.addHandler(console_handler)
    logger.setLevel(logging.DEBUG)

    return logger
