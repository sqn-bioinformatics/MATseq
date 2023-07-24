# 20220804 -> added line to create a log fold if it doesnt exist yet.
import datetime
import os


def log_entry(text, show, log_filename):
    if "logs" not in os.listdir("./"):
        os.mkdir("logs")

    text = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S: ") + text

    with open("logs/" + log_filename, "a") as file_object:
        file_object.write(text + "\n")
    if show:
        print(text)
