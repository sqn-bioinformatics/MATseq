import time


def human_readable_eta_1(seconds):
    if seconds == 0:
        return
    time_units = [("d", 86400), ("h", 3600), ("m", 60), ("s", 1)]
    result = ""
    for unit, divisor in time_units:
        value = seconds // divisor
        if value > 0 or unit == "s" and not result:
            result += f"{value}{unit}"
            seconds %= divisor
    return result


def human_readable_eta_2(seconds):
    days = seconds // 86400
    hours = seconds // 3600 % 24
    minutes = seconds // 60 % 60
    seconds = seconds % 60
    ret = str(round(days)) + "d" if days > 0 else ""
    ret += str(round(hours)) + "h" if hours > 0 else ""
    ret += str(round(minutes)) + "m" if minutes > 0 else ""
    ret += str(round(seconds)) + "s" if seconds > 0 and minutes < 1 else ""

    return ret


def main():
    time_in_seconds = 10245

    start_time = time.time()
    formatted_time_1 = human_readable_eta_1(time_in_seconds)
    end_time = time.time()
    elapsed_time_1 = end_time - start_time
    print(f"Formatted time: {formatted_time_1} in {elapsed_time_1*10000000:.10f}")

    start_time = time.time()
    formatted_time_2 = human_readable_eta_2(time_in_seconds)
    end_time = time.time()
    elapsed_time_2 = end_time - start_time
    print(f"Formatted size: {formatted_time_2} in {elapsed_time_2*10000000:.10f}")


if __name__ == "__main__":
    main()

    #  !! do this with snamekmake
    # # checks if a folder for downloaded files is present:
    # if not "raw_data" in os.listdir("experiment/"):
    #     os.mkdir("experiment/raw_data")
    # if not "fastq" in os.listdir("experiment/raw_data"):
    #     os.mkdir("experiment/raw_data/fastq")

    # # checks if a folder for downloaded files is present:
    # if not "raw_data" in os.listdir("experiment/"):
    #     os.mkdir("experiment/raw_data")
    # if not "fastq" in os.listdir("experiment/raw_data"):
    #     os.mkdir("experiment/raw_data/fastq")
