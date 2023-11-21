import http.cookiejar
import json
from typing import List
import requests
import os
import re
import time

from getpass import getpass
from mylogger import get_logger
from parses_experiment_info import get_experiment_info


# /home/t.afanasyeva/mambaforge/envs/genomescan_works/lib/python3.10/site-packages/gsport.py

# Sets path to ~/MATseq
path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

with open(os.path.join(path, "library/config_genomescan.json")) as json_data_file:
    config = json.load(json_data_file)
    HOST = config.get("HOST")
    GSPORT_VERSION = config.get("GSPORT_VERSION")
    USER_EMAIL = config.get("USER_EMAIL")
    USER_PASSWORD = config.get("USER_PASSWORD")
    TOKEN = config.get("TOKEN")


# Gets logger instance
logger = get_logger(__name__)


class UserLoggerIn:
    def __init__(self, session, logger):
        self.logger = logger
        self.session = session

    # loads cookies
    def check_for_cookies(self) -> str:
        try:
            cookies = http.cookiejar.MozillaCookieJar("library/gs_cookies.txt")
            cookies.load()
            response = requests.get(
                "https://portal.genomescan.nl/logged_in_api/", cookies=cookies
            )
            if json.loads(response.text)["logged_in"]:
                self.logger.info("you are already logged in")
            else:
                # Logs in using provoided user information
                self.logger.info("logging in...")
                cookies = self.log_in()

        except FileNotFoundError:
            self.logger.info("logging in...")
            cookies = self.log_in()

            return cookies

    def log_in(self):
        self.logger.info("opening session...")

        login_url = HOST + "/login/"
        otp_url = HOST + "/otp_ok/"

        # Performs the initial login
        response = self.session.get(login_url)

        csrf_token = re.search(
            'name="csrfmiddlewaretoken" value="(.+)"', response.text
        ).group(1)

        # Send login data
        login_data = {
            "username": USER_EMAIL,
            "password": USER_PASSWORD,
            "csrfmiddlewaretoken": csrf_token,
            "next": "/",
        }
        response = self.session.post(
            login_url, data=login_data, headers={"Referer": login_url}
        )
        token = "9FMZB"

        # Checks for 2FA requirement
        if 'name="csrfmiddlewaretoken" value="' in response.text:
            # token = input("Enter 2FA token: ")
            login_data = {
                "token": token,
                "username": USER_EMAIL,
                "csrfmiddlewaretoken": csrf_token,
                "next": "/",
            }

            response = self.session.post(
                otp_url,
                data=login_data,
                headers={
                    "Referer": login_url,
                    "User-Agent": "gsport " + GSPORT_VERSION,
                },
            )

        # Check for a successful login
        if "Welcome" in response.text:
            self.logger.info("username and password are correct")
        else:
            self.logger.error("login failed")

        return self.session.cookies


class FileHandler:
    def __init__(self, project_number, session, cookies, logger):
        self.project_number = project_number
        self.logger = logger
        self.session = session
        self.cookies = cookies

    # Downloads the list of file names
    def download_filenames(self) -> List:
        data_url = HOST + "data_api2/"

        response = self.session.get(
            data_url + str(self.project_number) + "/n",
            cookies=self.cookies,
            params={"cd": "."},
        )

        print(response.text)

        file_names = json.loads(response.text)  # check format of file_nmaes

        # except json.decoder.JSONDecodeError as e:
        #     logger.error(f"error reading response {response.text}: {e}")

        return file_names

    def format_data_size(size_in_bytes: int) -> str:
        """
        This function converts file size in bytes into units of powers of 10, following the International System of Units (SI)
        """
        if size_in_bytes == 0:
            return "0B"
        else:
            # Defines abbreviations of metric units
            units = ["B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB"]

            # Starts with the original size and no metric unit
            formatted_size = size_in_bytes
            unit = 0  # Index for the metric units

            # Chooses the appropriate metric unit
            while formatted_size >= 1000 and unit < len(units) - 1:
                formatted_size /= 1000.0
                unit += 1

            return f"{formatted_size:.1f} {units[unit]}"

    def format_time(seconds: float) -> str:
        if seconds == 0:
            return "0s"
        time_units = [("d", 86400), ("h", 3600), ("m", 60), ("s", 1)]
        result = ""
        for unit, divisor in time_units:
            value = seconds // divisor
            if value > 0 or unit == "s" and not result:
                result += f"{value}{unit}"
                seconds %= divisor
        return result

    def download_files(self):
        file_names = self.download_filenames()

        count = 0
        for name in file_names:
            count += 1
            file = self.download_file(name)
            logger.info(
                f"FileHandler.download_files: {count} files downloaded, {len(file_names) - count} files left to download"
            )

        return file

    def download_file(self, name):
        project = str(project)  # start to work on the code from here
        fsize = 0
        fname = ""
        fsize = file["size"]
        if fsize == 0:
            fsize = 1
        fname = file["name"]

        response = self.session.get(
            HOST + "/gen_session_file/",
            cookies=self.cookies,
            params={"project": project, "filename": "/" + "." + "/" + fname},
        )

        url = HOST + "/session_files2/" + project + "/" + response.text

        try:
            data_size = 0
            start = time.time()
            with self.session.get(url, stream=True, cookies=self.cookies) as r:
                save_file_name = ("/").join(["experiment/raw_data/fastq/", fname])

                with open(save_file_name, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:  # filter out keep-alive new chunks
                            f.write(chunk)
                            data_size += len(chunk)
                            rate = data_size // (time.time() - start)

                            logger.info(
                                f"file: {name}, ",
                                f"size: {self.format_data_size(fsize)}, {round(data_size / fsize * 100)}% of total size, ",
                                f"downloading speed: {self.format_data_size(rate)}/sec, time left to download: {self.format_time((fsize - data_size) / rate)}",
                            )

        except KeyboardInterrupt:
            logger.error(f"error downloading file {name}")

        logger.info(f"{name} downloaded")


def main():
    # Initializes an HTTP session
    session = requests.Session()

    experiment_info, _ = get_experiment_info()
    project_number = experiment_info.get("project_number")

    cuurent_user = UserLoggerIn(session, logger)
    cookies = cuurent_user.check_for_cookies()

    file_handler = FileHandler(project_number, session, cookies, logger)
    file_handler.download_files()
    logger.info("all files downloaded")


if __name__ == "__main__":
    main()
