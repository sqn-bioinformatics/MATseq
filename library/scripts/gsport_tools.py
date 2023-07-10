#from library.scripts.log_entry import log_entry
from log_entry import log_entry
import requests
import http.cookiejar
import re
from getpass import getpass
import json
import time


def gs_login(usr,psw,project,log_filename):    

    host='https://portal.genomescan.nl/'
    GSPORT_VERSION="1.6.2"



    log_entry("[GS login] Opening session...",True,log_filename)
    session = requests.Session()
    session.cookies = http.cookiejar.MozillaCookieJar('library/gs_cookies.txt')
    log_entry("[GS login] Get login page",True,log_filename)
    response = session.get('https://portal.genomescan.nl/' + "/login/")
    csrftoken = response.cookies['csrftoken']

    username = usr
    first_try = True


    while re.search('name="password"', response.text) is not None or first_try:
        if not first_try:
            log_entry("[GS login] Invalid credentials",True,log_filename)
            username=input("GS username: ")
            psw=getpass("Password: ")
        
        first_try = False
        login_data = dict(username=username,
                          #password=getpass("Password: "),
                          password=psw,
                          csrfmiddlewaretoken=csrftoken,
                          next='/')
        response = session.post(host + "/login/", data=login_data,
                                headers=dict(Referer=host + "/login/"))

    csrftoken = re.search('name="csrfmiddlewaretoken" value="(.+)"', response.text).group(1)

    first_try = True
    while re.search('name="csrfmiddlewaretoken" value="(.+)"', response.text) is not None or first_try:
        if not first_try:
            log_entry("[GS login] Invalid token",True,log_filename)
        first_try = False
        login_data = dict(token=input("Token: "), username=username, csrfmiddlewaretoken=csrftoken, next='/')
        response = session.post(host + "/otp_ok/", data=login_data,
                                headers={"Referer": host + "/login/",
                                         "User-Agent": "gsport " + GSPORT_VERSION
                                         })

    log_entry("[GS login] Success, saving cookies...",True,log_filename)
    session.cookies.save(ignore_discard=True)

    log_entry("[GS login] Done.",True,log_filename)
    #cookies = session.cookies
    #logged_in = True
    return(session.cookies, True)

def download_filenames(project,cookies, log_filename):
    # connects to genomescan and retrieves file names for the project 
    host='https://portal.genomescan.nl/'

    response = requests.get(host+'data_api2/'+str(project)+'/n',
                            cookies=cookies,
                            params={"cd": '.'})
    try:
        datafiles = json.loads(response.text)

    except json.decoder.JSONDecodeError:
        log_entry("[GS get filelist] Error reading response:"+response.text,True,log_filename)
        exit(1)
    return(datafiles)

 
def sizeofmetric_fmt(num, suffix='B'):
    for unit in ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z']:
        if abs(num) < 1000.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1000.0
    return "%.1f %s%s" % (num, 'Y', suffix)

def human_readable_eta(seconds):
    days = seconds // 86400
    hours = seconds // 3600 % 24
    minutes = seconds // 60 % 60
    seconds = seconds % 60
    ret = str(round(days))+'d' if days > 0 else ''
    ret += str(round(hours))+'h' if hours > 0 else ''
    ret += str(round(minutes))+'m' if minutes > 0 else ''
    ret += str(round(seconds))+'s' if seconds > 0 and minutes < 1 else ''
    return (ret)
 
def download_fastq(project,cookies,file,log_filename):
    project=str(project)
    host='https://portal.genomescan.nl/'
    fsize = 0
    fname = ''
    fsize = file['size']
    if fsize == 0:
        fsize = 1
    fname = file['name']

    response = requests.get(host + '/gen_session_file/', cookies=cookies,
                            params={"project": project,
                                    "filename": "/" + '.' + "/" +
                                    fname
                                    })
    url = host + '/session_files2/' + project + "/" + response.text


    try:
        dsize = 0
        start = time.time()
        with requests.get(url, stream=True, cookies=cookies) as r:
            save_filename = ('/').join(['experiment/raw_data/fastq/',fname])

            with open(save_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:  # filter out keep-alive new chunks
                        f.write(chunk)
                        dsize += len(chunk)
                        rate = dsize // (time.time() - start)

                        print("\r" + sizeofmetric_fmt(fsize) + " " +
                              str(round(dsize / fsize * 100)) + "% " +
                              str(sizeofmetric_fmt(rate)) + "/sec ",
                              "ETA:", human_readable_eta((fsize - dsize) / rate),
                              end='     ')

    except KeyboardInterrupt:
        log_entry('error downloading file'+file,True,log_filename) 
    
    log_entry(file['name']+' downloaded',True,log_filename) 
