import pandas as pd
import datetime
from getpass import getpass
import http.cookiejar
import json
import requests
import sys
import os

#from library.scripts.log_entry import log_entry
from log_entry import log_entry
#from library.scripts.load_exp_info import load_exp_info
from load_exp_info import load_exp_info
#from library.scripts.load_MATseq_config import load_MATseq_config
from load_config import load_config
#from library.scripts.gsport_tools import gs_login
from gsport_tools import gs_login
#from library.scripts.gsport_tools import download_filenames
from gsport_tools import download_filenames
#from library.scripts.gsport_tools import download_fastq
from gsport_tools import download_fastq

log_filename=log_filename=datetime.datetime.now().strftime("%Y%m%d_%H%M%S")+'.log'

exp_info,sample_dict=load_exp_info(log_filename)
MATseq_settings=load_config(log_filename)

if exp_info['usr_email']=='':
    log_entry('No GS user account found on exp properties.',True,log_filename)

if exp_info['usr_psw']=='':
    log_entry('No GS psw found on exp properties.',True,log_filename)

if exp_info['GS_project_number']=='':
    log_entry('No GS project found on exp properties.',True,log_filename)
    sys.exit(). sys.exit()
    
    
#check if session is open or cookies are present
try:
    cookies = http.cookiejar.MozillaCookieJar('library//gs_cookies.txt')
    cookies.load()
    if json.loads(requests.get('https://portal.genomescan.nl/' + '/logged_in_api/',
                        cookies=cookies).text)['logged_in']:
        logged_in=True
        log_entry("Already logged in!",True,log_filename)
    else:
        cookies,logged_in=gs_login(exp_info['usr_email'],
                               exp_info['usr_psw'],
                               exp_info['GS_project_number'],
                               log_filename)

#if no cookies logs in       
except FileNotFoundError:
    log_entry("[session] No cookies found. Logging in...",True,log_filename)
    cookies,logged_in=gs_login(exp_info['usr_email'],
                               exp_info['usr_psw'],
                               exp_info['GS_project_number'],
                               log_filename)


#Downloads the list of project files in genomescan servers
datafiles=download_filenames(exp_info['GS_project_number'],cookies,log_filename)
log_entry('File list downloaded from Genomescan servers',True,log_filename)

# checks if a folder for downloaded files is present:
if not 'raw_data' in os.listdir('experiment/'):
    os.mkdir('experiment/raw_data')
if not 'fastq' in os.listdir('experiment/raw_data'):
    os.mkdir('experiment/raw_data/fastq')

# checks if a folder for downloaded files is present:
if not 'raw_data' in os.listdir('experiment/'):
    os.mkdir('experiment/raw_data')
if not 'fastq' in os.listdir('experiment/raw_data'):
    os.mkdir('experiment/raw_data/fastq')

#checks if there are any files already downloaded
previous_downloaded_files=os.listdir('experiment/raw_data/fastq')


#remove old files from new filelist
count=0
temp=[]
for filename in datafiles:
    if not filename['name'] in previous_downloaded_files:
        temp.append(filename)
    bashcount=count+1
datafiles=temp

log_entry(str(len(datafiles))+' new files. '+str(count)+' files already downloaded',True,log_filename)


count=1        
for file in datafiles:
    download_fastq(exp_info['GS_project_number'],cookies,file,log_filename)
    log_entry(str(count)+' files downloaded, '+ str(len(datafiles)-count)+' to go!',True,log_filename)
    count=count+1







