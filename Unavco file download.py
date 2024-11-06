# -*- coding: utf-8 -*-
"""
Created on Tue May 14 10:33:02 2024

@author: Samy

This script downloads the desired RINEX files from the Unavco servers on 
the desired days. These files are downloaded as hatanaka zip files which
can be uncompressed through the bash terminal by going to the directory where
the compressed unavco files are stored and using the command 
<< uncompress *.Z >> without the <<>> symbols.
"""

import requests
from pathlib import Path
from earthscope_sdk.auth.device_code_flow import DeviceCodeFlowSimple
from earthscope_sdk.auth.auth_flow import NoTokensError


def get_unavco_name(year=2022,day=80,station='drao'):
    year=repr(year)
    assert(len(year)==4)
    day=f'{day:03d}'
    dir='data.unavco.org/archive/gnss/rinex/obs/'+year+'/'+day+'/'
    fname=station+day+'0.'+year[-2:]+'d.Z'
    url='https://'+dir+fname
    return url,fname

# Receiver station
station='drao'
year=2024
# Specify which days of the year you want to download
days=[279, 280, 281] 


# choose where you want the token saved - the default file name is sso_tokens.json
# if you want to keep the default name, set the path to a directory. 
# Include a file name to rename. 

for day in days:
    url,fname = get_unavco_name(year,day,station)

    # token_path should specify the same location the current script is
    # being run in, where the file name is token_dir
    token_path = 'C:/Users/Samy/Documents/GPS research (Sievers)/token_dir'

    # instantiate the device code flow subclass
    device_flow = DeviceCodeFlowSimple(Path(token_path))
    try:
    # get access token from local path
        device_flow.get_access_token_refresh_if_necessary()
    except NoTokensError:
    # if no token was found locally, do the device code flow
        device_flow.do_flow()
    token = device_flow.access_token

    # request a file and provide the token in the Authorization header
    file_name = Path(url).name
    #print(file_name)
    directory_to_save_file = Path.cwd() # where you want to save the downloaded file 
    #print(directory_to_save_file )

    r = requests.get(url, headers={"authorization": f"Bearer {token}"})

    if r.status_code == requests.codes.ok:
        # save the file
        with open(Path(directory_to_save_file) / ('unavco_files')/ (file_name), 'wb') as f: #b for binary
            print(f)        
            for data in r:                                              
                f.write(data)
    else:
        #problem occured
        print(f"failure: {r.status_code}, {r.reason}")
    





