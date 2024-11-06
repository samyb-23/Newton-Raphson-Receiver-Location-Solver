# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 10:45:40 2024

@author: Samy

This script downloads the TLE files without overwriting previous ones in case
comparison between TLEs is needed, otherwise no need to change filename after
each download. 
"""

from skyfield.api import load
import datetime

# CHANGE FILENAME TO NOT OVERWRITE PREVIOUS FILE, which I do here by specifying
# the date and time of download
cur_time = datetime.datetime.now().strftime("%Y-%m-%d %H-%M-%S")

fname = 'gnss' + ' ' + cur_time + '.tle'

stations_url='http://celestrak.org/NORAD/elements/gp.php?GROUP=gnss&FORMAT=tle'

print('Downloading new TLE', fname)

load.download(stations_url, filename='TLEs' + '/' + fname)