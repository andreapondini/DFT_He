# -*- coding: utf-8 -*-
"""
Created on Fri May 21 23:05:23 2021

@author: Andrea
"""
import configparser
from He import He
from sys import argv

config = configparser.ConfigParser()
config.read(argv[1])

SAMPLES = config.get('settings', 'SAMPLES')
R_MAX = config.get('settings', 'R_MAX')

SAMPLES = int(SAMPLES)
R_MAX = float(R_MAX)    

#after creating the object and launching the main method      
atom = He(R_MAX,SAMPLES)
atom.hdft(SAMPLES)
print("Total energy ",round(atom.total_energy,3)," a.u")
#the data is plotted and saved into a file
atom.save_data()
atom.plot_density()
atom.plot_potentials()
