# -*- coding: utf-8 -*-
"""
Created on Fri May 21 23:05:23 2021

@author: Andrea
"""
import configparser
from He import He

config = configparser.ConfigParser()
config.read("configuration.txt")

SAMPLES = config.get('settings', 'SAMPLES')
R_MAX = config.get('settings', 'R_MAX')
PREC_DFT = config.get('settings', 'PREC_DFT')
PREC_HSE = config.get('settings', 'PREC_HSE')
HSE_E_MIN = config.get('settings', 'HSE_E_MIN')

SAMPLES = int(SAMPLES)
R_MAX = float(R_MAX)    
PREC_DFT = float(PREC_DFT)
PREC_HSE = float(PREC_HSE)
HSE_E_MIN = int(HSE_E_MIN)

#after creating the object and launching the main method      
atom = He(R_MAX,SAMPLES)
atom.hdft(SAMPLES,PREC_DFT,PREC_HSE,HSE_E_MIN)
print("Total energy ",round(atom.total_energy,3)," a.u")
#the data is plotted and saved into a file
atom.save_data()
atom.plot_density()
atom.plot_potentials()
