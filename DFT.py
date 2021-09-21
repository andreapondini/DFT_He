# -*- coding: utf-8 -*-
"""
Created on Fri May 21 23:05:23 2021

@author: Andrea
"""
import configparser
from He import He
from plotting import save_data, plot_potentials, plot_density, save_plots
from sys import argv

config = configparser.ConfigParser()
#if the configuration file is not specified, use "configuration.txt"
if len(argv)>=2:
    config.read(argv[1])
else:
    config.read("configuration.txt")

SAMPLES = config.get('settings', 'SAMPLES')
R_MAX = config.get('settings', 'R_MAX')
PREC_DFT = config.get('settings', 'PREC_DFT')
PREC_HSE = config.get('settings', 'PREC_HSE')
HSE_E_MIN = config.get('settings', 'HSE_E_MIN')
density_plot_path = config.get('paths', 'density_plot')
potential_plot_path = config.get('paths', 'potentials_plot')
data_path = config.get('paths','data')

SAMPLES = int(SAMPLES)
R_MAX = float(R_MAX)    
PREC_DFT = float(PREC_DFT)
PREC_HSE = float(PREC_HSE)
HSE_E_MIN = float(HSE_E_MIN)

#after creating the object and launching the main method      
atom = He(R_MAX, SAMPLES)
atom.hdft(PREC_DFT, PREC_HSE, HSE_E_MIN)
print("Total energy ",round(atom.total_energy, 3)," a.u")
#the data is plotted and saved into files
save_data(atom,data_path)
fig1 = plot_density(atom)
fig2 = plot_potentials(atom)
save_plots(fig1, fig2, density_plot_path, potential_plot_path)
