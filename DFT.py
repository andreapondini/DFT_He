# -*- coding: utf-8 -*-
"""
Created on Fri May 21 23:05:23 2021

@author: Andrea
"""
import numpy as np
SAMPLES, R_MAX, = 4049, 50
NUCLEAR_CHARGE = N_ELECTRONS = 2 #for He
#NUCLEAR_CHARGE = N_ELECTRONS = 4 #for Be
pi=np.pi


class atom:
    def __init__(self):
        self.r = np.linspace(0,R_MAX,SAMPLES)
        self.E = 0
        self.rho = self.V_H = self.V_X = self.V_C = self.V_N = np.zeros(SAMPLES)
        self.u = np.zeros((SAMPLES,int(N_ELECTRONS/2))) #radial function for each point and each orbital
        self.V_N[1:] = NUCLEAR_CHARGE / self.r[1:]
        self.V_N[0]=0
