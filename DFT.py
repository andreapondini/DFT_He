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
        self.r = np.linspace(0,R_MAX,SAMPLES) #radius vector
        self.E = 0  #system energy
        self.rho = self.V_H = self.V_X = self.V_C = self.V_N = np.zeros(SAMPLES)
        self.u = np.zeros((SAMPLES,int(N_ELECTRONS/2))) #radial function for each point and each orbital
        self.V_N[1:] = NUCLEAR_CHARGE / self.r[1:]
        self.V_N[0]=0 #otherwise it diverges at 0

    def __compute_hartree_potential(self):
        #   Compute Hartree potential from solving the Poisson equation
        #   U_H''(r) = -rho(r) / r
        #   with the boundary conditions U_H(0) = 0, U_H(r_max) = n_electrons.
        U_H = np.zeros(SAMPLES)
        # First boundary condition: no charge within radius 0
        U_H[0] = 0
        # alpha 0 assumed. r_max boundary condition matching done after integration
        U_H[1] = 0
        # outwards integration using Verlet algorithm
        #DO I = 2, SAMPLES-1, 1 
        #U_H(I+1) = 2*U_H(I) - U_H(I-1) - STEP**2 * RHO(I)/R(I)
        #END DO
        step = (self.r[-1] - self.r[0]) / (SAMPLES-1)
        for i in range(1,SAMPLES-1):
            U_H[i+1] = 2*U_H[i] - U_H[i-1] - step**2 * self.rho[i]/self.r[i]
        # match boundary condition at r_max:
        # full charge of all electron within r_max
        alpha = (N_ELECTRONS - U_H[-1]) / self.r[-1]
        U_H = U_H + alpha * self.r[i]
        #finally, get the Hartree potential from U_H
        self.V_H[1:] = U_H[1:] / self.r[1:]
        self.V_H[0] = 0
    def __compute_exchange_potential(self):
        #numpy does
        self.V_X[1:] = -np.cbrt((3*self.rho[1:]/(4*pi**2*self.r[1:]**2))) #cube root
        self.V_X[0] = 0
    def __compute_correlation_potential(self):
        #Compute the correlation potential according to Ceperly-Alder
        #parameterization for the spin unpolarized system.
        A , B, C, D, GAM, BETA1, BETA2 = 0.0311, -0.048, 0.002, -0.0116, -0.1423, 1.0529, 0.3334
        rs = np.where(self.rho>1e-10,np.cbrt(3 * self.r**2 / self.rho),0)
        self.V_C = np.where(np.logical_and(self.rho>1e-10,rs<1),A*np.log(rs) + B - A/3 + C*2/3*rs*np.log(rs) + (2*D-C)*rs/3, 0)
        self.V_C = np.where(np.logical_and(self.rho>1e-10,rs<1e10),GAM / (1 + BETA1*rs**0.5 + BETA2*rs) * (1+BETA1*7/6*rs**0.5+BETA2*4/3*rs) / (1+BETA1*rs**0.5+BETA2*rs),0)

    def __compute_density(self):
        for N in range(0, int(N_ELECTRONS/2)):
            #every orbital is occipied with 2 electrons
            self.rho = self.rho + 2*self.u[:,N]**2
        
    def potential_energy(self, V):
        step = (self.r[-1] - self.r[0]) / (SAMPLES-1)
        E = np.sum(V*self.rho/2)
        return E * step    