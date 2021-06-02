# -*- coding: utf-8 -*-
"""
Created on Fri May 21 23:05:23 2021

@author: Andrea
"""
import numpy as np
import matplotlib.pyplot as plt
SAMPLES, R_MAX, = 4049, 50
NUCLEAR_CHARGE = N_ELECTRONS = 2
PREC = 1e-2
HSE_E_MIN = -20
pi=np.pi


class atom:
    def __init__(self):
        self.r = np.linspace(0,R_MAX,SAMPLES) #radius vector
        self.E = 0  #system energy
        self.rho = self.V_H = self.V_X = self.V_C = self.V_N = np.zeros(SAMPLES)
        self.u = np.zeros(SAMPLES) #radial function for each point and each orbital
        self.V_N[1:] = NUCLEAR_CHARGE / self.r[1:]
        self.V_N[0]=0 #otherwise it diverges at 0


    def __compute_hartree_potential(self):
        #   Compute Hartree potential from solving the Poisson equation
        #   U_H''(r) = -rho(r) / r
        #   with the boundary conditions U_H(0) = 0, U_H(r_max) = n_electrons.
        U_H = np.zeros(SAMPLES)
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
        #every orbital is occipied with 2 electrons
        self.rho = self.rho + 2*self.u**2

    def __compute_orbitals(self):
        #L is always 0 because s orbital, N = 1
        E = self.__hse_solve(1, 0) #N,L
        return E

    def __hse_normalize(self): #to normalize radial u wavefunction
        R=self.r*self.u
        area = np.trapz(R,self.r)
        self.u = self.u / area

    def __hse_integrate(self, L, E_N):
          step = (self.r[-1] - self.r[0]) / (SAMPLES-1)
          #inward integration
          self.u[-1] = self.r[-1]*np.exp(-self.r[-1])
          self.u[-2]= self.r[-2]*np.exp(-self.r[-2])
          #integrate inward using Verlet algorithm
          for i in range(SAMPLES-2,0,-1):
              self.u[i-1] = 2*self.u[i] - self.u[i+1] + step**2*(-2*E_N + 2*self.V[i] + L*(L+1)/self.r[i]**2)*self.u[i]
    
    def __hse_solve(self,N,L): #solve SE
        E_N=0
        E_max =  0
        E_min = HSE_E_MIN
        while np.abs(E_max-E_min) > PREC:
            E_N = (E_min+E_max)/2
            print("E_N " , E_N)
            self.__hse_integrate(L,E_N)
            nodes = 0 #look for nodes
            for i in range(0,SAMPLES-1):
                if self.u[i]*self.u[i+1] < 0: nodes+=1
           # breakpoint()
            #continue search in the above or below half of energy range
            if (nodes > N-L-1): E_max = E_N
            else: E_min = E_N 
        self.__hse_normalize()
        return E_N
        

    def potential_energy(self, V): #computes the energy of given potential
        step = (self.r[-1] - self.r[0]) / (SAMPLES-1)
        E = np.sum(V*self.rho/2)
        return E * step    

    def __hdft(self):
        last_total_energy = 1
        total_energy = 0
        while(np.abs(last_total_energy-total_energy)>PREC):
            last_total_energy = total_energy
            self.__compute_hartree_potential()
            print("V_H = ", self.V_H)
            self.__compute_exchange_potential()
            self.__compute_correlation_potential()
            self.V = self.V_N + self.V_H + self.V_X + self.V_C  
            E = self.__compute_orbitals()
            self.__compute_density()
            #total energy of the 2 electrons + potential energy
            total_energy = 2*E - self.potential_energy(self.V_H) - self.potential_energy(self.V_X)/2 - self.potential_energy(self.V_C)/2
            print("total energy ",round(total_energy,2),"E ",round(E,2), "VH ",round(self.potential_energy(self.V_H),2),"VX ",round(self.potential_energy(self.V_X),2),"V_C ",round(self.potential_energy(self.V_C),2))
            plt.plot(self.r[0:100],self.u[0:100])
        return total_energy

prova = atom()
prova._atom__hdft()
plt.plot(prova.r[0:100],prova.u[0,100])