# -*- coding: utf-8 -*-
"""
Created on Fri May 21 23:05:23 2021

@author: Andrea
"""
import numpy as np
import matplotlib.pyplot as plt
import configparser

config = configparser.ConfigParser()
config.read("configuration.txt")

NUCLEAR_CHARGE = 2
SAMPLES = config.get('settings', 'samples')
R_MAX = config.get('settings', 'r_max')
PREC_DFT = config.get('settings', 'prec_dft')
PREC_HSE = config.get('settings', 'prec_hse')
HSE_E_MIN = config.get('settings', 'hse_e_min')
plot1_path = config.get('paths', 'density_plot')
plot2_path = config.get('paths', 'potentials_plot')
pi = np.pi

SAMPLES = int(SAMPLES)
R_MAX = int(R_MAX)
PREC_DFT = float(PREC_DFT)
PREC_HSE = float(PREC_HSE)
HSE_E_MIN = int(HSE_E_MIN)

def hydrogen_like_wavefunc(x):
    """
    returns the wavefunction of He in the case
    where there's no interaction between the electrons
    """
    return x*np.exp(-NUCLEAR_CHARGE*x)/np.trapz((x*np.exp(-NUCLEAR_CHARGE*x))**2,x)**0.5
    
class He:
    def __init__(self):
        self.r = np.linspace(0,R_MAX,SAMPLES) #radial vector space
        self.E_k = 0  #kinetic energy
        self.total_energy = 0
        self.rho = np.zeros(SAMPLES)
        self.V_H = np.zeros(SAMPLES)
        self.V_X = np.zeros(SAMPLES)
        self.V_C = np.zeros(SAMPLES)
        self.V_N = np.zeros(SAMPLES)
        self.u = np.zeros(SAMPLES) #radial function
        self.V_N[1:] = - NUCLEAR_CHARGE / self.r[1:] #initializing nuclear potential
        self.V_N[0]=0 #otherwise it diverges at 0


    def __compute_hartree_potential(self):
        #   Compute Hartree potential from solving the Poisson equation
        #   U_H''(r) = -rho(r) / r
        #   with the boundary conditions U_H(0) = 0, U_H(r_max) = NUCLEAR_CHARGE.
        U_H = np.zeros(SAMPLES)
        U_H[0] = 0
        U_H[1] = 0
        # outwards integration using Verlet algorithm
        step = (self.r[-1] - self.r[0]) / (SAMPLES-1)
        for i in range(1,SAMPLES-1):
            U_H[i+1] = 2*U_H[i] - U_H[i-1] - step**2 * self.rho[i]/self.r[i]
        # match boundary condition at r_max:
        # full charge of all electron within r_max
        alpha = (NUCLEAR_CHARGE - U_H[-1]) / self.r[-1]
        U_H = U_H + alpha * self.r
        #get the Hartree potential from U_H
        self.V_H[1:] = U_H[1:] / self.r[1:]
        self.V_H[0] = 0
    def __compute_exchange_potential(self):
        self.V_X[1:] = -np.cbrt(3*self.rho[1:]/(4*pi**2*self.r[1:]**2)) #Slater Potential
        self.V_X[0] = 0 #otherwise it diverges at 0
    def __compute_correlation_potential(self):
        #Compute the correlation potential according to Ceperly-Alder
        #parameterization for the spin unpolarized system.
        A , B, C, D, GAM, BETA1, BETA2 = 0.0311, -0.048, 0.002, -0.0116, -0.1423, 1.0529, 0.3334
        for i in range(1,SAMPLES):
            if self.rho[i] < 1e-10: self.V_C[i] = 0
            else:
                rs = np.cbrt(3 * self.r[i]**2 / self.rho[i])
                if rs <1 : self.V_C[i] = A*np.log(rs) + B - A/3 + C*2/3*rs*np.log(rs) + (2*D-C)*rs/3
                elif rs < 1e10 : self.V_C[i] = GAM / (1 + BETA1*rs**0.5 + BETA2*rs) * (1+BETA1*7/6*rs**0.5+BETA2*4/3*rs) / (1+BETA1*rs**0.5+BETA2*rs)
                else: self.V_C[i]=0

    def __hse_normalize(self): #to normalize radial u wavefunction
        probability_density=self.u**2
        probability = np.trapz(probability_density,self.r) 
        #the probability of finding an electron has to be = 1
        probability_density = probability_density/probability 
        self.u = probability_density**0.5

    def __hse_integrate(self, L, E_N):
          step = (self.r[-1] - self.r[0]) / (SAMPLES-1)
          #setting boundary codition
          self.u[-1] = self.r[-1]*np.exp(-self.r[-1])
          self.u[-2]= self.r[-2]*np.exp(-self.r[-2])
          #integrate inward using Verlet algorithm
          for i in range(SAMPLES-2,0,-1):
              self.u[i-1] = 2*self.u[i] - self.u[i+1] + step**2*(-2*E_N + 2*self.V[i] + L*(L+1)/self.r[i]**2)*self.u[i]
    
    def __hse_solve(self,N,L): #solve SE
        E_N=0
        E_max =  0
        E_min = HSE_E_MIN
        while np.abs(E_max-E_min) > PREC_HSE:
            E_N = (E_min+E_max)/2
            self.__hse_integrate(L,E_N)
            nodes = 0 #look for nodes
            for i in range(0,SAMPLES-1):
                if self.u[i]*self.u[i+1] < 0: nodes+=1
            #continue search in the above or below half of energy range
            if (nodes > N-L-1): E_max = E_N
            else: E_min = E_N 
        self.__hse_normalize()
        return E_N
        

    def potential_energy(self, V): #computes the energy of given potential
        return np.trapz(V*self.u**2,self.r)

    def hdft(self):
        last_total_energy = 1
        self.total_energy = 0
        while(np.abs(last_total_energy-self.total_energy)>PREC_DFT):
            last_total_energy = self.total_energy
            #breakpoint()
            self.__compute_hartree_potential()
            self.__compute_exchange_potential()
            self.__compute_correlation_potential()
            self.V = self.V_N + self.V_H  + self.V_X + self.V_C 
            #L is always 0 because s orbital, N = 1
            self.E_k = self.__hse_solve(1, 0) #N,L
            self.rho = 2*self.u**2 #computes density, 2e- in 1s
            #total energy of the 2 electrons + potential energy
            self.total_energy = 2*self.E_k - self.potential_energy(self.V_H)  - self.potential_energy(self.V_X)/2 - self.potential_energy(self.V_C)/2 
        return self.total_energy
    
    def print_energy(self):
        print("Each electron kinetic energy ",round(self.E_k,3)," a.u")
        print("Total energy ",round(atom.total_energy,3)," a.u")
        
    def plot_density(self):
        """
        plots in the range [0:4] the electronic density computed
        the hydrogen-like result is plotted for comparison
        """
        fig, ax = plt.subplots(figsize=(7,3))
        ax.plot(atom.r[self.r<4],atom.rho[self.r<4],label='density')
        ax.plot(atom.r[self.r<4],2*hydrogen_like_wavefunc(atom.r[self.r<4])**2,label="hydrogen like density")
        ax.set(title='Wavefunction',xlabel='r [Å]')
        ax.legend(loc = 'upper right')
        ax.grid()
        fig.tight_layout()
        fig.savefig(plot1_path,dpi=100)
        
    def plot_potentials(self):
        fig, (ax1,ax2) = plt.subplots(2,1,figsize=(7,5))
        ax1.plot(atom.r[self.r<12],atom.V_H[self.r<12],label="V_H")
        ax1.plot(atom.r[self.r<12],atom.V_X[self.r<12],label="V_X")
        ax1.plot(atom.r[self.r<12],atom.V_C[self.r<12],label="V_C")
        ax1.set(title='Potentials',ylabel= 'Energy [a.u]')
        ax1.legend(loc = 'upper right')
        ax1.grid()
        ax2.set(xlabel='r [Å]',ylabel= 'Energy [a.u]')
        ax2.plot(atom.r[self.r<12],atom.V_C[self.r<12],label="V_C")
        ax2.legend(loc = 'lower right')
        ax2.grid()
        fig.tight_layout()
        fig.savefig(plot2_path,dpi=100)
        
atom = He()
atom.hdft()
atom.print_energy()
atom.plot_density()
atom.plot_potentials()

