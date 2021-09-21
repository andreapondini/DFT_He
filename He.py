import numpy as np

pi = np.pi
NUCLEAR_CHARGE = 2

class He:
    def __init__(self,R_MAX,SAMPLES):
        """
        Parameters
        ----------
        R_MAX : float
            maximum value of r.
        SAMPLES : integer
            number of divisions of the radial space.
        -----------
        
        Initializes all the attributes and computes nuclear potential
        
        """
        self.SAMPLES = SAMPLES
        self.r = np.linspace(0,R_MAX,SAMPLES) #radial vector space
        self.E_k = 0  #kinetic energy
        self.total_energy = 0
        self.rho = np.zeros(SAMPLES) #density
        self.V_H = np.zeros(SAMPLES) #hartree potential
        self.V_X = np.zeros(SAMPLES) #exchange potential
        self.V_C = np.zeros(SAMPLES) #correlation potential
        self.V_N = np.zeros(SAMPLES) #nuclear potential
        self.V = np.zeros(SAMPLES) #total potential
        self.u = np.zeros(SAMPLES) #radial wavefunction
        self.V_N[1:] = - NUCLEAR_CHARGE / self.r[1:] #initializing nuclear potential
        self.V_N[0]=0 #otherwise it diverges at 0


    def compute_hartree_potential(self):
        """
        Computes Hartree potential from solving the Poisson equation
        U_H''(r) = -rho(r) / r
        with the boundary conditions U_H(0) = 0, U_H(r_max) = NUCLEAR_CHARGE    
        """
        
        U_H = np.zeros(self.SAMPLES)
        U_H[0] = 0
        U_H[1] = 0
        # outwards integration using Verlet algorithm
        step = (self.r[-1] - self.r[0]) / (self.SAMPLES-1)
        for i in range(1,self.SAMPLES-1):
            U_H[i+1] = 2*U_H[i] - U_H[i-1] - step**2 * self.rho[i]/self.r[i]
        # match boundary condition at r_max:
        # full charge of all electron within r_max
        alpha = (NUCLEAR_CHARGE - U_H[-1]) / self.r[-1]
        U_H = U_H + alpha * self.r
        #get the Hartree potential from U_H
        self.V_H[1:] = U_H[1:] / self.r[1:]
        self.V_H[0] = 0
        
    def compute_exchange_potential(self):
        """
        Computes exchange potential using Slater approx
        """
        self.V_X[1:] = -np.cbrt(3*self.rho[1:]/(4*pi**2*self.r[1:]**2)) #Slater Potential
        self.V_X[0] = 0 #otherwise it diverges at 0
        
    def compute_correlation_potential(self):
        """        
        Computes the correlation potential according to Ceperly-Alder
        parameterization for the spin unpolarized system.
        """
        A , B, C, D, GAM, BETA1, BETA2 = 0.0311, -0.048, 0.002, -0.0116, -0.1423, 1.0529, 0.3334
        for i in range(1,self.SAMPLES):
            if self.rho[i] < 1e-10: self.V_C[i] = 0
            else:
                rs = np.cbrt(3 * self.r[i]**2 / self.rho[i])
                if rs <1 : self.V_C[i] = A*np.log(rs) + B - A/3 + C*2/3*rs*np.log(rs) + (2*D-C)*rs/3
                elif rs < 1e10 : 
                    e_c = GAM / (1 + BETA1*rs**0.5 + BETA2*rs)
                    self.V_C[i] = e_c * (1+BETA1*7/6*rs**0.5+BETA2*4/3*rs) / (1+BETA1*rs**0.5+BETA2*rs)
                else: self.V_C[i]=0
    
    def compute_total_potential(self):
        """        
        Computes the potential using the dedicated methods and then sums them
        all together into the total potential V
        """
        self.compute_exchange_potential()
        self.compute_hartree_potential()
        self.compute_correlation_potential()
        self.V = self.V_N + self.V_H  + self.V_X + self.V_C 
        
    def hse_normalize(self): 
        """
        Normalizes radial u wavefunction
        in order to have the probability of finding a single electron = 1
        """
        probability_density=self.u**2
        probability = np.trapz(probability_density,self.r) 
        #the probability of finding an electron has to be = 1
        probability_density = probability_density/probability 
        self.u = probability_density**0.5

    def hse_integrate(self, angular_momentum, E_N):
        """
        Parameters
        ----------
        angular_momentum : integer
            angular momentum
            always = 0 for He
            
        E_N : float
            energy to use in the solution of SE
        -----------------
            
        Solves the SE and updates u using Verlet's algorithm
        
        """
        L = angular_momentum
        step = (self.r[-1] - self.r[0]) / (self.SAMPLES-1)
        #setting boundary codition
        self.u[-1] = self.r[-1]*np.exp(-self.r[-1])
        self.u[-2]= self.r[-2]*np.exp(-self.r[-2])
        #integrate inward using Verlet algorithm
        for i in range(self.SAMPLES-2,0,-1):
            self.u[i-1] = 2*self.u[i] - self.u[i+1] + step**2*(-2*E_N + 2*self.V[i] + L*(L+1)/self.r[i]**2)*self.u[i]
    
    def hse_solve(self,N,angular_momentum,PREC_HSE,HSE_E_MIN): 
        """
        Parameters
        ----------
        N : integer
            principal quantum number
            always = 1 for He
            
        angular_momentum : integer
            angular momentum
            always = 0 for He
            
        PREC_HSE : float
            precision on the energy eignvalue of the SE solution
            
        HSE_E_MIN = float
            minimum energy value above which the energy eignvalues
            that are solution of the SE are searched
        --------------
            
        Solves Schrodinger problem with precision PREC_HSE on energy eignvalue,
        the bisection method is used to find the eignvalue.
        N,L were left as variables as the formulas are more clear.
        
        Returns
        ---------
        E_N : float
            kinetic energy of single electron
        """
        L = angular_momentum
        E_N=0
        E_max =  0
        E_min = HSE_E_MIN
        while np.abs(E_max-E_min) > PREC_HSE:
            E_N = (E_min+E_max)/2  #bisection method
            self.hse_integrate(L,E_N) #solve SE
            nodes = 0 #look for nodes
            for i in range(0,self.SAMPLES-1):
                if self.u[i]*self.u[i+1] < 0: nodes+=1
            #continue search in the above or below half of energy range
            if (nodes > N-L-1): E_max = E_N
            else: E_min = E_N 
        self.hse_normalize() #normalize WF
        return E_N
        
    def potential_energy(self, V): 
        """
        Parameters
        ----------
        V : 1d numpy.array
            a general potential
        ---------------
            
        Computes the energy of given potential V

        Returns
        ---------
            np.float64
            potential energy of the given potential
        """
        return np.trapz(V*self.u**2,self.r)

    def hdft(self,PREC_DFT,PREC_HSE,HSE_E_MIN):
        """
        Parameters
        PREC_DFT : float
            precision on the total energy at which the simulation stops
            
        PREC_HSE : float
            precision on the energy eignvalue of the SE solution
            
        HSE_E_MIN = float
            minimum energy value above which the energy eignvalues
            that are solution of the SE are searched
        ----------
        
        Reiterates the solution of the SE,
        updating the potentials each time that a new denisty is obtained
        """
        last_total_energy = 1
        self.total_energy = 0
        while(np.abs(last_total_energy-self.total_energy)>PREC_DFT):
            last_total_energy = self.total_energy
            self.compute_total_potential()
            #L is always 0 because s orbital, N = 1
            self.E_k = self.hse_solve(1, 0,PREC_HSE,HSE_E_MIN) #N,L
            self.rho = 2*self.u**2 #computes density, 2e- in 1s
            #total energy of the 2 electrons + potential energy
            self.total_energy = 2*self.E_k - self.potential_energy(self.V_H)  - self.potential_energy(self.V_X)/2 - self.potential_energy(self.V_C)/2    
        
