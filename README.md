# DFT simulation of the Helium atom 
This is a density functional theory program for the calculation of the ground state energy of the helium atom, using the Slater X-alpha exchange potential and the parametrization of the Ceperley-Alder correlation potential (Parametrisation given by Perdew and Zunger). In the DFT, the electronic orbitals are solution to the Schrödinger equation (SE) which depends on the electron density rather than on the individual electron orbitals. As such the DFT orbitals have no individual meaning but are used to construct the charge density. The main result of DFT is that there exists a form of the exchange correlation potential, depending only on the electron density ***n(r)***, that yields the exact ground state energy and density. The form of this functional is however unknown and we should relay on approximation, such as the local density approximation (LDA, in this example we will consider the CA-LDA parameterization).

Walter Kohn and Lu Jeu Scham have proposed to write the total functional ***E[n]*** as
***E[n] = T<sub>s</sub>[n] + E<sub>H</sub>[n] + E<sub>XC</sub>[n]***,

where the terms are respectively the non interacting system kinetic energy functional, the Hartree functional and the exchange functional.From ***δE[n]=0*** we obtain the total potential

<a href="https://imgbb.com/"><img src="https://i.ibb.co/6tGrpCd/Effective-potential.png" alt="Effective-potential" border="0" width = 400px height = auto></a>,

where ***V<sub>N</sub>*** is the culombic potential from the nucleus, the second term is the Hartree potential and the last term is the exchange potential. If ***E<sub>XC</sub>*** would be well known, then the self-consistent solution of the equation would give the exact ground-state density. Using the local density approximation (LDA) the potential is given by the sum of a exchange term and a correlation term: 
***V<sub>XC</sub><sup>LDA</sup> = V<sub>X</sub> + V<sub>C</sub>***. The Slater X potential was used for the Exchange potential while for the correlation interaction the Quantum Monte Carlo Parametrization from Ceperley-Adler was used.

## Structure of the program
Considering that the two electrons occupy the 1s-orbital, both the density and the Hartree potential are radially symmetric; thus, we can exploit this radial symmetry and solve the radial Schrödinger equation. The SE is solved directly using Verlet's algorithm. The program is organized in three steps:
1. Solution of the hydrogen-like radial SE using a simple integration algorithm combined with an interpolation routine in order to find the stationary states.
2. Calculation of the Hartree potential from the radial electron density.
3. Incorporation of the exchange-correlation potential (within the local density approximation, LDA).

<a href="https://ibb.co/474LnFd"><img src="https://i.ibb.co/zrPgTmJ/Immagine1.png" alt="Immagine1" border="0"></a>

## Solution of the radial SE
The radial symmetry is exploited in order to solve the radial Schrödinger equation on a radial vector ***r*** using the Verlet algorithm in order to integrate inward and find the electron energy ***E***

<a href="https://ibb.co/xjwWVT8"><img src="https://i.ibb.co/KqBdCvz/Verlet-SE.png" alt="Verlet-SE" border="0" width = 400px height = auto></a>

using the boundary conditions <a href="https://imgbb.com/"><img src="https://i.ibb.co/G5bp8Z9/boundary-SE.png" alt="boundary-SE" border="0" width = 260px height = auto></a>. The bisection method is used until the given precision on E is reached, checking the number of nodes in order to know if the correct energy is higher or lower with respect to the current one. 

## Computing the potentials
The nuclear potential is ***V<sub>N</sub> = n<sub>e</sub>/r***. If we define ***U''<sub>H</sub>(r) = -u(r)/r<sup>2</sup>***. We find ***U<sub>H</sub>(r)*** using the Verlet algorithm with the boundary conditions <a href="https://imgbb.com/"><img src="https://i.ibb.co/YdVQR26/boundary-Hartree.png" alt="boundary-Hartree" border="0" width = 260px height = auto></a>.

Within the LDA ***V<sub>X</sub> + V<sub>C</sub>*** depend locally on ***n(r)***. The Slater potential was used fo the exchange interaction:

<a href="https://imgbb.com/"><img src="https://i.ibb.co/HdDpL6K/Slater.png" alt="Slater" border="0" width = 180px height = auto></a>

The correlation term is implemented using the Ceperley Alder parametrization:

<a href="https://ibb.co/PjkSpNq"><img src="https://i.ibb.co/kBZdnxN/ec.png" alt="ec" border="0" width = 400px height = auto></a>

<a href="https://ibb.co/fpYyQyV"><img src="https://i.ibb.co/LtzcvcD/Vc.png" alt="Vc" border="0" width = 400px height = auto></a>

with ***r<sub>S</sub> = (4πn/3)<sup>1/3</sup>***.

The potential enegy for each potential is given by ***E<sub>i</sub> = ∫dr V<sub>i</sub>(r)u<sup>2</sup>(r)*** and the total energy of the system is therefore

***E<sub>tot</sub> = 2E - E<sub>H</sub> - E<sub>X</sub>/2 - E<sub>C</sub>/2***.

## Structure of the project
The main files of the program are:
1. [configuration.txt](https://github.com/andreapondini/DFT_He/blob/M/configuration.txt) : This file contains all the major parameters used in the simulation such as the precision  of the SE eignvalues (PREC_HSE) and the ground energy (PREC_DFT). Also the range (R_MAX) and the number of points (SAMPLES) in the radial space are specified, as well as the energy range (HSE_E_MIN) in which to look for the SE eignvalues. In addition to these parameters, the paths where the output files are saved are specified. If the user wants to modify these parameters, he has to change them in the text file.
2. [He.py](https://github.com/andreapondini/DFT_He/blob/M/He.py) : This file contains the main class of the programm. The He class has as attributes the radial WF (u) and density (rho), as well as the kinetic nergy of electrons (E_K), total energy (total_energy) the potentials (V_N,V_X,V_C,V_H,V). There are methods in order to compute the potentials as well as solving the SE and normalzing the WF. The main method is the hdft method which solves iteratively the SE updating the potentials as a new density is computed. Some methods in order to plot an save the results to files are present. The parameters used are taken from the [configuration](https://github.com/andreapondini/DFT_He/blob/M/configuration.txt) file.
3. [DFT.py](https://github.com/andreapondini/DFT_He/blob/M/DFT.py) : This is the script that actually runs the simulation. After instancing an object of the class [He](https://github.com/andreapondini/DFT_He/blob/M/He.py), the methods to perform the simulation and to save the results are called.
4. [testing.py](https://github.com/andreapondini/DFT_He/blob/M/testing.py) : contains all the tests for the methods of [He.py](https://github.com/andreapondini/DFT_He/blob/M/He.py) using hypothesis testing.

## Main results
<a href="https://ibb.co/nrYSXPV"><img src="https://i.ibb.co/K6fHtKQ/Density.png" alt="Density" border="0" width = 400px height = auto></a>
<a href="https://ibb.co/dPcM2sJ"><img src="https://i.ibb.co/Cbt87ps/Potentials.png" alt="Potentials" border="0" width = 400px height = auto></a>

The ground energy of the He atom using the  Ceperley-Alder LDA approximation and the Slater exchange potential is:

***E = -2.817 a.u***
