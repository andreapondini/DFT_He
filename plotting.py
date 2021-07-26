import matplotlib.pyplot as plt
import configparser
import numpy as np

NUCLEAR_CHARGE = 2

config = configparser.ConfigParser()
config.read("configuration.txt")

density_plot_path = config.get('paths', 'density_plot')
potential_plot_path = config.get('paths', 'potentials_plot')
data_path = config.get('paths','data')

def hydrogen_like_wavefunc(x):
    """
    returns the wavefunction of He in the case
    where there's no interaction between the electrons
    """
    return x*np.exp(-NUCLEAR_CHARGE*x)/np.trapz((x*np.exp(-NUCLEAR_CHARGE*x))**2,x)**0.5


def plot_density(atom):
    """
   Plots in the range [0:4] the electronic density computed
   the hydrogen-like result is plotted for comparison,
   the save path is taken from the config file
   """
    fig, ax = plt.subplots(figsize=(7,3))
    ax.plot(atom.r[atom.r<4],atom.rho[atom.r<4],label='density')
    ax.plot(atom.r[atom.r<4],2*hydrogen_like_wavefunc(atom.r[atom.r<4])**2,label="hydrogen like density")
    ax.set(title='Density',xlabel='r [Å]')
    ax.legend(loc = 'upper right')
    ax.grid()
    fig.tight_layout()
    return fig
    
def plot_potentials(atom):
    """
    Plots in the range [0:12] the different potentials computed,
    the save path is taken from the config file
    """
    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(7,4))
    ax1.plot(atom.r[atom.r<12],atom.V_H[atom.r<12],label="V_H")
    ax1.plot(atom.r[atom.r<12],atom.V_X[atom.r<12],label="V_X")
    ax1.plot(atom.r[atom.r<12],atom.V_C[atom.r<12],label="V_C")
    ax1.set(title='Potentials',ylabel= 'Energy [a.u]')
    ax1.legend(loc = 'upper right')
    ax1.grid()
    ax2.set(xlabel='r [Å]',ylabel= 'Energy [a.u]')
    ax2.plot(atom.r[atom.r<12],atom.V_C[atom.r<12],label="V_C")
    ax2.legend(loc = 'lower right')
    ax2.grid()
    fig.tight_layout()
    return fig
    
        
def save_data(atom, path = data_path):
    """
    Saves the main results into a txt file, 
    the save path is taken from the config file
    """
    zipped = zip(atom.r,atom.u,atom.rho,atom.V)
    np.savetxt(data_path,list(zipped),fmt='%.5e',delimiter='\t',header="R [Å] \t single electron WF \t density \t V [a.u.]")

def save_plots(fig1, fig2, path1 = density_plot_path, path2 = potential_plot_path):
        fig1.savefig(path1,dpi=100)
        fig2.savefig(path2,dpi=100)