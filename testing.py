from DFT import He
import configparser
import numpy as np
#import hypothesis
from hypothesis import strategies as st
from hypothesis import given
from hypothesis import settings

NUCLEAR_CHARGE = 2
config = configparser.ConfigParser()
config.read("configuration.txt")
@given(SAMPLES = st.integers(10,int(config.get('settings', 'SAMPLES'))), R_MAX = st.floats(4,float(config.get('settings', 'R_MAX'))))
def test_He_init(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    #tests edge values for r
    assert(test_atom.r[-1]==R_MAX)
    assert(test_atom.r[0]==0)
    #tests correct array dimension of attributes
    assert(len(test_atom.r)==SAMPLES)
    assert(len(test_atom.u)==SAMPLES)
    assert(len(test_atom.rho)==SAMPLES)
    assert(len(test_atom.V_H)==SAMPLES)
    assert(len(test_atom.V_N)==SAMPLES)
    assert(len(test_atom.V_C)==SAMPLES)
    assert(len(test_atom.V_X)==SAMPLES)
    #tests correct sign of nuclear potential
    assert(np.all(test_atom.V_N<=0))
    
@given(SAMPLES = st.integers(10,int(config.get('settings', 'SAMPLES'))),
       R_MAX = st.floats(4,float(config.get('settings', 'R_MAX'))))
def test_compute_hartee_potential(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    #using a random WF
    u=np.random.rand(SAMPLES)
    probability_density=u**2
    probability = np.trapz(probability_density,test_atom.r) 
    #the probability of finding an electron has to be = 1
    probability_density = probability_density/probability 
    test_atom.u = probability_density**0.5
    test_atom.rho=2*probability_density
    test_atom.compute_hartree_potential(SAMPLES)
    assert(np.all(test_atom.V_H >= 0))
    #tests that the boundary condition is matched
    np.testing.assert_almost_equal(test_atom.V_H[-1]*test_atom.r[-1],NUCLEAR_CHARGE,decimal=7)
    
@given(SAMPLES = st.integers(10,int(config.get('settings', 'SAMPLES'))),
       R_MAX = st.floats(4,float(config.get('settings', 'R_MAX'))))
def test_compute_correlation_potential(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    #using a random WF
    u=np.random.rand(SAMPLES)
    probability_density=u**2
    probability = np.trapz(probability_density,test_atom.r) 
    #the probability of finding an electron has to be = 1
    probability_density = probability_density/probability 
    test_atom.u = probability_density**0.5
    test_atom.rho=2*probability_density
    test_atom.compute_correlation_potential(SAMPLES)
    #tests correct sign and no divergance
    assert(np.all(test_atom.V_C<=0))
    assert(test_atom.V_C[0]==0)
    
@given(SAMPLES = st.integers(10,int(config.get('settings', 'SAMPLES'))),
       R_MAX = st.floats(4,float(config.get('settings', 'R_MAX'))))
def test_compute_exchange_potential(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    #using a random WF
    u=np.random.rand(SAMPLES)
    probability_density=u**2
    probability = np.trapz(probability_density,test_atom.r) 
    #the probability of finding an electron has to be = 1
    probability_density = probability_density/probability 
    test_atom.u = probability_density**0.5
    test_atom.rho=2*probability_density
    test_atom.compute_exchange_potential()
    #tests correct sign and no divergance
    assert(np.all(test_atom.V_X<=0))
    assert(test_atom.V_X[0]==0)



    
    
    
