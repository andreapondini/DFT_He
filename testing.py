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
def test_initialization(R_MAX,SAMPLES):
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
def test_compute_hartee(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    u=np.random.rand(SAMPLES)
    probability_density=u**2
    probability = np.trapz(probability_density,test_atom.r) 
    #the probability of finding an electron has to be = 1
    probability_density = probability_density/probability 
    test_atom.u = probability_density**0.5
    test_atom.rho=2*probability_density
    U_H = np.zeros(SAMPLES)
    U_H[0] = 0
    U_H[1] = 0
    # outwards integration using Verlet algorithm
    step = (test_atom.r[-1] - test_atom.r[0]) / (SAMPLES-1)
    for i in range(1,SAMPLES-1):
        U_H[i+1] = 2*U_H[i] - U_H[i-1] - step**2 * test_atom.rho[i]/test_atom.r[i]
    # match boundary condition at r_max:
    # full charge of all electron within r_max
    alpha = (NUCLEAR_CHARGE - U_H[-1]) / test_atom.r[-1]
    U_H = U_H + alpha * test_atom.r
    #tests that the potential is positive
    assert(np.all(U_H >= 0))
    #tests that the boundary condition is matched
    np.testing.assert_almost_equal(U_H[-1],NUCLEAR_CHARGE,decimal=7)
    
#def test_compute_exchange(R_MAX,SAMPLES):
 #   test_atom = He(R_MAX,SAMPLES)
  #  assert(np.all(test_atom.V_X>=0))
    
    
    
