from DFT import He
import configparser
import numpy as np
import hypothesis
from hypothesis import strategies as st
from hypothesis import given
from hypothesis import settings
from datetime import timedelta

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
    assert(len(test_atom.V)==SAMPLES)
    #tests correct sign of nuclear potential
    assert(np.all(test_atom.V_N<=0))
    
@given(SAMPLES = st.integers(10,int(config.get('settings', 'SAMPLES'))),
       R_MAX = st.floats(4,float(config.get('settings', 'R_MAX'))))
def test_compute_hartee_potential(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    #using a random WF
    test_atom.u=np.random.rand(SAMPLES)
    test_atom.hse_normalize()
    test_atom.rho = 2*test_atom.u**2
    test_atom.compute_hartree_potential(SAMPLES)
    assert(np.all(test_atom.V_H >= 0))
    #tests that the boundary condition is matched
    np.testing.assert_almost_equal(test_atom.V_H[-1]*test_atom.r[-1],NUCLEAR_CHARGE,decimal=7)
    
@given(SAMPLES = st.integers(10,int(config.get('settings', 'SAMPLES'))),
       R_MAX = st.floats(4,float(config.get('settings', 'R_MAX'))))
def test_compute_correlation_potential(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    #using a random WF
    test_atom.u=np.random.rand(SAMPLES)
    test_atom.hse_normalize()
    test_atom.rho = 2*test_atom.u**2
    test_atom.compute_correlation_potential(SAMPLES)
    #tests correct sign and no divergance
    assert(np.all(test_atom.V_C<=0))
    assert(test_atom.V_C[0]==0)
    
@given(SAMPLES = st.integers(10,int(config.get('settings', 'SAMPLES'))),
       R_MAX = st.floats(4,float(config.get('settings', 'R_MAX'))))
def test_compute_exchange_potential(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    #using a random WF
    test_atom.u=np.random.rand(SAMPLES)
    test_atom.hse_normalize()
    test_atom.rho = 2*test_atom.u**2
    test_atom.compute_exchange_potential()
    #tests correct sign and no divergance
    assert(np.all(test_atom.V_X<=0))
    assert(test_atom.V_X[0]==0)
    
@given(SAMPLES = st.integers(10,int(config.get('settings', 'SAMPLES'))),
       R_MAX = st.floats(4,float(config.get('settings', 'R_MAX'))))
def test_hse_normalize(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    #using a random WF
    test_atom.u=np.random.rand(SAMPLES)
    test_atom.hse_normalize()
    test_atom.rho = 2*test_atom.u**2
    probability = np.trapz(test_atom.u**2,test_atom.r)
    #tests positivity of wf and density
    assert(np.all(test_atom.rho >= 0))
    assert(np.all(test_atom.u >= 0))
    #tests good normalization 
    np.testing.assert_almost_equal(probability,1,decimal=7)

@given(SAMPLES = st.integers(10,int(config.get('settings', 'SAMPLES'))),
       R_MAX = st.floats(4,float(config.get('settings', 'R_MAX'))),
       E_N = st.floats(float(config.get('settings', 'HSE_E_MIN')),0))
def test_hse_integrate(R_MAX,SAMPLES,E_N):
    test_atom = He(R_MAX,SAMPLES)
    #using a random WF
    test_atom.u=np.random.rand(SAMPLES)
    test_atom.hse_normalize()
    test_atom.rho = 2*test_atom.u**2
    #computing potentials for random WF
    test_atom.compute_correlation_potential(SAMPLES)
    test_atom.compute_exchange_potential()
    test_atom.compute_hartree_potential(SAMPLES)
    test_atom.V = test_atom.V_N + test_atom.V_H + test_atom.V_X + test_atom.V_C
    #obtain new wavefunction integrating using a random energy
    test_atom.hse_integrate(1,SAMPLES,E_N)
    #this wavefunction is not an eignvalue but still needs to have basic properties
    #verify boundary condition at R_MAX
    assert(test_atom.u[-2] == test_atom.r[-2]*np.exp(-test_atom.r[-2]))
    assert(test_atom.u[-1] == test_atom.r[-1]*np.exp(-test_atom.r[-1]))
    #check continuity, using maximum difference from consecutive values of R_MAX/SAMPLES
    assert(abs(test_atom.u[i]-test_atom.u[i+1])< R_MAX/SAMPLES for i in range(0,SAMPLES-1))

@given(SAMPLES = st.integers(10,int(config.get('settings', 'SAMPLES'))),
       R_MAX = st.floats(4,float(config.get('settings', 'R_MAX'))))
@settings(deadline=timedelta(seconds=1))
def test_hse_solve(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    #using a random WF
    test_atom.u=np.random.rand(SAMPLES)
    #test_atom.u = np.zeros(SAMPLES)
    test_atom.hse_normalize()
    test_atom.rho = 2*test_atom.u**2
    #computing potentials using random WF
    test_atom.compute_correlation_potential(SAMPLES)
    test_atom.compute_exchange_potential()
    test_atom.compute_hartree_potential(SAMPLES)
    test_atom.V = test_atom.V_N + test_atom.V_H + test_atom.V_X + test_atom.V_C
    #compute new wavefunction using potentials
    test_atom.E_k = test_atom.hse_solve(1,0,SAMPLES)
    #1s orbitals have no nodes, wf needs to be positive
    #(checked before normalizing)
    assert(np.all(test_atom.u[1:] > 0))
    #energy eignvalue needs to be negative and within range
    assert(test_atom.E_k < 0 and test_atom.E_k > float(config.get('settings','HSE_E_MIN')))

#V is a list with too many values in order to perform hypothesis testing
def test_potential_energy():
    R_MAX = float(config.get('settings', 'R_MAX'))
    SAMPLES = int(config.get('settings', 'SAMPLES'))
    test_atom = He(R_MAX,SAMPLES)
    test_atom.V=(np.random.rand(SAMPLES)-0.5)*10
    pot_energy = test_atom.potential_energy(test_atom.V)
    assert(type(pot_energy) == np.float64)
    
    
    
    
