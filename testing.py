from DFT import He
import configparser
import numpy as np
#import hypothesis
from hypothesis import strategies as st
from hypothesis import given
from hypothesis import settings

config = configparser.ConfigParser()
config.read("configuration.txt")
#@settings(max_examples = 1)
@given(SAMPLES = st.integers(1,int(config.get('settings', 'SAMPLES'))), R_MAX = st.floats(1,float(config.get('settings', 'R_MAX'))))
def test_initialization(R_MAX,SAMPLES):
    test_atom = He(R_MAX,SAMPLES)
    assert(test_atom.r[0]==0)
    assert(len(test_atom.r)==SAMPLES)
    assert(len(test_atom.u)==SAMPLES)
    assert(len(test_atom.rho)==SAMPLES)
    assert(len(test_atom.V_H)==SAMPLES)
    assert(len(test_atom.V_N)==SAMPLES)
    assert(len(test_atom.V_C)==SAMPLES)
    assert(len(test_atom.V_X)==SAMPLES)
    assert(np.all(test_atom.V_N<=0))
    
    
