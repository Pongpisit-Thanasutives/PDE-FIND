import numpy as np
from scipy.io import loadmat
import matlab.engine

eng = matlab.engine.start_matlab()
npar = matlab.double
int32 = matlab.int32
