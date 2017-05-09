##
## Joel's fisher matrix
##

from __future__ import division, print_function
import numpy as np
import scipy.linalg as la
from FishMat import FishMat

class S4Fisher(FishMat):
    def __init__ (self):
        m=np.loadtxt('FishData/s4fisher05042017.txt')
        na=['tau','omegac','As','theta','Neff','mnu','omegab','ns']
        FishMat.__init__(self,na,m)

