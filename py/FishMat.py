##
## Super simple Fisher matrix class
##
from __future__ import division, print_function
import numpy as np
import scipy.linalg as la

class FishMat(object):

    def __init__ (self, plist, mat):
        self.plist=plist
        self.F=mat
        self.N=len(self.F)
        self.calcC()



    def calcC(self):
        self.C=la.inv(self.F+np.diag([1e-30]*self.N))
    
    def addF(self,F):
        if self.plist==F.plist:
            self.F+=F.F
            self.calcC()
        else:
            print("Params don't match")
            stop()

    def Error(self,pname):
        i=self.plist.index(pname)
        return np.sqrt(self.C[i,i])
    
