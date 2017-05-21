##
## Matrix for a tracer measuring Pk
##

from __future__ import division, print_function
import sys
import numpy as np
import scipy.linalg as la
from scipy.interpolate import interp1d
from FishMat import FishMat
from ParameterVec import DefaultParamList, Parameter
from ClassWrap import PkDiffer



class TracerPk(FishMat):

    def Pnoise(self,kpar,kperp,z):
        return 0

    def kNL(self,z):
        return 0.04+0.016*(1+z)**2.2

    def bias(self,z):
        return 1.

    def biaseta(self,z):
        return 1.

    def getInverseErrors(self):
        PkC=self.PkDiffer.cube0
        PkEI=[] #inverse errors
        for z,P,nm in zip(self.zvals,PkC,self.nmodes):
            Pn=self.Pnoise(self.kpar,self.kperp,z)
            PkE=(P+Pn)**2/nm
            knl=self.kNL(z)
            PkE[np.where(self.kt>knl)]=1e30
            PkEI.append(1/PkE)

        return PkEI
    
    def __init__ (self, zmin=2, zmax=6, dz=0.1, kmax=0.5, kzscal=0.1,dk=0.01,fsky=0.15):
        pl=DefaultParamList()
        ignorelist=['tau','As']
        N=len(pl)

        self.kvals=np.arange(dk/2,kmax,dk)
        self.Nk=len(self.kvals)
        self.dk=dk
        self.dz=dz
        self.fsky=fsky
        self.kpar=np.outer(self.kvals,np.ones(self.Nk))
        self.kperp=self.kpar.T
        self.kt=np.sqrt(self.kpar**2+self.kperp**2)
        self.mu=self.kpar/self.kt
        self.zvals=np.arange(zmin+dz/2,zmax,dz)
        for i,z in enumerate(self.zvals):
            pl.append (Parameter('b_delta_'+str(i),self.bias(z)))
            pl.append (Parameter('b_eta_'+str(i),self.biaseta(z)))
        Nwb=len(pl) # with bias parameters
        print("Setting up class...")
        self.PkDiffer=PkDiffer(pl,self.zvals, self.kvals, self.kperp, self.kpar)
        self.calcNModes()
        PkEI=self.getInverseErrors()

        
        eps=0.005
        Pderivs=[]
        print("Calculating derivatives... ", end='')
        for i1,p in enumerate(pl):
            print (" %s"%p.name, end='')
            sys.stdout.flush()
            if p.name not in ignorelist:
                Ders=self.PkDiffer.getDerivative(p,eps)
            else:
                Ders=None
            Pderivs.append(Ders)
        print("")
        F1=np.zeros((Nwb,Nwb))
        print ("Getting fisher: ",end='')
        for i1,D1 in enumerate(Pderivs):
            print (" %i"%i1,end='')
            sys.stdout.flush()
            if D1 is None:
                continue
            for i2,D2 in enumerate(Pderivs):
                if (i2<i1):
                    continue
                if D2 is None:
                    continue
                for zi,z in enumerate(self.zvals):
                    v=(D1[zi]*PkEI[zi]*D2[zi]).sum()
                    F1[i1,i2]+=v
                    F1[i2,i1]+=v
        print("")
        F1+=np.diag([1e-30]*Nwb)
        C=la.inv(F1)[:N,:N]
        F=la.inv(C)
        print (F1[:N,:N])
        print (F[:N,:N])
        FishMat.__init__(self,pl[:N],F)
        
    
    def calcNModes(self):
        self.nmodes=[]
        Da=self.PkDiffer.Da_fid
        for z in self.zvals:
            zlow=z-self.dz/2
            zhigh=z+self.dz/2
            V=self.fsky*4*np.pi/3*(Da(zhigh)**3-Da(zlow)**3)
            Vk=2*np.pi*self.kperp*self.dk*self.dk
            cnm=V*Vk/(2*(2*np.pi)**3)
            self.nmodes.append(cnm)
    
        
