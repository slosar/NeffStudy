#
# Matrix for a tracer measuring Pk.
#

from __future__ import division, print_function
import sys
import numpy as np
import scipy.linalg as la
from scipy.interpolate import interp1d
from FishMat import FishMat
from ParameterVec import DefaultParamList, Parameter
from ClassWrap import PkDiffer
import matplotlib.pyplot as plt

class TracerPk(FishMat):

    def Pnoise(self,kpar,kperp,z):
        return 0

    def kNL(self,z):
        """Returns the comoving wavenumber on the threshold
        of the linear and non-linear regimes.
        """
        return min(0.5,0.04+0.016*(1+z)**2.2)

    def bias(self,z): 
        return 1.
        
    def biaseta(self,z):
        return 1.

    def kmu2(self):
        """Returns a Nkmu2 x Nkmu2 matrix of (kt^i)(mu^2)^j values.
        """
        kmu2=[]
        for i in range(self.Nkmu2):
            for j in range(self.Nkmu2):
                km=(self.kt**i)+(self.mu**2.)**j
            kmu2.append(km)
        return kmu2

    def getInverseErrors(self):
        """Returns the inverse errors to be used in the Fisher matrix calculation.
        If SNR as a function of k_par and k_perp is specified for a given experiment,
        PkEI=1/\Delta P=SNR/P is returned.  
        If a specific experiment is not specified and hence no SNR is given,
        PkEI is calculated using the normal modes from calcNModes().
        """
        if (not hasattr(self,"nmodes")):
            self.calcNModes()
        PkC=self.PkDiffer.cube0
        PkEI=[] 
        if self.SNR=='None':
            for z,P,nm, in zip(self.zvals,PkC,self.nmodes):
                Pn=self.Pnoise(self.kpar,self.kperp,z)
                PkE=(P+Pn)**2/nm
                knl=self.kNL(z)
                PkE[np.where(self.kt>knl)]=1e30
                PkEI.append(1/PkE)
        else:
            for z,P,snr in zip(self.zvals, PkC, self.SNR):
                knl=self.kNL(z)
                PkE=P/snr
                PkE[np.where(self.kt>knl)]=1e30
                PkEI.append(1/PkE)
        return PkEI

    def calcNModes(self):
        self.nmodes=[]
        Da=self.PkDiffer.Da_fid
        for zlow,z,zhigh in zip(self.zmin,self.zvals,self.zmax):
            V=self.fsky*4*np.pi/3*(Da(zhigh)**3-Da(zlow)**3)
            Vk=2*np.pi*self.kperp*self.dk*self.dk
            cnm=V*Vk/(2*(2*np.pi)**3)
            self.nmodes.append(cnm)

    def __init__ (self, zvals, zmin, zmax, exp_name='None', SNR='None', kmax=0.5, dk=0.01, fsky=0.5, Nkmu2=3):
        pl=DefaultParamList()
        #ignorelist=['tau','As']
        ignorelist=[]
        N=len(pl)

        self.kvals=np.arange(dk,kmax+dk,dk)  
        self.Nk=len(self.kvals)
        self.dk=dk
        self.fsky=fsky
        self.kpar=np.outer(self.kvals,np.ones(self.Nk))
        self.kperp=self.kpar.T
        #plt.imshow(self.kperp)
        #plt.show()
        self.kt=np.sqrt(self.kpar**2+self.kperp**2)
        self.mu=self.kpar/self.kt
        self.Nkmu2=Nkmu2

        # Find the redshift values to be used in normal mode calculation, calcNModes()
        zmax=[]
        zmin=[]
        for i in range(len(zvals)):
            if i == 0:
                zlow=0.
                zhigh=2.*zvals[i]
            else:
                zlow=zmax[i-1]
                zhigh=zvals[i]+zlow-zmin[i-1]
            zmin.append(zlow)
            zmax.append(zhigh)
        self.zvals=zvals
        self.zmax=zmax
        self.zmin=zmin
        assert(zmax[1] == zmin[2])

        self.SNR=SNR[:,:self.Nk,:self.Nk]

        # Add bias parameters to parameter list
        for i,z in enumerate(self.zvals):
            pl.append(Parameter('b_delta_'+str(i), self.bias(z), ''))
            pl.append(Parameter('b_eta_'+str(i), self.biaseta(z), ''))

        # Add (kt^i)(mu^2)^j parameters to parameter list
        for i in range(self.Nkmu2):
            for j in range(self.Nkmu2):
                pl.append(Parameter('kmu2_'+str(i)+str(j), self.kmu2()[i][j],''))

        Nwb=len(pl) # with additional parameters

        # Calculate Fisher matrix
        print("Setting up class...")        
        self.PkDiffer=PkDiffer(pl,self.zvals, self.kvals, self.kperp, self.kpar, self.Nkmu2)
        PkEI=self.getInverseErrors()
        self.PkEI=PkEI
        eps=0.002
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
        self.Pderivs=Pderivs
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
                F1[i2,i1]=F1[i1,i2]
        F1+=np.diag([1e-30]*Nwb)
        C=la.inv(F1)[:N,:N]
#        plt.figure()
#        plt.imshow(np.log(F1),interpolation='nearest')
#        plt.colorbar()
#        plt.show()
        F=la.inv(C)
        print ('\n')
        print (F[:N,:N])
        FishMat.__init__(self, pl[:N], F)
        if exp_name != 'None':
            FishMat.saveF(self, F, exp_name)
